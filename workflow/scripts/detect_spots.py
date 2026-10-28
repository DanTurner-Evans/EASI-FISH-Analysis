"""Propose mRNA puncta in one channel of one sample with Spotiflow.

Writes every candidate the model accepts, with features for the amplitude
filter in filter_spots.py. Keeping the two apart means re-tuning the filter
never re-runs this step, which takes minutes per channel on CPU (Spotiflow's
3D models do not run on Apple MPS).

The stack is processed in overlapping slabs of planes: on 1920 x 1920 planes
Spotiflow peaks at ~12.5 GB + 0.125 GB per plane, so a whole 365-plane stack
would need ~58 GB. Each slab keeps only the spots whose centre lies in its
core, so spots in the overlap are neither lost nor counted twice. Runtime is
~10 s per plane on CPU.
"""
import math

import numpy as np
import pandas as pd
import zarr
from scipy import ndimage as ndi
from skimage.measure import block_reduce
from spotiflow.model import Spotiflow

from common import channel_index, setup_logging

log = setup_logging(snakemake)  # noqa: F821

p = snakemake.params.detection
store = zarr.open(snakemake.input[0], mode="r")
idx = channel_index(store, snakemake.params.emission_nm)
voxel_um = np.array(store.attrs["physical_pixel_sizes"])
volume = store[idx]
n_planes = volume.shape[0]
log.info("channel %s -> index %d, shape %s, voxel %s um",
         snakemake.wildcards.channel, idx, volume.shape, voxel_um)

model = Spotiflow.from_pretrained(p["model"])
block = (1, 4, 4)  # background is estimated on a coarse grid since it varies slowly
window = np.maximum(1, np.round(np.array(p["background_window_um"]) / (voxel_um * block))).astype(int)
# Overlap must cover a spot's axial extent and the background window.
pad = max(math.ceil(p["slab_overlap_um"] / voxel_um[0]), int(window[0]))
core = p["slab_planes"]

slabs = []
for start in range(0, n_planes, core):
    lo, hi = max(start - pad, 0), min(start + core + pad, n_planes)
    slab = np.asarray(volume[lo:hi]).astype(np.float32)
    n_tiles = (1, *(max(1, math.ceil(n / 480)) for n in slab.shape[1:]))
    points, details = model.predict(
        slab, prob_thresh=p["prob_thresh"], n_tiles=n_tiles, device="cpu", verbose=False,
    )
    in_core = (points[:, 0] + lo >= start) & (points[:, 0] + lo < start + core)
    points, prob = points[in_core], details.prob[in_core]

    # Amplitude = peak minus local background, in raw camera counts, so one
    # threshold means the same thing in every sample.
    centre = np.clip(np.round(points).astype(int), 0, np.array(slab.shape) - 1)
    peak = ndi.uniform_filter(slab, size=3)[tuple(centre.T)]  # 3x3x3 mean: ~ spot core
    coarse_bg = ndi.percentile_filter(block_reduce(slab, block, func=np.median), 30, size=tuple(window))
    background = coarse_bg[tuple(np.minimum(centre // block, np.array(coarse_bg.shape) - 1).T)]

    points[:, 0] += lo
    slabs.append(pd.DataFrame(points, columns=["z", "y", "x"]).assign(
        prob=prob, peak=peak, background=background, amplitude=peak - background))
    log.info("planes %d-%d: %d candidates", start, min(start + core, n_planes) - 1, len(points))

spots = pd.concat(slabs, ignore_index=True)
for axis, size in zip("zyx", voxel_um):
    spots.insert(spots.columns.get_loc("prob"), f"{axis}_um", spots[axis] * size)
spots.to_csv(snakemake.output[0], index=False)
log.info("Spotiflow (%s): %d candidates; amplitude percentiles 5/50/95: %s", p["model"], len(spots),
         np.percentile(spots.amplitude, [5, 50, 95]).round() if len(spots) else "n/a")
