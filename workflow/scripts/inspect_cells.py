"""Inspection stacks for checking the segmentation and spot filter in Fiji.

Written to data/processed/<sample>/inspect/:
  overview.tif         whole field on the segmentation grid (raw channels
                       block-averaged onto it; full resolution is ~2.5 GB)
  cell<id>.tif         full-resolution crop around each nucleus-free cell, the
                       labels upsampled nearest-neighbour
  cell<id>.png         z-planes through that cell: membrane, DAPI and a merge,
                       with its outline (yellow), other cells (cyan) and
                       Cellpose nuclei (red)

Both stacks hold, in order: membrane, DAPI, each probe, soma labels, nuclei
labels, soma outlines, then kept and rejected spot candidates for each probe
(rejected = below the sample's amplitude threshold). Nucleus-free cells get crops because, with
no nucleus to confirm them, they are the ones the shape filter decides.
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
import zarr
from scipy import ndimage as ndi
from skimage.measure import block_reduce
from skimage.segmentation import find_boundaries

from common import channel_index, min_amplitude, setup_logging

log = setup_logging(snakemake)  # noqa: F821

p = snakemake.params
seg = p.segmentation
block = np.array(seg["downsample"])
out = snakemake.output[0]
os.makedirs(out, exist_ok=True)

store = zarr.open(snakemake.input.image, mode="r")
vox_full = np.array(store.attrs["physical_pixel_sizes"])
vox = vox_full * block
raw = {"membrane": store[channel_index(store, seg["membrane_emission_nm"])],
       "DAPI": store[channel_index(store, seg["nuclei_emission_nm"])]}
for channel in p.channels:
    raw[channel] = store[channel_index(store, p.emission_nm[channel])]
soma = tifffile.imread(snakemake.input.soma)
nuclei = tifffile.imread(snakemake.input.nuclei)
cells = pd.read_csv(snakemake.input.table)
candidates = {channel: pd.read_csv(path) for channel, path in zip(p.channels, snakemake.input.candidates)}
threshold = {channel: min_amplitude(c.amplitude.to_numpy(), p.spot_filter, p.overrides.get(channel))[0]
             for channel, c in candidates.items()}
names = list(raw) + ["soma_labels", "nuclei_labels", "soma_outline"] + \
        [f"{channel}_{kind}" for channel in p.channels for kind in ("kept", "rejected")]


def spot_marks(shape, lo, scale):
    """One channel per probe and kind: a small disk at each spot centre that
    falls in the box starting at lo (full-res voxels) on a grid scaled by scale."""
    marks = []
    footprint = np.zeros((1, 5, 5), bool)
    footprint[0] = np.hypot(*np.mgrid[-2:3, -2:3]) <= 2
    for channel in p.channels:
        c = candidates[channel]
        kept = c.amplitude >= threshold[channel]
        for mask in (kept, ~kept):
            idx = np.floor((c.loc[mask, ["z", "y", "x"]].to_numpy() - lo) / scale).astype(int)
            idx = idx[np.all((idx >= 0) & (idx < shape), axis=1)]
            m = np.zeros(shape, bool)
            m[tuple(idx.T)] = True
            marks.append(ndi.binary_dilation(m, footprint) * 255)
    return marks


def save(path, channels, spacing):
    data = np.stack([c.astype(np.uint16) for c in channels], axis=1)  # ZCYX
    tifffile.imwrite(path, data, imagej=True, resolution=(1 / spacing[2], 1 / spacing[1]),
                     metadata={"axes": "ZCYX", "spacing": spacing[0], "unit": "um",
                               "mode": "composite", "Labels": names * len(data)})  # per image


# --- whole field on the segmentation grid ------------------------------------------
reduced = [block_reduce(np.asarray(a), tuple(block), np.mean) for a in raw.values()]
outline = find_boundaries(soma, mode="inner") * 255
save(f"{out}/overview.tif", reduced + [soma, nuclei, outline] + spot_marks(soma.shape, 0, block), vox)
log.info("overview %s, %d channels", soma.shape, len(names))


def stretch(img):
    lo, hi = np.percentile(img, [0.5, 99.8])
    return np.clip((img - lo) / (hi - lo + 1e-9), 0, 1)


# --- full-resolution crops around nucleus-free cells --------------------------------
half = np.array(p.crop_um) / 2
for _, c in cells[~cells.has_nucleus].iterrows():
    k = int(c.cell_id)
    centre = np.array([c.z_um, c.y_um, c.x_um])
    lo = np.maximum(((centre - half) / vox_full).astype(int), 0)
    hi = np.minimum(((centre + half) / vox_full).astype(int), raw["membrane"].shape)
    box = tuple(slice(a, b) for a, b in zip(lo, hi))
    crops = [np.asarray(a[box]).astype(np.float32) for a in raw.values()]
    grid = np.ix_(*[np.arange(a, b) // s for a, b, s in zip(lo, hi, block)])
    s, n = soma[grid], nuclei[grid]
    own = find_boundaries(s == k, mode="inner")
    others = find_boundaries(s, mode="inner") & ~own
    save(f"{out}/cell{k}.tif", crops + [s, n, own * 255 + others * 128] + spot_marks(s.shape, lo, 1),
         vox_full)

    planes = np.flatnonzero((s == k).any(axis=(1, 2)))
    planes = np.linspace(planes.min(), planes.max(), 5).round().astype(int)
    membrane, dapi = stretch(crops[0]), stretch(crops[1])
    fig, axes = plt.subplots(3, len(planes), figsize=(3 * len(planes), 9.4))
    for j, z in enumerate(planes):
        merge = np.dstack([dapi[z], membrane[z], dapi[z]])
        for i, (img, label) in enumerate([(membrane[z], "membrane"), (dapi[z], "DAPI"), (merge, "merge")]):
            ax = axes[i, j]
            ax.imshow(img, cmap="gray" if img.ndim == 2 else None)
            for mask, colour, width in [(own[z], "yellow", 0.8), (others[z], "cyan", 0.5),
                                        (n[z] > 0, "red", 0.6)]:
                if mask.any():
                    ax.contour(mask, [0.5], colors=colour, linewidths=width)
            ax.set(xticks=[], yticks=[])
            if i == 0:
                ax.set_title(f"z = {(lo[0] + z) * vox_full[0]:.1f} µm")
            if j == 0:
                ax.set_ylabel(label)
    axes[-1, 0].plot([5, 5 + 5 / vox_full[2]], [s.shape[1] - 8] * 2, color="white", lw=2)
    fig.suptitle(f"{snakemake.wildcards.sample} cell {k} (no nucleus), {c.soma_volume_um3:.0f} µm³; "
                 "yellow: this cell, cyan: other cells, red: Cellpose nuclei; bar 5 µm", fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{out}/cell{k}.png", dpi=110)
    plt.close(fig)
    log.info("cell %d: crop %s", k, s.shape)
