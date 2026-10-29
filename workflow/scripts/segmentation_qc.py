"""Segmentation QC figure, drawn from segment.py's saved labels.

Left: membrane MIP on the segmentation grid with somata with a nucleus
(magenta), nucleus-free cells (cyan) and soma-sized outline components that
the nucleus-free filter rejected (grey). The MIP saturates at
saturation_percentile so dim somata stay visible. Right: the fraction of each
nucleus's core inside the membrane outline, with the min_nucleus_inside cut.

The outline is recomputed as in segment.py (common.membrane_outline); the
rejected components are the soma-sized ones holding no cell.
"""
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

from common import channel_index, membrane_outline, setup_logging

log = setup_logging(snakemake)  # noqa: F821

p = snakemake.params.segmentation
store = zarr.open(snakemake.input.image, mode="r")
block = tuple(p["downsample"])
vox = np.array(store.attrs["physical_pixel_sizes"]) * block
membrane = block_reduce(np.asarray(store[channel_index(store, p["membrane_emission_nm"])]),
                        block, np.mean).astype(np.float32)
soma = tifffile.imread(snakemake.input.soma)
nuclei = tifffile.imread(snakemake.input.nuclei)
cells = pd.read_csv(snakemake.input.table)
with_nucleus = cells.cell_id[cells.has_nucleus].to_numpy()
without = cells.cell_id[~cells.has_nucleus].to_numpy()

_, _, outline = membrane_outline(membrane, p, vox)
components, n = ndi.label(outline)
volume = np.bincount(components.ravel(), minlength=n + 1) * vox.prod()
min_soma, max_soma = p["soma_volume_um3"]
holds_cell = np.zeros(n + 1, bool)
holds_cell[np.unique(components[soma > 0])] = True
rejected_ids = np.flatnonzero((volume >= min_soma) & (volume <= max_soma) & ~holds_cell)
rejected_ids = rejected_ids[rejected_ids > 0]
rejected = np.isin(components, rejected_ids)

ids = np.unique(nuclei[nuclei > 0])
depth = ndi.distance_transform_edt(nuclei > 0, sampling=vox)
core = (depth >= p["nucleus_core_um"]) & (nuclei > 0)
inside = np.nan_to_num(ndi.mean(outline[core], nuclei[core], ids))

fig, (ax_img, ax_hist) = plt.subplots(1, 2, figsize=(15, 7), gridspec_kw={"width_ratios": [1.6, 1]})
mip = membrane.max(axis=0)
lo, hi = np.percentile(mip, [1, snakemake.params.saturation_percentile])
g = np.clip((mip - lo) / (hi - lo + 1e-9), 0, 1)
rgb = np.dstack([g * 0.25, g, g * 0.25])
rgb[find_boundaries(rejected.max(axis=0))] = 0.6
rgb[find_boundaries(np.where(np.isin(soma, with_nucleus), soma, 0).max(axis=0))] = [1, 0, 1]
rgb[find_boundaries(np.where(np.isin(soma, without), soma, 0).max(axis=0))] = [0, 0.8, 1]
ax_img.imshow(rgb, interpolation="nearest")
for _, row in cells.iterrows():
    ax_img.annotate(str(row.cell_id), (row.x_um / vox[2], row.y_um / vox[1]), color="w", fontsize=8,
                    xytext=(6, -6), textcoords="offset points")
ax_img.set_title(f"membrane MIP (saturated at {snakemake.params.saturation_percentile:g}th pct)\n"
                 f"{len(with_nucleus)} somata with nucleus (magenta), {len(without)} without (cyan), "
                 f"{len(rejected_ids)} rejected (grey)", fontsize=10)
ax_img.set(xticks=[], yticks=[])
ax_hist.hist(inside, bins=np.linspace(0, 1, 33), color="0.6")
ax_hist.axvline(p["min_nucleus_inside"], color="r")
ax_hist.set(yscale="log", xlabel="fraction of nucleus core inside membrane outline", ylabel="nuclei",
            title=f"{len(ids)} nuclei near membrane signal")
fig.suptitle(snakemake.wildcards.sample)
fig.tight_layout()
fig.savefig(snakemake.output[0], dpi=90)
log.info("%d with nucleus, %d without, %d rejected", len(with_nucleus), len(without), len(rejected_ids))
