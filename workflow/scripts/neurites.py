"""Spot density in the driver line's neurites, against the tissue around them.

1. Neurite mask, at full resolution (processes are 0.2-1 um, below the
   segmentation grid): the membrane outline test from segment.py (global
   threshold and local contrast, same parameters), minus the counted somata
   (plus a margin) and minus fragments smaller than min_fragment_um3. Counted
   somata are those with a nucleus not listed in exclude_cells; nucleus-free
   cells stay in the mask, as on review most were dendrites. Membrane within
   nucleus_margin_um of any Cellpose nucleus is removed too: neuropil has no
   cell bodies, so this drops labelled somata the segmentation missed (dim
   ones in the soma layer), at the cost of the first micron or two of each
   primary neurite. The volume is processed in z slabs so deep stacks fit in
   memory.
2. Chance level: in neuropil, labelled processes are packed among unlabelled
   ones, so a spot on a neurite voxel could belong to a neighbour. The control
   is the density in a shell of tissue shell_um around the mask (no somata,
   nothing near a nucleus, no other membrane signal), computed on the
   segmentation grid, so neurites are compared with neuropil. Enrichment is
   neurite density / shell density; ~1 means neurites hold no more spots than
   the tissue they run through.

One neurite compartment per sample: neurites are not assigned to cells or
split by region.
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

from common import channel_index, setup_logging

log = setup_logging(snakemake)  # noqa: F821

p = snakemake.params
seg, np_ = p.segmentation, p.neurites
block = np.array(seg["downsample"])
store = zarr.open(snakemake.input.image, mode="r")
membrane = store[channel_index(store, seg["membrane_emission_nm"])]
vox = np.array(store.attrs["physical_pixel_sizes"])
vox_grid = vox * block
shape = np.array(membrane.shape)
sigma = np.array(np_["smooth_um"]) / vox
window = (2 * np.round(seg["local_window_um"] / vox) + 1).astype(int)
overlap = int(np.ceil(max(3 * sigma[0], window[0] // 2)))
slab = np_["slab_planes"]


def slabs():
    """(core slice, padded slice, offset of the core within the padded slab)."""
    for z in range(0, shape[0], slab):
        lo, hi = max(z - overlap, 0), min(z + slab + overlap, shape[0])
        yield slice(z, min(z + slab, shape[0])), slice(lo, hi), z - lo


def smoothed(padded):
    return ndi.gaussian_filter(np.asarray(membrane[padded]).astype(np.float32), sigma)


# --- 1. neurite mask ------------------------------------------------------------------
# Global statistics from every 7th smoothed voxel: the full float stack would be ~5 GB.
sample = np.concatenate([smoothed(padded)[off:off + core.stop - core.start].ravel()[::7]
                         for core, padded, off in slabs()])
median = np.median(sample)
mad = 1.4826 * np.median(np.abs(sample - median))
del sample

cells = pd.read_csv(snakemake.input.table)
counted = cells.cell_id[cells.has_nucleus & ~cells.cell_id.isin(p.exclude_cells)]
soma = np.isin(tifffile.imread(snakemake.input.soma), counted)
log.info("removing %d counted somata of %d cells", len(counted), len(cells))
r = np.floor(np_["soma_margin_um"] / vox_grid).astype(int)
zz, yy, xx = np.ogrid[-r[0]:r[0] + 1, -r[1]:r[1] + 1, -r[2]:r[2] + 1]
ball = (zz * vox_grid[0]) ** 2 + (yy * vox_grid[1]) ** 2 + (xx * vox_grid[2]) ** 2 <= np_["soma_margin_um"] ** 2
soma_zone = ndi.binary_dilation(soma, ball)
nuclei = tifffile.imread(snakemake.input.nuclei) > 0
near_nucleus = ndi.distance_transform_edt(~nuclei, sampling=vox_grid) <= np_["nucleus_margin_um"]
soma_zone |= near_nucleus
log.info("%.1f%% of the grid is within %s um of a nucleus", 100 * near_nucleus.mean(), np_["nucleus_margin_um"])
del nuclei, near_nucleus

foreground = np.zeros(shape, bool)
for core, padded, off in slabs():
    sm = smoothed(padded)
    local = ndi.maximum_filter(sm, size=tuple(window))
    fg = (sm > median + seg["membrane_threshold_mads"] * mad) & \
         (sm - median > seg["local_contrast"] * (local - median))
    foreground[core] = fg[off:off + core.stop - core.start]
in_soma = soma_zone[:, np.arange(shape[1]) // block[1]][:, :, np.arange(shape[2]) // block[2]]
neurite = foreground & ~in_soma
del in_soma
labels, n = ndi.label(neurite)
size = np.bincount(labels.ravel(), minlength=n + 1) * vox.prod()
keep = size >= np_["min_fragment_um3"]
keep[0] = False
neurite = keep[labels]
del labels
neurite_volume = neurite.sum() * vox.prod()
log.info("membrane median %.1f, MAD %.1f; neurite mask %.0f um3 in %d fragments >= %s um3",
         median, mad, neurite_volume, keep.sum(), np_["min_fragment_um3"])

# --- 2. control shell, on the segmentation grid --------------------------------------
neurite_grid = block_reduce(neurite, tuple(block), np.max)
foreground_grid = block_reduce(foreground, tuple(block), np.max)
del foreground
distance = ndi.distance_transform_edt(~neurite_grid, sampling=vox_grid)
inner, outer = np_["shell_um"]
shell = (distance > inner) & (distance <= outer) & ~soma_zone & ~foreground_grid
shell_volume = shell.sum() * vox_grid.prod()
tifffile.imwrite(snakemake.output.mask, (neurite_grid * 255 + shell * 80).astype(np.uint8), imagej=True,
                 resolution=(1 / vox_grid[2], 1 / vox_grid[1]),
                 metadata={"axes": "ZYX", "spacing": vox_grid[0], "unit": "um"})
log.info("shell %s-%s um: %.0f um3", inner, outer, shell_volume)

# --- counts ----------------------------------------------------------------------------
meta = p.meta
rows, spots_in = [], {}
for channel, path in zip(p.channels, snakemake.input.spots):
    spots = pd.read_csv(path)
    full = np.clip(np.round(spots[["z", "y", "x"]].to_numpy()).astype(int), 0, shape - 1)
    grid = np.minimum(full // block, np.array(shell.shape) - 1)
    on_neurite = neurite[tuple(full.T)]
    in_shell = shell[tuple(grid.T)] & ~on_neurite
    spots_in[channel] = spots[on_neurite]
    density = 100 * on_neurite.sum() / max(neurite_volume, 1e-9)
    background = 100 * in_shell.sum() / max(shell_volume, 1e-9)
    rows.append({
        "sample": snakemake.wildcards.sample,
        "line_label": meta["line_label"],
        "driver": meta["driver"],
        "neuron_type": meta["neuron_type"],
        "brain": meta["brain"],
        "receptor_set": meta["receptor_set"],
        "channel": channel,
        "receptor": meta["receptors"][channel],
        "neurite_volume_um3": neurite_volume,
        "n_spots_neurite": int(on_neurite.sum()),
        "neurite_per_100um3": density,
        "shell_volume_um3": shell_volume,
        "n_spots_shell": int(in_shell.sum()),
        "shell_per_100um3": background,
        "enrichment": density / background if background > 0 else np.nan,
    })
    log.info("%s: %d spots on neurites (%.3f per 100 um3), %d in shell (%.3f); enrichment %.2f",
             channel, on_neurite.sum(), density, in_shell.sum(), background, rows[-1]["enrichment"])
pd.DataFrame(rows).to_csv(snakemake.output.table, index=False)

# --- QC: membrane MIP with the mask, shell and the spots on neurites ------------------
mip = np.max([np.asarray(membrane[core]).max(axis=0) for core, _, _ in slabs()], axis=0)
mip = block_reduce(mip, (2, 2), np.max).astype(np.float32)
rgb = np.clip(np.dstack([mip] * 3) / np.percentile(mip, 99.8), 0, 1) * 0.8
up = lambda m: np.kron(m, np.ones((2, 2), bool))[:mip.shape[0], :mip.shape[1]]  # grid (4x) -> MIP (2x)
rgb[up(find_boundaries(shell.max(axis=0)))] = [0, 0.6, 0.8]
rgb[up(find_boundaries(neurite_grid.max(axis=0)))] = [1, 0, 1]
fig, ax = plt.subplots(figsize=(10, 10))
ax.imshow(rgb, interpolation="nearest")
for (channel, s), colour in zip(spots_in.items(), ["yellow", "orange"]):
    ax.scatter(s.x / 2, s.y / 2, s=14, facecolors="none", edgecolors=colour, linewidths=0.8,
               label=f"{channel} on neurites ({len(s)})")
ax.legend(loc="lower left")
ax.set(xticks=[], yticks=[],
       title=f"{snakemake.wildcards.sample}: neurite mask (magenta, {neurite_volume:.0f} µm³), "
             f"control shell (blue, {shell_volume:.0f} µm³)")
fig.tight_layout()
fig.savefig(snakemake.output.qc, dpi=80)
