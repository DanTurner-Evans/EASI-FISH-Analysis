"""Example figure of spot assignment to somata, in the example_plane region.

Left: membrane (green) and DAPI (blue) with the outlines of counted somata
(white) and their nuclei (cyan). Then one panel per probe: the raw probe with
its kept spots coloured by compartment, nuclear or cytoplasmic, and spots
outside any counted soma in grey. Outlines are drawn from the segmentation
grid in the middle projected plane (a maximum over planes merges somata stacked
in z); a spot is drawn if its centre lies in the projected planes.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
import zarr
from skimage.segmentation import find_boundaries

from common import channel_index, setup_logging

log = setup_logging(snakemake)  # noqa: F821

p = snakemake.params
ex = p.example
store = zarr.open(snakemake.input.image, mode="r")
vox = np.array(store.attrs["physical_pixel_sizes"])
z0, z1 = ex["planes"]
y0, x0 = ex["origin_px"]
size = ex["size_px"]
box = (slice(z0, z1 + 1), slice(y0, y0 + size), slice(x0, x0 + size))
block = np.array(p.downsample)

COMPARTMENT_COLOURS = {"nuclear": "#e87ba4", "cytoplasmic": "#eda100"}


def project(emission_nm):
    return np.asarray(store[channel_index(store, emission_nm)][box]).max(axis=0).astype(np.float32)


def stretch(img):
    lo, hi = np.percentile(img, [1, 99.8])
    return np.clip((img - lo) / (hi - lo + 1e-9), 0, 1)


def outline(labels):
    """Outline of a label image on the segmentation grid, cropped to the region
    at full resolution."""
    edge = find_boundaries(labels[(z0 + z1) // 2 // block[0]], mode="inner")
    full = np.repeat(np.repeat(edge, block[1], axis=0), block[2], axis=1)
    return full[y0:y0 + size, x0:x0 + size]


cells = pd.read_csv(snakemake.input.table)
counted = cells.cell_id[cells.has_nucleus & ~cells.cell_id.isin(p.exclude_cells)].to_numpy()
soma = tifffile.imread(snakemake.input.soma)
nuclei = tifffile.imread(snakemake.input.nuclei)
soma_edge = outline(np.where(np.isin(soma, counted), soma, 0))
nucleus_edge = outline(np.where(np.isin(nuclei, counted), nuclei, 0))

spots = pd.read_csv(snakemake.input.spots)
spots = spots[(spots.y >= y0) & (spots.y < y0 + size) & (spots.x >= x0) & (spots.x < x0 + size)
              & spots.z.round().between(z0, z1)]
spots.loc[~spots.cell_id.isin(counted), "compartment"] = "outside"

membrane = stretch(project(p.membrane_emission_nm))
dapi = stretch(project(p.nuclei_emission_nm))
overview = np.clip(membrane[..., None] * (0, 1, 0) + dapi[..., None] * (0.25, 0.25, 0.8), 0, 1)
overview[soma_edge] = 1
overview[nucleus_edge] = (0, 1, 1)

n = len(p.channels) + 1
fig, axes = plt.subplots(1, n, figsize=(5 * n, 5.6))
axes[0].imshow(overview)
axes[0].set_title("Membrane (green), DAPI (blue)\ncounted somata (white), their nuclei (cyan)")
bar = 10 / vox[2]
axes[0].plot([0.04 * size, 0.04 * size + bar], [0.94 * size] * 2, color="white", lw=6)
axes[0].text(0.04 * size + bar / 2, 0.91 * size, "10 µm", color="white", ha="center")
for ax, channel in zip(axes[1:], p.channels):
    img = stretch(project(p.emission_nm[channel]))
    rgb = np.dstack([img] * 3) * 0.8
    rgb[soma_edge] = 1
    rgb[nucleus_edge] = (0, 1, 1)
    ax.imshow(rgb)
    s = spots[spots.channel.astype(str) == channel]
    out = s.compartment == "outside"
    ax.scatter(s.x[out] - x0, s.y[out] - y0, s=40, facecolors="none", edgecolors="0.55", linewidths=0.8,
               label=f"outside counted somata ({out.sum()})")
    for compartment, colour in COMPARTMENT_COLOURS.items():
        m = s.compartment == compartment
        ax.scatter(s.x[m] - x0, s.y[m] - y0, s=70, facecolors="none", edgecolors=colour, linewidths=1.6,
                   label=f"{compartment} ({m.sum()})")
    ax.set_title(f"{channel} probe ({p.receptors[channel]}) with kept spots")
    ax.legend(loc="lower left", fontsize=9)
for ax in axes:
    ax.set(xticks=[], yticks=[])
fig.suptitle(f"{snakemake.wildcards.sample}: spots assigned to somata, max projection of planes {z0}-{z1} "
             f"({z0 * vox[0]:.1f}-{(z1 + 1) * vox[0]:.1f} µm deep)")
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(snakemake.output[0], dpi=110)
log.info("%d spots in region: %s", len(spots), spots.groupby(["channel", "compartment"]).size().to_dict())
