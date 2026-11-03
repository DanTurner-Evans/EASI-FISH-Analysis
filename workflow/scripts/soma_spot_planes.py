"""Plane-by-plane view of the spots assigned to somata, in the example_plane region.

One row per probe, one column per plane: the raw probe signal in that single
plane (contrast set over the whole slab, so planes compare), the outlines of the
counted somata in that plane (white), and a circle round each kept spot inside a
counted soma, coloured nuclear or cytoplasmic. The circle is solid in the plane
nearest the spot's centre and dashed in the other planes within
spot_half_depth_um of it, where the punctum is still in focus. Spots outside
counted somata are not drawn. Panels are cropped to the counted somata plus a
margin of crop_margin_um.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
import zarr
from matplotlib.lines import Line2D
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
planes = range(z0, z1 + 1)
block = np.array(p.downsample)
radius_px = ex["circle_radius_um"] / vox[2]

COMPARTMENT_COLOURS = {"nuclear": "#e87ba4", "cytoplasmic": "#eda100"}

cells = pd.read_csv(snakemake.input.table)
counted = cells.cell_id[cells.has_nucleus & ~cells.cell_id.isin(p.exclude_cells)].to_numpy()
soma = tifffile.imread(snakemake.input.soma)
soma = np.where(np.isin(soma, counted), soma, 0)


def outline(z):
    """Outline of the counted somata in plane z, cropped to the region at full resolution."""
    edge = find_boundaries(soma[z // block[0]], mode="inner")
    full = np.repeat(np.repeat(edge, block[1], axis=0), block[2], axis=1)
    return full[y0:y0 + size, x0:x0 + size]


spots = pd.read_csv(snakemake.input.spots)
half_depth = ex["spot_half_depth_um"] / vox[0]
spots = spots[(spots.y >= y0) & (spots.y < y0 + size) & (spots.x >= x0) & (spots.x < x0 + size)
              & spots.z.between(z0 - half_depth, z1 + half_depth) & spots.cell_id.isin(counted)]

edges = {z: outline(z) for z in planes}
ys, xs = np.nonzero(np.any(list(edges.values()), axis=0))
margin = ex["crop_margin_um"] / vox[2]
cy0, cy1 = int(max(ys.min() - margin, 0)), int(min(ys.max() + margin, size))
cx0, cx1 = int(max(xs.min() - margin, 0)), int(min(xs.max() + margin, size))
aspect = (cy1 - cy0) / (cx1 - cx0)
fig, axes = plt.subplots(len(p.channels), len(planes), figsize=(2.6 * len(planes), 2.6 * aspect * len(p.channels) + 1.2),
                         squeeze=False)
for row, channel in zip(axes, p.channels):
    stack = np.asarray(store[channel_index(store, p.emission_nm[channel])][z0:z1 + 1, y0:y0 + size,
                                                                            x0:x0 + size]).astype(np.float32)
    lo, hi = np.percentile(stack, [1, 99.8])
    s = spots[spots.channel.astype(str) == channel]
    for ax, z in zip(row, planes):
        rgb = np.dstack([np.clip((stack[z - z0] - lo) / (hi - lo + 1e-9), 0, 1)] * 3)
        rgb[edges[z]] = 1
        ax.imshow(rgb[cy0:cy1, cx0:cx1], extent=(cx0 - 0.5, cx1 - 0.5, cy1 - 0.5, cy0 - 0.5))
        here = s[(s.z - z).abs() <= half_depth]
        for spot in here.itertuples():
            centre = round(spot.z) == z
            ax.add_patch(plt.Circle((spot.x - x0, spot.y - y0), radius_px, fill=False,
                                    edgecolor=COMPARTMENT_COLOURS[spot.compartment],
                                    linewidth=1.4 if centre else 1.0, linestyle="-" if centre else "--"))
        ax.text(0.03, 0.97, f"{(here.z.round() == z).sum()} centred", transform=ax.transAxes, color="white",
                va="top", fontsize=8)
        if channel == p.channels[0]:
            ax.set_title(f"plane {z} ({z * vox[0]:.1f} µm)", fontsize=9)
        ax.set(xticks=[], yticks=[])
    row[0].set_ylabel(f"{channel} probe ({p.receptors[channel]})\n"
                      f"{len(s[s.z.round().between(z0, z1)])} spots in somata", fontsize=10)

bar = 5 / vox[2]
bx, by = cx0 + 0.05 * (cx1 - cx0), cy0 + 0.93 * (cy1 - cy0)
axes[0, 0].plot([bx, bx + bar], [by] * 2, color="white", lw=3)
axes[0, 0].text(bx + bar / 2, by - 0.03 * (cy1 - cy0), "5 µm", color="white", ha="center", va="bottom", fontsize=8)
handles = [Line2D([], [], color=c, marker="o", linestyle="", markerfacecolor="none", markersize=9) for c in
           COMPARTMENT_COLOURS.values()]
handles += [Line2D([], [], color="0.3", linestyle="-"), Line2D([], [], color="0.3", linestyle="--"),
            Line2D([], [], color="white", marker="s", linestyle="", markeredgecolor="0.3")]
fig.legend(handles, ["nuclear", "cytoplasmic", "spot centred in this plane",
                     f"spot centre within {ex['spot_half_depth_um']} µm", "counted soma outline (white)"],
           loc="lower center", ncol=5, frameon=False, fontsize=9)
fig.suptitle(f"{snakemake.wildcards.sample}: raw probe signal in each plane, with the spots assigned to somata")
fig.tight_layout(rect=(0, 0.06, 1, 0.95))
fig.savefig(snakemake.output[0], dpi=110)
log.info("%d spots in somata near the slab: %s", len(spots), spots.groupby(["channel", "compartment"]).size().to_dict())
