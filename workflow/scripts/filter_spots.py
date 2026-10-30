"""Keep candidate puncta above the sample's amplitude threshold and draw a QC figure.

The threshold is the valley between the noise and puncta peaks of this
sample's amplitude histogram (common.min_amplitude), unless overridden in the
config. For channels listed under spot_filter.ghosts, kept spots sitting just
above a brighter kept spot (a defocus ghost of it) are then dropped.

The figure shows the amplitude histogram with the threshold, and kept (green)
versus rejected (red) spots on a max projection of the densest region. The
noise population should sit left of the line and the puncta peak to its right.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import zarr
from scipy.spatial import cKDTree

from common import channel_index, min_amplitude, setup_logging

log = setup_logging(snakemake)  # noqa: F821

candidates = pd.read_csv(snakemake.input.candidates)
threshold, method = min_amplitude(candidates.amplitude.to_numpy(), snakemake.params.spot_filter,
                                  snakemake.params.override)
if method == "fallback":
    log.warning("no valley between noise and puncta peaks; using the fallback threshold %s", threshold)
keep = candidates.amplitude >= threshold
log.info("%d/%d candidates with amplitude >= %.1f (%s)", keep.sum(), len(candidates), threshold, method)


def ghosts_of(spots, xy_um, dz_um):
    """Spots with a brighter spot within xy_um laterally and dz_um deeper."""
    xy = spots[["y_um", "x_um"]].to_numpy()
    z, amplitude = spots.z_um.to_numpy(), spots.amplitude.to_numpy()
    lo, hi = dz_um
    ghost = np.zeros(len(spots), bool)
    for i, neighbours in enumerate(cKDTree(xy).query_ball_point(xy, xy_um)):
        ghost[i] = any(lo <= z[j] - z[i] <= hi and amplitude[j] > amplitude[i] for j in neighbours)
    return ghost


ghost = np.zeros(len(candidates), bool)
if snakemake.params.ghosts and keep.any():
    g = snakemake.params.ghosts
    ghost[np.flatnonzero(keep)] = ghosts_of(candidates[keep], g["xy_um"], g["dz_um"])
    keep &= ~ghost
    log.info("dropped %d defocus ghosts (brighter spot within %s um, %s um deeper)",
             ghost.sum(), g["xy_um"], g["dz_um"])
candidates[keep].to_csv(snakemake.output.spots, index=False)
log.info("kept %d/%d candidates", keep.sum(), len(candidates))

fig, (ax_hist, ax_img) = plt.subplots(1, 2, figsize=(14, 6.5))
bins = np.logspace(0, 5, 80)
ax_hist.hist(candidates.amplitude.clip(1, None), bins, color="0.6")
ax_hist.axvline(threshold, color="r")
ax_hist.set(xscale="log", xlabel="amplitude (counts above local background)", ylabel="candidates",
            title=f"{keep.sum()} of {len(candidates)} kept; threshold {threshold:.0f} ({method})"
                  + (f", {ghost.sum()} ghosts dropped" if ghost.any() else "")
                  + (": NO VALLEY, CHECK" if method == "fallback" else ""))

# Densest 400 x 400 px region, so the overlay shows the filter where it matters.
store = zarr.open(snakemake.input.image, mode="r")
volume = store[channel_index(store, snakemake.params.emission_nm)]
half = 200
if len(candidates):
    counts, ye, xe = np.histogram2d(candidates.y, candidates.x, bins=[volume.shape[1] // half] * 2,
                                    range=[[0, volume.shape[1]], [0, volume.shape[2]]])
    iy, ix = np.unravel_index(counts.argmax(), counts.shape)
    cy, cx = int((ye[iy] + ye[iy + 1]) / 2), int((xe[ix] + xe[ix + 1]) / 2)
else:
    cy, cx = volume.shape[1] // 2, volume.shape[2] // 2
y0, x0 = max(cy - half, 0), max(cx - half, 0)
mip = np.asarray(volume[:, y0:y0 + 2 * half, x0:x0 + 2 * half]).max(axis=0)
ax_img.imshow(mip, cmap="gray", vmin=np.percentile(mip, 1), vmax=np.percentile(mip, 99.8))
inside = candidates.y.between(y0, y0 + 2 * half) & candidates.x.between(x0, x0 + 2 * half)
for mask, colour in [(inside & keep, "lime"), (inside & ~keep & ~ghost, "red"), (inside & ghost, "cyan")]:
    ax_img.scatter(candidates.x[mask] - x0, candidates.y[mask] - y0, s=40,
                   facecolors="none", edgecolors=colour, linewidths=0.8)
ax_img.set(title=f"max projection, y {y0}, x {x0}: kept (green) / rejected (red) / ghost (cyan)",
           xticks=[], yticks=[])
fig.suptitle(f"{snakemake.wildcards.sample}  {snakemake.wildcards.channel}")
fig.tight_layout()
fig.savefig(snakemake.output.qc, dpi=80)
