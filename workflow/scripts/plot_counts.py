"""Per-cell receptor expression, one panel per receptor.

Driver line is encoded by x position, so color carries no extra information —
a single hue, no legend. Each point is one cell; the bar is the median.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import setup_logging

log = setup_logging(snakemake)  # noqa: F821

SERIES = "#2a78d6"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
SURFACE = "#fcfcfb"

df = pd.read_csv(snakemake.input[0])
receptors = sorted(df["receptor"].unique())
lines = sorted(df["line_label"].unique())
log.info("%d receptors x %d driver lines", len(receptors), len(lines))

fig, axes = plt.subplots(
    1, len(receptors), figsize=(2.1 * len(receptors), 3.4),
    sharey=True, facecolor=SURFACE,
)
axes = np.atleast_1d(axes)
rng = np.random.default_rng(0)

for ax, receptor in zip(axes, receptors):
    sub = df[df["receptor"] == receptor]
    ax.set_facecolor(SURFACE)
    for x, line_label in enumerate(lines):
        vals = sub.loc[sub["line_label"] == line_label, "spots_per_1k_vox"].to_numpy()
        if not len(vals):
            continue
        # Jitter so overlapping cells stay countable.
        ax.scatter(
            x + rng.uniform(-0.18, 0.18, len(vals)), vals,
            s=9, color=SERIES, alpha=0.55, linewidths=0,
        )
        ax.hlines(np.median(vals), x - 0.3, x + 0.3,
                  color=TEXT_PRIMARY, linewidth=2, zorder=3)

    ax.set_title(receptor, fontsize=9, color=TEXT_PRIMARY)
    ax.set_xticks(range(len(lines)))
    ax.set_xticklabels(lines, fontsize=8, color=TEXT_SECONDARY)
    ax.set_xlim(-0.6, len(lines) - 0.4)
    ax.tick_params(axis="y", labelsize=8, colors=TEXT_SECONDARY)
    # Recessive axes: keep the reader on the data.
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#d8d7d2")
    ax.grid(axis="y", color="#eceae5", linewidth=0.8)
    ax.set_axisbelow(True)

axes[0].set_ylabel("spots per 1000 voxels", fontsize=9, color=TEXT_SECONDARY)
fig.supxlabel("driver line", fontsize=9, color=TEXT_SECONDARY)
fig.tight_layout()
fig.savefig(snakemake.output[0], bbox_inches="tight", facecolor=SURFACE)
log.info("wrote %s", snakemake.output[0])
