"""Scatter plot summarising the soma counts for each receptor and driver line:
spots per cell, or spots per soma volume (params.value names the column).

Each dot is one brain, summarised by its median cell so brains with many cells
do not dominate; the short line is the mean of those brain medians. Receptors are ordered by receptor set, then channel, as in config.yaml.
The table holds only counted somata (see aggregate.py).
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import setup_logging

log = setup_logging(snakemake)  # noqa: F821

# Reference categorical slots 1-2 (dataviz skill palette), fixed order by line.
SERIES = ["#2a78d6", "#eb6834"]
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID = "#eceae5"
AXIS = "#d8d7d2"
SURFACE = "#fcfcfb"

df = pd.read_csv(snakemake.input[0])
receptors = [r for r in snakemake.params.receptor_order if r in set(df.receptor)]
set_of = snakemake.params.receptor_set_of
lines = sorted(df.line_label.unique())
label = {r.line_label: f"{r.line_label}: {r.driver}, {r.neuron_type}"
         for r in df.drop_duplicates("line_label").itertuples()}

brains = df.groupby(["line_label", "receptor", "sample"])[snakemake.params.value].median().rename("median").reset_index()
log.info("%d cells, %d brains", len(df) // df.channel.nunique(), brains["sample"].nunique())

fig, ax = plt.subplots(figsize=(9, 4.2), facecolor=SURFACE)
ax.set_facecolor(SURFACE)
rng = np.random.default_rng(0)
for i, line in enumerate(lines):
    offset = (i - (len(lines) - 1) / 2) * 0.3
    for x, receptor in enumerate(receptors):
        vals = brains.loc[(brains.line_label == line) & (brains.receptor == receptor), "median"].to_numpy()
        if not len(vals):
            continue
        ax.scatter(x + offset + rng.uniform(-0.05, 0.05, len(vals)), vals, s=36, color=SERIES[i],
                   edgecolors=SURFACE, linewidths=1.5, zorder=3)
        ax.hlines(vals.mean(), x + offset - 0.11, x + offset + 0.11, color=TEXT_PRIMARY, linewidth=1.5, zorder=2)

ax.set_xticks(range(len(receptors)))
ax.set_xticklabels([f"{r}\nset {set_of[r]}" for r in receptors], fontsize=8, color=TEXT_PRIMARY)
ax.set_xlim(-0.6, len(receptors) - 0.4)
ax.set_ylim(0, None)
ax.set_ylabel(snakemake.params.ylabel, fontsize=9, color=TEXT_SECONDARY)
ax.tick_params(axis="y", labelsize=8, colors=TEXT_SECONDARY)
ax.tick_params(axis="x", length=0)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
for side in ("left", "bottom"):
    ax.spines[side].set_color(AXIS)
ax.grid(axis="y", color=GRID, linewidth=0.8)
ax.set_axisbelow(True)

handles = [plt.Line2D([], [], marker="o", linestyle="", color=SERIES[i], markeredgecolor=SURFACE,
                      markersize=7) for i in range(len(lines))]
names = [f"{label[line]} ({df[df.line_label == line]['sample'].nunique()} brains, "
         f"{len(df[df.line_label == line]) // df.channel.nunique()} somata)" for line in lines]
handles.append(plt.Line2D([], [], color=TEXT_PRIMARY, linewidth=2))
names.append("mean of brains")
ax.legend(handles, names, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, frameon=False, fontsize=8,
          labelcolor=TEXT_PRIMARY, handlelength=1.2, columnspacing=1.5, borderaxespad=0.2)
ax.set_title(f"{snakemake.params.title}: each dot is one brain's median soma", fontsize=10, color=TEXT_PRIMARY,
             loc="left", pad=24)
fig.tight_layout()
for path in snakemake.output:
    fig.savefig(path, dpi=200, bbox_inches="tight", facecolor=SURFACE)
log.info("wrote %s", ", ".join(snakemake.output))
