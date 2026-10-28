"""Assign puncta to segmented somata; emit one tidy row per cell/receptor.

Labels are on the downsampled segmentation grid, so spot coordinates are scaled
onto it. A spot inside a soma is nuclear if it also falls in that soma's
nucleus (same ID in nuclei_labels.tif), cytoplasmic otherwise; in a soma with
no detected nucleus it is unknown. Every spot is also written out with its
assignment, for inspection.
"""
import numpy as np
import pandas as pd
import tifffile

from common import setup_logging

log = setup_logging(snakemake)  # noqa: F821

soma = tifffile.imread(snakemake.input.soma)
nuclei = tifffile.imread(snakemake.input.nuclei)
cells = pd.read_csv(snakemake.input.table)
no_nucleus = cells.cell_id[~cells.has_nucleus].to_numpy()
block = np.array(snakemake.params.downsample)
meta = snakemake.params.meta
compartments = ["nuclear", "cytoplasmic", "unknown"]

assigned, rows = [], []
for channel, spots_path in zip(snakemake.params.channels, snakemake.input.spots):
    spots = pd.read_csv(spots_path)
    coords = np.floor(spots[["z", "y", "x"]].to_numpy() / block).astype(int)
    inside = np.all((coords >= 0) & (coords < np.array(soma.shape)), axis=1)
    if (~inside).any():
        log.warning("%s: %d spots outside the label grid", channel, (~inside).sum())
    cell = np.zeros(len(spots), int)
    cell[inside] = soma[tuple(coords[inside].T)]
    nuclear = np.zeros(len(spots), bool)
    nuclear[inside] = nuclei[tuple(coords[inside].T)] == cell[inside]
    spots["channel"] = channel
    spots["cell_id"] = cell
    spots["compartment"] = np.select([cell == 0, np.isin(cell, no_nucleus), nuclear],
                                     ["outside", "unknown", "nuclear"], "cytoplasmic")
    assigned.append(spots)
    log.info("%s: %d/%d spots in a soma (%s)", channel, (cell > 0).sum(), len(spots),
             spots.compartment.value_counts().to_dict())

    counts = (spots[spots.cell_id > 0].groupby(["cell_id", "compartment"]).size()
              .unstack(fill_value=0).reindex(columns=compartments, fill_value=0))
    for _, c in cells.iterrows():
        n = counts.loc[c.cell_id] if c.cell_id in counts.index else pd.Series(0, index=compartments)
        rows.append({
            "sample": snakemake.wildcards.sample,
            "line_label": meta["line_label"],
            "driver": meta["driver"],
            "neuron_type": meta["neuron_type"],
            "brain": meta["brain"],
            "receptor_set": meta["receptor_set"],
            "channel": channel,
            "receptor": meta["receptors"][channel],
            "cell_id": int(c.cell_id),
            "has_nucleus": bool(c.has_nucleus),
            "soma_volume_um3": c.soma_volume_um3,
            "nucleus_volume_um3": c.nucleus_volume_um3,
            "touches_edge": bool(c.touches_edge),
            "n_spots": int(n.sum()),
            "n_nuclear": int(n["nuclear"]),
            "n_cytoplasmic": int(n["cytoplasmic"]),
            "n_unknown": int(n["unknown"]),
        })

pd.concat(assigned, ignore_index=True).to_csv(snakemake.output.spots, index=False)
pd.DataFrame(rows).to_csv(snakemake.output.counts, index=False)
log.info("wrote %d rows for %d cells (%d without a nucleus)", len(rows), len(cells), len(no_nucleus))
