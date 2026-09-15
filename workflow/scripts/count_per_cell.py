"""Assign detected puncta to labelled cells; emit one tidy row per cell/receptor."""
import numpy as np
import pandas as pd
import tifffile

from common import setup_logging

log = setup_logging(snakemake)  # noqa: F821

labels = tifffile.imread(snakemake.input.labels)
meta = snakemake.params.meta
rows = []

for channel, spots_path in zip(snakemake.params.channels, snakemake.input.spots):
    spots = pd.read_csv(spots_path)
    if spots.empty:
        log.warning("no spots in %s", spots_path)
        assigned = pd.Series(dtype=int)
    else:
        coords = spots[["z", "y", "x"]].round().astype(int).to_numpy()
        # Drop anything outside the volume before indexing.
        inside = np.all((coords >= 0) & (coords < np.array(labels.shape)), axis=1)
        dropped = (~inside).sum()
        if dropped:
            log.warning("%d spots outside volume bounds", dropped)
        cell_ids = labels[tuple(coords[inside].T)]
        # label 0 is background: puncta not inside any segmented soma.
        assigned = pd.Series(cell_ids[cell_ids > 0]).value_counts()
        log.info("%s: %d/%d spots fell inside a cell",
                 channel, int(assigned.sum()), len(spots))

    for cell_id in np.unique(labels[labels > 0]):
        rows.append({
            "sample": snakemake.wildcards.sample,
            "line_label": meta["line_label"],
            "driver": meta["driver"],
            "neuron_type": meta["neuron_type"],
            "brain": meta["brain"],
            "receptor_set": meta["receptor_set"],
            "channel": channel,
            "receptor": meta["receptors"][channel],
            "cell_id": int(cell_id),
            "cell_volume_vox": int((labels == cell_id).sum()),
            "n_spots": int(assigned.get(cell_id, 0)),
        })

pd.DataFrame(rows).to_csv(snakemake.output[0], index=False)
log.info("wrote %d rows", len(rows))
