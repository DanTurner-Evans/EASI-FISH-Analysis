"""Concatenate per-sample tables into the single analysis-ready table.

Only somata with a detected nucleus are counted. Nucleus-free cells stay in
each sample's cells.csv and counts.csv (has_nucleus = False) for inspection,
but most turned out to be dendrites, axons, trachea or fat globules. Cells
listed under exclude_cells in the config (found not to be neurons on review)
are dropped too.
"""
import pandas as pd

from common import setup_logging

log = setup_logging(snakemake)  # noqa: F821

frames = [pd.read_csv(path) for path in snakemake.input]
combined = pd.concat(frames, ignore_index=True)
combined = combined[combined["has_nucleus"]].copy()
for sample, ids in snakemake.params.exclude_cells.items():
    present = set(combined.loc[combined["sample"] == sample, "cell_id"])
    missing = set(ids) - present
    if missing:  # a re-segmentation may have renumbered the cells
        log.warning("%s: excluded cells %s not among its counted cells", sample, sorted(missing))
    combined = combined[~((combined["sample"] == sample) & combined["cell_id"].isin(ids))]
    log.info("%s: excluded cells %s", sample, sorted(set(ids) & present))
combined["spots_per_100um3"] = (
    100 * combined["n_spots"] / combined["soma_volume_um3"].clip(lower=1)
)
combined.to_csv(snakemake.output[0], index=False)
log.info("%d cells across %d samples", len(combined), combined["sample"].nunique())
