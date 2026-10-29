"""Concatenate per-sample tables into the single analysis-ready table.

Only somata with a detected nucleus are counted. Nucleus-free cells stay in
each sample's cells.csv and counts.csv (has_nucleus = False) for inspection,
but most turned out to be dendrites, axons, trachea or fat globules.
"""
import pandas as pd

from common import setup_logging

log = setup_logging(snakemake)  # noqa: F821

frames = [pd.read_csv(path) for path in snakemake.input]
combined = pd.concat(frames, ignore_index=True)
combined = combined[combined["has_nucleus"]].copy()
combined["spots_per_100um3"] = (
    100 * combined["n_spots"] / combined["soma_volume_um3"].clip(lower=1)
)
combined.to_csv(snakemake.output[0], index=False)
log.info("%d cells across %d samples", len(combined), combined["sample"].nunique())
