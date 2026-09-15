"""Concatenate per-sample tables into the single analysis-ready table."""
import pandas as pd

from common import setup_logging

log = setup_logging(snakemake)  # noqa: F821

frames = [pd.read_csv(path) for path in snakemake.input]
combined = pd.concat(frames, ignore_index=True)
combined["spots_per_1k_vox"] = (
    1000 * combined["n_spots"] / combined["cell_volume_vox"].clip(lower=1)
)
combined.to_csv(snakemake.output[0], index=False)
log.info("%d cells across %d samples", len(combined), combined["sample"].nunique())
