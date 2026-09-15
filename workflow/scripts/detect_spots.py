"""Detect mRNA puncta in one channel of one sample.

PLACEHOLDER THRESHOLDS. Laplacian-of-Gaussian blob detection is the right family
of method for EASI-FISH puncta, but min_sigma/max_sigma/threshold in
config.yaml are guesses. Tune them on a sample with known expression and check
the detected count by eye before running the full set.
"""
import numpy as np
import pandas as pd
import zarr
from skimage.feature import blob_log

from common import setup_logging

log = setup_logging(snakemake)  # noqa: F821

store = zarr.open(snakemake.input[0], mode="r")
channel_names = list(store.attrs["channel_names"])
wavelength = snakemake.wildcards.channel

# Match the requested wavelength to a channel index by name.
matches = [i for i, name in enumerate(channel_names) if wavelength in str(name)]
if len(matches) != 1:
    raise ValueError(
        f"Expected exactly one channel matching {wavelength!r}, "
        f"found {matches} in {channel_names}"
    )
volume = np.asarray(store[matches[0]]).astype(np.float32)
log.info("channel %r -> index %d, shape %s", wavelength, matches[0], volume.shape)

# Normalise so the threshold means the same thing across samples.
lo, hi = np.percentile(volume, [1, 99.9])
volume = np.clip((volume - lo) / max(hi - lo, 1e-6), 0, 1)

p = snakemake.params.detection
blobs = blob_log(
    volume,
    min_sigma=p["min_sigma"],
    max_sigma=p["max_sigma"],
    threshold=p["threshold"],
)
log.info("detected %d spots", len(blobs))

pd.DataFrame(blobs, columns=["z", "y", "x", "sigma_z", "sigma_y", "sigma_x"]).to_csv(
    snakemake.output[0], index=False
)
