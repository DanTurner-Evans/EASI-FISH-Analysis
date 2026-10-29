"""Helpers shared by the workflow scripts.

Snakemake puts the script's directory on sys.path, so these import normally.
"""
import logging
import sys

import numpy as np
from scipy.signal import find_peaks
from scipy.stats import gaussian_kde


def setup_logging(snakemake):
    """Send stdout/stderr and logging output to the rule's log file."""
    handler = logging.FileHandler(snakemake.log[0], mode="w")
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    root = logging.getLogger()
    root.setLevel(logging.INFO)
    root.addHandler(handler)
    sys.stderr = sys.stdout = handler.stream
    return root


def channel_index(store, emission_nm, tolerance_nm=15):
    """Index of the zarr channel whose recorded emission wavelength matches."""
    recorded = store.attrs["emission_nm"]
    matches = [i for i, em in enumerate(recorded) if abs(em - emission_nm) <= tolerance_nm]
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one channel with emission {emission_nm} nm "
            f"(+/- {tolerance_nm}), found {matches} in {recorded}"
        )
    return matches[0]


def min_amplitude(amplitudes, spot_filter, override=None):
    """Amplitude threshold for one sample and channel, and how it was set.

    The lowest point of the smoothed log10-amplitude histogram between the
    noise peak (the first peak) and the puncta peak (the tallest after it).
    Returns (threshold, method) with method "valley", "override" or "fallback"
    (no second peak, so no valley).
    """
    if override is not None:
        return float(override), "override"
    x = np.log10(amplitudes[amplitudes >= 2])  # a spike of candidates sits at ~1
    if len(x) < 10:
        return float(spot_filter["fallback"]), "fallback"
    grid = np.linspace(np.log10(2), x.max(), 500)
    density = gaussian_kde(x, bw_method=spot_filter["smoothing"])(grid)
    peaks, _ = find_peaks(density, prominence=spot_filter["min_prominence"] * density.max())
    if len(peaks) < 2:
        return float(spot_filter["fallback"]), "fallback"
    noise, puncta = peaks[0], peaks[1:][np.argmax(density[peaks[1:]])]
    valley = noise + np.argmin(density[noise:puncta + 1])
    return float(10 ** grid[valley]), "valley"
