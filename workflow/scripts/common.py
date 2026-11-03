"""Helpers shared by the workflow scripts.

Snakemake puts the script's directory on sys.path, so these import normally.
"""
import logging
import sys

import numpy as np
from scipy import ndimage as ndi
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

    Noise is modelled as a Gaussian on the smoothed log10-amplitude histogram:
    centred on the noise peak (the first peak), with the width of the peak's
    left half at half height (the right half is overlaid by puncta). The
    threshold is the first amplitude above the peak where that Gaussian is
    below max_noise_fraction of the histogram. Puncta often form only a
    shoulder on the noise peak, with no valley between them.
    Returns (threshold, method) with method "gaussian", "override" or
    "fallback" (too few candidates, or the noise never falls that low).
    """
    if override is not None:
        return float(override), "override"
    x = np.log10(amplitudes[amplitudes >= 2])  # a spike of candidates sits at ~1
    if len(x) < 10:
        return float(spot_filter["fallback"]), "fallback"
    grid = np.linspace(np.log10(2), x.max(), 500)
    density = gaussian_kde(x, bw_method=spot_filter["smoothing"])(grid)
    peaks, _ = find_peaks(density, prominence=spot_filter["min_prominence"] * density.max())
    noise = peaks[0]
    half = np.flatnonzero(density[:noise + 1] >= density[noise] / 2)[0]
    sigma = (grid[noise] - grid[half]) / np.sqrt(2 * np.log(2))
    model = density[noise] * np.exp(-0.5 * ((grid - grid[noise]) / sigma) ** 2)
    above = np.flatnonzero(model[noise:] < spot_filter["max_noise_fraction"] * density[noise:])
    if not len(above):
        return float(spot_filter["fallback"]), "fallback"
    return float(10 ** grid[noise + above[0]]), "gaussian"


def grid_ball(radius_um, vox):
    """Ellipsoidal footprint of the given physical radius on an anisotropic grid."""
    r = np.floor(radius_um / vox).astype(int)
    z, y, x = np.ogrid[-r[0]:r[0] + 1, -r[1]:r[1] + 1, -r[2]:r[2] + 1]
    return (z * vox[0]) ** 2 + (y * vox[1]) ** 2 + (x * vox[2]) ** 2 <= radius_um ** 2


def membrane_outline(membrane, p, vox):
    """Step 1 of segment.py, shared with its QC figure: (smoothed membrane,
    foreground, outline). Foreground is above a global threshold and bright
    relative to the brightest membrane nearby; the outline is the foreground
    closed, filled and opened to strip neurites."""
    smooth = ndi.gaussian_filter(membrane, 0.5)
    median = np.median(smooth)
    mad = 1.4826 * np.median(np.abs(smooth - median))
    local_max = ndi.maximum_filter(smooth, footprint=grid_ball(p["local_window_um"], vox))
    fg = (smooth > median + p["membrane_threshold_mads"] * mad) & \
         (smooth - median > p["local_contrast"] * (local_max - median))
    closed = ndi.binary_closing(fg, grid_ball(p["close_um"], vox), border_value=0)
    filled = ndi.binary_fill_holes(closed) | np.stack([ndi.binary_fill_holes(plane) for plane in closed])
    return smooth, fg, ndi.binary_opening(filled, grid_ball(p["neurite_um"], vox))
