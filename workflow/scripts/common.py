"""Helpers shared by the workflow scripts.

Snakemake puts the script's directory on sys.path, so these import normally.
"""
import logging
import sys


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
