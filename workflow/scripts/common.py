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
