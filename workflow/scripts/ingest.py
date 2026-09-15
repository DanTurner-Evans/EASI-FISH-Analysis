"""Convert one CZI to chunked zarr.

Reading a CZI is slow and single-threaded; every downstream rule reads the zarr
instead so that cost is paid exactly once per sample.
"""
import numpy as np
import zarr
from bioio import BioImage

from common import setup_logging

log = setup_logging(snakemake)  # noqa: F821  (injected by Snakemake)

img = BioImage(snakemake.input[0])
log.info("dims=%s  channels=%s", img.dims, img.channel_names)

# bioio returns TCZYX; this data is a single timepoint z-stack.
data = img.get_image_data("CZYX", T=0)
log.info("array shape=%s dtype=%s", data.shape, data.dtype)

store = zarr.open(
    snakemake.output[0],
    mode="w",
    shape=data.shape,
    dtype=data.dtype,
    # Chunk whole z-planes per channel: matches how segmentation and spot
    # detection walk the stack.
    chunks=(1, 1, *data.shape[-2:]),
)
store[:] = data
store.attrs["channel_names"] = list(img.channel_names)
store.attrs["physical_pixel_sizes"] = [
    getattr(img.physical_pixel_sizes, ax) for ax in ("Z", "Y", "X")
]
store.attrs["source"] = str(snakemake.input[0])
log.info("wrote %s", snakemake.output[0])
