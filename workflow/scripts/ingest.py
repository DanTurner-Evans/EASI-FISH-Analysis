"""Convert one CZI to chunked zarr.

Reading a CZI is slow and single-threaded; every downstream rule reads the zarr
instead so that cost is paid exactly once per sample.
"""
import xml.etree.ElementTree as ET

import zarr
from bioio import BioImage

from common import setup_logging

log = setup_logging(snakemake)  # noqa: F821  (injected by Snakemake)

img = BioImage(snakemake.input[0])
log.info("dims=%s  channels=%s", img.dims, img.channel_names)

# Channel names (Cam1-T1, ...) say nothing about the fluorophore; the emission
# wavelength in the Zeiss metadata does. The first <Channel> element carrying it
# wins, as the same channel Id recurs elsewhere in the XML.
md = img.metadata
root = md if isinstance(md, ET.Element) else ET.fromstring(md)
emission = {}
for ch in root.iter("Channel"):
    em = ch.find("EmissionWavelength")
    if em is not None and em.text and ch.get("Name") not in emission:
        emission[ch.get("Name")] = float(em.text)
emission_nm = [emission[str(name)] for name in img.channel_names]
log.info("emission wavelengths (nm): %s", emission_nm)

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
store.attrs["channel_names"] = [str(n) for n in img.channel_names]
store.attrs["emission_nm"] = emission_nm
store.attrs["physical_pixel_sizes"] = [
    getattr(img.physical_pixel_sizes, ax) for ax in ("Z", "Y", "X")
]
store.attrs["source"] = str(snakemake.input[0])
log.info("wrote %s", snakemake.output[0])
