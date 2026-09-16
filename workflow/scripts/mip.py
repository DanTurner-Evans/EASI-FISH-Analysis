"""Maximum intensity projection of the membrane-tag channel, for visual QC.

Contrast is a percentile stretch, so brightness is comparable in shape but not
in absolute intensity across samples. A scale bar is burned into the corner.

Reads the CZI directly rather than the ingested zarr: only one channel is
needed, and this avoids writing a full zarr copy of every sample just for QC.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from bioio import BioImage

from common import setup_logging

log = setup_logging(snakemake)  # noqa: F821  (injected by Snakemake)

img = BioImage(snakemake.input[0])
channel = snakemake.params.channel_index
opts = snakemake.params.mip
log.info("dims=%s  channel %d of %s", img.dims, channel, img.channel_names)

# Only the requested channel is read; the full four-channel stack never is.
mip = img.get_image_data("ZYX", T=0, C=channel).max(axis=0)

lo, hi = np.percentile(mip, opts["percentiles"])
log.info("raw range %d-%d, display range %.0f-%.0f", mip.min(), mip.max(), lo, hi)
disp = np.clip((mip.astype(float) - lo) / (hi - lo), 0, 1)

rgb = np.zeros((*disp.shape, 3))
rgb[..., 1] = disp

# Scale bar, bottom right.
um_per_px = img.physical_pixel_sizes.X
margin, thickness = 60, 15
length = int(round(opts["scale_bar_um"] / um_per_px))
y0 = rgb.shape[0] - margin
x1 = rgb.shape[1] - margin
rgb[y0:y0 + thickness, x1 - length:x1] = 1
log.info("scale bar %g µm = %d px", opts["scale_bar_um"], length)

plt.imsave(snakemake.output[0], rgb)
log.info("wrote %s", snakemake.output[0])
