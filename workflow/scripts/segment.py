"""Label ring-neuron cell bodies from the reporter channel.

PLACEHOLDER SEGMENTATION. The Otsu + watershed approach below is a structurally
correct starting point, not a validated one — inspect the labels against a known
sample and replace this with Cellpose or a tuned classical pipeline before
trusting any counts downstream.
"""
import numpy as np
import tifffile
import zarr
from scipy import ndimage as ndi
from skimage.feature import peak_local_max
from skimage.filters import threshold_otsu
from skimage.segmentation import watershed

from common import setup_logging

log = setup_logging(snakemake)  # noqa: F821

store = zarr.open(snakemake.input[0], mode="r")
reporter = np.asarray(store[snakemake.params.channel_index])
log.info("segmenting channel %d, shape %s", snakemake.params.channel_index, reporter.shape)

smoothed = ndi.gaussian_filter(reporter.astype(np.float32), sigma=(1, 2, 2))
mask = smoothed > threshold_otsu(smoothed)
mask = ndi.binary_opening(mask, iterations=2)

# Split touching somata by distance-transform watershed.
distance = ndi.distance_transform_edt(mask)
peaks = peak_local_max(distance, labels=mask, min_distance=10, exclude_border=False)
markers = np.zeros(distance.shape, dtype=np.int32)
markers[tuple(peaks.T)] = np.arange(1, len(peaks) + 1)
labels = watershed(-distance, markers, mask=mask)

log.info("labelled %d objects", labels.max())
tifffile.imwrite(snakemake.output[0], labels.astype(np.uint16))
