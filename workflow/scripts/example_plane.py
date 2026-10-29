"""Example figure: raw channels and detected spots in one region of one sample.

Top row: membrane, each probe and DAPI, raw (max projection over a few planes,
about one spot deep). Bottom row: composite, each probe with its detected
(yellow) and rejected (red, below the sample's amplitude threshold) candidates, and the composite
with detected spots. A spot is drawn if its centre lies in the projected planes.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import zarr

from common import channel_index, min_amplitude, setup_logging

log = setup_logging(snakemake)  # noqa: F821

p = snakemake.params
ex = p.example
store = zarr.open(snakemake.input.image, mode="r")
vox = np.array(store.attrs["physical_pixel_sizes"])
z0, z1 = ex["planes"]
y0, x0 = ex["origin_px"]
size = ex["size_px"]
box = (slice(z0, z1 + 1), slice(y0, y0 + size), slice(x0, x0 + size))


def project(emission_nm):
    return np.asarray(store[channel_index(store, emission_nm)][box]).max(axis=0).astype(np.float32)


def stretch(img):
    lo, hi = np.percentile(img, [1, 99.8])
    return np.clip((img - lo) / (hi - lo + 1e-9), 0, 1)


membrane = stretch(project(p.membrane_emission_nm))
dapi = stretch(project(p.nuclei_emission_nm))
probes = {channel: stretch(project(p.emission_nm[channel])) for channel in p.channels}
spots, threshold = {}, {}
for channel, path in zip(p.channels, snakemake.input.candidates):
    c = pd.read_csv(path)
    threshold[channel] = min_amplitude(c.amplitude.to_numpy(), p.spot_filter, p.overrides.get(channel))[0]
    in_box = (c.y >= y0) & (c.y < y0 + size) & (c.x >= x0) & (c.x < x0 + size)
    c = c[in_box & c.z.round().between(z0, z1)]  # spot centred in a projected plane
    spots[channel] = (c, c.amplitude >= threshold[channel])

# Composite: green membrane, blue DAPI, then magenta and cyan for the probes.
probe_colours = [(1, 0, 1), (0, 1, 1)]
composite = np.zeros((size, size, 3))
composite += membrane[..., None] * (0, 1, 0)
composite += dapi[..., None] * (0.2, 0.2, 0.6)
for img, colour in zip(probes.values(), probe_colours):
    composite += img[..., None] * colour
composite = np.clip(composite, 0, 1)

n = len(probes) + 2
fig, axes = plt.subplots(2, n, figsize=(4.5 * n, 10.2), gridspec_kw={"hspace": 0.12})
top = [("Membrane (driver line), raw", membrane)] + \
      [(f"{channel} probe, raw", img) for channel, img in probes.items()] + [("DAPI, raw", dapi)]
for ax, (title, img) in zip(axes[0], top):
    ax.imshow(img, cmap="gray")
    ax.set_title(title)
bar = 10 / vox[2]
axes[0, 0].plot([0.04 * size, 0.04 * size + bar], [0.94 * size] * 2, color="white", lw=6)
axes[0, 0].text(0.04 * size + bar / 2, 0.91 * size, "10 µm", color="white", ha="center")

axes[1, 0].imshow(composite)
axes[1, 0].set_title("Composite")
for ax, (channel, img) in zip(axes[1, 1:-1], probes.items()):
    c, kept = spots[channel]
    ax.imshow(img, cmap="gray")
    for mask, colour, label in [(kept, "yellow", f"detected ({kept.sum()})"),
                                (~kept, "red", f"rejected: amplitude < {threshold[channel]:.0f} "
                                               f"({(~kept).sum()})")]:
        ax.scatter(c.x[mask] - x0, c.y[mask] - y0, s=90, facecolors="none", edgecolors=colour,
                   linewidths=1.2, label=label)
    ax.legend(loc="lower left", fontsize=9)
    ax.set_title(f"{channel} probe + spots")
axes[1, -1].imshow(composite)
for (channel, (c, kept)), colour in zip(spots.items(), probe_colours):
    axes[1, -1].scatter(c.x[kept] - x0, c.y[kept] - y0, s=90, facecolors="none", edgecolors=[colour],
                        linewidths=1.2, label=f"{channel} spots ({kept.sum()})")
axes[1, -1].legend(loc="lower left", fontsize=9)
axes[1, -1].set_title("Composite + detected spots")
for ax in axes.flat:
    ax.set(xticks=[], yticks=[], xlim=(-0.5, size - 0.5), ylim=(size - 0.5, -0.5))

fig.suptitle(f"{snakemake.wildcards.sample}: max projection of planes {z0}-{z1} "
             f"({z0 * vox[0]:.1f}-{z1 * vox[0]:.1f} µm deep, one spot depth); "
             "spots shown if centred in that range", fontsize=13)
fig.text(0.5, 0.01, "Composite: green membrane, " +
         ", ".join(f"{colour} {channel}" for channel, colour in zip(probes, ["magenta", "cyan"])) +
         ", blue DAPI. Each panel contrast-stretched independently (1-99.8 percentile), "
         "so brightness is not comparable between channels.", ha="center", fontsize=10)
fig.subplots_adjust(left=0.01, right=0.99, top=0.93, bottom=0.04, wspace=0.04)
fig.savefig(snakemake.output[0], dpi=80)
log.info("wrote %s", snakemake.output[0])
