"""Segment the somata of driver-line neurons from the membrane channel.

1. Membrane outline: voxels above a global threshold that are also bright
   relative to the brightest membrane within local_window_um. The relative test
   cuts the out-of-focus glow that a bright soma throws over its neighbours.
   Small gaps are closed, enclosed holes filled (in 3D and per plane), and
   anything thinner than neurite_um (neurites, fibre bundles) is opened away.
2. Nuclei: Cellpose on DAPI, in the box around the outline. A
   nucleus belongs to a labelled neuron if most of its core (deeper than
   nucleus_core_um inside it) lies within the outline; unlabelled neighbours
   only touch the outline with their edges. Nuclei in neuropil (surrounded by
   solid outline) are rejected.
3. Somata: each outline voxel within max_soma_radius_um of a labelled nucleus
   goes to the nearest one. (A watershed on membrane intensity lets a soma
   flood along the bright shell of the soma stacked below it.)
4. Outline components with no labelled nucleus at all are kept as
   nucleus-free cells (has_nucleus = False), for inspection only: they are not
   counted, because on review most were dendrites, axons, trachea or fat
   globules. They must be soma-sized and soma-shaped:
   a thick body (some voxel at least orphan_min_core_um from its surface) that
   is evenly filled with membrane signal (coefficient of variation at most
   orphan_max_cv). Neuropil that survives the opening (D1 dendrites, D3
   arbours) is patchy, and fibre bundles are thin.

The QC figure is drawn by segmentation_qc.py from the saved labels, so it can
be redrawn without re-running Cellpose.

Everything runs on a grid downsampled by `downsample` (near-isotropic ~0.5 um
voxels); the label images are written on that grid. A cell with a nucleus
takes the nucleus's ID, so nuclei_labels.tif identifies each soma's nucleus.
"""
import numpy as np
import pandas as pd
import tifffile
import torch
import zarr
from cellpose import models
from scipy import ndimage as ndi
from skimage.measure import block_reduce

from common import channel_index, membrane_outline, setup_logging

log = setup_logging(snakemake)  # noqa: F821

p = snakemake.params.segmentation
store = zarr.open(snakemake.input[0], mode="r")
block = tuple(p["downsample"])
vox = np.array(store.attrs["physical_pixel_sizes"]) * block


def load(emission_nm):
    return block_reduce(np.asarray(store[channel_index(store, emission_nm)]),
                        block, np.mean).astype(np.float32)


membrane = load(p["membrane_emission_nm"])
dapi = load(p["nuclei_emission_nm"])
log.info("grid %s, voxel %s um", membrane.shape, vox.round(3))

# --- 1. membrane outline ------------------------------------------------------------
smooth, fg, outline = membrane_outline(membrane, p, vox)
components, n_components = ndi.label(outline)
log.info("membrane foreground %.3f%% of voxels; outline %.3f%% in %d components",
         100 * fg.mean(), 100 * outline.mean(), n_components)

# --- 2. nuclei, in one box around the outline --------------------------------------
# One box, not one per component: Cellpose-SAM misses nuclei in small crops
# (a D3_brain2a nucleus found in the whole volume was lost in crops of 96-192 px).
min_soma, max_soma = p["soma_volume_um3"]
component_volume = ndi.sum_labels(np.ones_like(components), components,
                                  np.arange(1, n_components + 1)) * vox.prod()
big = np.isin(components, np.flatnonzero(component_volume >= min_soma) + 1)
nuclei = np.zeros(dapi.shape, np.int32)
if big.any():
    pad = np.ceil(8.0 / vox).astype(int)
    lo = np.maximum(np.argwhere(big).min(0) - pad, 0)
    hi = np.minimum(np.argwhere(big).max(0) + pad + 1, dapi.shape)
    box = tuple(slice(a, b) for a, b in zip(lo, hi))
    device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    model = models.CellposeModel(gpu=device.type != "cpu", device=device)
    masks, _, _ = model.eval(dapi[box], do_3D=True, z_axis=0, anisotropy=vox[0] / vox[1])
    nuclei[box] = masks
    log.info("Cellpose on box %s: %d nuclei", dapi[box].shape, masks.max())

ids = np.unique(nuclei[nuclei > 0])
volume = ndi.sum_labels(np.ones_like(nuclei), nuclei, ids) * vox.prod()
min_nucleus, max_nucleus = p["nucleus_volume_um3"]  # outside: fragments, glia, merged nuclei
sized = (volume >= min_nucleus) & (volume <= max_nucleus)
nuclei[np.isin(nuclei, ids[~sized])] = 0
ids, volume = ids[sized], volume[sized]
depth = ndi.distance_transform_edt(nuclei > 0, sampling=vox)
core = (depth >= p["nucleus_core_um"]) & (nuclei > 0)
inside = np.nan_to_num(ndi.mean(outline[core], nuclei[core], ids))
labelled = ids[inside >= p["min_nucleus_inside"]]
# Neuropil is solid outline, so any nucleus within it (glia) counts as inside.
# Share out the space within neuropil_check_um of labelled nuclei, each voxel
# to the nearest: a soma's share is mostly empty beyond its own outline (or
# claimed by its neighbours), while a nucleus in neuropil is surrounded by it.
dist, nearest = ndi.distance_transform_edt(~np.isin(nuclei, labelled), sampling=vox, return_indices=True)
share = np.where(dist <= p["neuropil_check_um"], nuclei[tuple(nearest)], 0)
crowding = ndi.mean(outline, share, labelled)
log.info("%d nuclei of %s-%s um3; %d with core inside outline >= %s, %d of them in neuropil "
         "(outline fills > %s of their share)", len(ids), min_nucleus, max_nucleus, len(labelled),
         p["min_nucleus_inside"], (crowding > p["max_crowding"]).sum(), p["max_crowding"])
labelled = labelled[crowding <= p["max_crowding"]]

# --- 3. somata grown from labelled nuclei -----------------------------------------
is_labelled = np.isin(nuclei, labelled)
dist, nearest = ndi.distance_transform_edt(~is_labelled, sampling=vox, return_indices=True)
soma = np.where((outline & (dist <= p["max_soma_radius_um"])) | is_labelled,
                nuclei[tuple(nearest)], 0)
# Fill each soma's enclosed holes, without taking voxels from a neighbour.
for k, box in enumerate(ndi.find_objects(soma), 1):
    if box is None:
        continue
    own = ndi.binary_fill_holes(soma[box] == k)
    own = np.stack([ndi.binary_fill_holes(plane) for plane in own])
    soma[box][own & (soma[box] == 0)] = k

# --- 4. nucleus-free cells ---------------------------------------------------------
# Only whole components of soma size: cutting larger ones up turns neuropil
# (ellipsoid body, bulb) into soma-sized pieces.
has_labelled = np.unique(components[is_labelled & outline])
candidates = [c for c in range(1, n_components + 1)
              if c not in has_labelled and min_soma <= component_volume[c - 1] <= max_soma]
orphans = []
component_boxes = ndi.find_objects(components)
for c in candidates:
    box = component_boxes[c - 1]
    own = components[box] == c
    core = ndi.distance_transform_edt(np.pad(own, 1), sampling=vox).max()
    values = smooth[box][own]
    cv = values.std() / values.mean()
    keep = core >= p["orphan_min_core_um"] and cv <= p["orphan_max_cv"]
    log.info("nucleus-free component %d: %.0f um3, core %.2f um, cv %.2f -> %s",
             c, component_volume[c - 1], core, cv, "kept" if keep else "rejected")
    if keep:
        orphans.append(c)
orphan_ids = list(range(int(nuclei.max()) + 1, int(nuclei.max()) + 1 + len(orphans)))
for c, k in zip(orphans, orphan_ids):
    soma[components == c] = k
log.info("%d nucleus-free cells, %d soma-sized components rejected as neuropil or fibre",
         len(orphan_ids), len(candidates) - len(orphans))

# --- table and label images ------------------------------------------------------
cell_ids = np.concatenate([labelled, orphan_ids]).astype(int)
centres = np.array(ndi.center_of_mass(soma > 0, soma, cell_ids)).reshape(-1, 3) * vox
edge = np.zeros(soma.shape, bool)
for axis in range(3):
    edge[(slice(None),) * axis + (0,)] = edge[(slice(None),) * axis + (-1,)] = True
nucleus_volume = dict(zip(ids, volume))
nucleus_inside = dict(zip(ids, inside))
table = pd.DataFrame({
    "cell_id": cell_ids,
    "has_nucleus": np.isin(cell_ids, labelled),
    "z_um": centres[:, 0], "y_um": centres[:, 1], "x_um": centres[:, 2],
    "soma_volume_um3": ndi.sum_labels(np.ones_like(soma), soma, cell_ids) * vox.prod(),
    "nucleus_volume_um3": [nucleus_volume.get(k, np.nan) for k in cell_ids],
    "nucleus_inside": [nucleus_inside.get(k, np.nan) for k in cell_ids],
    # A soma cut by the edge of the volume is only partly imaged, so undercounted.
    "touches_edge": np.isin(cell_ids, np.unique(soma[edge])),
})
table.to_csv(snakemake.output.table, index=False)
for path, labels in [(snakemake.output.soma, soma), (snakemake.output.nuclei, nuclei)]:
    tifffile.imwrite(path, labels.astype(np.uint16), imagej=True, resolution=(1 / vox[2], 1 / vox[1]),
                     metadata={"axes": "ZYX", "spacing": vox[0], "unit": "um"})
log.info("soma volumes (um3): %s", table.soma_volume_um3.round().tolist())

