# EASI-FISH-Analysis

## Overview
Analysis of [EASI-FISH](https://www.cell.com/cell/fulltext/S0092-8674(21)01339-8) data of octopamine receptor-related mRNAs in the ring neurons of the *Drosophila* central complex

## Data naming convention

| Label | Line name | Neuron type |
|-------|-----------|-------------|
| A     | SS67611   | ER3p_a/b    |
| B     | SS04766   | ER3p_a/b    |
| C     | SS00238   | ER4d        |
| D     | SS17020   | ER1_a       |

| Label | receptor 1, 546 nm | receptor 2, 647 nm |
|-------|--------------------|--------------------|
| 1     | Oamb-B1            | Oct-TyR_B3         |
| 2     | Oct $\beta$ 1R-B3  | Oct $\beta$ 2R-B4  |
| 3     | Oct $\alpha$ 2R-B2 | Oct $\beta$ 3R-B5  |

The data consists of czi files with the naming convention {letter}{number}_brain{sample_number}{view}.
The letter corresponds to the genetic line used, as outlined in the first table above.
The number corresponds to the subset of receptors-related mRNAs tagged, as outlined in the second table above.
The view can be blank, a, or b, depending on how many different configurations were used to image the brain.

## Data structure

Each czi file has four color channels, recorded in a z stack throughout the brain.
The green channel is a membrane tag for neurons in the given genetic line.
The 546 and 647 nm channels are the mRNA probes, as described in the second table in the *Data naming convention* section
The fourth channel is DAPI staining.

## Analysis approach
1. Determine which voxels map to the neuron of interest (somata only so far; neurites are not yet included). The membrane channel is thresholded (globally, and relative to the brightest membrane nearby, to cut out-of-focus glow), filled, and stripped of neurites by a morphological opening. [Cellpose](https://github.com/MouseLand/cellpose) segments nuclei in DAPI, and a nucleus whose core lies inside the filled membrane belongs to a labelled neuron; each soma is the filled membrane nearest its nucleus. Soma-sized filled objects without a nucleus are kept but flagged (`has_nucleus = False` in `cells.csv`): in `C1_brain1a` they are bulb neuropil, not somata, so they are left out of the plots. Per-sample QC figures are written to `results/figures/segmentation/`.
2. Localize the mRNAs. [Spotiflow](https://github.com/weigertlab/spotiflow)'s pretrained 3D smFISH model proposes puncta, and candidates whose peak is less than `min_amplitude` counts above the local background are dropped (`config/config.yaml`). Spotiflow was chosen over a Laplacian-of-Gaussian detector (Big-FISH) because it ignores the defocus rings around bright puncta, which Big-FISH picks up. RS-FISH was not benchmarked. Per-sample QC plots are written to `results/figures/spots/`.
3. Colocalize the mRNAs with the neuron voxels to determine how many mRNAs are in each cell type. Each spot in a soma is tagged nuclear or cytoplasmic, or unknown in a soma without a detected nucleus.