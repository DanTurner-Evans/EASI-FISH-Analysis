# EASI-FISH-Analysis

## Overview
Analysis of [EASI-FISH](https://www.cell.com/cell/fulltext/S0092-8674(21)01339-8) data of octopamine receptor-related mRNAs in the ring neurons of the *Drosophila* central complex

## Data naming convention

| Label | Line name | Neuron type |
|-------|-----------|-------------|
| A     | SS67611   | ER3p_a/b    |
| B     | SS17020   | ER1_a       |
| C     | SS00238   | ER4d        |
| D     | SS04766   | ER3p_a/b    |

| Label | receptor 1, 546 nm | receptor 2, 647 nm |
|-------|--------------------|--------------------|
| 1     | Oamb-B1            | Oct-TyR_B3         |
| 2     | Oct $\beta$ 1R-B3  | Oct $\beta$ 2R-B4  |
| 3     | Oct $\alpha$ 2R-B2 | Oct $\beta$ 3R-B5  |

The data consists of czi files with the naming convention {letter}{number}_brain{sample_number}{view}.
The letter corresponds to the genetic line used, as outlined in the first table above.
The number corresponds to the subset of receptors-related mRNAs tagged, as outlies in the second table above.
The view can be blank, a, or b, depending on how many different configurations were used to image the brain.

## Data structure

Each czi file has four color channels, recorded in a z stack throughout the brain.
The green channel is a membrane tag for neurons in the given genetic line.
The 546 and 647 nm channels are the mRNA probes, as desrobed in the second table in the *Data naming convention* section
The

## Analysis approach
1. Determine which voxels map to the neuron of interest
2. Localize the mRNAs
3. Coloclaize the mRNAs with the neuron voxels to determine how many mRNAs are in each cell type.