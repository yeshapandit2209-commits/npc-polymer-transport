## Overview

We use coarse-grained molecular dynamics (LAMMPS) to study how charge patterning, local chain stiffness and hydrophobic cohesion control passive cargo transport through a polymer-grafted cylindrical pore, a model of the nuclear pore complex (NPC).

The repository has two parts:

1. **Single-chain simulations** (`blocklength/`, `single_chain_k_theta/`): how charge patterning and bending stiffness set the size of one polymer chain (Fig. 1c, ESI Fig. S2).
2. **Pore simulations** (`no_hydrophobic/`, `fA_0.2/`, `fA_0.4/`): cargo transport and free energy in the grafted pore without hydrophobic beads and with hydrophobic fractions f_A = 0.2 and 0.4 (Figs. 2 to 6, ESI Figs. S3 to S9).

Raw trajectories (DCD files) are not included because of their size. They are available from the corresponding author on request.

## Requirements

- LAMMPS (29 Aug 2024 version), built with the MOLECULE and KSPACE packages
- VMD (for building data files with the .tcl scripts and for snapshots)
- Python 3.8 or newer with numpy, scipy, pandas, matplotlib
- MDAnalysis (for the density profiles)

```
pip install numpy scipy pandas matplotlib MDAnalysis
```

## Model in brief

| Quantity | Value |
|---|---|
| Units | reduced Lennard-Jones units (σ, k_BT) |
| Single chain | N = 64 beads, charge blocklength λ = 1, 2, 4, 8, 16, 32 (k_θ = 0) |
| Single chain, stiffness | uncharged, k_θ = 0 to 50 |
| Pore | radius 10σ, length 25σ |
| Brush | 88 chains (11 rings × 8 chains), 40 beads each |
| Rigid fraction, f_rigid | 0, 0.25, 0.5, 0.75, 1.0 (k_θ = 5 for rigid segments) |
| Hydrophobic fraction, f_A | 0, 0.2, 0.4 (ε_AA = 1.5 k_BT, cutoff 2.5σ) |
| Cargo size, σ_c | 0.5, 1, 2, ..., 8σ (purely repulsive) |
| Transport runs | 100 trajectories per brush and cargo size, up to 2 × 10⁶ steps |
| Umbrella sampling | σ_c = 4σ, windows z₀ = −50 to 50σ (1σ apart), k = 10, 3 replicates |
