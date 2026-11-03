# npc-polymer-transport
Overview

We use coarse-grained molecular dynamics (LAMMPS) to study how charge patterning, local chain stiffness and hydrophobic cohesion control passive cargo transport through a polymer-grafted cylindrical pore, a model of the nuclear pore complex (NPC).

The repository has two parts:

Single-chain simulations (blocklength/, single_chain_k_theta/): how charge patterning and bending stiffness set the size of one polymer chain (Fig. 1c, ESI Fig. S2).

Pore simulations (no_hydrophobic/, fA_0.2/, fA_0.4/): cargo transport and free energy in the grafted pore without hydrophobic beads and with hydrophobic fractions f_A = 0.2 and 0.4 (Figs. 2 to 6, ESI Figs. S3 to S9).

Raw trajectories (DCD files) are not included because of their size. They are available from the corresponding author on request.

Requirements
LAMMPS (29 Aug 2024 version), built with the MOLECULE and KSPACE packages
VMD (for building data files with the .tcl scripts and for snapshots)
Python 3.8 or newer with numpy, scipy, pandas, matplotlib
MDAnalysis (for the density profiles)
pip install numpy scipy pandas matplotlib MDAnalysis

Model in brief

Units	- reduced Lennard-Jones units (sigma, k_B T)
Single chain - N = 64 beads, charge blocklength lambda = 1, 2, 4, 8, 16, 32 (k_theta = 0)
Single chain, stiffness	- uncharged, k_theta = 0 to 50

Pore	- radius 10 sigma, length 25 sigma
Brush	- 88 chains (11 rings x 8 chains), 40 beads each
Rigid fraction, f_rigid	- 0, 0.25, 0.5, 0.75, 1.0 (k_theta = 5 for rigid segments)

Hydrophobic fraction, f_A	- 0, 0.2, 0.4 (epsilon_AA = 1.5 k_B T, cutoff 2.5 sigma)
Cargo size, sigma_c	- 0.5, 1, 2, ..., 8 sigma (purely repulsive)
Transport runs	- 100 trajectories per brush and cargo size, up to 2 x 10^6 steps

Umbrella sampling	- sigma_c = 4 sigma, windows z0 = -50 to 50 sigma (1 sigma apart), k = 10, 3 replicates

