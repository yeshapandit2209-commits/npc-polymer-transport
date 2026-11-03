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
Quantity	Value
Units	reduced Lennard-Jones units (sigma, k_B T)
Single chain	N = 64 beads, charge blocklength lambda = 1, 2, 4, 8, 16, 32 (k_theta = 0)
Single chain, stiffness	uncharged, k_theta = 0 to 50
Pore	radius 10 sigma, length 25 sigma
Brush	88 chains (11 rings x 8 chains), 40 beads each
Rigid fraction, f_rigid	0, 0.25, 0.5, 0.75, 1.0 (k_theta = 5 for rigid segments)
Hydrophobic fraction, f_A	0, 0.2, 0.4 (epsilon_AA = 1.5 k_B T, cutoff 2.5 sigma)
Cargo size, sigma_c	0.5, 1, 2, ..., 8 sigma (purely repulsive)
Transport runs	100 trajectories per brush and cargo size, up to 2 x 10^6 steps
Umbrella sampling	sigma_c = 4 sigma, windows z0 = -50 to 50 sigma (1 sigma apart), k = 10, 3 replicates

Full details are in the Methods section of the paper.

How to reproduce
Single chain: charge blocklength (blocklength/)
01_system_generation/: generate_saw.py writes a self-avoiding chain; build_data.tcl (VMD) assigns charges in blocks of length lambda and writes the LAMMPS data file.
02_simulation/: run lmp -in in.blocklength.lammps for each lambda. It writes the radius of gyration over time.
03_analysis/: python3 plot_rg.py <folder with rg_bl*.dat> gives R_g/sqrt(N), averaged over the second half (Fig. 1c).
Single chain: bending stiffness (single_chain_k_theta/)
01_system_generation/: generate_saw.py and build_data_simple.tcl build an uncharged chain.
02_simulation/: run lmp -in in.angle_k.lammps for each k_theta.
Analyse with blocklength/03_analysis/plot_rg.py <folder with rg_k*.dat> (ESI Fig. S2).
Pore brushes (no_hydrophobic/, fA_0.2/, fA_0.4/)
01_system_generation/: build and equilibrate the grafted pore for each f_rigid.
02_simulation/: set up the transport runs (100 trajectories per cargo size) and the umbrella windows.
03_analysis/: assign the outcome of each trajectory, count passed and trapped cargo, fit the size cutoff sigma_50 (eqn 6 of the paper, with 95% bootstrap intervals) and compute density profiles.
04_pmf/: solve WHAM for the umbrella windows to get F(z) and F(r, z).

Which folder gives which figure
Figure	Content	Folder
Fig. 1c	R_g vs time for each lambda	blocklength/
ESI Fig. S2	R_g vs k_theta	single_chain_k_theta/
Fig. 2	brush structure, density profiles	no_hydrophobic/03_analysis/
Fig. 3a,b	P_pass and P_trap vs cargo size	no_hydrophobic/03_analysis/
Fig. 3c,d	PMF, entry barrier, well depth	no_hydrophobic/04_pmf/
Fig. 4	effect of hydrophobic beads, sigma_50	fA_0.2/, fA_0.4/ (with no_hydrophobic/)
Fig. 5	PMF with hydrophobic beads, density vs P_pass	fA_0.2/, fA_0.4/
Fig. 6	transport regimes	all three brush folders
ESI Figs. S6, S8, S9	PMFs, density profiles and transport with hydrophobic beads	fA_0.2/, fA_0.4/
File naming
0, 25, 50, 75, 100: rigid fraction f_rigid in percent.
k5: bending constant k_theta = 5 for the rigid segments.
bl<lambda>: charge blocklength; k<value>: bending constant of the single chain.

Acknowledgements

This work was supported by the Core Research Grant CRG/2022/008633 (SERB/ANRF, Government of India). Simulations were run on the HPC facility of IIT Gandhinagar and the PARAM Ananta supercomputer (National Supercomputing Mission).
