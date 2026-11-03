# npc-polymer-transport

Coarse-grained molecular dynamics of cargo transport through a model nuclear pore complex.

Contact: Mithun Radhakrishna (mithunr@iitgn.ac.in), Yesha Pandit (yesha.pandit@iitgn.ac.in)
Department of Chemical Engineering, Indian Institute of Technology Gandhinagar, India.

---

## Overview

We use coarse-grained molecular dynamics (LAMMPS) to study how charge patterning, local chain stiffness and hydrophobic cohesion control passive cargo transport through a polymer-grafted cylindrical pore, a model of the nuclear pore complex (NPC).

The repository has two parts:

1. **Single-chain simulations** (`blocklength/`, `single_chain_k_theta/`): how charge patterning and bending stiffness set the size of one polymer chain.
2. **Pore simulations** (`no_hydrophobic/`, `fA_0.2/`, `fA_0.4/`): cargo transport and free energy in the grafted pore without hydrophobic beads and with hydrophobic fractions f_A = 0.2 and 0.4.

Raw trajectories (DCD files) are not included because of their size. They are available from the authors on request.

---

## Repository structure

```
npc-polymer-transport/
├── blocklength/                  single charged chain, charge blocklength λ = 1 to 32
│   ├── 01_system_generation/     generate_saw.py, build_data.tcl
│   ├── 02_simulation/            in.blocklength.lammps
│   └── 03_analysis/              plot_rg.py
├── single_chain_k_theta/         single uncharged chain, bending constant k_θ = 0 to 50
│   ├── 01_system_generation/     generate_saw.py, build_data_simple.tcl
│   └── 02_simulation/            in.angle_k.lammps
├── no_hydrophobic/               pore brushes without hydrophobic beads
│   ├── 01_system_generation/     npc.py
│   ├── 02_simulation/            setup_runs.py, run_all.sh, lammps.config
│   ├── 03_analysis/              dwell_hydro.py, counts_from_runs.py, fit_sigma50.py, compute_density_profiles.py
│   └── 04_pmf/                   setup_umbrella.py, setup_hydro.py, wham_half.py
├── fA_0.2/                       pore brushes with hydrophobic fraction f_A = 0.2
│   ├── 01_system_generation/     build_cyl_system_tunable.py, chain.tcl
│   ├── 02_simulation/            setup_equil.py, setup_hydro.py
│   ├── 03_analysis/              dwell_hydro.py, counts_from_runs.py, fit_sigma50.py, compute_density_profiles.py
│   └── 04_pmf/                   setup_umbrella.py, setup_hydro.py, wham_half.py, diffusion_z.py
├── fA_0.4/                       same as fA_0.2/, with f_A = 0.4
└── README.md
```

---

## Requirements

- LAMMPS (29 Aug 2024 version), built with the MOLECULE and KSPACE packages
- VMD (for building data files with the .tcl scripts and for snapshots)
- Python 3.8 or newer with numpy, scipy, pandas, matplotlib
- MDAnalysis (for the density profiles)

```
pip install numpy scipy pandas matplotlib MDAnalysis
```

---

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

---

## 1. Single chain: charge blocklength (`blocklength/`)

| Folder | File | Purpose |
|---|---|---|
| `01_system_generation/` | `generate_saw.py` | generates a self-avoiding random walk of 64 beads |
| | `build_data.tcl` | VMD script; assigns charges in blocks of length λ and writes the LAMMPS data file |
| `02_simulation/` | `in.blocklength.lammps` | LAMMPS input: WCA + Coulomb (PPPM), NVT; writes R_g over time |
| `03_analysis/` | `plot_rg.py` | R_g/√N over time and its average over the second half of each run |

```
cd blocklength/03_analysis
python3 plot_rg.py <folder with rg_bl*.dat>
```

## 2. Single chain: bending stiffness (`single_chain_k_theta/`)

| Folder | File | Purpose |
|---|---|---|
| `01_system_generation/` | `generate_saw.py` | generates a self-avoiding random walk of 64 beads |
| | `build_data_simple.tcl` | VMD script; writes the data file of an uncharged chain |
| `02_simulation/` | `in.angle_k.lammps` | LAMMPS input with a harmonic angle potential of strength k_θ |

Analyse with the same script:
```
python3 blocklength/03_analysis/plot_rg.py <folder with rg_k*.dat>
```

---

## 3. Pore simulations (`no_hydrophobic/`, `fA_0.2/`, `fA_0.4/`)

`fA_0.2/` and `fA_0.4/` contain the same files; they differ only in the hydrophobic fraction f_A.

### 01_system_generation: build the grafted pore

| File | Purpose |
|---|---|
| `npc.py` | builds the brush without hydrophobic beads for each f_rigid |
| `build_cyl_system_tunable.py` | builds the brush with a fraction f_A of hydrophobic (A) beads placed at random along each chain |
| `chain.tcl` | VMD script that writes the LAMMPS data file |

### 02_simulation: equilibration and transport runs

| File | Purpose |
|---|---|
| `setup_equil.py` | equilibrates the brush with hydrophobic beads |
| `setup_runs.py` | sets up the transport runs without hydrophobic beads |
| `setup_hydro.py` | sets up the transport runs with hydrophobic beads (ε_AA = 1.5 k_BT, cutoff 2.5σ) |
| `run_all.sh`, `lammps.config` | job script and LAMMPS settings for running all trajectories |

Each set-up script writes 100 independent trajectories for every cargo size (σ_c = 0.5 to 8σ). A trajectory ends when the cargo leaves through the opposite end of the pore, moves more than 60σ away, or reaches 2 × 10⁶ steps.

### 03_analysis: transport

Run in this order:

| Step | File | Purpose |
|---|---|---|
| 1 | `dwell_hydro.py` | reads each trajectory and assigns its outcome: passed, trapped, escaped or outside |
| 2 | `counts_from_runs.py` | counts passed and trapped trajectories per brush and cargo size; writes `counts.csv` |
| 3 | `fit_sigma50.py` | fits P_pass(σ_c) = P₀ / (1 + exp[(σ_c − σ₅₀)/w]) by binomial maximum likelihood; size cutoff σ₅₀ with 95% bootstrap interval and goodness-of-fit p |
| 4 | `compute_density_profiles.py` | polymer density along the pore axis, inside and outside the pore |

```
cd no_hydrophobic/03_analysis
python3 dwell_hydro.py --flat 0 25 50 75 100 --name nh
python3 counts_from_runs.py none dwell_runs.csv
python3 fit_sigma50.py 16          # number of CPU cores
```
For `fA_0.2/` and `fA_0.4/`, use the f_A value as the label, for example `python3 counts_from_runs.py 0.2 dwell_runs.csv`.

### 04_pmf: free energy of a 4σ cargo

| File | Purpose |
|---|---|
| `setup_umbrella.py` | writes the umbrella windows (z₀ = −50 to 50σ, k = 10, 3 replicates). It imports `setup_hydro.py` for the force field, so keep both in this folder |
| `wham_half.py` | solves WHAM with Newton's method; PMF F(z), s.d. over replicates, half-to-half convergence check and 2D map F(r, z) |
| `diffusion_z.py` | position-dependent diffusion coefficient D(z) of the cargo |

```
cd no_hydrophobic/04_pmf
python3 setup_umbrella.py 0k5_05.data 50k5_05.data 100k5_05.data --sigmas 4 --reps 3 --zmin -50 --zmax 50 --no-start
# run the windows, then
python3 wham_half.py 0k5/umbrella_s4.0 50k5/umbrella_s4.0 100k5/umbrella_s4.0
```

---

## File naming

- `0`, `25`, `50`, `75`, `100`: rigid fraction f_rigid in percent.
- `k5`: bending constant k_θ = 5 for the rigid segments.
- `bl<λ>`: charge blocklength of the single chain; `k<value>`: bending constant of the single chain.
- `s20`, `4`: brushes with f_A = 0.2 and 0.4.

---


## Acknowledgements

 Simulations were run on the High-Performance Computing Facility of IIT Gandhinagar and the PARAM Ananta supercomputer (National Supercomputing Mission).
