"""
Compute rho_z(z) and rho_r(r) [inside / outside pore] density profiles
for a series of FG-Nup rigidity fractions (f = 0, 25, 50, 75, 100),
and write the results to CSV files for later plotting.

Requires: MDAnalysis, numpy, pandas
    pip install MDAnalysis numpy pandas

Expected files (edit FRACTIONS / filename patterns below if different):
    cyl_0.data,   dump_0.dcd
    cyl_25.data,  dump_25.dcd
    cyl_50.data,  dump_50.dcd
    cyl_75.data,  dump_75.dcd
    cyl_100.data, dump_100.dcd

Each .data file must contain a "Pair Coeffs" section like:
    # 1  1     <- local type 1 -> canonical 1 (anchor)
    # 2  3     <- local type 2 -> canonical 3 (flexible)
    # 3  4     <- local type 3 -> canonical 4 (wall)
This mapping is parsed automatically, so the script does not care
which local type number means what in a given file.
"""

import os
import re
import numpy as np
import pandas as pd
import MDAnalysis as mda

# ----------------------------- CONFIG -----------------------------------

FRACTIONS = [0, 25, 50, 75, 100]

DATA_PATTERN = "cyl_{f}.data"     # topology file per fraction
DUMP_PATTERN = "dump_{f}k5.dcd"     # trajectory file per fraction

# atom_style used when the .data files were written (VMD/TopoTools "full")
ATOM_STYLE = "id resid type charge x y z"

# canonical type ids (from the Pair Coeffs comment column)
CANON_ANCHOR = 1
CANON_RIGID = 2
CANON_FLEXIBLE = 3
CANON_WALL = 4

# include the anchor bead (chain root) as part of the FG-Nup density?
INCLUDE_ANCHOR = True

# geometry: nuclear-envelope wall sits at z = +/-12.5 in your system
# (checked directly from the data files) -> pore channel = |z| < Z_PORE_HALF
Z_PORE_HALF = 12.5

# outer radius used to normalize rho_z (choose something that safely
# encloses all FG-Nup beads at all z, e.g. half the box width)
R_MAX = 75.0

# binning
NBINS_Z = 100
NBINS_R = 60
R_MAX_HIST = 40.0   # max r (in same units as coordinates) for rho_r histograms

# trajectory handling
EQUIL_FRACTION = 0.5   # skip this fraction of frames as equilibration
STRIDE = 1              # use every STRIDE-th frame after equilibration

# --------------------------------------------------------------------------


def parse_pair_coeff_map(data_file):
    """Return {local_type(int): canonical_type(int)} from the Pair Coeffs
    comment block of a LAMMPS data file."""
    mapping = {}
    with open(data_file) as fh:
        lines = fh.readlines()

    in_block = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("# Pair Coeffs"):
            in_block = True
            continue
        if in_block:
            if stripped.startswith("# Bond Coeffs") or (
                stripped == "" and mapping
            ):
                break
            m = re.match(r"#\s*(\d+)\s+(\d+)\s*$", stripped)
            if m:
                local_t, canon_t = int(m.group(1)), int(m.group(2))
                mapping[local_t] = canon_t
    if not mapping:
        raise RuntimeError(f"Could not parse Pair Coeffs mapping from {data_file}")
    return mapping


def get_role_selection(universe, type_map):
    """Return an AtomGroup of the FG-Nup beads to use for density calcs."""
    types = universe.atoms.types.astype(int)
    canon = np.vectorize(type_map.get)(types)

    wanted = [CANON_RIGID, CANON_FLEXIBLE]
    if INCLUDE_ANCHOR:
        wanted.append(CANON_ANCHOR)

    mask = np.isin(canon, wanted)
    return universe.atoms[mask]


def compute_profiles_for_fraction(f):
    data_file = DATA_PATTERN.format(f=f)
    dump_file = DUMP_PATTERN.format(f=f)

    if not os.path.exists(data_file):
        raise FileNotFoundError(
            f"Missing topology file '{data_file}' for f={f}. "
            f"You need a data file with the same Pair-Coeffs convention "
            f"as cyl_50.data for every f value."
        )
    if not os.path.exists(dump_file):
        raise FileNotFoundError(f"Missing trajectory file '{dump_file}' for f={f}")

    type_map = parse_pair_coeff_map(data_file)

    u = mda.Universe(
        data_file, dump_file,
        topology_format="DATA", format="DCD",
        atom_style=ATOM_STYLE,
    )

    fg_atoms = get_role_selection(u, type_map)
    print(f"[f={f}] {len(fg_atoms)} FG-Nup beads selected "
          f"(anchor included: {INCLUDE_ANCHOR})")

    z_bins = np.linspace(u.dimensions[2] * -0.5, u.dimensions[2] * 0.5, NBINS_Z + 1)
    r_bins = np.linspace(0, R_MAX_HIST, NBINS_R + 1)

    z_hist_sum = np.zeros(NBINS_Z)
    r_in_hist_sum = np.zeros(NBINS_R)
    r_out_hist_sum = np.zeros(NBINS_R)

    n_frames_total = len(u.trajectory)
    start = int(n_frames_total * EQUIL_FRACTION)
    frame_indices = range(start, n_frames_total, STRIDE)
    n_used = 0

    for i in frame_indices:
        u.trajectory[i]
        pos = fg_atoms.positions
        x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]
        r = np.sqrt(x ** 2 + y ** 2)

        # rho_z: all FG-Nup beads
        hz, _ = np.histogram(z, bins=z_bins)
        z_hist_sum += hz

        # rho_r inside pore (|z| < Z_PORE_HALF)
        inside = np.abs(z) < Z_PORE_HALF
        hr_in, _ = np.histogram(r[inside], bins=r_bins)
        r_in_hist_sum += hr_in

        # rho_r outside pore (|z| >= Z_PORE_HALF)
        outside = ~inside
        hr_out, _ = np.histogram(r[outside], bins=r_bins)
        r_out_hist_sum += hr_out

        n_used += 1

    if n_used == 0:
        raise RuntimeError(f"No frames used for f={f}; check STRIDE/EQUIL_FRACTION")

    dz = z_bins[1] - z_bins[0]
    dr = r_bins[1] - r_bins[0]
    z_centers = 0.5 * (z_bins[:-1] + z_bins[1:])
    r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])

    # number density: <count> / (frames * bin volume)
    vol_z = dz * np.pi * R_MAX ** 2
    rho_z = z_hist_sum / (n_used * vol_z)

    z_range_in = 2 * Z_PORE_HALF
    vol_r_in = 2 * np.pi * r_centers * dr * z_range_in
    rho_r_in = r_in_hist_sum / (n_used * vol_r_in)

    # "outside" spans both sides of the envelope combined; use full box
    # minus the pore slab as the effective z-extent
    z_range_out = u.dimensions[2] - z_range_in
    vol_r_out = 2 * np.pi * r_centers * dr * z_range_out
    rho_r_out = r_out_hist_sum / (n_used * vol_r_out)

    return (
        pd.DataFrame({"f": f, "z": z_centers, "rho_z": rho_z}),
        pd.DataFrame({"f": f, "r": r_centers, "rho_r": rho_r_in}),
        pd.DataFrame({"f": f, "r": r_centers, "rho_r": rho_r_out}),
    )


def main():
    z_all, r_in_all, r_out_all = [], [], []

    for f in FRACTIONS:
        try:
            dz, dr_in, dr_out = compute_profiles_for_fraction(f)
        except FileNotFoundError as e:
            print(f"SKIPPING f={f}: {e}")
            continue
        z_all.append(dz)
        r_in_all.append(dr_in)
        r_out_all.append(dr_out)

    if not z_all:
        print("No fractions processed - check your file names/paths.")
        return

    pd.concat(z_all, ignore_index=True).to_csv("density_z.csv", index=False)
    pd.concat(r_in_all, ignore_index=True).to_csv("density_r_inside.csv", index=False)
    pd.concat(r_out_all, ignore_index=True).to_csv("density_r_outside.csv", index=False)

    print("\nWrote: density_z.csv, density_r_inside.csv, density_r_outside.csv")


if __name__ == "__main__":
    main()
