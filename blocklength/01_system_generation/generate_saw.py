#!/usr/bin/env python3
"""
Self-avoiding random walk (SAW) generator for a single coarse-grained
polymer chain, one bead per residue.

Generates ONE backbone geometry (PDB). This same backbone is reused
for every blocklength -- only the charge pattern (assigned later in
the VMD/TCL step) differs between blocklengths, so all conditions
start from the same initial chain conformation.
"""

import numpy as np
import argparse


def generate_saw(n_beads, bond_length=1.0, min_dist=1.0, box=100.0,
                  max_attempts_per_bead=500, max_restarts=200, seed=None):
    rng = np.random.default_rng(seed)

    for restart in range(max_restarts):
        coords = [np.array([box / 2, box / 2, box / 2])]
        success = True
        for i in range(1, n_beads):
            placed = False
            for attempt in range(max_attempts_per_bead):
                v = rng.normal(size=3)
                v /= np.linalg.norm(v)
                candidate = coords[-1] + bond_length * v

                dists = np.linalg.norm(np.array(coords) - candidate, axis=1)
                if np.all(dists > min_dist):
                    coords.append(candidate)
                    placed = True
                    break
            if not placed:
                success = False
                break
        if success:
            print(f"SAW generated on restart {restart}, n={n_beads}")
            return np.array(coords)

    raise RuntimeError(
        "Failed to generate SAW after max_restarts; "
        "loosen min_dist or increase max_attempts_per_bead"
    )


def write_pdb(coords, filename, resname="BEA"):
    with open(filename, "w") as f:
        for i, xyz in enumerate(coords, start=1):
            f.write(
                f"ATOM  {i:5d}  CA  {resname} A{i:4d}    "
                f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00\n"
            )
        f.write("END\n")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=64, help="chain length")
    ap.add_argument("--bond", type=float, default=1.0, help="bond length (sigma units)")
    ap.add_argument("--mindist", type=float, default=1.0, help="min bead-bead distance")
    ap.add_argument("--box", type=float, default=100.0, help="cubic box side, chain centered")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", type=str, default="chain.pdb")
    args = ap.parse_args()

    coords = generate_saw(args.n, args.bond, args.mindist, args.box, seed=args.seed)
    write_pdb(coords, args.out)
    print(f"Wrote {args.out}")
