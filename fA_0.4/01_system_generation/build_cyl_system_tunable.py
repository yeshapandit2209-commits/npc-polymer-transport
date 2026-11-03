"""
Cylindrical FG-Nup brush -- initial configuration + LAMMPS data file builder.

Same geometry/growth engine as npc.py (unchanged). The only thing replaced
is how "A" beads (the hydrophobic monomer, atom name "A") get placed along
each chain. You now set EXACT counts, independently per domain, via:

    A_PLAN = {
        "rigid": {"singlet": 2, "doublet": 1, "triplet": 1},
        "flex":  {"singlet": 1, "doublet": 1, "triplet": 2},
    }

These are per-chain counts -- every chain gets exactly this many
singlet/doublet/triplet A-clusters in its rigid domain, and exactly this
many in its flex domain. Only the POSITION of each cluster along the chain
(and the order they're laid down in) is randomized -- different every run
unless you seed rng below. Nothing about the counts themselves is random,
so you have full independent control over singlet/doublet/triplet in the
rigid domain versus the flex domain.

Everything else (ring anchor placement, self-avoiding chain growth, wall
pseudoatom grid, LAMMPS/xyz writers) is untouched from npc.py.
"""

from collections import defaultdict

import numpy as np

# ================= Geometry (measured, unchanged from npc.py) =================

R_cyl = 10.0
L_cyl = 25.0
BOX = 150.0

n_rings = 11
chains_per_ring = 8
n_chains = n_rings * chains_per_ring          # 88

n_monomers = 40
n_atoms_chain = n_chains * n_monomers          # 3520

# ================= Bonded parameters (measured, unchanged) =================

sigma = 1.0
sigma_wall = 1.0
r0 = 1.0
bond_tolerance = 0.05
bond_len_min = (1 - bond_tolerance) * r0
bond_len_max = (1 + bond_tolerance) * r0

intra_sigma = 0.95 * sigma
inter_sigma = 0.90 * sigma

# Rigidity: fraction of each chain (after the anchor) that is C_rigid,
# the rest is C_flex. rigid domain = positions [1, rigid_length],
# flex domain = positions [rigid_length+1, n_monomers-1].
rigidity_percent = 100
rigid_length = min(int((rigidity_percent / 100) * n_monomers), n_monomers - 1)
flex_length = (n_monomers - 1) - rigid_length

# ================= TUNABLE: A-bead (hydrophobic monomer) placement =================
# STEP 1 -- pick a scope:
#   "whole"  -> ONE percentage applies to the whole chain (both domains
#               together). Set A_PERCENT and A_RUN_SPLIT below.
#   "domain" -> tune rigid and flex separately, each with its own
#               percentage and its own singlet/doublet/triplet split.
#               Set A_PERCENT_RIGID/A_RUN_SPLIT_RIGID and
#               A_PERCENT_FLEX/A_RUN_SPLIT_FLEX below.

A_SCOPE = "whole"

# ---- used when A_SCOPE == "whole" ----
A_PERCENT = 80             # % of the whole chain (positions 1..n_monomers-1) that is A
A_RUN_SPLIT = {            # % of A-clusters that are singlets/doublets/triplets (must sum to 100)
    "singlet": 100,
    "doublet": 0,
    "triplet": 0,
}

# ---- used when A_SCOPE == "domain" ----
A_PERCENT_RIGID = 20
A_RUN_SPLIT_RIGID = {"singlet": 100, "doublet": 0, "triplet": 0}

A_PERCENT_FLEX = 20
A_RUN_SPLIT_FLEX = {"singlet": 40, "doublet": 20, "triplet": 40}

# ---- nothing below this line needs editing ----

assert A_SCOPE in ("whole", "domain")
for _split in ([A_RUN_SPLIT] if A_SCOPE == "whole" else [A_RUN_SPLIT_RIGID, A_RUN_SPLIT_FLEX]):
    assert abs(sum(_split.values()) - 100) < 1e-6, f"run split must sum to 100: {_split}"


def _percent_to_counts(percent, run_split, domain_len):
    """Converts a percentage + singlet/doublet/triplet split into exact
    integer counts (deterministic rounding -- same result every time for
    a given input, so every chain gets identical counts; only cluster
    POSITION is randomized)."""
    ps, pd, pt = run_split["singlet"] / 100, run_split["doublet"] / 100, run_split["triplet"] / 100
    target_residues = int(round(percent / 100 * domain_len))
    avg_len = ps * 1 + pd * 2 + pt * 3
    if target_residues <= 0 or avg_len <= 0:
        return {"singlet": 0, "doublet": 0, "triplet": 0}
    N = max(1, round(target_residues / avg_len))
    n_s = round(ps * N)
    n_d = round(pd * N)
    n_t = max(0, N - n_s - n_d)
    # nudge to hit target_residues exactly, adjusting the largest bucket first
    residues = n_s * 1 + n_d * 2 + n_t * 3
    diff = target_residues - residues
    for key, length in sorted([("triplet", 3), ("doublet", 2), ("singlet", 1)], key=lambda x: -x[1]):
        counts = {"singlet": n_s, "doublet": n_d, "triplet": n_t}
        while diff >= length and length > 0:
            counts[key] += 1
            diff -= length
        n_s, n_d, n_t = counts["singlet"], counts["doublet"], counts["triplet"]
    if diff != 0:
        n_s += diff  # absorb any tiny leftover (+/-1 or 2) into singlets
        n_s = max(0, n_s)
    return {"singlet": n_s, "doublet": n_d, "triplet": n_t}


if A_SCOPE == "whole":
    _whole_len = n_monomers - 1
    _whole_counts = _percent_to_counts(A_PERCENT, A_RUN_SPLIT, _whole_len)
    A_PLAN = {"whole": _whole_counts}
else:
    A_PLAN = {
        "rigid": _percent_to_counts(A_PERCENT_RIGID, A_RUN_SPLIT_RIGID, rigid_length),
        "flex": _percent_to_counts(A_PERCENT_FLEX, A_RUN_SPLIT_FLEX, flex_length),
    }

_domains_to_check = [("whole", n_monomers - 1)] if A_SCOPE == "whole" else [("rigid", rigid_length), ("flex", flex_length)]
for _dom, _dom_len in _domains_to_check:
    _needed = (A_PLAN[_dom]["singlet"] * 1 + A_PLAN[_dom]["doublet"] * 2
               + A_PLAN[_dom]["triplet"] * 3)
    assert _needed <= _dom_len, (
        f"A_PLAN['{_dom}'] needs {_needed} residues but that domain is only "
        f"{_dom_len} residues long"
    )
    _pct_str = f"{100*_needed/_dom_len:.1f}%" if _dom_len > 0 else "n/a (domain is 0 residues long)"
    print(f"[A_PLAN] {_dom}: singlet={A_PLAN[_dom]['singlet']} "
          f"doublet={A_PLAN[_dom]['doublet']} triplet={A_PLAN[_dom]['triplet']}  "
          f"({_needed}/{_dom_len} = {_pct_str})")

# ================= Wall / pseudoatom grid (unchanged from npc.py) =================

wall_dz = 1.0
wall_arc_spacing = 1.0
wall_min_anchor_dist = 1.0 * sigma_wall

output_xyz = "100_80.xyz"
output_data = "100_80.data"

rng = np.random.default_rng()  # swap for np.random.default_rng(seed) for reproducibility


# ================= Ring anchor placement (unchanged) =================

def ring_anchor(chain_id):
    ring = chain_id // chains_per_ring
    slot = chain_id % chains_per_ring
    z_vals = np.linspace(-L_cyl / 2, L_cyl / 2, n_rings)
    z = z_vals[ring]
    theta = 2 * np.pi * slot / chains_per_ring
    return np.array([R_cyl * np.cos(theta), R_cyl * np.sin(theta), z])


# ================= NEW: clustered A-bead pattern generator =================


def _place_runs_in_range(eligible_positions, run_lengths, rng):
    """Place run_lengths (each a block of consecutive A positions) inside
    eligible_positions, spread across the WHOLE eligible range (random
    partition of leftover space into gaps, same idea as before) -- BUT
    now guarantees every INTERNAL gap (the one between two consecutive
    blocks) is at least 1 residue, so two distinct blocks (e.g. a singlet
    next to a doublet) never end up physically adjacent and merge into
    what would incorrectly be read back as one longer run (e.g. a
    'phantom triplet'). Only the leading gap (before the first block) and
    trailing gap (after the last block) are allowed to be 0. If there
    isn't even enough slack to guarantee 1 residue between every block,
    places as many blocks as fit with best-effort separation."""
    a_positions = set()
    L = len(eligible_positions)
    total_block_len = sum(run_lengths)
    if total_block_len > L:
        run_lengths = list(run_lengths)
        while run_lengths and sum(run_lengths) > L:
            run_lengths.pop()
        total_block_len = sum(run_lengths)

    if not run_lengths:
        return a_positions

    n_blocks = len(run_lengths)
    internal_gaps = max(0, n_blocks - 1)
    slack = L - total_block_len

    guaranteed = min(internal_gaps, slack)   # how many internal gaps we can force to >=1
    remaining_slack = slack - guaranteed

    k = n_blocks + 1  # gap slots: [before, internal x (n_blocks-1), after]
    extra = (rng.multinomial(remaining_slack, [1.0 / k] * k)
             if remaining_slack > 0 else np.zeros(k, dtype=int))

    gaps = extra.copy()
    for gi in range(guaranteed):
        gaps[1 + gi] += 1  # add the reserved +1 to internal gap slots (indices 1..internal_gaps)

    cursor = gaps[0]
    for i, run_len in enumerate(run_lengths):
        for j in range(run_len):
            a_positions.add(eligible_positions[cursor + j])
        cursor += run_len + gaps[i + 1]

    return a_positions


def _eligible_positions(domain):
    if domain == "whole":
        return list(range(1, n_monomers))
    elif domain == "rigid":
        return list(range(1, rigid_length + 1))
    else:  # "flex"
        return list(range(rigid_length + 1, n_monomers))


def build_chain_type_sequence(rng):
    """Builds ONE chain's type list (length n_monomers, index 0 = C_anchor)
    using the exact A_PLAN counts (derived from your percentages above).
    Called fresh per chain so the random placement differs chain to
    chain, but the counts are identical for every chain -- exactly what
    the percentages resolve to."""
    types = ["C_anchor"] + [
        "C_rigid" if i <= rigid_length else "C_flex"
        for i in range(1, n_monomers)
    ]

    a_positions = set()
    for domain in A_PLAN:
        plan = A_PLAN[domain]
        blocks = [1] * plan["singlet"] + [2] * plan["doublet"] + [3] * plan["triplet"]
        rng.shuffle(blocks)
        a_positions |= _place_runs_in_range(_eligible_positions(domain), blocks, rng)

    for idx in a_positions:
        types[idx] = "A"

    return types


# ================= Cell list / spatial hash (unchanged) =================

class CellList:
    def __init__(self, cell_size):
        self.cell_size = cell_size
        self.cells = defaultdict(list)

    def _key(self, p):
        return (
            int(np.floor(p[0] / self.cell_size)),
            int(np.floor(p[1] / self.cell_size)),
            int(np.floor(p[2] / self.cell_size)),
        )

    def add(self, p):
        self.cells[self._key(p)].append(p)

    def _neighbor_points(self, p):
        ix, iy, iz = self._key(p)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    bucket = self.cells.get((ix + dx, iy + dy, iz + dz))
                    if bucket:
                        yield from bucket

    def too_close(self, p, min_dist):
        min_dist2 = min_dist * min_dist
        for q in self._neighbor_points(p):
            if np.sum((p - q) ** 2) < min_dist2:
                return True
        return False


# ================= Self-avoiding chain growth (unchanged except A-pattern call) =================

def initialize_all_chains(wall_positions):
    positions = []
    types = []
    global_cells = CellList(cell_size=sigma)
    wall_cells = CellList(cell_size=sigma_wall)
    for w in wall_positions:
        wall_cells.add(w)

    # FIX (bug 2): every anchor position is deterministic and known up
    # front (ring_anchor() doesn't depend on growth at all), so register
    # ALL of them in global_cells before any chain grows -- not just the
    # ones whose full chain has already finished. Otherwise an earlier
    # chain can freely grow a bead right on top of a later chain's
    # not-yet-placed anchor (verified: measured a 0.73 violation from
    # exactly this before the fix).
    for cid in range(n_chains):
        global_cells.add(ring_anchor(cid))

    for chain_id in range(n_chains):

        chain_types = build_chain_type_sequence(rng)

        start = ring_anchor(chain_id)
        success = False

        for chain_attempt in range(20):

            chain = [start]

            for i in range(1, n_monomers):
                placed = False

                for _ in range(500):
                    v = rng.normal(0, 1, 3)
                    v /= np.linalg.norm(v)
                    d = rng.uniform(0.95, 1.05)
                    new = chain[-1] + v * d

                    if new[0] ** 2 + new[1] ** 2 > (R_cyl - sigma / 2) ** 2:
                        continue

                    if any(np.sum((new - p) ** 2) < (0.90 * sigma) ** 2 for p in chain):
                        continue

                    if global_cells.too_close(new, 0.90 * sigma):
                        continue

                    # FIX (bug 1): reject growth positions that land too
                    # close to a wall pseudoatom. Previously ONLY anchor
                    # beads were kept clear of the wall (via the holes
                    # punched in generate_cylinder_shell) -- every other
                    # bead had no wall-exclusion check at all (verified:
                    # measured beads as close as 0.54 to a wall atom
                    # before this fix, 215+ beads within 0.85).
                    if wall_cells.too_close(new, 0.90 * sigma_wall):
                        continue

                    chain.append(new)
                    placed = True
                    break

                if not placed:
                    break

            if len(chain) == n_monomers:
                success = True
                break

        if not success:
            raise RuntimeError(
                f"Chain {chain_id} could not be placed (system too dense)"
            )

        for i, p in enumerate(chain):
            positions.append(p)
            types.append(chain_types[i])
            global_cells.add(p)

    return np.array(positions), types


def initialize_all_chains_with_restarts(wall_positions, max_restarts=30):
    for restart in range(1, max_restarts + 1):
        try:
            return initialize_all_chains(wall_positions)
        except RuntimeError as e:
            print(f"[restart {restart}/{max_restarts}] {e} -- "
                  f"restarting the whole packing from chain 0")
    raise RuntimeError(
        f"Could not pack all {n_chains} chains in {max_restarts} full "
        f"restarts. The geometry (R_cyl={R_cyl}, n_monomers={n_monomers}, "
        f"n_chains={n_chains}) may simply be too dense for this "
        f"random-walk packer -- consider growing chains from both ends, "
        f"or a lattice-based initial packing instead."
    )


# ================= Wall pseudoatom grid (unchanged) =================

def generate_cylinder_shell(anchor_positions):
    n_z = int(L_cyl / wall_dz) + 1
    circumference = 2 * np.pi * R_cyl
    n_theta = int(circumference / wall_arc_spacing)

    z_vals = np.linspace(-L_cyl / 2, L_cyl / 2, n_z)
    anchor_positions = np.asarray(anchor_positions)

    pseudo = []
    for z in z_vals:
        for k in range(n_theta):
            theta = 2 * np.pi * k / n_theta
            p = np.array([R_cyl * np.cos(theta), R_cyl * np.sin(theta), z])
            d2 = np.sum((anchor_positions - p) ** 2, axis=1)
            if d2.min() >= wall_min_anchor_dist ** 2:
                pseudo.append(p)

    return np.asarray(pseudo, dtype=np.float64)


# ================= Writers (unchanged except type_order already dynamic) =================

def write_xyz(fname, pos, pseudo, types):
    with open(fname, "w") as f:
        f.write(f"{len(pos) + len(pseudo)}\n")
        f.write("Cycle 0\n")
        for p, t in zip(pos, types):
            f.write(f"{t} {p[0]:.3f} {p[1]:.3f} {p[2]:.3f}\n")
        for p in pseudo:
            f.write(f"X {p[0]:.3f} {p[1]:.3f} {p[2]:.3f}\n")


def write_lammps_data(fname, pos, pseudo, types):
    n_poly = len(pos)
    n_wall = len(pseudo)
    n_atoms = n_poly + n_wall

    # Atom type numbering built from whatever bead types actually show up,
    # in this fixed preference order, then "wall" last.
    type_order = ["C_anchor", "C_rigid", "C_flex", "A"]
    present = [t for t in type_order if t in set(types)]
    type_id = {t: i + 1 for i, t in enumerate(present)}
    wall_type_id = len(present) + 1
    type_id["wall"] = wall_type_id
    mass_labels = {"C_anchor": "anchor", "C_rigid": "rigid",
                   "C_flex": "flex", "A": "A bead"}

    print("Atom type numbering for this run:")
    for t in present:
        print(f"    type {type_id[t]} = {t}")
    print(f"    type {wall_type_id} = wall")

    bonds = []
    for c in range(n_chains):
        off = c * n_monomers
        for i in range(n_monomers - 1):
            bonds.append((off + i + 1, off + i + 2))

    angles = []
    for c in range(n_chains):
        off = c * n_monomers
        for i in range(n_monomers - 2):
            a, b, cc = off + i + 1, off + i + 2, off + i + 3
            atype = 1 if i == 0 else 2
            angles.append((atype, a, b, cc))

    half = BOX / 2

    with open(fname, "w") as f:
        f.write("LAMMPS data file. Generated by build_cyl_system_tunable.py\n")
        f.write(f" {n_atoms} atoms\n")
        f.write(f" {len(bonds)} bonds\n")
        f.write(f" {len(angles)} angles\n")
        f.write(" 0 dihedrals\n 0 impropers\n")
        f.write(f" {len(present) + 1} atom types\n 1 bond types\n 2 angle types\n")
        f.write(" 0 dihedral types\n 0 improper types\n\n")
        f.write(f" {-half:.6f} {half:.6f}  xlo xhi\n")
        f.write(f" {-half:.6f} {half:.6f}  ylo yhi\n")
        f.write(f" {-half:.6f} {half:.6f}  zlo zhi\n\n")

        f.write(" Masses\n\n")
        for t in present:
            f.write(f" {type_id[t]} 1.000000 # {mass_labels[t]}\n")
        f.write(f" {wall_type_id} 1.000000 # wall\n\n")

        f.write(" Atoms # full\n\n")
        aid = 0
        for p, t in zip(pos, types):
            aid += 1
            f.write(f"{aid} 1 {type_id[t]} 0.000000 {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
        for p in pseudo:
            aid += 1
            f.write(f"{aid} 1 {wall_type_id} 0.000000 {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
        f.write("\n")

        f.write(" Bonds\n\n")
        for i, (a, b) in enumerate(bonds, 1):
            f.write(f"{i} 1 {a} {b}\n")
        f.write("\n")

        f.write(" Angles\n\n")
        for i, (atype, a, b, c) in enumerate(angles, 1):
            f.write(f"{i} {atype} {a} {b} {c}\n")


# ================= Verification helper: check achieved A fractions =================

def report_A_stats(types):
    """Prints achieved singlet/doublet/triplet COUNTS per domain against
    A_PLAN's targets, plus the resulting overall/rigid/flex A percentages,
    so you can confirm the output matches what you set before using the
    files. Since A_PLAN counts are exact and applied identically to every
    chain, per-chain counts should match A_PLAN exactly (any mismatch
    means a block didn't fit its domain -- check the fit-check output at
    startup)."""
    per_chain = [types[c * n_monomers:(c + 1) * n_monomers] for c in range(n_chains)]

    def run_lengths_in(seq_slice):
        runs = []
        i, n = 0, len(seq_slice)
        while i < n:
            if seq_slice[i] == "A":
                j = i
                while j < n and seq_slice[j] == "A":
                    j += 1
                runs.append(j - i)
                i = j
            else:
                i += 1
        return runs

    n_rigid_tot = n_flex_tot = 0
    a_rigid_tot = a_flex_tot = 0
    rigid_runs, flex_runs = [], []
    for chain in per_chain:
        rigid_part = chain[1:rigid_length + 1]
        flex_part = chain[rigid_length + 1:]
        n_rigid_tot += len(rigid_part)
        n_flex_tot += len(flex_part)
        a_rigid_tot += rigid_part.count("A")
        a_flex_tot += flex_part.count("A")
        rigid_runs.extend(run_lengths_in(rigid_part))
        flex_runs.extend(run_lengths_in(flex_part))

    n_a_total = a_rigid_tot + a_flex_tot
    n_total = n_rigid_tot + n_flex_tot

    print("\n=== A-bead placement check (per chain, averaged over all chains) ===")
    for dom_name, runs in (("Rigid", rigid_runs), ("Flex", flex_runs)):
        s = sum(1 for r in runs if r == 1)
        d = sum(1 for r in runs if r == 2)
        t = sum(1 for r in runs if r >= 3)
        plan_key = "whole" if A_SCOPE == "whole" else dom_name.lower()
        plan = A_PLAN.get(plan_key)
        if plan is not None and A_SCOPE == "domain":
            print(f"{dom_name} domain: singlet {s}/{n_chains*plan['singlet']}  "
                  f"doublet {d}/{n_chains*plan['doublet']}  "
                  f"triplet {t}/{n_chains*plan['triplet']}  "
                  f"(target per chain: S{plan['singlet']}/D{plan['doublet']}/T{plan['triplet']})")
        else:
            print(f"{dom_name} domain: singlet {s}  doublet {d}  triplet {t}  "
                  f"(A_SCOPE='whole' -- these split naturally between domains, not separately targeted)")

    print(f"\nOverall A fraction:  {100 * n_a_total / n_total:.1f}%  ({n_a_total}/{n_total})")
    if n_rigid_tot > 0:
        print(f"Rigid-domain A frac: {100 * a_rigid_tot / n_rigid_tot:.1f}%  ({a_rigid_tot}/{n_rigid_tot})")
    else:
        print(f"Rigid-domain A frac: n/a (rigid domain is 0 residues long -- rigidity_percent=0)")
    if n_flex_tot > 0:
        print(f"Flex-domain A frac:  {100 * a_flex_tot / n_flex_tot:.1f}%  ({a_flex_tot}/{n_flex_tot})")
    else:
        print(f"Flex-domain A frac:  n/a (flex domain is 0 residues long)")


# ================= Run =================

def run():
    anchor_positions = [ring_anchor(c) for c in range(n_chains)]
    pseudo = generate_cylinder_shell(anchor_positions)
    print(f"wall pseudoatoms: {len(pseudo)}  (file has 2706)")

    pos, types = initialize_all_chains_with_restarts(pseudo)

    report_A_stats(types)

    write_xyz(output_xyz, pos, pseudo, types)
    write_lammps_data(output_data, pos, pseudo, types)
    print(f"wrote {output_xyz} and {output_data}: "
          f"{len(pos) + len(pseudo)} atoms total")


if __name__ == "__main__":
    run()
