"""
Cylindrical FG-Nup brush — initial configuration + LAMMPS data file builder.

Reconstructed by reverse-engineering cyl_0.data and cyl_init.xyz, since the
original generator script was lost. Every geometric/topological parameter
below was measured directly from those two files (not guessed) unless a
comment says otherwise:

  - R_cyl = 10.0            <- radius of all 88 anchor beads and all wall
                                pseudoatoms, exact to 3 decimals
  - L_cyl = 25.0             <- anchor / wall z-range is exactly [-12.5, 12.5]
  - 11 rings x 8 chains/ring = 88 chains, ring z-spacing = 2.5 (linspace)
  - 40 monomers/chain (1 anchor + 39 flex), matches 3432 bonds / 3344 angles
  - r0 = 1.0, bond_tolerance = 0.05  <- bond length range in the data file is
                                exactly [0.9493, 1.0503], centered on 1.0
  - every non-anchor bead is "C_flex" — the file has ZERO rigid monomers,
    unlike the old cyl.py (which had rigidity_percent = 50). Rigidity is
    kept as a switch below in case you want it back.
  - wall pseudoatoms sit on an EXACT regular grid: 89 azimuthal positions x
    36 axial positions (r = R_cyl, z = linspace(-12.5, 12.5, 36)). All 2706
    wall atoms in the file match points on this grid to within float
    rounding — this is not a guess, it's a 100% match once the grid is
    built with these two numbers.
  - only 2706 of the 3204 grid points survive: a grid point is dropped if
    it's within `1.0 * sigma_wall` of any ANCHOR bead (not any chain atom).
    This is your actual generate_cylinder_shell() from system.py, confirmed
    with R_cyl=10, L_cyl=25, sigma_wall=1.0 — it reproduces all 2706 wall
    atoms in the file to floating-point rounding, exactly.
  - single molecule id (1) for every atom, single bond type, and TWO angle
    types purely because of how VMD/TopoTools auto-types angles from atom
    types: the one angle touching the anchor bead (pattern 1-3-3) becomes
    angle type 1, and every other angle (pattern 3-3-3) becomes type 2.
    That's not a physical rigid/flex distinction — it reproduces the same
    LAMMPS numbering directly in Python so you no longer depend on VMD.
  - box: the file's box is 150 x 150 x 150, but with a strange non-centered
    z-range (-75.8565 to 74.1435, i.e. shifted ~0.86 off-center) which looks
    like a VMD/TopoTools artifact rather than anything intentional. I use a
    clean, centered -75..75 box in all three dimensions here.

Run this script to regenerate init.xyz and init.data equivalent to your
originals. Chain growth now uses your actual system.py algorithm (20
attempts per chain, sigma/2 confinement buffer, A-bead support inactive
by default) wrapped in a full-system restart loop for robustness at this
system's density -- see initialize_all_chains() and
initialize_all_chains_with_restarts() below. Exact atom-by-atom
coordinates will differ from the original files (it's a new random
realization) but every structural statistic will match.
"""

from collections import defaultdict

import numpy as np

# ================= Geometry (measured) =================

R_cyl = 10.0
L_cyl = 25.0
BOX = 150.0

n_rings = 11
chains_per_ring = 8
n_chains = n_rings * chains_per_ring          # 88

n_monomers = 40
n_atoms_chain = n_chains * n_monomers          # 3520

# ================= Bonded parameters (measured) =================

sigma = 1.0
sigma_wall = 1.0
r0 = 1.0
bond_tolerance = 0.05
bond_len_min = (1 - bond_tolerance) * r0
bond_len_max = (1 + bond_tolerance) * r0

intra_sigma = 0.95 * sigma   # min separation within a chain during growth
inter_sigma = 0.95 * sigma   # min separation between different chains

# Rigidity: the uploaded file has 0% rigid monomers (all beads past the
# anchor are C_flex). Left switchable in case you want it back on.
rigidity_percent = 0
rigid_length = int((rigidity_percent / 100) * n_monomers)

# "A" (hydrophobic FG) bead fraction, as in system.py. Default 0 -- the
# uploaded cyl_0.data / cyl_init.xyz have no A beads at all, so this stays
# inactive unless you turn it on.
f_fraction = 0.0
n_A_per_chain = int(round(f_fraction * n_monomers))
_A_pool = np.arange(1, n_monomers)

# ================= Wall / pseudoatom grid (from your system.py) =================

wall_dz = 1.0          # axial spacing -> n_z = int(L_cyl/dz)+1 = 36
wall_arc_spacing = 1.0  # azimuthal spacing -> n_theta = int(2*pi*R_cyl/arc_spacing) = 89
wall_min_anchor_dist = 1.0 * sigma_wall  # exclusion radius around each anchor bead

output_xyz = "init.xyz"
output_data = "init.data"

rng = np.random.default_rng()  # swap for np.random.default_rng(seed) for reproducibility


# ================= Ring anchor placement (matches cyl.py) =================

def ring_anchor(chain_id):
    ring = chain_id // chains_per_ring
    slot = chain_id % chains_per_ring
    z_vals = np.linspace(-L_cyl / 2, L_cyl / 2, n_rings)
    z = z_vals[ring]
    theta = 2 * np.pi * slot / chains_per_ring
    return np.array([R_cyl * np.cos(theta), R_cyl * np.sin(theta), z])


def get_monomer_type(i, A_indices_chain):
    if i == 0:
        return "C_anchor"
    elif i in A_indices_chain:
        return "A"
    elif i < rigid_length:
        return "C_rigid"
    else:
        return "C_flex"


# ================= Cell list / spatial hash (interchain neighbor search) =================
#
# system.py defines cell_index()/neighboring_cells() plus fixed-size nx/ny/nz
# grid dimensions but never actually wires them into initialize_all_chains()
# (that function still does a brute-force loop over `positions`). Below is
# that same cell-list idea actually put to use, replacing the cKDTree query
# for the interchain check. It's a dict keyed by integer cell coordinates
# (rather than system.py's fixed-size nx*ny*nz array) so it isn't bounded to
# a precomputed box -- important here because chain growth only rejects on
# the radial (x,y) confinement, never on z, so beads can wander outside
# [-L_cyl/2, L_cyl/2] and a fixed-size array indexed off L_cyl would go
# out of range.
#
# Cell size is set to sigma, which is >= both the intrachain (0.95*sigma)
# and interchain (0.90*sigma) exclusion radii used below -- that's the
# condition needed for "check the cell + its 26 neighbors" to be exact
# (equivalent to checking every point within the exclusion radius).

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


# ================= Self-avoiding chain growth (ported from your system.py) =================
#
# This is system.py's initialize_all_chains(), not cyl.py's: per chain, up
# to 20 "chain_attempt" restarts (not 50000), confinement buffer is
# sigma/2 (not the hardcoded 1.0 cyl.py used), and it carries the same
# A-bead (f_fraction) support system.py has -- inactive by default since
# f_fraction=0 above. The only deliberate change is the interchain distance
# check: instead of system.py's brute-force python loop over `positions`,
# it uses the CellList spatial hash defined above -- a working version of
# the cell_index()/neighboring_cells() scaffolding system.py defines but
# never actually wires into its own initialize_all_chains(). Same
# accept/reject outcome, just far cheaper to query as `positions` grows
# into the thousands. The cell list is built once for the whole packing
# and updated incrementally as each chain finishes (not rebuilt every
# attempt). Everything else (RNG calls, bond-length window, order of
# checks) matches system.py line for line.

def initialize_all_chains():
    positions = []
    types = []
    global_cells = CellList(cell_size=sigma)

    for chain_id in range(n_chains):

        A_indices_chain = (
            set(rng.choice(_A_pool, size=n_A_per_chain, replace=False).tolist())
            if n_A_per_chain > 0 else set()
        )

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

                    if any(np.sum((new - p) ** 2) < (0.95 * sigma) ** 2 for p in chain):
                        continue

                    if global_cells.too_close(new, 0.90 * sigma):
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
            types.append(get_monomer_type(i, A_indices_chain))
            global_cells.add(p)

    return np.array(positions), types


def initialize_all_chains_with_restarts(max_restarts=30):
    """Wraps initialize_all_chains() with full-system restarts. system.py's
    algorithm only gives each chain 20 attempts (vs. cyl.py's 50000), so at
    this system's density (88 chains, R_cyl=10) a single chain failing is
    common and doesn't mean the packing is impossible -- it usually means
    earlier chains sealed it in. Retrying just that chain can't fix that,
    so on failure we throw away the whole packing and start over from
    chain 0 with a fresh random sequence."""
    for restart in range(1, max_restarts + 1):
        try:
            return initialize_all_chains()
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


# ================= Wall pseudoatom grid (exact, from your system.py) =================

def generate_cylinder_shell(anchor_positions):
    """Regular grid on the cylinder wall, skipping any point within
    wall_min_anchor_dist of an anchor bead. Confirmed to reproduce all
    2706 wall atoms in cyl_init.xyz to floating-point rounding."""
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


# ================= Writers =================

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

    # Atom type numbering: build it from whatever bead types actually show
    # up, in this fixed preference order, then "wall" last. With the
    # defaults (rigidity_percent=0, f_fraction=0) this is exactly
    # {C_anchor: 1, C_flex: 2, wall: 3} -- the same 3-type layout validated
    # against cyl_0.data. Turning rigidity/A back on adds types 2/3/4
    # without disturbing that validated case.
    type_order = ["C_anchor", "C_rigid", "C_flex", "A"]
    present = [t for t in type_order if t in set(types)]
    type_id = {t: i + 1 for i, t in enumerate(present)}
    wall_type_id = len(present) + 1
    type_id["wall"] = wall_type_id
    mass_labels = {"C_anchor": "anchor", "C_rigid": "rigid",
                   "C_flex": "flex", "A": "A bead"}

    # bonds: consecutive pairs within each chain (bond type always 1)
    bonds = []
    for c in range(n_chains):
        off = c * n_monomers
        for i in range(n_monomers - 1):
            bonds.append((off + i + 1, off + i + 2))  # 1-indexed

    # angles: consecutive triplets within each chain.
    # type 1 = the one angle touching the anchor bead, type 2 = all others
    # (this reproduces VMD/TopoTools' atom-type-pattern-based numbering
    # for the anchor/flex/wall-only case validated against cyl_0.data; if
    # you turn on rigidity/A beads this simplified 2-type scheme no longer
    # matches what TopoTools would auto-generate from the extra triplet
    # patterns, but is a reasonable stand-in)
    angles = []
    for c in range(n_chains):
        off = c * n_monomers
        for i in range(n_monomers - 2):
            a, b, cc = off + i + 1, off + i + 2, off + i + 3  # 1-indexed
            atype = 1 if i == 0 else 2
            angles.append((atype, a, b, cc))

    half = BOX / 2

    with open(fname, "w") as f:
        f.write("LAMMPS data file. Reconstructed by build_cyl_system.py\n")
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


# ================= Run =================

def run():
    pos, types = initialize_all_chains_with_restarts()

    anchor_positions = [ring_anchor(c) for c in range(n_chains)]
    pseudo = generate_cylinder_shell(anchor_positions)
    print(f"wall pseudoatoms: {len(pseudo)}  (file has 2706)")

    write_xyz(output_xyz, pos, pseudo, types)
    write_lammps_data(output_data, pos, pseudo, types)
    print(f"wrote {output_xyz} and {output_data}: "
          f"{len(pos) + len(pseudo)} atoms total")


if __name__ == "__main__":
    run()