#!/usr/bin/env python3
"""
npc_2d_accessibility.py

True 2D cross-sectional accessibility analysis for NPC (nuclear pore complex)
FG-Nup transport channels, replacing the crude 1D radial-gap approximation
with a proper occupancy-grid / distance-transform approach.

For every trajectory frame and every axial (z) bin:
    1. FG-Nup / polymer beads within the slice are collected and projected
       onto the XY (cross-sectional) plane.
    2. Each bead is inflated by (bead_radius + cargo_radius) -- the Minkowski
       sum radius that a spherical cargo center must avoid.
    3. A high-resolution occupancy grid over the pore cross-section
       (r <= pore_radius) is built and split into accessible / blocked cells.
    4. scipy.ndimage distance transforms are used to compute, for every
       accessible cell, the distance to the nearest blocked/out-of-pore cell.
       The maximum of this distance transform is the radius of the largest
       circle that can be inscribed anywhere in the slice -- i.e. the largest
       cargo that could physically pass through that cross-section.

Per-frame-per-zbin results are saved to CSV, then averaged over frames to
give P_open(z), mean accessible area fraction, mean/variance of the largest
inscribed radius, plus the full distribution.

Works with both atom-numbering schemes used in the cyl_*.data files:
    3-type systems: type 1 = anchor, type 2 = backbone,        type 3 = FG
    4-type systems: type 1 = anchor, type 2 = rigid, type 3 = flexible, type 4 = FG
All polymer bead types are treated as pore-occupying obstacles (the anchor
and backbone/rigid/flexible beads are physically part of the Nup chains and
occlude the channel just as much as the FG beads do); this is controlled by
--occupancy-types if a narrower selection is ever wanted.

Author: analysis pipeline for Yesha's NPC transport project.
"""

import argparse
import os
import sys
import json

import numpy as np
import pandas as pd
import MDAnalysis as mda
from scipy import ndimage
from scipy.spatial import cKDTree

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Publication-quality plot defaults
# ----------------------------------------------------------------------
plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "font.size": 12,
    "axes.labelsize": 13,
    "axes.titlesize": 13,
    "legend.fontsize": 10,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "font.family": "sans-serif",
    "axes.grid": True,
    "grid.alpha": 0.3,
    "savefig.bbox": "tight",
})


# ----------------------------------------------------------------------
# Atom-type scheme detection
# ----------------------------------------------------------------------
def detect_type_scheme(universe):
    """Return a dict describing which numeric LAMMPS type is which role,
    based on the number of distinct atom types (3-type vs 4-type system)."""
    ntypes = len(np.unique(universe.atoms.types))
    if ntypes <= 3:
        scheme = {"1": "anchor", "2": "backbone", "3": "FG"}
        label = "3-type (anchor/backbone/FG)"
    else:
        scheme = {"1": "anchor", "2": "rigid", "3": "flexible", "4": "FG"}
        label = "4-type (anchor/rigid/flexible/FG)"
    return scheme, label, ntypes


# ----------------------------------------------------------------------
# Grid construction (built once, reused for every frame/zbin)
# ----------------------------------------------------------------------
class PoreGrid:
    """High-resolution XY occupancy grid over the pore cross-section."""

    def __init__(self, pore_radius, grid_res):
        self.pore_radius = pore_radius
        self.grid_res = grid_res

        n = int(np.ceil(2 * pore_radius / grid_res)) + 1
        edges = np.linspace(-pore_radius, pore_radius, n)
        # cell-center coordinates
        xc = 0.5 * (edges[:-1] + edges[1:])
        yc = 0.5 * (edges[:-1] + edges[1:])
        self.nx, self.ny = len(xc), len(yc)
        self.xc, self.yc = xc, yc

        XX, YY = np.meshgrid(xc, yc, indexing="xy")
        self.XX, self.YY = XX, YY
        RR = np.sqrt(XX ** 2 + YY ** 2)
        # domain = cells whose center lies within the pore cross-section
        self.domain_mask = RR <= pore_radius          # shape (ny, nx)
        self.n_domain = int(self.domain_mask.sum())
        self.cross_section_area = np.pi * pore_radius ** 2

        # flat coordinates of domain cells, and a KD-tree over them so that
        # for each bead we can cheaply find every grid cell within its
        # inflated radius (much faster than the reverse per-cell query).
        self.domain_xy = np.column_stack(
            (XX[self.domain_mask], YY[self.domain_mask])
        )
        self.domain_tree = cKDTree(self.domain_xy)
        # map from flat-domain index back into (row, col) of the full grid
        self.domain_rc = np.argwhere(self.domain_mask)  # (n_domain, 2) -> row, col

    def occupancy_and_distance(self, bead_xy, effective_radius):
        """
        Given bead centers (Nbeads,2) within this slice and a single
        inflation radius, return:
            accessible_fraction : float
            largest_radius      : float (sigma units)
            accessible_full     : (ny,nx) bool array (True = accessible),
                                   NaN-outside-domain style not applied here;
                                   cells outside the pore are marked False.
            dist_full            : (ny,nx) float array, EDT distance in sigma
                                    (0 outside domain / blocked cells)
        """
        blocked_domain = np.zeros(self.n_domain, dtype=bool)

        if len(bead_xy) > 0 and effective_radius > 0:
            hits = self.domain_tree.query_ball_point(bead_xy, r=effective_radius)
            for idx_list in hits:
                if idx_list:
                    blocked_domain[idx_list] = True

        accessible_domain = ~blocked_domain
        accessible_fraction = accessible_domain.sum() / self.n_domain

        # build full-grid boolean array (outside-pore cells count as blocked
        # for the purposes of the distance transform, which is exactly what
        # we want: the pore wall itself limits how large a circle can be)
        accessible_full = np.zeros((self.ny, self.nx), dtype=bool)
        rows, cols = self.domain_rc[:, 0], self.domain_rc[:, 1]
        accessible_full[rows, cols] = accessible_domain

        # scipy.ndimage distance transform: distance (in grid cells) from
        # every True cell to the nearest False cell
        dist_cells = ndimage.distance_transform_edt(accessible_full)
        dist_sigma = dist_cells * self.grid_res

        largest_radius = float(dist_sigma.max()) if accessible_domain.any() else 0.0

        return accessible_fraction, largest_radius, accessible_full, dist_sigma


# ----------------------------------------------------------------------
# Main per-system analysis
# ----------------------------------------------------------------------
def run_analysis(args):
    os.makedirs(args.outdir, exist_ok=True)

    # ---- load trajectory -------------------------------------------------
    if args.dcd:
        u = mda.Universe(args.data, args.dcd, topology_format="DATA", format="DCD")
    else:
        u = mda.Universe(args.data, topology_format="DATA")
        print("[warn] No DCD supplied -- analyzing the single frame in the "
              "LAMMPS data file only.", file=sys.stderr)

    scheme, scheme_label, ntypes = detect_type_scheme(u)
    print(f"[info] Detected {scheme_label}, {ntypes} atom types, "
          f"{len(u.atoms)} polymer beads, {len(u.trajectory)} frame(s).")

    if args.occupancy_types:
        wanted = set(args.occupancy_types.split(","))
        polymer = u.select_atoms(" or ".join(f"type {t}" for t in wanted))
        print(f"[info] Restricting occupancy to types {wanted} "
              f"({len(polymer)} beads).")
    else:
        polymer = u.atoms  # all beads occupy pore volume

    bead_radius = args.fg_diameter / 2.0
    cargo_radius = args.cargo_diameter / 2.0
    effective_radius = bead_radius + cargo_radius
    print(f"[info] bead_radius={bead_radius:.3f} sigma, "
          f"cargo_radius={cargo_radius:.3f} sigma, "
          f"effective_radius={effective_radius:.3f} sigma")

    # ---- axial bins --------------------------------------------------
    z_lo, z_hi = -args.pore_length / 2.0, args.pore_length / 2.0
    z_edges = np.linspace(z_lo, z_hi, args.nz + 1)
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])

    # ---- grid ----------------------------------------------------------
    grid = PoreGrid(args.pore_radius, args.grid_res)
    print(f"[info] Grid: {grid.nx}x{grid.ny} cells, res={args.grid_res} sigma, "
          f"{grid.n_domain} cells inside pore cross-section.")

    # z-bins selected for XY accessibility-probability heatmaps
    n_select = min(args.n_heatmap_z, args.nz)
    select_idx = np.linspace(0, args.nz - 1, n_select).round().astype(int)
    select_idx = sorted(set(select_idx.tolist()))
    heatmap_accum = {i: np.zeros((grid.ny, grid.nx)) for i in select_idx}
    heatmap_counts = {i: 0 for i in select_idx}

    nframes = len(u.trajectory)
    largest_radius_mat = np.full((args.nz, nframes), np.nan)  # (z, time) heatmap

    records = []
    for fidx, ts in enumerate(u.trajectory):
        pos = polymer.positions  # (N,3), Angstrom-equivalent = sigma units here
        for zi in range(args.nz):
            zlo, zhi = z_edges[zi], z_edges[zi + 1]
            sel = (pos[:, 2] >= zlo) & (pos[:, 2] < zhi)
            bead_xy = pos[sel, :2]

            acc_frac, largest_r, acc_full, dist_full = grid.occupancy_and_distance(
                bead_xy, effective_radius
            )
            is_open = int(largest_r >= cargo_radius)

            records.append({
                "frame": fidx,
                "z_bin": zi,
                "z_center": z_centers[zi],
                "n_beads": int(sel.sum()),
                "accessible_fraction": acc_frac,
                "largest_inscribed_radius": largest_r,
                "channel_open": is_open,
            })

            largest_radius_mat[zi, fidx] = largest_r

            if zi in heatmap_accum:
                heatmap_accum[zi] += acc_full.astype(float)
                heatmap_counts[zi] += 1

        if (fidx + 1) % max(1, nframes // 10) == 0 or fidx == nframes - 1:
            print(f"[info] processed frame {fidx + 1}/{nframes}")

    df = pd.DataFrame.from_records(records)
    per_frame_csv = os.path.join(args.outdir, f"{args.tag}_per_frame_zbin.csv")
    df.to_csv(per_frame_csv, index=False)
    print(f"[info] wrote {per_frame_csv} ({len(df)} rows)")

    # ---- averages over frames per z-bin --------------------------------
    summary = df.groupby("z_bin").agg(
        z_center=("z_center", "first"),
        P_open=("channel_open", "mean"),
        mean_accessible_fraction=("accessible_fraction", "mean"),
        std_accessible_fraction=("accessible_fraction", "std"),
        mean_largest_radius=("largest_inscribed_radius", "mean"),
        var_largest_radius=("largest_inscribed_radius", "var"),
        std_largest_radius=("largest_inscribed_radius", "std"),
        n_frames=("frame", "nunique"),
    ).reset_index()
    summary["tag"] = args.tag
    summary["cargo_diameter"] = args.cargo_diameter
    summary["fg_diameter"] = args.fg_diameter
    summary_csv = os.path.join(args.outdir, f"{args.tag}_zbin_summary.csv")
    summary.to_csv(summary_csv, index=False)
    print(f"[info] wrote {summary_csv}")

    # ---- distribution CSV (pooled largest_inscribed_radius) ------------
    dist_csv = os.path.join(args.outdir, f"{args.tag}_largest_radius_distribution.csv")
    df[["frame", "z_bin", "z_center", "largest_inscribed_radius"]].to_csv(
        dist_csv, index=False
    )
    print(f"[info] wrote {dist_csv}")

    # =====================================================================
    # PLOTS
    # =====================================================================
    tag = args.tag

    # 1. P_open(z)
    fig, ax = plt.subplots(figsize=(6, 4.2))
    ax.plot(summary["z_center"], summary["P_open"], marker="o", lw=2, color="#1f6f8b")
    ax.axhline(1.0, color="gray", lw=0.8, ls="--")
    ax.set_xlabel(r"$z$ ($\sigma$)")
    ax.set_ylabel(r"$P_{\mathrm{open}}(z)$")
    ax.set_title(f"Channel-open probability along pore axis ({tag})")
    ax.set_ylim(-0.02, 1.02)
    fig.savefig(os.path.join(args.outdir, f"{tag}_1_Popen_vs_z.png"))
    plt.close(fig)

    # 2. mean accessible area fraction vs z
    fig, ax = plt.subplots(figsize=(6, 4.2))
    ax.plot(summary["z_center"], summary["mean_accessible_fraction"],
            marker="o", lw=2, color="#c1440e")
    ax.fill_between(
        summary["z_center"],
        summary["mean_accessible_fraction"] - summary["std_accessible_fraction"].fillna(0),
        summary["mean_accessible_fraction"] + summary["std_accessible_fraction"].fillna(0),
        color="#c1440e", alpha=0.2,
    )
    ax.set_xlabel(r"$z$ ($\sigma$)")
    ax.set_ylabel("Mean accessible area fraction")
    ax.set_title(f"Cross-sectional accessible area fraction ({tag})")
    ax.set_ylim(0, 1.02)
    fig.savefig(os.path.join(args.outdir, f"{tag}_2_accessible_fraction_vs_z.png"))
    plt.close(fig)

    # 3. mean largest inscribed radius vs z (with variance band)
    fig, ax = plt.subplots(figsize=(6, 4.2))
    std = summary["std_largest_radius"].fillna(0)
    ax.plot(summary["z_center"], summary["mean_largest_radius"],
            marker="o", lw=2, color="#2e7d32")
    ax.fill_between(summary["z_center"],
                     summary["mean_largest_radius"] - std,
                     summary["mean_largest_radius"] + std,
                     color="#2e7d32", alpha=0.2, label="±1 s.d.")
    ax.axhline(cargo_radius, color="k", lw=1.2, ls="--", label="cargo radius")
    ax.set_xlabel(r"$z$ ($\sigma$)")
    ax.set_ylabel(r"Largest inscribed radius ($\sigma$)")
    ax.set_title(f"Largest passable circle along pore axis ({tag})")
    ax.legend()
    fig.savefig(os.path.join(args.outdir, f"{tag}_3_largest_radius_vs_z.png"))
    plt.close(fig)

    # 4. heatmaps: accessibility probability in XY at selected z
    n_panels = len(select_idx)
    ncols = min(5, n_panels)
    nrows = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 3.0 * nrows),
                              squeeze=False)
    extent = [-args.pore_radius, args.pore_radius, -args.pore_radius, args.pore_radius]
    for k, zi in enumerate(select_idx):
        r, c = divmod(k, ncols)
        ax = axes[r][c]
        prob = heatmap_accum[zi] / max(1, heatmap_counts[zi])
        prob_masked = np.where(grid.domain_mask, prob, np.nan)
        im = ax.imshow(prob_masked, extent=extent, origin="lower",
                        cmap="viridis", vmin=0, vmax=1)
        circ = plt.Circle((0, 0), args.pore_radius, fill=False, color="white", lw=1)
        ax.add_patch(circ)
        ax.set_title(f"z={z_centers[zi]:.1f}$\\sigma$", fontsize=10)
        ax.set_xticks([]); ax.set_yticks([])
    for k in range(n_panels, nrows * ncols):
        r, c = divmod(k, ncols)
        axes[r][c].axis("off")
    fig.suptitle(f"XY accessibility probability at selected z ({tag})", y=1.02)
    cbar = fig.colorbar(im, ax=axes, shrink=0.8, label="Accessibility probability")
    fig.savefig(os.path.join(args.outdir, f"{tag}_4_xy_accessibility_heatmaps.png"),
                bbox_inches="tight")
    plt.close(fig)

    # 5. heatmap: largest free-space radius as function of (z, time/frame)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    im = ax.imshow(
        largest_radius_mat,
        aspect="auto", origin="lower", cmap="magma",
        extent=[0, nframes, z_lo, z_hi],
    )
    fig.colorbar(im, ax=ax, label=r"Largest inscribed radius ($\sigma$)")
    ax.set_xlabel("Frame")
    ax.set_ylabel(r"$z$ ($\sigma$)")
    ax.set_title(f"Largest free-space radius: z vs time ({tag})")
    fig.savefig(os.path.join(args.outdir, f"{tag}_5_radius_z_time_heatmap.png"))
    plt.close(fig)

    # 6. histogram of largest inscribed radius (pooled)
    fig, ax = plt.subplots(figsize=(6, 4.2))
    ax.hist(df["largest_inscribed_radius"], bins=40, color="#5b3a9e", alpha=0.85)
    ax.axvline(cargo_radius, color="k", lw=1.5, ls="--", label="cargo radius")
    ax.set_xlabel(r"Largest inscribed radius ($\sigma$)")
    ax.set_ylabel("Count (all frames × z-bins)")
    ax.set_title(f"Distribution of largest inscribed radius ({tag})")
    ax.legend()
    fig.savefig(os.path.join(args.outdir, f"{tag}_6_largest_radius_histogram.png"))
    plt.close(fig)

    # ---- also stash a machine-readable run-config record ----------------
    with open(os.path.join(args.outdir, f"{tag}_run_config.json"), "w") as fh:
        json.dump(vars(args), fh, indent=2)

    print(f"[done] {tag}: mean P_open = {summary['P_open'].mean():.3f}, "
          f"mean accessible fraction = {summary['mean_accessible_fraction'].mean():.3f}")

    return summary_csv, per_frame_csv


def build_argparser():
    p = argparse.ArgumentParser(
        description="True 2D cross-sectional accessibility analysis of NPC "
                    "FG-Nup transport channels (occupancy grid + distance "
                    "transform, no 1D radial-gap approximation)."
    )
    p.add_argument("--data", required=True, help="LAMMPS data file (topology + coords)")
    p.add_argument("--dcd", default=None, help="DCD trajectory (optional; if omitted, "
                    "only the single frame in --data is analyzed)")
    p.add_argument("--cargo-diameter", type=float, required=True,
                    help="Cargo diameter, sigma units")
    p.add_argument("--fg-diameter", type=float, required=True,
                    help="FG/polymer bead diameter, sigma units")
    p.add_argument("--pore-radius", type=float, default=10.0, help="Pore radius, sigma")
    p.add_argument("--pore-length", type=float, default=25.0, help="Pore length, sigma "
                    "(bins span -length/2 to +length/2)")
    p.add_argument("--nz", type=int, default=25, help="Number of axial bins")
    p.add_argument("--grid-res", type=float, default=0.25,
                    help="XY occupancy grid resolution, sigma per cell")
    p.add_argument("--n-heatmap-z", type=int, default=5,
                    help="Number of z-slices shown in the XY accessibility heatmap panel")
    p.add_argument("--occupancy-types", default=None,
                    help="Comma-separated LAMMPS atom type numbers to treat as "
                    "pore-occupying obstacles (default: all types in the data file)")
    p.add_argument("--outdir", default="./npc_accessibility_out", help="Output directory")
    p.add_argument("--tag", default="system",
                    help="Label used in filenames/titles/legends, e.g. rigidity percentage "
                    "such as 'rigidity_50'")
    return p


if __name__ == "__main__":
    args = build_argparser().parse_args()
    run_analysis(args)