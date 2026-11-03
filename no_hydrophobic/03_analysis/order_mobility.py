#!/usr/bin/env python3
"""
Nematic order and bead mobility vs distance r from the pore axis, inside the pore (|z| <= 12.5), cargo-free brush.
Tests the idea: f_rigid = 0.5 = ordered stiff scaffold near the wall + disordered, mobile chain ends in the core.
python3 + numpy + scipy only.

  python3 order_mobility.py --dt 0.2 --out om_noA \
     --sys 0=DATA:DCD[:REF] 25=DATA:DCD[:REF] 50=... 75=... 100=...

  DATA  data file with the same atom order as the DCD (write_data output of the run that made the DCD)
  DCD   cargo-free trajectory. For MOBILITY use frames <= 1 tau apart (fast dump, e.g. every 0.2 tau).
        Nematic order works with any DCD (your long equilibrium DCD is best).
  REF   original data file with '# name' comments in Masses (1 anchor, 2 rigid, 3 flexible, 4 shell, 5 A),
        needed if DATA has no names.
  --dt  tau between DCD frames.  --lags  mobility lag times in tau (default 1 5).

What it computes, per radial bin (default 1 sigma wide, r = 0..10), averaged over frames:
  ORDER (from bond vectors u = unit vector along each polymer bond, bond midpoint decides the bin)
    S        nematic order = largest eigenvalue of Q = <(3 u u - I)/2>.  0 = random, 1 = all bonds parallel.
    P2_rad   <(3 cos^2 a - 1)/2>, a = angle between bond and the radial direction (axis -> bead).
             +1 = bonds point straight out from the wall (stretched brush), -0.5 = bonds lie perpendicular to r.
    P2_ax    the same with the pore axis (z). +1 = bonds along the pore.
    Bonds are split by type: 'rigid' (both beads rigid), 'flex' (any flexible or A bead), 'all'.
  MOBILITY (beads binned by their r at the start of each lag)
    D_eff(lag) = <|r(t+lag) - r(t)|^2> / (6 lag)   in sigma^2/tau. Larger = more mobile, easier to push aside.
    Split by bead name: rigid, flex, A, all.  Also split by place along the chain: 'base' = first half from the
    anchor, 'tip' = second half (free end).
Outputs
  OUT_order.csv     system, r, group, S, P2_rad, P2_ax, n_bonds
  OUT_mobility.csv  system, r, group, lag, D_eff, n
  OUT.png           quick look: S(r), P2_rad(r), D_eff(r) for each system
"""
import argparse, struct, sys
from collections import defaultdict
import numpy as np

p = argparse.ArgumentParser()
p.add_argument("--sys", nargs="+", required=True)
p.add_argument("--ref", default=None)
p.add_argument("--dt", type=float, required=True)
p.add_argument("--lags", type=float, nargs="+", default=[1.0, 5.0])
p.add_argument("--skip", type=float, default=0.0, help="fraction of frames dropped at the start")
p.add_argument("--every", type=int, default=1, help="use every n-th frame for ORDER (mobility uses all)")
p.add_argument("--zpore", type=float, default=12.5)
p.add_argument("--rmax", type=float, default=10.0)
p.add_argument("--dr", type=float, default=1.0)
p.add_argument("--out", default="om")
a = p.parse_args()


def sections(L):
    idx = {}
    for k, l in enumerate(L):
        s = l.split("#")[0].strip()
        if s in ("Masses", "Atoms", "Bonds", "Angles", "Velocities"): idx[s] = k
    return idx


def block(L, k):
    i = k + 1
    while not L[i].strip(): i += 1
    j = i
    while j < len(L) and L[j].strip(): j += 1
    return i, j


def names_from(fn):
    L = open(fn).read().split("\n"); i, j = block(L, sections(L)["Masses"]); nm = {}
    for l in L[i:j]:
        b, _, c = l.partition("#"); t = int(b.split()[0])
        nm[t] = int(c.split()[0]) if c.strip() else None
    return nm


def topology(data, ref):
    nm = names_from(ref)
    if any(v is None for v in nm.values()): sys.exit("%s has no '# name' comments in Masses - add :REF" % ref)
    L = open(data).read().split("\n"); idx = sections(L)
    i, j = block(L, idx["Atoms"]); ids, ty = [], []
    for l in L[i:j]:
        q = l.split("#")[0].split(); ids.append(int(q[0])); ty.append(int(q[2]))
    ids = np.array(ids); o = np.argsort(ids); ids, ty = ids[o], np.array(ty)[o]
    pos = {d: k for k, d in enumerate(ids)}; name = np.array([nm[t] for t in ty])
    cnt = defaultdict(int)
    for t in ty: cnt[t] += 1
    single = np.array([cnt[t] == 1 for t in ty]); name[single] = -1          # cargo (if any) ignored
    i, j = block(L, idx["Bonds"]); bonds = []; adj = defaultdict(list)
    for l in L[i:j]:
        q = l.split(); u, v = pos[int(q[2])], pos[int(q[3])]; bonds.append((u, v)); adj[u].append(v); adj[v].append(u)
    bonds = np.array(bonds)
    seq = np.full(len(ty), -1); clen = np.zeros(len(ty), int)
    for a0 in np.where(name == 1)[0]:                                          # walk each chain from its anchor
        path, prev, cur = [a0], None, a0
        while True:
            nx = [n for n in adj[cur] if n != prev]
            if not nx: break
            prev, cur = cur, nx[0]; path.append(cur)
        seq[path] = np.arange(len(path)); clen[path] = len(path) - 1
    poly = np.isin(name, [2, 3, 5]) & (seq > 0)
    bonds = bonds[poly[bonds[:, 0]] & poly[bonds[:, 1]]]
    return name, seq, clen, poly, bonds


def dcd(fn):
    try:
        from MDAnalysis.lib.formats.libdcd import DCDFile
        f = DCDFile(fn)
        def read(k):
            f.seek(k); return f.read().xyz.astype(np.float32)
        return f.n_frames, f.header["natoms"], read
    except ImportError:
        pass
    fh = open(fn, "rb"); head = fh.read(100)
    e = "<" if struct.unpack("<i", head[:4])[0] == 84 else ">"
    ic = struct.unpack(e + "20i", head[8:88]); cell = ic[10] != 0
    fh.seek(92); n = struct.unpack(e + "i", fh.read(4))[0]; pos = 92 + 4 + n + 4
    fh.seek(pos + 4); N = struct.unpack(e + "i", fh.read(4))[0]; pos += 12
    blk = 4 * N + 8; fs = (56 if cell else 0) + 3 * blk
    fh.seek(0, 2); nfr = (fh.tell() - pos) // fs
    def read(k):
        X = np.empty((N, 3), dtype=np.float32); b = pos + k * fs + (56 if cell else 0)
        for c in range(3):
            fh.seek(b + c * blk + 4); X[:, c] = np.frombuffer(fh.read(4 * N), dtype=e + "f4")
        return X
    return nfr, N, read


rb = np.arange(0, a.rmax + 1e-9, a.dr); nb = len(rb) - 1; rc = 0.5 * (rb[1:] + rb[:-1])
order_rows, mob_rows = [], []
for item in a.sys:
    lab, rest = item.split("=", 1); parts = rest.split(":")
    data, dc = parts[0], parts[1]; ref = parts[2] if len(parts) > 2 else (a.ref or data)
    name, seq, clen, poly, bonds = topology(data, ref)
    nfr, N, rd = dcd(dc)
    if N != len(name): sys.exit("%s: DCD has %d atoms, data file %d" % (lab, N, len(name)))
    k0 = int(a.skip * nfr); frames = range(k0, nfr)
    print("%s: %d frames, %d polymer beads, %d polymer bonds" % (lab, nfr - k0, poly.sum(), len(bonds)), flush=True)
    rig_b = (name[bonds[:, 0]] == 2) & (name[bonds[:, 1]] == 2)
    groups = {"all": np.ones(len(bonds), bool), "rigid": rig_b, "flex": ~rig_b}
    Q = {g: np.zeros((nb, 3, 3)) for g in groups}; PR = {g: np.zeros(nb) for g in groups}
    PA = {g: np.zeros(nb) for g in groups}; NB = {g: np.zeros(nb) for g in groups}
    for k in list(frames)[::a.every]:
        X = rd(k); u = X[bonds[:, 1]] - X[bonds[:, 0]]; u /= np.linalg.norm(u, axis=1)[:, None]
        m = 0.5 * (X[bonds[:, 1]] + X[bonds[:, 0]]); r = np.hypot(m[:, 0], m[:, 1])
        ok = (np.abs(m[:, 2]) <= a.zpore) & (r < a.rmax); bi = np.digitize(r, rb) - 1
        rhat = np.c_[m[:, 0], m[:, 1], np.zeros(len(m))] / np.maximum(r, 1e-6)[:, None]
        c_rad = np.einsum("ij,ij->i", u, rhat); c_ax = u[:, 2]
        uu = 1.5 * u[:, :, None] * u[:, None, :] - 0.5 * np.eye(3)[None]
        for g, gm in groups.items():
            s = ok & gm
            np.add.at(Q[g], bi[s], uu[s]); np.add.at(PR[g], bi[s], 1.5 * c_rad[s] ** 2 - 0.5)
            np.add.at(PA[g], bi[s], 1.5 * c_ax[s] ** 2 - 0.5); np.add.at(NB[g], bi[s], 1)
    for g in groups:
        for b in range(nb):
            n = NB[g][b]
            if n < 10: continue
            S = np.linalg.eigvalsh(Q[g][b] / n).max()
            order_rows.append(dict(system=lab, r=rc[b], group=g, S=S, P2_rad=PR[g][b] / n, P2_ax=PA[g][b] / n, n_bonds=int(n)))
    # mobility
    half = (seq > 0) & (seq <= clen / 2)
    mgroups = {"all": poly, "rigid": poly & (name == 2), "flex": poly & (name == 3), "A": poly & (name == 5),
               "base": poly & half, "tip": poly & ~half}
    for lag in a.lags:
        L = int(round(lag / a.dt))
        if L < 1 or L >= nfr - k0:
            print("   lag %g tau: not possible with dt %g and %d frames - skipped" % (lag, a.dt, nfr - k0)); continue
        S2 = {g: np.zeros(nb) for g in mgroups}; NN = {g: np.zeros(nb) for g in mgroups}
        starts = list(range(k0, nfr - L, max(1, L // 2)))
        for k in starts:
            X0 = rd(k); X1 = rd(k + L); r0 = np.hypot(X0[:, 0], X0[:, 1])
            ok = (np.abs(X0[:, 2]) <= a.zpore) & (r0 < a.rmax); bi = np.digitize(r0, rb) - 1; d2 = ((X1 - X0) ** 2).sum(1)
            for g, gm in mgroups.items():
                s = ok & gm; np.add.at(S2[g], bi[s], d2[s]); np.add.at(NN[g], bi[s], 1)
        for g in mgroups:
            for b in range(nb):
                if NN[g][b] >= 10:
                    mob_rows.append(dict(system=lab, r=rc[b], group=g, lag=lag, D_eff=S2[g][b] / NN[g][b] / (6 * lag), n=int(NN[g][b])))
    print("   done", flush=True)

import csv
for nm_, rows in (("order", order_rows), ("mobility", mob_rows)):
    if rows:
        with open("%s_%s.csv" % (a.out, nm_), "w") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); [w.writerow(r) for r in rows]
try:
    import pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    O = pd.DataFrame(order_rows); M = pd.DataFrame(mob_rows)
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.2))
    for s, g in O[O.group == "all"].groupby("system", sort=False):
        ax[0].plot(g.r, g.S, marker="o", ms=3, label=s); ax[1].plot(g.r, g.P2_rad, marker="o", ms=3, label=s)
    if len(M):
        lag = M.lag.min()
        for s, g in M[(M.group == "all") & (M.lag == lag)].groupby("system", sort=False):
            ax[2].plot(g.r, g.D_eff, marker="o", ms=3, label=s)
        ax[2].set(ylabel="D_eff (lag %g tau)" % lag)
    ax[0].set(ylabel="nematic order S"); ax[1].set(ylabel="P2 (bond vs radial)"); ax[1].axhline(0, color="0.6", lw=0.6)
    for x in ax: x.set_xlabel("r (sigma)")
    ax[0].legend(fontsize=7, title="f_rigid"); fig.tight_layout(); fig.savefig(a.out + ".png", dpi=200)
except Exception as ex:
    print("plot skipped:", ex)
print("saved %s_order.csv, %s_mobility.csv, %s.png" % (a.out, a.out, a.out))
