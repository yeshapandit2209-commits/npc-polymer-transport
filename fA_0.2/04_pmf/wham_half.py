#!/usr/bin/env python3
"""
WHAM for umbrella windows written by setup_umbrella.py.
Half pore (windows z0 >= 0): mirrored about z = 0.  Full pore (windows on both sides, e.g. --zmin -50):
no mirroring; reference = upper bulk, and the lower bulk is printed as a convergence check.

    python wham_half.py 50_tp/umbrella_s4.0                 # one system, one sigma
    python wham_half.py 0/umbrella_s4.0 50/umbrella_s4.0 100/umbrella_s4.0   # overlay several

Reads   <dir>/rep*/z*/colvar.dat   (columns: step x y z; header holds z0 and k)
Does    WHAM per replicate (Kumar et al. 1992), U_i = 0.5 k (z - z0_i)^2, kT = 1
        PMF = 0 in BULK: mean over the far end z >= zmax - REF (default 5 sigma)
        (never the midplane end - that is the barrier top)
        mirror: PMF(-z) = PMF(z)
        error bars = SD over replicates; convergence = first half vs second half of every window
Writes  pmf_<system>_s<sigma>.csv  (z, pmf_mean, pmf_sd, pmf_pooled, pmf_half1, pmf_half2, n_reps)
        pmf_<system>_s<sigma>_rz.csv/.png : 2D map W(r, z) from the same windows (r is not biased,
        so the radial profile at any height comes free; W includes the 1/r area correction)
        pmf_<system>_s<sigma>.png, and pmf_compare.png when several dirs are given
Checks  prints windows with poor overlap (neighbouring histograms that barely touch)
"""
import argparse, glob, os, re, warnings
import numpy as np
warnings.filterwarnings("ignore", category=RuntimeWarning)


def read_window(path):
    z0 = k = None
    with open(path) as f:
        head = f.readline()
    m = re.search(r"z0 = ([-\d.eE+]+), k = ([-\d.eE+]+)", head)
    if m:
        z0, k = float(m.group(1)), float(m.group(2))
    else:                                  # fall back on the folder name z<z0>
        z0 = float(re.search(r"z(-?[\d.]+)$", os.path.dirname(path)).group(1)); k = 10.0
    d = np.loadtxt(path, comments="#", ndmin=2)
    if not len(d):
        return z0, k, np.array([])
    read_window.last_r = np.hypot(d[:, 1], d[:, 2])        # radial coordinate of the same samples
    return z0, k, d[:, 3]


def wham(windows, bins, tol=1e-10, max_iter=200000, return_f=False):
    """windows: list of (z0, k, z_samples). Returns (bin centres, -ln P(z)) (nan where empty).
    Solved as a convex minimisation (Zhu & Hummer 2012) with Newton's method: exact, fast, and
    stable for 100+ windows (the plain self-consistent iteration drifts / stalls there)."""
    c = 0.5 * (bins[1:] + bins[:-1])
    H = np.array([np.histogram(z, bins)[0] for _, _, z in windows], float)
    N = H.sum(axis=1)
    U = np.array([0.5 * k * (c - z0) ** 2 for z0, k, _ in windows])
    n = H.sum(axis=0)
    use = n > 0
    Uu, nu = U[:, use], n[use]
    logN = np.log(N)

    def logden(f):                         # log sum_i N_i exp(f_i - U_ij), per bin
        a = logN[:, None] + f[:, None] - Uu
        m = a.max(axis=0)
        return m + np.log(np.exp(a - m).sum(axis=0))

    def obj(g):                            # g = f_1..f_{K-1}, f_0 = 0
        f = np.concatenate([[0.0], g])
        ld = logden(f)
        val = (nu * ld).sum() - (N * f).sum()
        w = np.exp(logN[:, None] + f[:, None] - Uu - ld[None, :])     # N_i e^{f_i-U_ij} / den_j
        grad = (w * nu[None, :]).sum(axis=1) - N
        return val, grad[1:]

    # Newton's method on the convex objective (exact Hessian, K x K): converges to machine
    # precision in ~10-30 steps even for 100+ windows, no scipy needed.
    K = len(windows)
    g = np.zeros(K - 1)
    val, grad = obj(g)
    for it in range(200):
        f = np.concatenate([[0.0], g])
        ld = logden(f)
        w = np.exp(logN[:, None] + f[:, None] - Uu - ld[None, :])       # K x bins
        Hs = np.diag((w * nu).sum(axis=1)) - (w * nu) @ w.T              # Hessian in f
        Hs = Hs[1:, 1:] + 1e-12 * np.eye(K - 1)
        step = np.linalg.solve(Hs, grad)
        t = 1.0
        while True:                                                      # backtracking line search
            gn = g - t * step
            vn, grn = obj(gn)
            if vn <= val - 1e-4 * t * grad @ step or t < 1e-8:
                break
            t *= 0.5
        g, val, grad = gn, vn, grn
        if np.abs(grad).max() < 1e-7 * max(1.0, N.max()) * 1e-3:
            break
    else:
        print(f"  WARNING: WHAM Newton did not converge (gradient {np.abs(grad).max():.1e})")
    f = np.concatenate([[0.0], g])
    out = np.full(len(c), np.nan)
    out[use] = np.log(nu) - logden(f)
    out[use] = -out[use]
    return (c, out, f) if return_f else (c, out)


def pmf_rz(windows, rlist, f, zbins, ref, dr=0.25):
    """2D free energy W(r, z) from the AXIAL windows (bias acts on z only, r is sampled freely):
    P(r,z) = sum_i n_i(r,z) / sum_i N_i exp(f_i - U_i(z));  W = -ln[P / r]  (r = area Jacobian),
    set to 0 in the upper bulk (z >= zmax - ref)."""
    rmax = max(r.max() for r in rlist)
    if rmax < 2 * dr:
        return np.array([]), np.array([]), np.zeros((0, 0))
    rbins = np.arange(0.0, rmax + dr, dr)
    zc = 0.5 * (zbins[1:] + zbins[:-1]); rc = 0.5 * (rbins[1:] + rbins[:-1])
    n = np.zeros((len(rc), len(zc)))
    N = np.array([len(z) for _, _, z in windows], float)
    for (z0, k, z), r in zip(windows, rlist):
        n += np.histogram2d(r, z, bins=[rbins, zbins])[0]
    U = np.array([0.5 * k * (zc - z0) ** 2 for z0, k, _ in windows])
    a = np.log(N)[:, None] + f[:, None] - U
    m = a.max(axis=0); logden = m + np.log(np.exp(a - m).sum(axis=0))
    W = np.full(n.shape, np.nan)
    ok = n > 0
    W[ok] = -(np.log(n[ok]) - np.broadcast_to(logden, n.shape)[ok] - np.log(np.broadcast_to(rc[:, None], n.shape)[ok]))
    bulk = (zc >= zc.max() - ref)
    W -= np.nanmean(W[:, bulk][rc < 8.0])
    return rc, zc, W


def overlap_report(windows, bins):
    ws = sorted(windows, key=lambda w: w[0])
    bad = []
    for (z0a, _, za), (z0b, _, zb) in zip(ws, ws[1:]):
        ha = np.histogram(za, bins)[0] / max(len(za), 1)
        hb = np.histogram(zb, bins)[0] / max(len(zb), 1)
        ov = np.minimum(ha, hb).sum()
        if ov < 0.03:
            bad.append((z0a, z0b, ov))
    return bad


def analyse(umb_dir, bw, ref, zmax_user, skip=0.0):
    reps = sorted(glob.glob(os.path.join(umb_dir, "rep*")))
    per_rep, pooled, halves = [], [], ([], [])
    rs = []                                              # r samples, same order as the pooled windows
    for r in reps:
        wins = []
        for p in sorted(glob.glob(os.path.join(r, "z*", "colvar.dat"))):
            z0, k, z = read_window(p)
            if len(z) > 10:
                cut = int(skip * len(z))                 # drop the first part of production (slow chain relaxation)
                wins.append((z0, k, z[cut:]))
                rs.append(read_window.last_r[cut:])
        if wins:
            per_rep.append((os.path.basename(r), wins))
    if not per_rep:
        raise SystemExit(f"{umb_dir}: no colvar.dat found (runs not finished?)")
    allz = np.concatenate([z for _, w in per_rep for _, _, z in w])
    z0min = min(z0 for _, w in per_rep for z0, _, _ in w)
    full = z0min < -ref                      # windows on both sides -> no mirroring
    zhi = allz.max() if zmax_user is None else zmax_user
    zlo = allz.min() if full else 0.0
    # edges offset by an odd amount so they never coincide with the 1e-4 precision of colvar.dat
    # (samples sitting exactly on edges shift every histogram the same way -> drift over many windows)
    bins = np.arange(zlo - 3.7e-5, zhi + bw, bw)
    print(f"  {'FULL pore (no mirroring)' if full else 'HALF pore (mirrored about z = 0)'}")

    def ref_shift(c, pmf):
        m = (c >= c.max() - ref) & np.isfinite(pmf)          # bulk above the pore
        return pmf - np.nanmean(pmf[m])

    curves = []
    for name, wins in per_rep:
        c, p = wham(wins, bins)
        curves.append(ref_shift(c, p))
        bad = overlap_report(wins, bins)
        print(f"  {name}: {len(wins)} windows, {sum(len(z) for _, _, z in wins)} samples"
              + (f"; poor overlap between windows {[(a, b) for a, b, _ in bad]}" if bad else "; overlap OK"))
        for z0, k, z in wins:
            pooled.append((z0, k, z))
            h = len(z) // 2
            halves[0].append((z0, k, z[:h])); halves[1].append((z0, k, z[h:]))
    curves = np.array(curves)
    c, pp, fpool = wham(pooled, bins, return_f=True); pp = ref_shift(c, pp)
    map2d = pmf_rz(pooled, rs, fpool, bins, ref)
    _, h1 = wham(halves[0], bins); h1 = ref_shift(c, h1)
    _, h2 = wham(halves[1], bins); h2 = ref_shift(c, h2)
    z0max = max(z0 for z0, _, _ in pooled)
    edge = (c > z0max + 0.5) | (c < z0min - 0.5)          # beyond the outermost windows: unreliable
    for arr in (pp, h1, h2):
        arr[edge] = np.nan
    curves[:, edge] = np.nan
    mean = np.nanmean(curves, axis=0)
    sd = np.nanstd(curves, axis=0, ddof=1) if len(curves) > 1 else np.full(len(c), np.nan)
    if full:
        lo = (c <= c.min() + ref) & np.isfinite(pp)
        print(f"  check: PMF in the lower bulk = {np.nanmean(pp[lo]):+.2f} kT "
              f"(should be ~0: both bulks are the same state; a large value = not converged)")
        return dict(z=c, mean=mean, sd=sd, pooled=pp, h1=h1, h2=h2, n=len(curves), map2d=map2d, full=True)
    keep = c >= 0                                        # half pore, then mirror
    c, mean, sd, pp, h1, h2 = (x[keep] for x in (c, mean, sd, pp, h1, h2))
    mir = lambda x: np.concatenate([x[::-1], x])
    zfull = np.concatenate([-c[::-1], c])
    return dict(z=zfull, mean=mir(mean), sd=mir(sd), pooled=mir(pp), h1=mir(h1), h2=mir(h2),
                n=len(curves), map2d=map2d, full=False)


def label_of(umb_dir):
    d = os.path.normpath(umb_dir)
    sysname = os.path.basename(os.path.dirname(d)) or "sys"
    s = re.search(r"umbrella_s([\d.]+)", d)
    return sysname, (s.group(1) if s else "")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dirs", nargs="+", help="<SYS>/umbrella_s<sigma> folders")
    ap.add_argument("--bin", type=float, default=0.2)
    ap.add_argument("--ref", type=float, default=5.0, help="bulk reference = last REF sigma of the z range")
    ap.add_argument("--zmax", type=float, default=None)
    ap.add_argument("--skip", type=float, default=0.0,
                    help="fraction of each window's production to discard first (e.g. 0.3); use if half1 and half2 differ")
    a = ap.parse_args()
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    res = {}
    for d in a.dirs:
        sysname, s = label_of(d)
        print(f"{sysname}, sigma {s}:")
        R = analyse(d, a.bin, a.ref, a.zmax, a.skip)
        res[(sysname, s)] = R
        out = f"pmf_{sysname}_s{s}"
        np.savetxt(out + ".csv", np.column_stack([R["z"], R["mean"], R["sd"], R["pooled"], R["h1"], R["h2"],
                                                  np.full(len(R["z"]), R["n"])]),
                   delimiter=",", header="z,pmf_mean,pmf_sd,pmf_pooled,pmf_half1,pmf_half2,n_reps", comments="",
                   fmt="%.5g")
        mid = np.nanmax(R["pooled"][np.abs(R["z"]) <= 0.5]) if np.any(np.abs(R["z"]) <= 0.5) else np.nan
        print(f"  PMF at midplane {mid:.2f} kT, max {np.nanmax(R['pooled']):.2f} kT, "
              f"min {np.nanmin(R['pooled']):.2f} kT;  half1 vs half2 max diff "
              f"{np.nanmax(np.abs(R['h1'] - R['h2'])):.2f} kT  -> {out}.csv")
        rc, zc, W = R["map2d"]
        if W.size == 0 or not np.isfinite(W).any():
            print("  (no radial data for a 2D map)")
        else:
          if not R["full"]:                                   # half pore: mirror the map too
              keep = zc >= 0
              zc = np.concatenate([-zc[keep][::-1], zc[keep]]); W = np.concatenate([W[:, keep][:, ::-1], W[:, keep]], axis=1)
          RR, ZZ = np.meshgrid(rc, zc, indexing="ij")
          np.savetxt(out + "_rz.csv", np.column_stack([RR.ravel(), ZZ.ravel(), W.ravel()]), delimiter=",",
                     header="r,z,W_kT", comments="", fmt="%.4g")
          fig2, ax2 = plt.subplots(figsize=(4.2, 6))
          vmax = np.nanpercentile(W, 98) if np.isfinite(W).any() else 1
          im = ax2.pcolormesh(rc, zc, W.T, shading="auto", cmap="viridis", vmin=np.nanmin(W), vmax=vmax)
          ax2.axhline(12.5, color="w", lw=0.6, ls="--"); ax2.axhline(-12.5, color="w", lw=0.6, ls="--")
          ax2.set_xlabel("cargo r (σ)"); ax2.set_ylabel("cargo z (σ)"); ax2.set_title(f"W(r, z): {sysname}, σ {s}", fontsize=9)
          fig2.colorbar(im, ax=ax2, label="kT"); fig2.tight_layout(); fig2.savefig(out + "_rz.png", dpi=200); plt.close(fig2)
        fig, ax = plt.subplots(figsize=(6, 3.8))
        ax.axvspan(-12.5, 12.5, color="#eee", zorder=0, label="pore")
        if R["n"] > 1:
            ax.fill_between(R["z"], R["mean"] - R["sd"], R["mean"] + R["sd"], color="#9ecae1", label=f"mean ± SD ({R['n']} reps)")
        ax.plot(R["z"], R["pooled"], color="#08519c", lw=1.8, label="all reps pooled")
        ax.plot(R["z"], R["h1"], color="#999", lw=0.8, ls="--", label="1st / 2nd half")
        ax.plot(R["z"], R["h2"], color="#999", lw=0.8, ls="--")
        ax.axhline(0, color="k", lw=0.5)
        ax.set_xlabel("cargo z (σ)"); ax.set_ylabel("PMF (kT)"); ax.set_title(f"{sysname}, cargo σ = {s}")
        ax.legend(fontsize=7, frameon=False); fig.tight_layout(); fig.savefig(out + ".png", dpi=200); plt.close(fig)
    if len(res) > 1:
        fig, ax = plt.subplots(figsize=(6.5, 4))
        ax.axvspan(-12.5, 12.5, color="#eee", zorder=0)
        for (sysname, s), R in res.items():
            l, = ax.plot(R["z"], R["pooled"], lw=1.6, label=f"{sysname} σ{s}")
            if R["n"] > 1:
                ax.fill_between(R["z"], R["mean"] - R["sd"], R["mean"] + R["sd"], color=l.get_color(), alpha=0.2)
        ax.axhline(0, color="k", lw=0.5); ax.set_xlabel("cargo z (σ)"); ax.set_ylabel("PMF (kT)")
        ax.legend(fontsize=7, frameon=False); fig.tight_layout(); fig.savefig("pmf_compare.png", dpi=200)
        print("-> pmf_compare.png")


if __name__ == "__main__":
    main()
