#!/usr/bin/env python3
"""
Position-dependent diffusion D(z) of the cargo from umbrella windows (Hummer, New J. Phys. 2005).
For each window:  D = var(z)^2 / integral_0^inf C(t) dt,   C(t) = <dz(0) dz(t)>,  dz = z - <z>
This holds for underdamped (Langevin) motion too, as long as the integral is taken over the whole decay.
Reference: free cargo with 'fix langevin ... damp' and mass m has D0 = kT * damp / m  (= 1 for damp 1, m 1).

  python3 diffusion_z.py 0k5/umbrella_s4.0 25k5/umbrella_s4.0 50k5/umbrella_s4.0 75k5/umbrella_s4.0 100k5/umbrella_s4.0
  options: --dt 0.05 (time between colvar lines = print every 10 x 0.005) --skip 0.1 --tmax 50

Reads <dir>/rep*/z*/colvar.dat  (step x y z; header holds z0).
Writes  diffz_<SYS>.csv  (z0, mean z, D per replicate, D mean, D sd)  and  diffz.png
"""
import argparse, glob, os, re
import numpy as np

p = argparse.ArgumentParser()
p.add_argument("dirs", nargs="+")
p.add_argument("--dt", type=float, default=0.05)
p.add_argument("--skip", type=float, default=0.1)
p.add_argument("--tmax", type=float, default=50.0, help="longest lag used (tau)")
p.add_argument("--D0", type=float, default=1.0, help="bulk D for the reference line")
p.add_argument("--out", default="diffz")
a = p.parse_args()


def acf(x):
    x = x - x.mean(); n = len(x); f = np.fft.rfft(x, 2 * n)
    c = np.fft.irfft(f * np.conj(f))[:n] / np.arange(n, 0, -1)
    return c


def D_window(z, dt, tmax):
    var = z.var()
    c = acf(z)
    L = min(len(c) - 1, int(tmax / dt))
    c = c[:L]
    # integrate until the decay is complete: first lag after which |C| stays below 2% of C(0) for 20 lags
    small = np.abs(c) < 0.02 * c[0]
    stop = L
    for i in range(1, L - 20):
        if small[i:i + 20].all():
            stop = i; break
    I = np.trapezoid(c[:stop + 1], dx=dt) if hasattr(np, "trapezoid") else np.trapz(c[:stop + 1], dx=dt)
    return var ** 2 / I if I > 0 else np.nan, var, stop * dt


results = {}
for d in a.dirs:
    sysname = d.rstrip("/").split("/")[-2] if "/" in d.rstrip("/") else d
    rows = {}
    for f in sorted(glob.glob(os.path.join(d, "rep*", "z*", "colvar.dat"))):
        head = open(f).readline()
        m = re.search(r"z0\s*=\s*(-?[\d.]+)", head)
        if not m: continue
        z0 = float(m.group(1)); rep = f.split(os.sep)[-3]
        try:
            z = np.loadtxt(f, comments="#", usecols=3)
        except Exception:
            continue
        z = z[int(a.skip * len(z)):]
        if len(z) < 2000: continue
        D, var, tc = D_window(z, a.dt, a.tmax)
        rows.setdefault(z0, {})[rep] = (D, z.mean(), var, tc)
    if not rows:
        print("%s: no colvar.dat found" % d); continue
    out = []
    reps = sorted({r for v in rows.values() for r in v})
    with open("%s_%s.csv" % (a.out, sysname), "w") as fh:
        fh.write("z0,z_mean," + ",".join("D_" + r for r in reps) + ",D_mean,D_sd,var_mean,tcorr_mean\n")
        for z0 in sorted(rows):
            v = rows[z0]; Ds = [v[r][0] if r in v else np.nan for r in reps]
            zm = np.nanmean([v[r][1] for r in v]); Dm = np.nanmean(Ds); Dsd = np.nanstd(Ds)
            vm = np.nanmean([v[r][2] for r in v]); tm = np.nanmean([v[r][3] for r in v])
            fh.write("%g,%.3f,%s,%.4f,%.4f,%.4f,%.3f\n" % (z0, zm, ",".join("%.4f" % x for x in Ds), Dm, Dsd, vm, tm))
            out.append((z0, Dm, Dsd))
    out = np.array(out); results[sysname] = out
    inside = out[np.abs(out[:, 0]) <= 12.5, 1]; bulk = out[np.abs(out[:, 0]) >= 40, 1]
    print("%-8s windows %3d | D bulk %.3f | D in pore %.3f (min %.3f) | ratio %.2f"
          % (sysname, len(out), np.nanmean(bulk), np.nanmean(inside), np.nanmin(inside), np.nanmean(inside) / np.nanmean(bulk)))

if results:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    C = ["#86b6ef", "#3987e5", "#256abf", "#184f95", "#0d366b", "#555555"]
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    ax.axvspan(-12.5, 12.5, color="#f0efec", zorder=0)
    for i, (s, o) in enumerate(results.items()):
        ax.plot(o[:, 0], o[:, 1], "-", color=C[i % len(C)], lw=2, label=s)
        ax.fill_between(o[:, 0], o[:, 1] - o[:, 2], o[:, 1] + o[:, 2], color=C[i % len(C)], alpha=0.15, lw=0)
    ax.axhline(a.D0, color="#6b6b6b", ls="--", lw=1); ax.text(ax.get_xlim()[1], a.D0, " free", va="center", color="#6b6b6b", fontsize=8)
    ax.set_xlabel(r"$z$ ($\sigma$)"); ax.set_ylabel(r"$D(z)$  ($\sigma^2/\tau$)"); ax.set_ylim(bottom=0)
    ax.legend(fontsize=9); ax.grid(axis="y", color="#e4e4e4")
    fig.tight_layout(); fig.savefig(a.out + ".png", dpi=250)
    print("wrote %s.png and %s_<SYS>.csv" % (a.out, a.out))
