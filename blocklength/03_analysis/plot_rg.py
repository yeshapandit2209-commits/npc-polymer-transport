"""
Radius of gyration of the single-chain simulations (N = 64).
Input : LAMMPS 'fix ave/time' output of compute gyration, one file per run:
        rg_bl<lambda>.dat  (charged chain, charge blocklength lambda, k_theta = 0)
        rg_k<k>.dat        (uncharged chain, bending constant k_theta)
        Format: '# ...' header lines, then 'TimeStep  Rg'.
Output: Rg_vs_time_lambda.pdf  (main text Fig. 1c)
        Rg_vs_ktheta.pdf       (ESI Fig. S2)
        rg_summary.csv         (Rg/sqrt(N), mean and s.d. over the second half)
Usage : python3 plot_rg.py <folder with the .dat files>
"""
import sys, glob, re, os
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

N = 64; SQN = np.sqrt(N); TAU_PER_STEP = 0.01          # time step 0.01 tau
folder = sys.argv[1] if len(sys.argv) > 1 else "."

def read(f):
    # keep only 'step Rg' lines (a restarted run may leave lines with a different number of columns)
    rows = [l.split() for l in open(f) if not l.startswith("#")]
    d = np.array([r for r in rows if len(r) == 2], float); return d[:, 0], d[:, 1]

def summary(rg):
    h = rg[len(rg) // 2:] / SQN                          # second half = equilibrated part
    return h.mean(), h.std(ddof=1)

rows = []
# ---------- Fig. 1c: Rg/sqrt(N) over time for each blocklength
BL = sorted(glob.glob(os.path.join(folder, "*rg_bl*.dat")), key=lambda f: int(re.search(r"bl(\d+)", f).group(1)))
fig, ax = plt.subplots(figsize=(3.3, 2.4))
for f in BL:
    lam = int(re.search(r"bl(\d+)", f).group(1)); t, rg = read(f); m, sd = summary(rg)
    line, = ax.plot(t * TAU_PER_STEP, rg / SQN, lw=0.8, label=rf"$\lambda$ = {lam}")
    ax.axhline(m, color=line.get_color(), ls="--", lw=0.8)
    rows.append(("lambda", lam, m, sd))
ax.set_xlabel(r"Time ($\tau$)"); ax.set_ylabel(r"$R_g/\sqrt{N}$"); ax.legend(frameon=False, fontsize=6)
fig.tight_layout(); fig.savefig("Rg_vs_time_lambda.pdf")

# ---------- ESI Fig. S2: Rg/sqrt(N) vs k_theta (uncharged chain)
K = sorted(glob.glob(os.path.join(folder, "*rg_k*.dat")), key=lambda f: int(re.search(r"rg_k(\d+)", f).group(1)))
ks, ms, sds = [], [], []
for f in K:
    k = int(re.search(r"rg_k(\d+)", f).group(1)); _, rg = read(f); m, sd = summary(rg)
    ks.append(k); ms.append(m); sds.append(sd); rows.append(("k_theta", k, m, sd))
fig, ax = plt.subplots(figsize=(3.3, 2.4))
ax.errorbar(ks, ms, yerr=sds, marker="o", ms=4, lw=1.2, capsize=2, color="k")
ax.set_xlabel(r"Bending constant, $k_\theta$"); ax.set_ylabel(r"$R_g/\sqrt{N}$")
fig.tight_layout(); fig.savefig("Rg_vs_ktheta.pdf")

pd.DataFrame(rows, columns=["series", "value", "Rg_over_sqrtN_mean", "Rg_over_sqrtN_sd"]).to_csv("rg_summary.csv", index=False)
print(pd.DataFrame(rows, columns=["series", "value", "mean", "sd"]).round(3).to_string(index=False))
