"""
Step 2: size cutoff sigma_50 of every brush (main text eqn 6, ESI Section S1.1).
  P_pass(sigma_c) = P0 / (1 + exp[(sigma_c - sigma_50)/w]),   w >= 0.15 sigma
  - binomial maximum likelihood on (passes, total) at each cargo size
  - 95% interval of sigma_50: parametric bootstrap, 300 resampled data sets
  - goodness of fit p: deviance of the data vs 400 bootstrap deviances
Input : counts.csv (from counts_from_runs.py). Denominator = n_total (100 per size).
Output: fits.csv
Usage : python3 fit_sigma50.py [n_cores]
"""
import sys, numpy as np, pandas as pd
from scipy.optimize import minimize
from multiprocessing import Pool

def model(p, x):
    P0 = 1 / (1 + np.exp(-p[0])); w = max(np.exp(p[2]), 0.15)    # P0 kept in (0,1), w >= 0.15
    return P0 / (1 + np.exp((x - p[1]) / w))

def nll(p, x, n, k):
    P = np.clip(model(p, x), 1e-12, 1 - 1e-12)
    return -(k * np.log(P) + (n - k) * np.log(1 - P)).sum()

def fit(x, n, k, starts=None):
    starts = starts or [[a, s, np.log(0.5)] for s in (0.5, 1, 2, 3, 4, 5, 6, 7, 8) for a in (-2, -1, 0)]
    best = None
    for s0 in starts:
        r = minimize(nll, s0, args=(x, n, k), method="Nelder-Mead", options=dict(maxiter=4000, xatol=1e-6, fatol=1e-9))
        if best is None or r.fun < best.fun: best = r
    return best.x

def deviance(n, k, P):
    P = np.clip(P, 1e-12, 1 - 1e-12)
    a = np.where(k > 0, k * np.log(np.where(k > 0, k, 1) / (n * P)), 0)
    b = np.where(n - k > 0, (n - k) * np.log(np.where(n - k > 0, n - k, 1) / (n * (1 - P))), 0)
    return 2 * (a + b).sum()

def job(arg):
    (fa, fr), q = arg
    q = q.sort_values("sigma"); x = q.sigma.values; n = q.n_total.values; k = q.n_pass.values
    p = fit(x, n, k); P = model(p, x); rng = np.random.default_rng(1)
    near = lambda p: [list(p), [p[0], p[1] - 1, p[2]], [p[0], p[1] + 1, p[2]], [-1, 3, np.log(0.5)], [-1, 6, np.log(0.5)]]
    boot = [fit(x, n, rng.binomial(n, P), near(p))[1] for _ in range(300)]
    D0 = deviance(n, k, P); worse = 0
    for _ in range(400):
        kk = rng.binomial(n, P); worse += deviance(n, kk, model(fit(x, n, kk, near(p)), x)) >= D0 - 1e-9
    lo, hi = np.percentile(boot, [2.5, 97.5])
    return dict(f_A=fa, f_rigid=fr, P0=1 / (1 + np.exp(-p[0])), sigma50=p[1], w=max(np.exp(p[2]), 0.15),
                ci_low=lo, ci_high=hi, gof_p=worse / 400)

if __name__ == "__main__":
    c = pd.read_csv("counts.csv", dtype={"f_A": str})
    groups = [g for g in c.groupby(["f_A", "f_rigid"]) if g[1].n_pass.sum() >= 10]   # f_A = 0.8: too few passes (12 of 4500)
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 4) as pool:
        res = pd.DataFrame(pool.map(job, groups))
    res.to_csv("fits.csv", index=False); print(res.round(3).to_string(index=False))
