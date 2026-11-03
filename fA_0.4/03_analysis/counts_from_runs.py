"""
Step 1: count the outcome of every transport trajectory.
Input : the per-run tables written by dwell_hydro.py (one row per trajectory):
        system, sigma, run, outcome, ...   outcome = passed / trapped / escaped / outside
        e.g.  dwell_runs.csv (runs 1-50) + dwell_..._runs51_100_runs.csv (runs 51-100)
Output: counts.csv with, for every brush (f_A, f_rigid) and cargo size:
        n_total (100), n_pass, n_trap, n_kept (= total minus 'outside')
Usage : python3 counts_from_runs.py <fA label> <csv> [<csv> ...]  -> appends to counts.csv
        e.g. python3 counts_from_runs.py none dwell_runs.csv dwell_nh_runs51_100_runs.csv
"""
import sys, os, pandas as pd
fa, files = sys.argv[1], sys.argv[2:]
d = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
d["f_rigid"] = d.system.str.split("_").str[0].astype(int) / 100      # '50_nh' -> 0.5
d = d.drop_duplicates(["f_rigid", "sigma", "run"])
g = d.groupby(["f_rigid", "sigma"]).outcome
out = pd.DataFrame({"n_total": g.size(), "n_pass": g.apply(lambda s: (s == "passed").sum()),
                    "n_trap": g.apply(lambda s: (s == "trapped").sum()),
                    "n_kept": g.apply(lambda s: (s != "outside").sum())}).reset_index()
out.insert(0, "f_A", fa)
out.to_csv("counts.csv", mode="a", header=not os.path.exists("counts.csv"), index=False)
print(out.to_string(index=False))
