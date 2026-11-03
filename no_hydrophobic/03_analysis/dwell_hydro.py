#!/usr/bin/env python3
"""
Dwell times of the cargo in the pore, from the DCD of every finished run.
Works with python 3.6 + numpy (no MDAnalysis needed).

Layouts
  hydrophobe runs (setup_hydro.py):   <SYS>/sigma_<s>/run_<n>/dump_run<n>.dcd
      python3 dwell_hydro.py 0_sk5 25_sk5 50_sk5 75_sk5 100_sk5
      python3 dwell_hydro.py                       -> every <SYS>/sigma_* folder here
  no-hydrophobe runs (setup_runs.py):  <pct>_run_<n>  (sigma 0.5)  and  <pct>_<s>_run_<n>
      python3 dwell_hydro.py --flat 0 25 50 75 100 --name nh

Definitions (pore = |z| <= 12.5 and r <= 10; frame spacing = 100 steps x 0.005 = 0.5 tau)
  (the cargo starts AT the top mouth, z = 12.5, so "first entry" is t = 0 for every run)
  t_pass       release -> exit through the bottom (PASSED runs) = total translocation time
  transit      last time the cargo was above the top mouth (z >= 12.5) -> exit
               = dwell of the successful attempt only (failed back-and-forth attempts removed)
  t_inside     total time spent inside the pore (any outcome)
  reached_mid  1 if the cargo ever got below z = 0 inside the pore
Outputs
  dwell_runs.csv      one line per run
  dwell_summary.csv   per system x sigma: n_pass, median/mean t_pass and transit of passed runs,
                      mean t_inside of passed runs, fraction reaching midpore, mean t_inside of escaped runs
"""
import argparse, glob, os, re, struct, sys
import numpy as np

p = argparse.ArgumentParser()
p.add_argument("systems", nargs="*")
p.add_argument("--flat", action="store_true", help="no-hydrophobe layout <pct>_<s>_run_<n>; give pct values")
p.add_argument("--name", default="nh", help="system suffix for --flat (e.g. 50 -> 50_nh)")
p.add_argument("--dt", type=float, default=0.005)
p.add_argument("--every", type=int, default=100, help="DCD every N steps")
p.add_argument("--zpore", type=float, default=12.5)
p.add_argument("--rpore", type=float, default=10.0)
p.add_argument("--out", default="dwell")
a = p.parse_args()
DTF = a.dt * a.every


# ---------------- minimal LAMMPS DCD reader: returns cargo (x,y,z) per frame ----------------
def read_dcd_atom(path, idx):
    with open(path, "rb") as f:
        buf = f.read()
    if len(buf) < 100:
        return None
    e = "<"
    if struct.unpack("<i", buf[:4])[0] != 84:
        e = ">"
    icntrl = struct.unpack(e + "20i", buf[8:88])
    has_cell = icntrl[10] != 0
    pos = 92
    n = struct.unpack(e + "i", buf[pos:pos + 4])[0]; pos += 4 + n + 4          # title block
    pos += 4; natoms = struct.unpack(e + "i", buf[pos:pos + 4])[0]; pos += 8   # natoms block
    blk = 4 * natoms + 8
    fsize = (56 if has_cell else 0) + 3 * blk
    nfr = (len(buf) - pos) // fsize
    out = np.empty((nfr, 3))
    for k in range(nfr):
        b = pos + k * fsize + (56 if has_cell else 0)
        for c in range(3):
            o = b + c * blk + 4 + 4 * idx
            out[k, c] = struct.unpack(e + "f", buf[o:o + 4])[0]
    return out


# ---------------- find the cargo atom index from the run's data file ----------------
_cargo_cache = {}
def cargo_index(rundir):
    inp = os.path.join(rundir, "in.lammps")
    df = None
    for line in open(inp):
        m = re.match(r"\s*read_data\s+(\S+)", line)
        if m: df = m.group(1); break
    if df is None: raise RuntimeError("no read_data in " + inp)
    df = os.path.normpath(os.path.join(rundir, df))
    if df in _cargo_cache: return _cargo_cache[df]
    lines = open(df).read().split("\n")
    i = next(k for k, l in enumerate(lines) if re.match(r"\s*Atoms", l))
    col = 1 if "atomic" in lines[i] else 2           # full / molecular / bond: id mol type ...
    ids, types = [], []
    for l in lines[i + 1:]:
        s = l.split("#")[0].split()
        if not s:
            if ids: break
            continue
        if not re.match(r"^-?\d+$", s[0]): break
        ids.append(int(s[0])); types.append(int(s[col]))
    ids = np.array(ids); types = np.array(types)
    ut, cnt = np.unique(types, return_counts=True)
    single = ut[cnt == 1]
    if len(single) != 1: raise RuntimeError("cannot identify the single cargo atom in " + df)
    cid = ids[types == single[0]][0]
    idx = int(np.searchsorted(np.sort(ids), cid))       # DCD is sorted by atom ID
    _cargo_cache[df] = idx
    return idx


def outcome(log):
    t = open(log).read()
    if "completed" not in t: return None
    if "cargo PASSED" in t: return "passed"
    if "OUTSIDE" in t: return "outside"
    if "ESCAPED" in t: return "escaped"
    return "trapped"


# ---------------- collect run folders ----------------
jobs = []   # (system, sigma, run, dir)
if a.flat:
    for pct in a.systems or ["0", "25", "50", "75", "100"]:
        for d in sorted(glob.glob("%s_run_*" % pct)):
            jobs.append(("%s_%s" % (pct, a.name), 0.5, int(d.split("_")[-1]), d))
        for d in sorted(glob.glob("%s_*_run_*" % pct)):
            m = re.match(r"%s_([0-9.]+)_run_(\d+)$" % pct, os.path.basename(d.rstrip("/")))
            if m: jobs.append(("%s_%s" % (pct, a.name), float(m.group(1)), int(m.group(2)), d))
else:
    systems = a.systems or sorted(set(os.path.dirname(x.rstrip("/"))
                                      for x in glob.glob("*/sigma_*/")))
    for S in systems:
        for d in glob.glob(os.path.join(S, "sigma_*", "run_*")):
            m = re.search(r"sigma_([0-9.]+)/run_(\d+)$", d)
            if m: jobs.append((S.rstrip("/"), float(m.group(1)), int(m.group(2)), d))
jobs.sort(key=lambda j: (j[0], j[1], j[2]))
if not jobs: sys.exit("no run folders found")

rows = []
for n, (S, s, r, d) in enumerate(jobs):
    log = os.path.join(d, "stdout.log")
    if not os.path.exists(log): continue
    oc = outcome(log)
    if oc is None: continue
    dcds = glob.glob(os.path.join(d, "*.dcd"))
    if not dcds: rows.append((S, s, r, oc, 0, np.nan, np.nan, np.nan, "no dcd")); continue
    try:
        xyz = read_dcd_atom(dcds[0], cargo_index(d))
    except Exception as ex:
        rows.append((S, s, r, oc, 0, np.nan, np.nan, np.nan, str(ex)[:60])); continue
    if xyz is None or len(xyz) == 0:
        rows.append((S, s, r, oc, 0, np.nan, np.nan, np.nan, "empty dcd")); continue
    z = xyz[:, 2]; rr = np.hypot(xyz[:, 0], xyz[:, 1])
    inside = (np.abs(z) < a.zpore) & (rr <= a.rpore)
    t_end = (len(z) - 1) * DTF
    t_in = inside.sum() * DTF
    mid = int(((z < 0) & inside).any())
    if oc == "passed":
        top = np.where(z[:-1] >= a.zpore)[0]              # last time above the top mouth before exit
        t_top = (top[-1] if len(top) else 0) * DTF
        rows.append((S, s, r, oc, mid, t_end, t_end - t_top, t_in, ""))
    else:
        rows.append((S, s, r, oc, mid, t_end, np.nan, t_in, ""))
    if (n + 1) % 200 == 0: print("  %d / %d runs read" % (n + 1, len(jobs)))

with open(a.out + "_runs.csv", "w") as f:
    f.write("system,sigma,run,outcome,reached_mid,t_end,transit,t_inside,note\n")
    for x in rows:
        f.write("%s,%g,%d,%s,%d,%s,%s,%s,%s\n" % (x[0], x[1], x[2], x[3], x[4],
                *["" if not np.isfinite(v) else "%.1f" % v for v in x[5:8]], x[8]))

def st(v):
    v = np.array([u for u in v if np.isfinite(u)])
    if not len(v): return (0, np.nan, np.nan, np.nan)
    return (len(v), np.median(v), v.mean(), v.std(ddof=1) if len(v) > 1 else np.nan)

hdr = ("system", "sigma", "N", "pass", "tpass_med", "tpass_mean", "transit_med", "transit_mean", "tin_pass", "mid_frac", "tin_esc")
print("\n%-9s %5s %4s %4s | %9s %10s | %11s %12s | %8s | %8s %7s" % hdr)
with open(a.out + "_summary.csv", "w") as f:
    f.write("system,sigma,N,n_pass,tpass_median,tpass_mean,tpass_sd,transit_median,transit_mean,transit_sd,"
            "tinside_pass_mean,frac_reached_mid,tinside_escaped_mean\n")
    pk = lambda k: (int(re.match(r"\d+", k[0]).group()) if re.match(r"\d+", k[0]) else 0, k[0], k[1])
    for S, s in sorted(set((x[0], x[1]) for x in rows), key=pk):
        g = [x for x in rows if x[0] == S and x[1] == s and x[3] != "outside"]
        P = [x for x in g if x[3] == "passed"]
        n, tpm, tpa, tps = st([x[5] for x in P])
        _, trm, tra, trs = st([x[6] for x in P])
        tip = st([x[7] for x in P])[2]
        mid = np.mean([x[4] for x in g]) if g else np.nan
        tie = st([x[7] for x in g if x[3] == "escaped"])[2]
        f.write("%s,%g,%d,%d,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.2f,%.1f\n" % (S, s, len(g), n, tpm, tpa, tps, trm, tra, trs, tip, mid, tie))
        print("%-9s %5g %4d %4d | %9.1f %10.1f | %11.1f %12.1f | %8.1f | %8.2f %7.1f" % (S, s, len(g), n, tpm, tpa, trm, tra, tip, mid, tie))
print("\nsaved %s_runs.csv and %s_summary.csv   (times in tau; OUTSIDE runs left out of the summary)" % (a.out, a.out))
