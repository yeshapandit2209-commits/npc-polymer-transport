#!/usr/bin/env python3
"""
Equilibrate brush data files (no cargo) with the k = 5 angle rule and A-A epsilon 1.5.
Python 3.6+, no numpy.

    python3 setup_equil.py --groups s,db  *.data        (cluster 1)
    python3 setup_equil.py --groups tp    *.data        (cluster 2: one group gets all 48 cores)

  * group of a file = letters after "<pct>_" in its name:  50_s20.data -> s,  25_db20.data -> db.
    25_cs20.data is treated as s (it is an all-singlet file); change with --alias.
  * one folder per system:  equil/<SYS>/   (SYS = file name without .data, e.g. 50_s20)
  * one SLURM job per cluster runs both groups at once, each on its own share of the cores:
        24 cores for group 1, 24 for group 2  (--cores 48)
        cores per run = share // number of systems in the group  (5 systems -> 4 each, 4 -> 6 each)
  * resumable: restart files every --restart-every steps; resubmit the same job and unfinished
    runs continue with "run N upto"; finished runs are skipped.
  * at the end each run writes <SYS>_eq.data (Masses keep their "# name" comments) for setup_hydro.py.

Force field (from lammps_*.config, with the requested changes)
    all pairs WCA (eps 1, sigma 1, rc 2^(1/6)), except A-A (name 5): eps --eps-aa (1.5), rc 2.5
    bond harmonic 100 1.0, special_bonds lj 0 1 1, angle harmonic theta0 180,
    K by sorted name triple (k = 5 system; A = name 5):
      1-2-2 5 | 2-2-2 5 | 2-2-3 3 | 2-3-3 2 | 1-2-5 3 | 2-2-5 3 | 2-5-5 2 | 2-3-5 2 | everything else 0
    anchors (name 1) and shell (name 4) frozen; polymer NVT at T = 1 (Tdamp 1.0), dt 0.005
"""
import argparse, os, re, sys, glob
from collections import defaultdict, OrderedDict

K_TABLE = {(1, 2, 2): 5.0, (2, 2, 2): 5.0, (2, 2, 3): 3.0, (2, 3, 3): 2.0,
           (1, 2, 5): 3.0, (2, 2, 5): 3.0, (2, 5, 5): 2.0, (2, 3, 5): 2.0}

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("datafiles", nargs="+")
p.add_argument("--groups", required=True, help="one or two groups for this cluster, comma separated: s,db   or   tp")
p.add_argument("--alias", nargs="*", default=["cs=s"], help="group renames, default cs=s")
p.add_argument("--eps-aa", type=float, default=1.5)
p.add_argument("--steps", type=int, default=100000000, help="equilibration steps (config: 1e8)")
p.add_argument("--dump-every", type=int, default=50000)
p.add_argument("--thermo", type=int, default=10000)
p.add_argument("--restart-every", type=int, default=1000000)
p.add_argument("--seed", type=int, default=12345)
p.add_argument("--cores", type=int, default=48)
p.add_argument("--np", type=int, nargs="+", default=None, help="override cores per run, one value per group")
p.add_argument("--outdir", default="equil")
p.add_argument("--partition", default="standard")
p.add_argument("--time", default="7-00:00:00")
p.add_argument("--modules", nargs="*", default=["gcc/12.4", "openmpi/4.1.8"])
p.add_argument("--purge", action="store_true", help="run module purge before loading (off by default: keeps the default MPI on PATH)")
p.add_argument("--lmp", default="", help="LAMMPS executable (else auto-detect in the job)")
a = p.parse_args()
a.groups = [g for g in re.split(r"[,\s]+", a.groups) if g]
alias = dict(x.split("=") for x in a.alias)


# ------------------------------------------------------------------ data handling
def sections(L):
    idx = {}
    for k, l in enumerate(L):
        s = l.split("#")[0].strip()
        if s in ("Masses", "Atoms", "Bonds", "Angles", "Velocities", "Dihedrals", "Impropers",
                 "Pair Coeffs", "Bond Coeffs", "Angle Coeffs"):
            idx[s] = k
    return idx


def block(L, k):
    i = k + 1
    while not L[i].strip():
        i += 1
    j = i
    while j < len(L) and L[j].strip():
        j += 1
    return i, j


def prepare(src, dst):
    L = open(src).read().split("\n")
    idx = sections(L)
    tname = {}
    i, j = block(L, idx["Masses"])
    for l in L[i:j]:
        body, _, com = l.partition("#")
        t = int(body.split()[0])
        try:
            tname[t] = int(com.split()[0])
        except (ValueError, IndexError):
            sys.exit("%s: Masses line without '# name' comment: %r" % (src, l))
    ai, aj = block(L, idx["Atoms"])
    name, count = {}, defaultdict(int)
    for l in L[ai:aj]:
        q = l.split("#")[0].split()
        name[int(q[0])] = tname[int(q[2])]; count[int(q[2])] += 1
    names = set(tname.values())
    if not {1, 4} <= names:
        sys.exit("%s: need names 1 (anchor) and 4 (shell); found %s" % (src, sorted(names)))
    if any(c == 1 for c in count.values()):
        print("   note: %s has a single-atom type (cargo?) -> it is treated as polymer and will move" % src)
    # angles: one type per sorted name triple
    ni, nj = block(L, idx["Angles"])
    combos, rows = OrderedDict(), []
    for l in L[ni:nj]:
        q = l.split("#")[0].split()
        tri = [int(x) for x in q[2:5]]
        key = tuple(sorted(name[x] for x in tri))
        if key not in combos:
            combos[key] = len(combos) + 1
        rows.append("%s %d %d %d %d" % (q[0], combos[key], tri[0], tri[1], tri[2]))
    L[ni:nj] = rows
    L = [re.sub(r"^\s*\d+(\s+angle types)", " %d\\1" % len(combos), l) if l.strip().endswith("angle types") else l
         for l in L]
    out, skip = [], False                     # drop VMD '# Angle Coeffs' comment block (types changed)
    for l in L:
        if l.strip().startswith("# Angle Coeffs"):
            skip = True; continue
        if skip and l.strip().startswith("#"):
            continue
        skip = False; out.append(l)
    open(dst, "w").write("\n".join(out))
    return tname, combos


def types_of(tname, n):
    return [t for t, v in sorted(tname.items()) if v == n]


# ------------------------------------------------------------------ LAMMPS input
def ff_block(tname, combos):
    ts = sorted(tname)
    A = types_of(tname, 5)
    has_A = bool(A)
    s = ["pair_style      lj/cut %s" % ("2.5" if has_A else "1.122462048")]
    for i in ts:
        for j in ts:
            if j < i: continue
            if tname[i] == 5 and tname[j] == 5:
                s.append("pair_coeff %d %d %.3f 1.0 2.5          # A-A attractive" % (i, j, a.eps_aa))
            else:
                s.append("pair_coeff %d %d 1.0 1.0 1.122462048   # WCA  %d-%d" % (i, j, tname[i], tname[j]))
    s += ["pair_modify shift yes", "",
          "bond_style harmonic", "bond_coeff * 100.0 1.0", "special_bonds lj 0.0 1.0 1.0", "",
          "angle_style harmonic"]
    for key, t in combos.items():
        s.append("angle_coeff %-3d %4.1f 180.0   # %s" % (t, K_TABLE.get(key, 0.0), "-".join(map(str, key))))
    return "\n".join(s)


def groups_block(tname):
    anc = types_of(tname, 1); sh = types_of(tname, 4)
    pol = [t for t in sorted(tname) if t not in anc + sh]
    return ("group anchor  type %s\ngroup shell   type %s\ngroup polymer type %s\n"
            "group fixed   union anchor shell" % (" ".join(map(str, anc)), " ".join(map(str, sh)), " ".join(map(str, pol))))


COMMON_TAIL = """
neighbor        0.3 bin
neigh_modify    delay 0 every 1 check yes
neigh_modify    exclude group fixed fixed

fix freeze fixed setforce 0.0 0.0 0.0
compute tpol polymer temp
compute rg polymer gyration
compute zmax polymer reduce max z
compute zmin polymer reduce min z
thermo_style custom step c_tpol pe ke etotal c_rg c_zmin c_zmax   # c_tpol = polymer temperature
thermo_modify norm no
"""

IN_FRESH = """# equilibration of {sys}  (k = 5 angles, A-A eps {eps})  -- fresh start
units           lj
dimension       3
boundary        f f f
atom_style      full
read_data       {sys}.data
comm_modify     cutoff 3.0

{groups}

{ff}
{tail}
velocity fixed set 0.0 0.0 0.0
thermo 100
minimize 1.0e-8 1.0e-8 10000 100000
reset_timestep 0

velocity polymer create 1.0 {seed} mom yes rot yes
fix nvt_polymer polymer nvt temp 1.0 1.0 1.0
timestep 0.005
thermo {thermo}
dump traj all dcd {dump} eq_{sys}_part0.dcd
restart {rst} restart.a restart.b
run {steps} upto
write_restart restart.final
write_data {sys}_eq.data nocoeff
print "Equilibration completed"
"""

IN_RESUME = """# equilibration of {sys} -- resume from the newest restart file (set by the job: -var rfile ...)
units           lj
dimension       3
boundary        f f f
atom_style      full
read_restart    ${{rfile}}
comm_modify     cutoff 3.0
{tail}
velocity fixed set 0.0 0.0 0.0
fix nvt_polymer polymer nvt temp 1.0 1.0 1.0
timestep 0.005
thermo {thermo}
variable s0 equal step
dump traj all dcd {dump} eq_{sys}_part${{s0}}.dcd
restart {rst} restart.a restart.b
run {steps} upto
write_restart restart.final
write_data {sys}_eq.data nocoeff
print "Equilibration completed"
"""

RESTORE = r'''#!/usr/bin/env python3
# put the "# name" comments back into Masses of <SYS>_eq.data (write_data drops them)
import sys, re
sysname = sys.argv[1]
names = {}
for l in open(sysname + ".data"):
    pass
L0 = open(sysname + ".data").read().split("\n")
i = next(k for k, l in enumerate(L0) if l.split("#")[0].strip() == "Masses") + 1
while not L0[i].strip(): i += 1
while i < len(L0) and L0[i].strip():
    b, _, c = L0[i].partition("#"); names[int(b.split()[0])] = c.strip(); i += 1
L = open(sysname + "_eq.data").read().split("\n")
i = next(k for k, l in enumerate(L) if l.split("#")[0].strip() == "Masses") + 1
while not L[i].strip(): i += 1
while i < len(L) and L[i].strip():
    t = int(L[i].split()[0]); L[i] = L[i].split("#")[0].rstrip() + " # " + names[t]; i += 1
open(sysname + "_eq.data", "w").write("\n".join(L))
print("names restored in %s_eq.data" % sysname)
'''

SLURM = """#!/bin/bash
#SBATCH --job-name=eq_{tag}
#SBATCH --partition={partition}
#SBATCH --nodes=1
#SBATCH --ntasks={cores}
#SBATCH --cpus-per-task=1
#SBATCH --time={time}
#SBATCH --output=eq_{tag}_%j.out
#SBATCH --error=eq_{tag}_%j.err
# Equilibration: {desc}, all at once.
# Resubmit the same file to continue unfinished runs:  sbatch {me}
{modules}
export OMP_NUM_THREADS=1
cd "${{SLURM_SUBMIT_DIR:-$(pwd)}}"

LMP_BIN=${{LMP_BIN:-{lmp}}}
if [ -z "$LMP_BIN" ] || [ ! -x "$LMP_BIN" ]; then
    for c in $(command -v lmp 2>/dev/null) $(command -v lmp_mpi 2>/dev/null) \\
             $HOME/lammps-29Aug2024/build/lmp $HOME/lammps/build/lmp; do
        [ -x "$c" ] && LMP_BIN=$c && break
    done
fi
[ -x "$LMP_BIN" ] || {{ echo "!!! LAMMPS not found - set LMP_BIN at the top or: sbatch --export=ALL,LMP_BIN=/path/lmp {me}"; exit 1; }}
echo "LAMMPS: $LMP_BIN"
export LMP_BIN

run_one() {{
    d=$1; np=$2; sys=$(basename $d)
    cd "$d" || return
    if grep -q "Equilibration completed" log.lammps 2>/dev/null && [ -f ${{sys}}_eq.data ]; then
        echo "skip  $sys (done)"; return
    fi
    rf=$(ls -t restart.a restart.b 2>/dev/null | head -1)
    if [ -n "$rf" ]; then
        n=$(ls log_part*.lammps 2>/dev/null | wc -l)
        echo "resume $sys from $rf  ($np cores)  $(date)"
        mpirun --bind-to none -np $np "$LMP_BIN" -in in.resume -var rfile $rf -log log_part$n.lammps > stdout_resume$n.log 2>&1
        cp log_part$n.lammps log.lammps
    else
        echo "start $sys  ($np cores)  $(date)"
        mpirun --bind-to none -np $np "$LMP_BIN" -in in.lammps > stdout.log 2>&1
    fi
    if grep -q "Equilibration completed" log.lammps; then
        python3 ../restore_names.py $sys > /dev/null 2>&1 || echo "  (restore names later: cd $d; python3 ../restore_names.py $sys)"
        echo "done  $sys  $(date)"
    else
        echo "STOP  $sys  (time limit or error; see $d/log.lammps)  $(date)"
    fi
}}
export -f run_one

run_group() {{
    g=$1; np=$2; n=$(wc -l < jobs_$g.txt)
    echo "group $g: $n systems, $np cores each"
    xargs -a jobs_$g.txt -P $n -I{{}} bash -c 'run_one "$@"' _ {{}} $np
}}

{launch}
wait
echo "all finished $(date)"
"""

# ------------------------------------------------------------------ main
os.makedirs(a.outdir, exist_ok=True)
groups = OrderedDict((g, []) for g in a.groups)
skipped = []
for f in a.datafiles:
    base = os.path.basename(f)
    m = re.match(r"(\d+)_([a-zA-Z]+)\d*\.data$", base)
    if not m:
        skipped.append(base); continue
    g = alias.get(m.group(2), m.group(2))
    if g in groups:
        groups[g].append(f)
    else:
        skipped.append(base)
if skipped:
    print("not in groups %s (ignored): %s" % ("/".join(a.groups), " ".join(sorted(skipped))))
for g in groups:
    if not groups[g]:
        sys.exit("no data files for group '%s'" % g)
    groups[g].sort(key=lambda f: int(re.match(r"\d+", os.path.basename(f)).group()))

if len(a.groups) > 2: sys.exit("give one or two groups")
share = a.cores // len(a.groups)
nps = a.np or [max(1, share // len(groups[g])) for g in groups]
open(os.path.join(a.outdir, "restore_names.py"), "w").write(RESTORE)
for (g, files), npg in zip(groups.items(), nps):
    print("\n=== group %s: %d systems, %d cores each (%d of %d)" % (g, len(files), npg, npg * len(files), share))
    dirs = []
    for f in files:
        sysname = os.path.basename(f)[:-5]
        d = os.path.join(a.outdir, sysname); os.makedirs(d, exist_ok=True)
        tname, combos = prepare(f, os.path.join(d, sysname + ".data"))
        kw = dict(sys=sysname, eps=a.eps_aa, groups=groups_block(tname), ff=ff_block(tname, combos),
                  tail=COMMON_TAIL, seed=a.seed, thermo=a.thermo, dump=a.dump_every, rst=a.restart_every,
                  steps=a.steps)
        open(os.path.join(d, "in.lammps"), "w").write(IN_FRESH.format(**kw))
        open(os.path.join(d, "in.resume"), "w").write(IN_RESUME.format(**kw))
        nA = len(types_of(tname, 5))
        ks = ", ".join("%s:%g" % ("-".join(map(str, k)), K_TABLE.get(k, 0.0)) for k in combos)
        print("  %-10s types %s  A types %s  angles %s" % (sysname, dict(sorted(tname.items())), types_of(tname, 5), ks))
        dirs.append(sysname)
    open(os.path.join(a.outdir, "jobs_%s.txt" % g), "w").write("\n".join(dirs) + "\n")

tag = "_".join(a.groups)
me = "job_equil_%s.slurm" % tag
mods = ("module purge\n" if a.purge else "") + "\n".join("module load %s" % m for m in a.modules) if a.modules else "# (no modules)"
desc = " and ".join("group %s (%d systems x %d cores)" % (g, len(groups[g]), n) for g, n in zip(groups, nps))
launch = "\n".join("run_group %s %d &" % (g, n) for g, n in zip(groups, nps))
open(os.path.join(a.outdir, me), "w").write(SLURM.format(
    tag=tag, desc=desc, launch=launch, partition=a.partition, cores=a.cores, time=a.time, modules=mods,
    lmp=a.lmp, me=me))
print("\nwritten: %s/  (one folder per system), %s/%s" % (a.outdir, a.outdir, me))
print("submit:  cd %s && sbatch %s" % (a.outdir, me))
