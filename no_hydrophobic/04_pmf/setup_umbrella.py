#!/usr/bin/env python3
"""
Axial umbrella sampling of the cargo PMF(z), HALF pore (z = 0 .. ZMAX), mirrored in analysis.
Same force field as the transport runs (imports setup_hydro.py: types, cargo, WCA + A-A pairs,
your angle table K_TABLE), so PMFs and transport runs describe the same pore.

    python setup_umbrella.py 50_tp05.data 75_tp05.data --sigmas 4 --reps 3
    python setup_umbrella.py 0_05.data 25_05.data 50_05.data 75_05.data 100_05.data --sigmas 4
    python setup_umbrella.py 50_tp05.data --sigmas 4 --zmin -50     # FULL pore (hydrophobes: not symmetric)

Window scheme (Patel, Soni & Panwar, PCCP 2026: spacing 1 l0, k = 10, mirrored about z = 0):
    z0 = 0, 1, ..., ZMAX (default 50: beyond the tallest brush ~40 + cargo radius)
    U = 0.5 k (z - z0)^2   (LAMMPS fix spring tether; checked: E = 5 at k = 10, dz = 1)
Per window (one folder each):
    cargo put at (0, 0, z0) and HELD -> minimize -> RELAX steps (chains settle around it)
    -> release: Langevin on cargo (as in transport runs) + z spring + flat-bottom cylinder
       (cargo r <= RC, default 10: keeps it over the pore instead of wandering sideways)
    -> EQUIL steps (discarded) -> PROD steps, cargo x y z written to colvar.dat every 10 steps.
Folders:  <SYS>/umbrella_s<sigma>/rep<r>/z<z0>/in.lammps  (+ colvar.dat, log, restart)
Job file: <SYS>/job_umb_<SYS>.slurm  (same cluster settings as setup_hydro.py), submitted unless --no-start.
Analyse with:  python wham_half.py <SYS>/umbrella_s4.0          (axial; also gives the 2D map W(r,z))

RADIAL mode (--radial z1 z2 ...): at each height the cargo z is held by the z spring (k) and r is
windowed: r0 = 0, dr, ..., rmax with U = 0.5 kr (r - r0)^2 (fix spring tether K 0 0 NULL r0).
    python setup_umbrella.py 50_tp05.data --sigmas 4 --radial 12.5 6 0 -6 -12.5
Folders: <SYS>/umbrella_r_s<sigma>/rep<r>/zs<z>_r<r0>/ ;  analyse with  python wham_radial.py <SYS>/umbrella_r_s4.0
"""
import argparse, os, random, re, shutil, subprocess, sys
import setup_hydro as H

TEMPLATE = """# ---- umbrella window: {sys}, cargo sigma {sigma}, rep {rep}, {title} ----
units           lj
dimension       3
boundary        f f f
atom_style      full
read_data       {datafile}
comm_modify     cutoff {comm_cut:.2f}
log             log.lammps

group anchor  type {t_anchor}
group polymer type {t_polymer}
group shell   type {t_shell}
group cargo   type {t_cargo}
group fixed   union anchor shell

set group cargo x {x0} y 0.0 z {z0}              # window centre
fix freeze fixed setforce 0.0 0.0 0.0
fix hold_cargo cargo setforce 0.0 0.0 0.0      # cargo held during minimize + relax
velocity fixed set 0.0 0.0 0.0

pair_style      lj/cut {pair_cut:.9f}
{pair_block}
pair_modify shift yes

bond_style harmonic
bond_coeff 1 100.0 1.0
special_bonds lj 0.0 1.0 1.0
angle_style harmonic
{angle_block}

neighbor 1.0 bin
neigh_modify delay 0 every 1 check yes
neigh_modify exclude group fixed fixed

compute cargox cargo reduce sum x
compute cargoy cargo reduce sum y
compute cargoz cargo reduce sum z
thermo_style custom step temp pe c_cargox c_cargoy c_cargoz
thermo_modify norm no

# ---- 1. minimize + relax chains around the held cargo ----
min_modify dmax 0.1
minimize 1.0e-8 1.0e-8 10000 100000
velocity polymer create 1.0 {seed} mom yes rot yes
velocity cargo set 0.0 0.0 0.0
fix nvt_polymer polymer nvt temp 1.0 1.0 0.5
fix nve_cargo cargo nve
timestep 0.005
thermo 5000
run {relax}

# ---- 2. release: Langevin cargo + umbrella spring on z + flat-bottom cylinder r <= {rc:g} ----
unfix hold_cargo
fix lang_cargo cargo langevin 1.0 1.0 {damp:g} {lseed}
fix umbrella cargo spring tether {k:g} NULL NULL {z0} 0.0
{rspring}region rcyl cylinder z 0.0 0.0 {rwall:g} INF INF side in units box
fix rwall cargo wall/region rcyl harmonic 50.0 1.0 2.0     # pushes back only for r > {rc:g}
run {equil}

# ---- 3. production: cargo position every {every} steps ----
fix colvar all print {every} "$(step) $(c_cargox:%.4f) $(c_cargoy:%.4f) $(c_cargoz:%.4f)" file colvar.dat screen no title "# step x y z   ({title})"
{dump_line}run {prod}
write_restart restart.umb
print "Umbrella window {sys} sigma {sigma} rep {rep} z0 {z0} completed"
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("datafiles", nargs="+")
    ap.add_argument("--out", default=".")
    ap.add_argument("--sigmas", nargs="+", default=["4"])
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--zmin", type=float, default=0.0)
    ap.add_argument("--zmax", type=float, default=50.0)
    ap.add_argument("--dz", type=float, default=1.0)
    ap.add_argument("--k", type=float, default=10.0)
    ap.add_argument("--radial", nargs="+", type=float, default=None, metavar="Z",
                    help="RADIAL mode: windows in r = sqrt(x^2+y^2) at these heights (z held by the z spring)")
    ap.add_argument("--rmax", type=float, default=9.0)
    ap.add_argument("--dr", type=float, default=0.5)
    ap.add_argument("--kr", type=float, default=10.0)
    ap.add_argument("--rc", type=float, default=10.0, help="flat-bottom cylinder radius for the cargo")
    ap.add_argument("--relax", type=int, default=5000)
    ap.add_argument("--equil", type=int, default=20000)
    ap.add_argument("--prod", type=int, default=500000)
    ap.add_argument("--every", type=int, default=10)
    ap.add_argument("--dump", type=int, default=0, help="also write a DCD every N steps (0 = none; DCDs are big)")
    ap.add_argument("--damp", type=float, default=1.0)
    ap.add_argument("--aa-eps", type=float, default=2.0)
    ap.add_argument("--aa-cut", type=float, default=2.5)
    ap.add_argument("--seed", type=int, default=24680)
    ap.add_argument("--no-start", action="store_true")
    ap.add_argument("--par", type=int, default=5)
    ap.add_argument("--np", type=int, default=4)
    ap.add_argument("--ntasks", type=int, default=20)
    ap.add_argument("--partition", default="small")
    ap.add_argument("--time", default="3-00:00:00")
    ap.add_argument("--lmp", default=H.LMP_BIN)
    a = ap.parse_args()
    if a.par * a.np > a.ntasks:
        sys.exit(f"--par x --np = {a.par * a.np} is more than --ntasks {a.ntasks}")
    rng = random.Random(a.seed)
    base = os.path.abspath(a.out); os.makedirs(base, exist_ok=True)
    nz = int(round((a.zmax - a.zmin) / a.dz)) + 1
    zs = [round(a.zmin + i * a.dz, 3) for i in range(nz)]
    radial = a.radial is not None
    if radial:
        nr = int(round(a.rmax / a.dr)) + 1
        wins = [(zh, round(i * a.dr, 3)) for zh in a.radial for i in range(nr)]
    else:
        wins = [(z0, 0.0) for z0 in zs]
    tag = "umbr" if radial else "umb"
    slurms = []
    for src in a.datafiles:
        fname = os.path.basename(src)
        sysname = re.sub(r"_?(05)?\.data$", "", fname)
        sysdir = os.path.join(base, sysname); os.makedirs(sysdir, exist_ok=True)
        run_data = os.path.join(sysdir, "umb_" + fname)
        T = H.prepare_data(src, run_data, H.k_rule, True)
        print(f"\n=== {sysname}  ({fname})")
        for m in T["msgs"]:
            print("   " + m)
        print(f"   types: anchor={T['anchor']} polymer={T['polymer']} shell={T['shell']} "
              f"A={T['A'] or 'none'} cargo={T['cargo']}")
        for at, key, kk in T["table"]:
            print(f"      angle type {at:2d}  {'-'.join(map(str, key)):6s} K = {kk:g}")
        jobs = []
        for s in a.sigmas:
            sigma = float(s)
            for rep in range(1, a.reps + 1):
                for z0, r0 in wins:
                    if radial:
                        d = os.path.join(sysdir, f"umbrella_r_s{sigma:.1f}", f"rep{rep}", f"zs{z0:g}_r{r0:g}")
                        title = f"zs = {z0:g}, r0 = {r0:g}, kz = {a.k:g}, kr = {a.kr:g}, U = 0.5 kr (r - r0)^2"
                        rspring = (f"fix umb_r cargo spring tether {a.kr:g} 0.0 0.0 NULL {r0:g}"
                                   f"     # radial window: U = 0.5 kr (r - r0)^2\n")
                    else:
                        d = os.path.join(sysdir, f"umbrella_s{sigma:.1f}", f"rep{rep}", f"z{z0:g}")
                        title = f"z0 = {z0:g}, k = {a.k:g}, U = 0.5 k (z - z0)^2"
                        rspring = ""
                    os.makedirs(d, exist_ok=True)
                    txt = TEMPLATE.format(
                        sys=sysname, sigma=s, rep=rep, z0=f"{z0:g}", x0=f"{r0:g}", title=title, rspring=rspring,
                        datafile=os.path.relpath(run_data, d),
                        comm_cut=max(4.0, H.WCA * (sigma + 1) / 2 + 1.5),
                        pair_cut=max(H.WCA, a.aa_cut if T["A"] else H.WCA),
                        t_anchor=T["anchor"], t_polymer=" ".join(map(str, T["polymer"])),
                        t_shell=T["shell"], t_cargo=T["cargo"],
                        pair_block=H.pair_block(T, sigma, a.aa_eps, a.aa_cut),
                        angle_block=T["angle_block"],
                        seed=rng.randint(10000, 9999999), lseed=rng.randint(10000, 9999999),
                        relax=a.relax, equil=a.equil, prod=a.prod, every=a.every,
                        k=a.k, rc=a.rc, rwall=a.rc + 2.0, damp=a.damp,
                        dump_line=(f"dump traj all dcd {a.dump} dump.dcd\n" if a.dump > 0 else ""))
                    open(os.path.join(d, "in.lammps"), "w").write(txt)
                    jobs.append(os.path.relpath(d, sysdir))
        random.Random(a.seed + 7).shuffle(jobs)
        open(os.path.join(sysdir, f"jobs_{tag}.txt"), "w").write("\n".join(jobs) + "\n")
        sl = os.path.join(sysdir, f"job_{tag}_{sysname}.slurm")
        open(sl, "w").write(H.SLURM.format(sys=f"{tag}_" + sysname, partition=a.partition, ntasks=a.ntasks,
                                           time=a.time, sysdir=sysdir, njobs=len(jobs), par=a.par,
                                           np=a.np, env=H.ENV_BLOCK, lmp=a.lmp)
                            .replace("jobs.txt", f"jobs_{tag}.txt")
                            .replace("slurm_%j", f"slurm_{tag}_%j")
                            .replace("Simulation .*: cargo", "Umbrella window .* completed"))
        os.chmod(sl, 0o755); slurms.append(sl)
        if radial:
            print(f"   {len(jobs)} RADIAL windows ({len(a.sigmas)} sigma x {a.reps} reps x heights {a.radial} x "
                  f"r0 = 0..{a.rmax:g} step {a.dr:g}) -> {sysdir}/umbrella_r_s*/rep*/zs*_r*   job {sl}")
        else:
            print(f"   {len(jobs)} windows ({len(a.sigmas)} sigma x {a.reps} reps x {nz} z0 from "
                  f"{a.zmin:g} to {a.zmax:g}) -> {sysdir}/umbrella_s*/rep*/z*   job {sl}")
    if a.no_start:
        print("\n--no-start: submit with:  " + "  ;  ".join(f"sbatch {s}" for s in slurms)); return
    if shutil.which("sbatch"):
        for s in slurms:
            r = subprocess.run(["sbatch", s], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            print(f"sbatch {os.path.basename(s)}: {(r.stdout or r.stderr).strip()}")
    else:
        for s in slurms:
            log = s.replace(".slurm", ".log")
            subprocess.Popen(["nohup", "bash", s], stdout=open(log, "w"), stderr=subprocess.STDOUT,
                             start_new_session=True)
            print(f"no sbatch here -> started locally: bash {s}  (log {log})")


if __name__ == "__main__":
    main()
