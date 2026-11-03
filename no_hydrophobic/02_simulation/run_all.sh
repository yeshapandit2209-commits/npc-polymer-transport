#!/bin/bash
# Runs every folder listed in jobs.txt.
#   PAR   = how many runs at the same time
#   NP    = MPI processes per run   (PAR x NP should not exceed the cores you have)
# Examples (20 cores):
#   bash run_all.sh            -> 5 runs at a time, 4 cores each   (default, best throughput)
#   PAR=1 NP=20 bash run_all.sh -> one run at a time on 20 cores     (like run_50_transport.sh)
# Finished runs (stdout.log contains 'completed') are skipped, so you can restart safely.
LMP_BIN=${LMP_BIN:-~/lammps-29Aug2024/build_full/lmp}
PAR=${PAR:-5}
NP=${NP:-4}
JOBS=/data/yesha/npc/md/passive/1/2/runs_clean/jobs.txt
export LMP_BIN NP
run_one() {
    d="$1"
    if grep -q "completed" "$d/stdout.log" 2>/dev/null; then echo "skip  $d"; return; fi
    if ( cd "$d" && mpirun -np $NP "$LMP_BIN" -in in.lammps > stdout.log 2>&1 ); then
        echo "done  $d"
    else
        echo "FAIL  $d (see $d/stdout.log)" | tee -a failed_runs.txt
    fi
}
export -f run_one
echo "running $(wc -l < $JOBS) runs: $PAR at a time, $NP MPI processes each"
xargs -a $JOBS -P $PAR -I{} bash -c 'run_one "$@"' _ {}
echo "all runs finished"; [ -f failed_runs.txt ] && echo "failed runs listed in failed_runs.txt"
