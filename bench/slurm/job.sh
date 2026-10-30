#!/bin/bash -l
# Run bench/run commands as a Slurm job on the BioHPC BSCB cluster, from a node-local work dir.
#
#   sbatch --dependency=singleton [-t HH:MM:SS] bench/slurm/job.sh '<bench/run args>' ['<bench/run args>' ...]
#   e.g. sbatch --dependency=singleton -t 12:00:00 bench/slurm/job.sh 'plan --tier standard --max-threads 64' \
#            'exec --suites equivalence --parallel 2' 'exec --suites paired --max-load 1000' report
#
# Every timed run so far ran on cbsubscb19 (2x AMD EPYC 9554, AVX-512) from /workdir/$USER/PINTS, in a 64-CPU
# allocation, so this job does too: nodes differ in CPU model, which changes speed and numpy's SIMD kernels. The
# node-local repo is brought to the shared repo's HEAD (kit and branches) before the commands run, and the work
# dir's results are copied back to the shared one afterwards, also when the job fails or times out. If the
# node-local copy is gone (scratch cleanup), it is rebuilt from the shared repo. The fixed job name together with
# --dependency=singleton keeps two bench jobs from sharing the work dir.
#SBATCH --job-name=pints-bench
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=256G
#SBATCH --nodelist=cbsubscb19
#SBATCH --partition=regular,long7
#SBATCH --time=24:00:00
#SBATCH --output=/home2/ayh8/PINTS/bench_work/slurm-logs/%x.%j.out
set -euo pipefail

SRC=${PINTS_BENCH_SRC:-/home2/ayh8/PINTS}   # shared repo: kit, branches, and the copy of the results
NODE=/workdir/$USER/PINTS                     # node-local repo and work dir
export PINTS_BENCH_WORK=$NODE/bench_work
unset PYTHONPATH                              # the login profile copies PATH into it

echo "job $SLURM_JOB_ID on $(hostname), cpus $(taskset -pc $$ | sed 's/.*: //'), kit $(git -C "$SRC" rev-parse --short HEAD)"
if [ ! -d "$NODE/.git" ]; then
  echo "no node-local copy; rebuilding from $SRC"
  git clone -q "$SRC" "$NODE"
  rsync -a --exclude impl/ --exclude venv/ --exclude slurm-logs/ "$SRC/bench_work/" "$PINTS_BENCH_WORK/"
fi

sync_back() {
  rsync -a --exclude data/ --exclude venv/ --exclude impl/ --exclude uv-cache/ --exclude tools/ \
    --exclude slurm-logs/ "$PINTS_BENCH_WORK/" "$SRC/bench_work/" && echo "results copied to $SRC/bench_work"
}
trap sync_back EXIT
trap 'exit 143' TERM

git -C "$NODE" fetch -q "$SRC" '+refs/heads/*:refs/remotes/shared/*'
git -C "$NODE" checkout -q -f --detach "$(git -C "$SRC" rev-parse HEAD)"
bash "$NODE/bench/setup.sh" >"$PINTS_BENCH_WORK/setup.log" 2>&1 || { tail -20 "$PINTS_BENCH_WORK/setup.log"; exit 1; }
"$NODE/bench/run" sysinfo >/dev/null

for cmd in "$@"; do
  echo "==> bench/run $cmd  ($(date '+%F %T'), load $(cut -d' ' -f1 /proc/loadavg))"
  # shellcheck disable=SC2086
  "$NODE/bench/run" $cmd
done
