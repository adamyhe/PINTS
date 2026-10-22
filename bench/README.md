# PINTS optimisation benchmark kit

This kit checks equivalence, runtime, peak memory and thread scaling. It compares the upstream reference against the optimisation branches:

| variant | checkout | what it is |
|---|---|---|
| `ref` | `main` (upstream 81b2dd9) | reference |
| `gate1` | `perf/gate1-bit-identical` | faster, must be byte-identical |
| `gate2` | `perf/gate2-fast-math` with `PINTS_FAST_MATH=1` | compiled ZIP fit, judged against a drift budget |
| `gate2_off` | `perf/gate2-fast-math`, flag unset | must be byte-identical |
| `trackb` | `fix/trackB-flags`, no new flags | must be byte-identical |
| `pr4` | `fix/empty-frame-dtype` (`main` + one fix) | must be byte-identical |
| `io` | `perf/io-memory` (gate1 + I/O and memory commits) | must be byte-identical |

Each variant is a git worktree pinned to the commit in `variants.json`. The kit runs that worktree's unmodified `scripts/pints_caller`.

## Requirements

- Linux x86_64. A 64+ core machine with ≥ 128 GB RAM covers every thread count in the plans.
- About 25 GB free disk: 5.8 GB of ENCODE data, plus prepared BAMs and run outputs.
- [uv](https://docs.astral.sh/uv/), `git`, `curl`, and for building samtools/htslib: `gcc`, `make` and the zlib headers (`zlib1g-dev` / `zlib-devel`). If `bgzip`, `tabix`, `samtools` and `bedtools` are already on `PATH`, the kit uses them and builds nothing. Otherwise it builds htslib/samtools 1.22.1 and downloads the static bedtools 2.31.0 binary (Biomni used 2.31.1, which has no static build; every variant uses the same bedtools, so this cannot affect equivalence).
- Outbound HTTPS to encodeproject.org and github.com.

## Steps

```bash
git clone https://github.com/adamyhe/PINTS.git && cd PINTS
git checkout bench/external
bash bench/setup.sh                    # ~10 min: venv, tools, implementation worktrees -> ./bench_work
bench/run sysinfo                      # check CPU count, cgroup limits, governor
bench/run fetch                        # downloads + verifies ENCODE data (md5), subsets BAMs
bench/run plan --tier standard         # prints the job list size and an estimated duration
tmux new -s pints                      # long run: keep it in tmux/screen
bench/run exec --parallel 2            # equivalence runs 2 at a time, then timing/scaling one at a time
bench/run status                       # from another shell: progress and failures
bench/run report                       # bench_work/report/report.md + bench_work/bench_results_<host>_<date>.tar.gz
```

Set `PINTS_BENCH_WORK=/fast/disk/pints_bench` before `setup.sh` (and keep it set) to put everything on another disk. Prefer local NVMe over NFS, because PINTS writes many intermediate files.

`exec` can be interrupted at any time. Rerunning it skips runs that already succeeded. `--retry-failed` reruns failures. `--suites timing` (or `equivalence`, `scaling`) restricts it to one suite. `--dry-run` lists the remaining runs.

**Bring back** `bench_work/bench_results_<host>_<date>.tar.gz`. It holds all run records, memory traces, phase logs, sysinfo and the report, but no data or peak files. The gate1/gate2 candidate tables needed for the drift check are summarised in the report, and the run directories stay under `bench_work/runs/` if you need to dig further.

## Re-running after new commits

Each record stores the commit it ran (`impl_sha`). A run counts as done only if that commit is the one its variant is pinned to now in `bench/variants.json`. After the branches gain commits, update the kit and rerun in the same work directory:

```bash
git pull                                      # new SHAs in bench/variants.json
bash bench/setup.sh                           # moves the implementation worktrees to the new commits
bench/run plan --tier standard --max-threads 64
bench/run exec --suites equivalence --parallel 2
bench/run exec                                # timing + scaling
bench/run report
```

Variants whose commit did not change, such as `ref`, keep their runs, so only the changed implementations rerun. `plan` estimates each job from its last measured run. Older records stay in `runs.jsonl`, and the report's section 6 compares each variant's new commit with the one it replaced on the same timed jobs.

`pr4` (`main` plus the empty-frame dtype fix) runs at reference speed, so it gets a reduced plan: equivalence on every dataset in the default config plus all three configs on H3 (the dataset where the fix triggers), and timing on H3.

## I/O and memory commits (`io`)

`perf/io-memory` adds three commits to gate1, from the profiling results:
- block-wise replicate merge, for H1's memory peak;
- early, background deletion of coverage dumps;
- the stratified-filtering table written in parallel chunks, without `bedtools sort`.

It must match the reference byte for byte like gate1. The `paired` suite times it against gate1 in the same session: alternating, at 1, 16 and 64 threads on the timing datasets, with fresh gate1 runs, so session-to-session noise does not enter the comparison. In the same work directory:

```bash
git pull && bash bench/setup.sh                               # adds the io checkout
bench/run plan --tier standard --suites equivalence,paired --max-threads 64
bench/run exec --suites equivalence --parallel 2              # only the io runs are new
bench/run exec --suites paired                                # ~2.5 h
bench/run report                                              # section 7: io vs gate1
bench/run profile --variant io                                # optional, profiles io like gate1
```

## Profiling

`bench/run profile` profiles one implementation (default `gate1`, default config) on H1 and H3. It runs in the same work directory, needs `bench/run fetch` to have been run, and takes about 1.5 h:

```bash
git pull && bash bench/setup.sh               # setup installs py-spy
bench/run profile --dry-run                   # lists the 9 runs
bench/run profile                             # -> bench_work/profile_results_<host>_<date>.tar.gz
```

Each mode is its own run, so the profilers don't distort each other:

| Mode | Runs | What it records |
|---|---|---|
| `trace` | H1 at 1, 16, 64 threads; H3 at 1, 64 | how long the main functions take in the parent and in every pool worker (including the bedtools and tabix calls), memory and I/O of each process every 0.25 s, and every file that appears in the run directory, with size and lifetime |
| `pyspy` | H1, H3 at 64 threads | sampled call stacks of all processes (py-spy) |
| `cprofile` | H1, H3 at 1 thread | exact call counts and times of the main process, which does all the work at 1 thread |

- **Output check:** each run's peak files are compared with the benchmark runs of the same commit; the progress line says if they match.
- **Clean-up:** the run directories are deleted afterwards; only the profiles are kept (tens of MB).
- **Options:**
  - `--modes`, `--datasets` and `--variant` restrict or change the runs;
  - finished runs are skipped, and `--force` reruns them.
- **If py-spy fails:** most likely the node doesn't allow tracing processes (`/proc/sys/kernel/yama/ptrace_scope` above 1). The other two modes still work.

## Tiers

| tier | equivalence | timing | scaling | rough duration |
|---|---|---|---|---|
| `smoke` | T1_bw, T1_bam, default config, all variants | ref/gate1/gate2 at 1 and 8 threads | T1_bam, 1–8 threads | ~15 min |
| `standard` | all 19 datasets × {default, sensitive, disable_small} × all variants at 32 threads, plus determinism checks for gate1 (1 thread, hash seed 1234, all CPUs) on T1–T5 | 8 datasets (T2–T5, H1–H3, H6) × {1, 16, 64} threads × 3 reps (ref at 1 thread: 1 rep) | T3, T4, H1, H3: gate1/gate2 at 1–64 threads and all CPUs, 2 reps; ref at 1/4/16/64 | ~67 h sequential at 64 CPUs, ~50 h with `exec --parallel 2` (about 2 days) |
| `full` | as standard, determinism also on H1–H3 | all 19 datasets, + all CPUs, 3 reps | 8 datasets, 10 thread counts, 3 reps | ~130 h sequential, ~115 h with `--parallel 2` (about 5 days) |

Run `bench/run plan --tier <t>` to see the job count and estimated hours for the machine you are on. Thread counts are capped at the CPUs available to the process. `--max-threads N` lowers the cap, and `--datasets a,b` restricts the plan to some datasets. Most of the time goes to the reference: it barely scales with threads, and Biomni measured 30+ minutes per run on genome-wide data even at 16 threads for the larger sets. To get results sooner, run `bench/run exec --suites equivalence --parallel 2` first, then `bench/run report` (a partial report is fine), then the timing and scaling suites.

## What is measured

Each run gets a fresh output directory under `bench_work/runs/`, and its record is appended to `bench_work/results/runs.jsonl`:

- **Equivalence:** SHA-256 of the final `*_bidirectional_peaks.bed`, `*_divergent_peaks.bed` and `*_unidirectional_peaks.bed`. Every successful gate1/gate2_off/trackb run is compared with the reference run of the same dataset and config, whatever the thread count or hash seed. The report also checks that the reference agrees with itself.
- **Gate 2 drift:** gate1 and gate2 equivalence runs keep their per-strand candidate tables. The report compares candidate coordinates, p/padj relative differences, significance flips and final-set Jaccard against Biomni's pre-registered budget.
- **Runtime:** wall time, plus CPU time from `wait4` rusage (`cpu_util` is the average number of busy cores). Timed runs run one at a time and wait for the 1-minute load average to fall below `--max-load` (default 2). Implementations alternate order between repetitions.
- **Peak memory**, sampled every 0.25 s over the whole process tree:
  - `peak_pss_gb`: summed PSS from `/proc/<pid>/smaps_rollup`, the primary metric. It splits pages shared between forked workers instead of counting them once per worker.
  - `peak_rss_sum_gb`: summed RSS, which overcounts shared pages.
  - `peak_sys_delta_gb`: MemTotal − MemAvailable above the pre-run baseline. It includes page cache and is noisy if anything else is running.
  - `max_single_rss_gb`: the largest single process.
  - `cgroup_peak_gb`: cgroup v2 `memory.peak`, only with `exec --systemd-scope`, which runs each job in its own `systemd-run --user --scope`. This is the most faithful number if your system supports it.
- **Phases:** with `PINTS_PHASE_TIMING=1` the optimised branches print `[phase]` markers. The runner timestamps them and records the duration and peak PSS of each phase. The reference prints none.

The report states whether the reference and each variant scale with threads. `ref` itself only parallelises parts of the pipeline, so its curve flattens early.

## Notes

- Runs pin `OMP/MKL/OPENBLAS/NUMBA_NUM_THREADS=1`, `PYTHONHASHSEED` (0 unless the run is a seed check), `LC_ALL=C`, and `TMPDIR` set to the run directory.
- For stable timings, use the `performance` CPU governor if you can, keep the machine otherwise idle, and avoid running on a node shared with other jobs. `sysinfo` records the governor, THP setting, cgroup limits and `numpy.show_runtime()` (the SIMD extensions in use, which matter for Gate 2 drift).
- `fetch` checks every ENCODE file against its md5 and records the SHA-256 of every prepared input. It also flags any prepared file whose SHA-256 differs from the Biomni run's. A different samtools version can recompress subset BAMs differently with identical records, so such a flag is informational.
