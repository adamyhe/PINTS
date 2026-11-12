## Keep coverage in an in-memory sparse store instead of per-base `.npy` dumps

Replaces the per-chromosome `.npy` coverage dumps with an in-memory store of the non-zero bases. **The output is
byte-identical to the base branch.** One behaviour changes: `.npy` dumps left in the output directory by an
earlier run are no longer reused (below). Stacked on B4 (`pr/04-io`); review the last ten commits.

### What changes for users

- **Old `.npy` dumps are no longer written or reused.** On `main`, `get_read_signal` and `get_coverage_bw` reuse
  per-chromosome dumps they find in the output directory under the same file prefix, and log "Loading cache from
  previous result". A successful run deletes its dumps, so this only helps when a run is repeated after a crash.
  It can also silently reuse dumps from a different input that was run with the same output directory and prefix,
  and call peaks on the wrong coverage. With the store, coverage is always built from the BAM or bigWig; dumps
  already in the output directory are left untouched. The public functions `get_read_signal` and
  `get_coverage_bw` keep `main`'s behaviour.
- **`PINTS_COVERAGE_STORE=0`** (read once, in `pints/calling_engine.py`) keeps the dump-based path, with `main`'s
  reuse of old dumps. It gives the same results and is kept for differential testing.
- **No temporary coverage files.** An H1 (two BAM replicates) or H3 (BAM with control) run writes 0.1–0.7 GB
  instead of about 37 GB.
- **Memory.** Lower on replicate and control inputs, slightly higher on some others (Costs, below).

The reuse change is its own commit (`5638567`), disjoint from the last one, so it can be dropped on its own if you
prefer to keep the reuse. Keeping it would mean loading old dumps into tracks.

**At the branch heads**, 2026-11-09 (setup and provenance under Timing and memory), after / before wall time:

| Dataset | Threads | Before | After | After / before | Peak PSS | Input parsing | Element definition |
|---|---|---|---|---|---|---|---|
| `H1_bam_reps_karpas422` | 16 | 58.1 s | 54.5 s | 0.939 | 4.90 → 3.30 GB | 6.8 → 3.1 s | 2.8 → 2.6 s |
| `H1_bam_reps_karpas422` | 64 | 34.7 s | 31.7 s | 0.911 | 5.96 → 4.32 GB | 6.6 → 2.8 s | 1.7 → 2.2 s |
| `H3_bam_ct_gm12878_grocap` | 16 | 78.0 s | 70.5 s | 0.904 | 5.17 → 3.94 GB | 4.1 → 3.3 s | 3.1 → 3.0 s |
| `H3_bam_ct_gm12878_grocap` | 64 | 49.7 s | 42.2 s | 0.849 | 5.16 → 5.16 GB | 3.3 → 2.3 s | 1.9 → 2.3 s |
| `T4_bam_reps_hct116` | 16 | 32.6 s | 30.4 s | 0.931 | 2.13 → 2.13 GB | 4.7 → 2.7 s | 1.6 → 1.3 s |
| `T4_bam_reps_hct116` | 64 | 19.9 s | 17.2 s | 0.869 | 2.85 → 2.89 GB | 5.0 → 2.6 s | 1.7 → 1.5 s |
| `T5_bam_ct_k562_grocap` | 16 | 13.6 s | 12.5 s | 0.923 | 1.19 → 1.27 GB | 1.7 → 1.4 s | 0.8 → 1.0 s |
| `T5_bam_ct_k562_grocap` | 64 | 9.4 s | 8.6 s | 0.913 | 1.67 → 1.84 GB | 2.0 → 1.5 s | 0.8 → 1.2 s |

**cgroup memory peak** (`memory.peak`, which counts the page cache as well as process memory), B4 → this branch,
measured on 2026-11-09:

| Dataset | 16 threads | 64 threads |
|---|---|---|
| `H1_bam_reps_karpas422` | 38.0 → 3.4 GB | 38.7 → 4.7 GB |
| `H3_bam_ct_gm12878_grocap` | 38.5 → 4.2 GB | 40.8 → 5.4 GB |
| `T4_bam_reps_hct116` | 5.5 → 2.2 GB | 5.8 → 3.1 GB |
| `T5_bam_ct_k562_grocap` | 3.7 → 1.2 GB | 3.8 → 1.8 GB |

On the genome-wide datasets, the base's peak is mostly the page cache of the per-base `.npy` dumps it writes. That
cache can be reclaimed, so a job given less memory does not fail for it, but the kernel has to write back and
evict it under the limit; on a machine without a limit it pushes out other cached data.

### Costs

- **Peak PSS rises by up to +9% at 16 and 64 threads on bigWig and small inputs**: for example H2 3.06 → 3.10 GB
  and T3 2.11 → 2.30 GB at 16 threads (eight-dataset comparison below), and T5 1.67 → 1.84 GB at 64 threads in the
  branch-head run above. It falls on replicate and control inputs (H1, H3, H6), and at 1 thread on every dataset.
- **Element definition is 0.4–0.5 s slower at 64 threads on H1, H3 and T5** in the branch-head run (not traced
  further); it is faster on T4.

### Why

PINTS writes the coverage of every chromosome strand as a dense per-base `.npy` file, at each stage: parsed
input, each replicate and their merge, control and normalised sample. On the base branch's code, an H1 or H3 run
writes 1176 dumps, 36.7–36.8 GB (profile details under Timing and memory).

More than 99.9% of those bases are zero. Measured on the benchmark's inputs, 0.02–0.07% of the bases of a strand
carry signal (section 2 of the design note,
[`bench/design/sparse-coverage.md`](https://github.com/adamyhe/PINTS/blob/bench/external/bench/design/sparse-coverage.md);
the commit bodies' "step N" refer to its steps). Writing the dumps took most of the input-parsing time at 64
threads, and all of it passes through the page cache, which cgroup v2 charges to the job.

### What changes

- **`pints/coverage_store.py` (new).** A `CoverageTrack` holds the non-zero bases of one chromosome strand as
  sorted runs, or as points for BAM 5' counts. It returns exactly what the same slice, sum or maximum of the dense
  array returns, dtype included. Readers use its methods (`window`, `fill`, `clear`, `nonzero`, `sum`, `max`,
  `is_zero`), never the arrays, so the layout can change behind them.
- **Pool workers reach tracks through a module-level registry, by key.** Tasks carry keys; the workers inherit
  the registry when the pool forks, and tracks are read-only after that. A track refuses to be pickled, so no task
  can carry one by mistake.
- **Every stage builds tracks:**
  - BAM parsing, from the contig- and region-parallel paths and the serial parser. Dump and track share one count
    builder, so the dtype rule and the uint32 wrap-around are in one place.
  - bigWig loading, from the intervals (same block queries, checks and element-wise transform as the block-wise
    loader; runs that become 0 are dropped). Also the dense `values()` fallback, the loader for pyBigWig without
    numpy, and the intersection of the two strands' chromosomes with its read-count adjustment.
  - The replicate merge sums the replicates on the union of their runs in int32 (wrapping, as the in-place int32
    additions do), then casts to the dtype chosen from the maximum. Chromosome order, zero tracks and read counts
    are the same.
  - Control normalisation evaluates `max(exp - input * sf, 0)` only where the sample has signal. For a finite
    `sf >= 0` and non-negative inputs, the dense result is 0 everywhere else, in every integer dtype, and the
    dtype chosen from the maximum cannot change (the argument is in `normalize_track`'s docstring). Other scale
    factors use the dense expression.
- **The pipeline passes registry keys where it passed dump paths.** Element definition reads block windows.
  Phases 1 and 2 fill a full-length zero buffer over the span they read (the subpeaks; the context plus
  `check_window`'s flanks and the 2.5 kb of the local environment). The buffers come from their own anonymous
  mappings: from the C heap, a freed buffer's pages stayed in the worker, about 150 MB per worker on T4.
- **Control parsing and normalisation run before the pool forks**, so that the workers see the normalised tracks.
  Element definition, which runs on the un-normalised sample, keeps reading those tracks.
- **Releasing coverage drops tracks**; nothing is left to delete.
- The fractional-value warning of the bigWig loaders is also logged by the store's loaders (still once per
  chromosome, as before).
- **Each element-definition task gets only its chromosome's entries** (last commit, `0e60047`). Before, every
  task received every sample's whole chromosome → coverage dict, although it reads one entry. There are two tasks
  per contig, so the pickled arguments grow with the square of the contig count. The registry keys are small, but
  T4's BAM header lists 2066 contigs: pickling 4132 tasks, each carrying dicts of 2066 keys that pickle through a
  Python-level `__reduce__`, took about 4.3 s at 64 threads, against about 1.2 s with `.npy` paths in the dicts,
  while the tasks' own work was the same (8.4 s against 8.6 s, summed over workers).

Public functions are unchanged: `get_read_signal` and `get_coverage_bw` still write dumps and return their paths,
and `pints_boundary_extender`, which calls them, still works from files.

### Risks

- **Memory grows with the number of distinct coverage runs**, not with genome length. The benchmark inputs carry
  signal on 0.02–0.07% of bases. Each run stores a start, an end and a value, so a dense input whose value
  changes at most bases approaches per-base storage, can exceed the dense array, and is held in memory rather
  than on disk.
- **The parent holds every track across the fork.** Workers read them through copy-on-write pages inherited from
  the parent, so the tracks count once, but they stay resident in the parent for the whole run.
- **Normalisation against the control no longer runs in parallel with the store.** Control parsing still uses
  its own short-lived pools, but `normalize_tracks_using_input` runs in the main process, before the shared pool
  exists; the dump path (`PINTS_COVERAGE_STORE=0`) still passes `threads` to `normalize_using_input`. On H3 the
  step is faster anyway (8.0 → 1.0 s at 16 threads, 7.7 → 1.1 s at 64), because it only visits the sample's
  non-zero runs; with dense inputs it would be serial work.

### Timing and memory

The benchmark kit is [`bench/`](https://github.com/adamyhe/PINTS/tree/bench/external/bench) on the fork's
`bench/external` branch: [`bench/variants.json`](https://github.com/adamyhe/PINTS/blob/bench/external/bench/variants.json)
pins each measured variant to a commit, and
[`bench/datasets.json`](https://github.com/adamyhe/PINTS/blob/bench/external/bench/datasets.json) describes the
datasets. PINTS runs on one thread unless `--thread` is given; the 16- and 64-thread rows use `--thread`.

**Branch-head table** (top of this description): `bench/run ab --fresh`, B4 against this branch, on an AMD EPYC
9554 server (64 CPUs) on 2026-11-09, three alternating repetitions per cell, median wall time, identical output in
every cell. Another user's single-core jobs ran on the server's other cores throughout, which makes both sides
slower than on an idle server (B4 on H1 at 16 threads: 58.1 s here, 46.2 s in B4's idle-server run); the ratios
compare runs made under the same load. On H3, control parsing and normalisation falls from 8.0 to 1.0 s
(16 threads) and 7.7 to 1.1 s (64 threads).

**cgroup table**: B4 against this branch's code. Each run ran in its own Slurm step, so it had its own cgroup. The
inputs were read in an earlier step, so their page cache was charged there. Two runs each on the same server; the
two runs agreed within 0.2 GB.

<details><summary>Eight datasets at 1, 16 and 64 threads, phases and temporary files (earlier commits)</summary>

`bench/run` paired suite on the same server, in two sessions on 2026-10-31. In each, a combined commit of B2 and
B3 (`f8110c7`) and one other commit ran alternately: three pairs at 16 and 64 threads, one at 1 thread; medians.
The two commits compared below are:

- `6d2056c`, an earlier version of B4: B4's code without A1's fixes (they add no work on these runs apart from the
  bigWig warning's check of the loaded values).
- `adada92`, an earlier version of this branch: this branch's code without A1's fixes, without the bigWig warning
  in the store's loaders, without the no-reuse commit (which does nothing in a fresh output directory) and without
  the last commit, **plus** an earlier version of B6's first commit (BED sort, merge and concatenation without
  bedtools). That commit changes only the filtering and pairing steps, so the phases below are not affected by it.

`f8110c7`'s wall times in the two sessions are within 4% of each other.

**Temporary files** (profiles of the same two commits, 1 and 64 threads):

| Dataset | Base | This branch |
|---|---|---|
| `H1_bam_reps_karpas422` | 36.8–37.4 GB written, 1176 `.npy` | 0.1–0.5 GB, no `.npy` |
| `H3_bam_ct_gm12878_grocap` | 37.0–37.8 GB written, 1176 `.npy` | 0.2–0.7 GB, no `.npy` |

**Peak PSS** (summed over all processes), base → this branch:

| Dataset | 1 thread | 16 threads | 64 threads |
|---|---|---|---|
| `H1_bam_reps_karpas422` | 2.29 → 1.86 GB | 4.81 → 3.34 GB | 5.44 → 4.39 GB |
| `H2_bw_caco2` | 2.48 → 1.94 GB | 3.06 → 3.10 GB | 4.11 → 4.45 GB |
| `H3_bam_ct_gm12878_grocap` | 2.68 → 2.14 GB | 4.60 → 3.92 GB | 5.32 → 5.13 GB |
| `H6_bam_reps_bw_a549r2_karpas` | 2.05 → 1.72 GB | 4.67 → 2.94 GB | 5.14 → 4.09 GB |
| `T2_bam_k562` | 0.68 → 0.53 GB | 0.97 → 0.97 GB | 1.25 → 1.21 GB |
| `T3_bw_gw_a549` | 1.64 → 1.34 GB | 2.11 → 2.30 GB | 2.89 → 3.13 GB |
| `T4_bam_reps_hct116` | 1.05 → 0.88 GB | 2.08 → 2.16 GB | 2.74 → 2.88 GB |
| `T5_bam_ct_k562_grocap` | 0.63 → 0.51 GB | 1.19 → 1.26 GB | 1.72 → 1.78 GB |

At 1 thread, peak PSS falls by 16–22% on every dataset. At 16 and 64 threads, it falls on the replicate and
control datasets (H1, H3, H6), and changes by at most +9% elsewhere. PSS leaves out the page cache, so it does not
show the 37 GB of dumps that no longer go through it.

**Phases** (seconds, base → this branch):

| Dataset | Threads | Input parsing | Control parsing + normalisation | Element definition | Phase 2 tests |
|---|---|---|---|---|---|
| `H1_bam_reps_karpas422` | 1 | 45.3 → 19.9 | – | 25.0 → 23.7 | 179.8 → 178.7 |
| `H1_bam_reps_karpas422` | 64 | 6.7 → 2.8 | – | 1.6 → 2.1 | 10.1 → 9.8 |
| `H3_bam_ct_gm12878_grocap` | 1 | 32.8 → 26.5 | 50.5 → 8.0 | 28.5 → 26.6 | 218.0 → 219.5 |
| `H3_bam_ct_gm12878_grocap` | 64 | 3.3 → 2.2 | 6.9 → 1.1 | 1.8 → 2.2 | 13.1 → 12.9 |
| `T3_bw_gw_a549` | 1 | 9.3 → 1.7 | – | 19.4 → 18.1 | 85.8 → 85.3 |
| `T4_bam_reps_hct116` | 16 | 4.7 → 2.6 | – | 1.3 → 4.3 | 14.0 → 14.5 |
| `T4_bam_reps_hct116` | 64 | 5.0 → 2.7 | – | 1.4 → 4.6 | 6.1 → 6.0 |

Input parsing and control normalisation get faster everywhere. In these measurements element definition got
slower at 64 threads: by 0.4–0.5 s on H1 and H3, and by about 3 s on T4. The last commit removes the T4 slowdown
(the rows above predate it); the 0.4–0.5 s on H1 and H3 remains (see the branch-head table).

Wall time relative to `f8110c7` in the same session, at 64 threads, goes from 0.62 to 0.55 on H1 and from 0.63 to
0.54 on H3; these totals include the bedtools-free commit.

</details>

<details><summary>The last commit on its own</summary>

`bench/run ab` on 2026-11-08, an earlier version of this branch without the last commit (`6774d93`) against one
with it (`80a791b`, the same code as `769beba` below), on the same server, three alternating repetitions per
cell, median wall time, identical output in every cell:

| Dataset | Threads | Before | After | After / before | Element definition |
|---|---|---|---|---|---|
| `T4_bam_reps_hct116` | 16 | 27.1 s | 23.9 s | 0.885 | 4.35 → 1.28 s |
| `T4_bam_reps_hct116` | 64 | 18.4 s | 15.4 s | 0.835 | 4.61 → 1.47 s |
| `H1_bam_reps_karpas422`, `H3_bam_ct_gm12878_grocap`, `T5_bam_ct_k562_grocap` | 16, 64 | | | 0.999–1.030 | unchanged |

</details>

<details><summary>Provenance of the measured commits</summary>

- Branch-head run and cgroup table: B4 at `0333d58` (now `03a9b20`) against this branch at `769beba` (now
  `0e60047`), which includes the last commit and none of B6. Each pair has the same code apart from comments,
  test names and B3's phase-timer fix, which affects only stderr when `PINTS_PHASE_TIMING` is set.
- `f8110c7`, `6d2056c`, `adada92` and `6774d93` live under `refs/archive/*` on the fork
  (github.com/adamyhe/PINTS).
- Benchmark equivalence of `adada92`: all 134 runs (19 datasets; default, `--disable-small` and sensitive
  configurations; 1 to 64 threads) are byte-identical to `main`.
- The `strace` check below ran on an earlier head of B10 (`3ac124d`), whose coverage code is this branch's.

</details>

### Testing

- `tests/unit/test_coverage_store.py`: `window`, `fill` and `clear` against the dense array, including empty
  tracks, negative bounds and windows past the chromosome ends; sums and maxima; `track_to_npy` writes
  `np.save`'s bytes; registry keys pickle small and tracks refuse to; float tracks are rejected.
- `tests/diff/test_differential_store.py` (Hypothesis): tracks equal `np.load` of the dumps, values and dtype,
  for BAM (library layouts, reverse complement, chromosome filters, fast and serial parsers, region splits) and
  bigWig (dtype limits, block sizes, chromosomes missing from one strand, the dense fallback, pyBigWig without
  numpy). The track merge and normalisation equal the file-based ones, including scale factors of 0, -0, NaN,
  infinity and below 0. Stale dumps with the same prefix in the output directory are neither read nor changed.
- `tests/test_bugfixes.py`: the store's bigWig loaders log the non-integer warning.
- Full suite: 265 passed.
- chr22 test data, with bigWig, BAM, and `--annotation-gtf --focused-chrom chr22` input, at 1 and 4 threads: the
  bidirectional, divergent and unidirectional peaks are byte-identical to `main`.
- `strace` of the top of the stack (B10, whose coverage code is this branch's; commit under provenance) on chr22
  BAM, bigWig, BAM with control, and two BAM replicates, at 1 and 4 threads: no `.npy` file is opened for reading
  or writing. With `PINTS_COVERAGE_STORE=0`, the BAM run at 1 thread writes 170 dumps.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
