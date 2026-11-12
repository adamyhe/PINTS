## Parallel engine: BAM and bigWig parsing, chunked window tests, one worker pool, bounded memory

Stacked on B2 (`pr/03a-gate1-kernels`); review only the last 17 commits (none of them is a merge).

**The output is byte-identical to the base branch.** This PR changes how the work is split across processes and how
coverage is held in memory, not what is computed. Each fast path checks that its inputs are within the cases it
reproduces exactly, and runs the original code otherwise. See Testing for what was checked at which commit.

**Peak memory rises with `--thread` on some inputs.** On a genome-wide BAM with two replicates (`H1`, below), peak
PSS goes from 12.3 GB in the release to 25.6 GB at 16 threads (B2 and B3 together; see Timing), and it grows with
the number of workers. B4 (I/O and memory) brings it down: in B4's paired runs, `H1` at 16 threads falls from 20.67
to 4.81 GB.

### What changes for users

- No new command-line options. `--thread` now also parallelises parsing, element definition, normalisation,
  replicate merging, filtering and pairing (see below); with `--thread 1` everything runs in the main process.
- Peak memory: see above and Timing.
- A new environment variable for diagnostics, `PINTS_PHASE_TIMING=1` (documented in the README; see Diagnostics).
- With an input control and `--thread` > 1, up to about twice `--thread` worker processes exist at once (see
  Risks).
- The non-integer bigWig warning from A1 is still logged once per chromosome.

### Why

In the current release, only the window tests run in parallel, with one task per chromosome and strand, and the two
strands one after the other. BAM and bigWig parsing, element definition, input normalisation, replicate merging,
filtering and pairing all run in the main process. So a run hardly gets faster beyond a few threads: on the
benchmark datasets below, the release is at most 20% faster at 64 threads than at 16. The largest chromosome sets the
length of the parallel step, and the serial steps account for the rest.

### What changes

**Parsing (`io_engine`)**
- BAM files with an index are parsed one contig per task (`parse_bam_by_contig`). Contigs with no reads in the
  index statistics get an all-zero track without being read. If there is no index, or the index statistics do not
  match the records fetched, the serial parser runs as before.
- bigWig files are loaded in parallel, one chromosome per task, and each track is built from `bw.intervals()`
  instead of a dense float64 `bw.values()` array. Overlapping or out-of-range intervals fall back to the dense
  path. The coverage dumps and read counts are byte-identical.
- Replicate merging and input normalisation run on a pool, chromosome by chromosome, in the same order and with
  the same read-count accumulation as before.

**Element definition**
- Window sums and non-zero regions are computed with vectorised NumPy (`window_sums`, `nonzero_window_regions`)
  instead of a Python loop over windows. Tests compare them with `sliding_window` and with the bedtools merge they
  replace.
- Chromosomes are defined in parallel on the pool, costliest (largest coverage dump) first.

**Window tests (`evaluate_strand_jobs`)**
- All samples and both strands are evaluated in one scheduled pass, in two phases:
  - phase 1 computes the per-subpeak statistics in blocks of subpeaks, then the per-chromosome thresholds;
  - phase 2 runs the local-environment tests in chunks.
- A chromosome can be split for phase 2 because the LER cache is keyed by subpeaks within 5 kb (half the largest
  flank) of the tested peak, so two runs of subpeaks more than 10 kb apart never share a key. Chunks are cut only at
  gaps of at least 15 kb. Each chunk is evaluated with a fresh cache, and the results are joined in file order.
- Each chunk reads its subpeaks from an in-memory index (`SubpeakIndex`) that returns what `pysam.TabixFile.fetch`
  returns. Files the index does not support (tied or unsorted starts, non-positive lengths) use tabix, unchunked.
- The likelihood-ratio test is skipped outside `--sensitive`, since only `--sensitive` uses its p-value.

**Filtering and pairing**
- Stratified filtering of every sample and strand runs on the pool.
- With several samples (`--dont-merge-reps`), each sample is paired on the pool.
- `merge_opposite_peaks` reads the opposite-strand candidate file once into an in-memory overlap index
  (`PeakBedIndex`), instead of making one tabix query per significant peak. The index returns the same fields in the
  same order as tabix, and raises the same errors for bad regions and unknown chromosomes. Files outside the subset it
  supports (header lines, CR line endings, unsorted starts, non-contiguous chromosomes) use tabix. Tabix is also
  kept when there are few queries relative to the file size, where building the index costs more than it saves.

**One pool and bounded memory**
- `peak_calling` creates one worker pool of `--thread` processes and keeps it from element definition to the end
  of the run. BAM parsing, bigWig loading, replicate merging and normalisation use short-lived pools of their own.
  With `--thread 1` there are no pools and everything runs in the main process. `_make_pool` uses the `fork` start
  method where it exists. The pool initializer sets the statistical model in each worker, so the workers do not rely
  on fork to inherit it. The Numba kernel from B2 is compiled before the fork (`numba_warmup`), so the workers do not
  compile it again.
- After phase 1, each subpeak bed is compressed and indexed (bgzip, tabix) on the pool asynchronously while the
  main process goes on to phase 2. The compressed file is read only by the tabix fallback, which waits for the
  pending indexing first, and by housekeeping; all pending indexing is awaited after phase 2.
- Element definition and input normalisation work on blocks of 2 Mb instead of whole chromosomes, and the bigWig
  loader writes each track in blocks of 4 Mi bases:
  - window sums never straddle a block, and a region that touches a block end is carried into the next block, so
    regions, subpeaks and summits are those of the whole-array code;
  - at 2 Mb, the float64 temporaries (16 MB) stay below glibc's 32 MiB mmap threshold and are reused from the heap
    (8 Mb blocks were mmap'd and unmapped each time);
  - freed buffers are returned to the OS with `malloc_trim` (glibc only).
- Element definition reads the coverage dumps with plain file reads instead of resident memory maps. Inputs outside
  the contract of the block-wise code (non-integer, object dtype, mismatched lengths) use the original code.

**Diagnostics**
- With `PINTS_PHASE_TIMING=1`, `peak_calling` writes one `[phase]` line per stage to stderr. The clock starts when
  `peak_calling` starts, and each line gives the wall time since the previous line. The labels are stage names:
  `input parsing`, `element definition`, `subpeak bookkeeping`, `control parsing+normalisation`, the window-test
  stages `phase1 setup`, `phase1 stats`, `phase1 thresholds`, `phase2 setup` and `phase2 tests` (once per pass; at
  one thread there is one pass per sample and strand), then `input normalisation/eval` and `filtering+pairing+qc`.
  Without the variable nothing is printed. It never changes the output. The README documents it as a debugging
  option.

**Fixes to earlier commits**
- `876108a` passes A1's non-integer bigWig warning and `--peak-width` on to the block-wise bigWig loader and
  element definition. The default path uses those since the earlier commits of this PR, so without it the warning
  and the option did not reach it. Output without `--peak-width` is unchanged.
- `c23eaf5` fixes a regression introduced earlier in this PR by `947a1ad` (per-sample pairing on the pool). With
  several samples (`--dont-merge-reps`) and `--thread` > 1, that path read `kwargs["stringent_only"]`, but the
  command line stores the option as `stringent_pairs_only`, so `--stringent-pairs-only` was ignored there. It now
  uses `_stringent_only(kwargs)`, like the other pairing calls.

**Review order** (`git log pr/03a-gate1-kernels..pr/03b-gate1-engine`), grouped by topic, oldest first in each:

- Parsing: `1f99ef3` BAM by contig; `225dc04` parallel bigWig loading.
- Element definition: `986bdcb` vectorised window sums (and the LRT skip outside `--sensitive`).
- Window tests: `36cde41` chunked tests and `SubpeakIndex`; `522238d` one persistent pool, all samples and strands
  in one pass, parallel element definition.
- Filtering and pairing: `7415cce` filtering on the pool (and bigWig tracks from intervals); `7b79977`
  `PeakBedIndex`; `fc93ee5` tabix kept for few queries.
- Pool and memory: `7b31929` subpeak compression, replicate merging and normalisation on the pool; `947a1ad` Numba
  warm-up and per-sample pairing on the pool; `e78f6ea` block-wise element definition, bigWig loading and
  normalisation; `8c03e5c` 2 Mb blocks.
- Diagnostics: `76ad22f` `PINTS_PHASE_TIMING`.
- Tests: `a3c48ab` parallel engine; `ad38a0b` concatenation dtypes and `PeakBedIndex`.
- Fixes: `876108a` bigWig warning and `--peak-width`; `c23eaf5` `--stringent-pairs-only`.

### Risks and limitations

- **Peak memory** rises with `--thread` on some inputs (see the lead and Timing). At one thread it is lower than the
  release's on every dataset measured except `H2` and `T3`, where it is about the same.
- **Process count.** The persistent pool is created before the control is parsed and the samples are normalised.
  Control BAM parsing and `normalize_using_input(..., threads=...)` open their own short-lived pools while the
  persistent pool's workers sit idle, so with an input control (`--ct-bam` or control bigWigs) and `--thread` > 1 up
  to about 2 × `--thread` worker processes exist at once. B5 (sparse coverage) moves control parsing and
  normalisation before the pool is created (with its default in-memory coverage store).
- **Platform.** `_make_pool` uses `fork` where available and falls back to the default start method, but the
  parallel bigWig loader, input normalisation and replicate merging call `multiprocessing.get_context("fork")` with
  no fallback. The parallel paths therefore assume a platform with `fork` (Linux, macOS); with `--thread` > 1 they
  would raise on Windows. All testing was on Linux. `malloc_trim` is used only where glibc provides it.

### Timing

The benchmark kit is [`bench/`](https://github.com/adamyhe/PINTS/tree/bench/external/bench) on the fork's
`bench/external` branch: [`bench/variants.json`](https://github.com/adamyhe/PINTS/blob/bench/external/bench/variants.json)
pins each measured variant to a commit, and
[`bench/datasets.json`](https://github.com/adamyhe/PINTS/blob/bench/external/bench/datasets.json) describes the
datasets. PINTS runs on one thread unless `--thread` is given. All runs were on an AMD EPYC 9554 server.

**This PR on its own**, against B2 (`pr/03a-gate1-kernels`), three alternating runs each, median wall time:

| Dataset | 16 threads | 64 threads |
|---|---|---|
| `H3_bam_ct_gm12878_grocap` | 455.5 → 86.9 s (0.191) | 456.7 → 69.4 s (0.152) |
| `T4_bam_reps_hct116` | 446.1 → 35.6 s (0.080) | 449.1 → 26.8 s (0.060) |
| `T5_bam_ct_k562_grocap` | 291.3 → 16.1 s (0.055) | 291.1 → 13.6 s (0.047) |

Outputs are identical. B2 (exact kernels) alone gives 0.56–0.74 at 16 threads on the same datasets (see its
description).

**B2 and B3 together against upstream `main`** (`81b2dd9`), measured on a combined commit of B2 and B3 (`f8110c7`):

| Dataset | Input | 1 thread | 16 threads | 64 threads | Peak PSS, 16 threads | Peak PSS, 64 threads |
|---|---|---|---|---|---|---|
| `H1_bam_reps_karpas422` | PRO-cap BAM, 2 replicates, genome-wide | 2975 → 530 s (5.6×) | 612 → 72.8 s (8.4×) | 551 → 54.7 s (10.1×) | 12.3 → 25.6 GB | 12.3 → 24.0 GB |
| `H2_bw_caco2` | PRO-cap bigWig, genome-wide | 3518 → 580 s (6.1×) | 566 → 67.6 s (8.4×) | 569 → 49.2 s (11.6×) | 6.0 → 5.5 GB | 7.1 → 8.5 GB |
| `H3_bam_ct_gm12878_grocap` | GRO-cap BAM + input, genome-wide | 2674 → 663 s (4.0×) | 717 → 90.7 s (7.9×) | 719 → 73.5 s (9.8×) | 12.6 → 6.5 GB | 12.6 → 9.6 GB |
| `H6_bam_reps_bw_a549r2_karpas` | PRO-cap BAM, genome-wide | 2129 → 382 s (5.6×) | 563 → 57.3 s (9.8×) | 450 → 40.1 s (11.2×) | 12.1 → 5.1 GB | 12.1 → 7.8 GB |
| `T2_bam_k562` | PRO-cap BAM, chr1, chr21, chrX | 321 → 40.4 s (7.9×) | 290 → 9.3 s (31×) | 285 → 8.4 s (34×) | 6.7 → 1.6 GB | 6.6 → 4.1 GB |
| `T3_bw_gw_a549` | PRO-cap bigWig, genome-wide | 1667 → 281 s (5.9×) | 389 → 34.7 s (11.2×) | 391 → 26.4 s (14.8×) | 5.7 → 3.9 GB | 6.9 → 6.6 GB |
| `T4_bam_reps_hct116` | PRO-cap BAM, 2 replicates, chr17–22 | 1886 → 269 s (7.0×) | 800 → 35.8 s (22×) | 807 → 27.4 s (29×) | 7.5 → 2.9 GB | 7.5 → 3.5 GB |
| `T5_bam_ct_k562_grocap` | GRO-cap BAM + input, chr19–22 | 583 → 90.0 s (6.5×) | 395 → 16.1 s (25×) | 395 → 13.7 s (29×) | 6.9 → 2.0 GB | 6.9 → 4.3 GB |

Median wall time of three runs (one run for the reference at 1 thread). Peak PSS is summed over all PINTS processes.
The reference does not have A1's empty-frame dtype fix, which alone accounts for part of the H3 gain (717 → 625 s at
16 threads, measured on H3 only).

Peak memory rises on `H1` at 16 and 64 threads (about twice the release's), and on `H2` at 64 threads. At one
thread it is lower everywhere except `H2` and `T3`, where it is about the same (+2% and +4%). On `H1` it grows with
the number of workers (3.5, 9.1, 15.6 and 25.6 GB at 1, 4, 8 and 16 threads). B4 (I/O and memory) works on this.

<details><summary>Provenance of these measurements</summary>

- The per-PR table was measured on an earlier version of this branch (`06d7b74`), whose `pints/` differs from
  `c23eaf5` only by the phase-timer fix (it affects only stderr with `PINTS_PHASE_TIMING` set), against an earlier
  version of B2 (`f45b706`) with the same `pints/`.
- The combined commit `f8110c7` (archived) is the head of the branch from which B2 and B3 were split. Its `pints/`
  code equals this branch's except for:
  - A1's changes, which it does not have: the `get_elbow` arithmetic, option passing, `Poisson.lrt`, one-strand
    pairing and the bigWig warning. Of these, only the bigWig warning runs on these datasets with default
    options, and it is one comparison per interval value;
  - where `_concat_result_frames` sits and its docstring;
  - setup.py formatting.
- The combined runs used the `performance` governor, a 64-CPU allocation and node-local disk. The reference ran on
  2026-10-16/17 and `f8110c7` on 2026-10-19, not interleaved. Other users' jobs were on the node: the load recorded
  before each run was about 150 for the 16- and 64-thread H6 runs of `f8110c7`, and about 186 for the 16- and
  64-thread T3 runs of the reference. The empty-frame fix alone was timed as its own commit (`de0febe`).
- Environment: Python 3.12.13, NumPy 2.1.0, pandas 2.3.1, SciPy 1.15.0, Numba 0.61.2. The `T` datasets were used
  while tuning; the `H` datasets were held out.

</details>

### Testing

At the branch head (`c23eaf5`):

- New differential tests (Hypothesis where inputs are random) compare the working tree with `pints_ref` bit for bit:
  - `tests/diff/test_differential_parallel.py`: vectorised window sums and regions; contig-parallel BAM parsing and
    its fallback without an index; all-zero tracks; parallel bigWig loading, replicate merging and normalisation;
  - `tests/diff/test_differential_engine.py`: `SubpeakIndex` against tabix; forced chunking on a pool (with checks
    that it really splits); several samples and strands in one pass; the tabix fallback;
  - `tests/diff/test_differential_memory.py`: block-wise element definition against the whole-array code and the
    reference, including int32 wrap-around and the fallbacks; block-wise bigWig loading and normalisation; the
    `.npy` headers written by hand;
  - `tests/diff/test_differential_pairing.py`: `PeakBedIndex` against `pysam.TabixFile`, its fallback outside the
    supported subset, the choice between the two, and `merge_opposite_peaks` output against `pints_ref`;
  - `tests/diff/test_differential_extra.py`: multi-chromosome `peaks_single_strand`, on a pool and without one;
    pairing; `cut_peaks_dry_run`.
- `tests/unit/test_concat_bins_frames.py`: the phase-2 result frames keep their dtypes when some are empty.
- `tests/test_bugfixes.py` gains two tests: the parallel bigWig loader warns about non-integer values, and
  `--peak-width` reaches the block-wise element definition.
- Full suite (`pytest tests`, Python 3.12.13, NumPy 2.1.0, pandas 2.3.1, SciPy 1.15.0, Numba 0.61.2): 211 passed.
- chr22 test data, with bigWig, BAM, and `--annotation-gtf --focused-chrom chr22` input, at 1 and 4 threads: the
  bidirectional, divergent and unidirectional peaks are byte-identical to `main`.
- The per-PR timing runs above: outputs identical to B2's.

The wider thread-count check was run on the combined commit `f8110c7`, not on this head:

<details><summary>Equivalence suite on <code>f8110c7</code></summary>

626 runs over 19 datasets, three configurations (default, `--disable-small`, `--sensitive`), 1 to 64 threads and
`PYTHONHASHSEED` 0, 99 and 1234. The bidirectional, divergent and unidirectional peaks of every run are
byte-identical to upstream `main`. `f8110c7` differs from this branch as listed under Timing.

</details>

🤖 Generated with [Claude Code](https://claude.com/claude-code)
