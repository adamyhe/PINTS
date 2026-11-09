## Parallel engine: BAM and bigWig parsing, chunked window tests, one worker pool, bounded memory

Stacked on the statistical-kernels PR (`pr/03a-gate1-kernels`); review only the last 17 commits (none of them is a
merge).

**The output is byte-identical to the base branch, at every thread count tested (1 to 64).** This PR changes how
the work is split across processes and how coverage is held in memory, not what is computed. Each fast path checks that its inputs are
within the cases it reproduces exactly, and runs the original code otherwise.

### Why

In the current release, only the window tests run in parallel, with one task per chromosome and strand, and the two
strands one after the other. BAM and bigWig parsing, element definition, input normalisation, replicate merging,
filtering and pairing all run in the main process. So a run hardly gets faster beyond a few threads: on the
benchmark datasets below, the release is at most 20% faster at 64 threads than at 16. The largest chromosome sets the
length of the parallel step, and the serial steps account for the rest.

### What changes

The commits fall into six groups, plus a one-line fix at the end.

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
  With `--thread 1` there are no pools and everything runs in the main process. The pool uses the `fork` start method where it exists.
  The pool initializer sets the statistical model in each worker, so the workers do not rely on fork to inherit it.
  The Numba kernel from the previous PR is compiled before the fork, so the workers do not compile it again.
- Element definition and input normalisation work on blocks of 2 Mb instead of whole chromosomes, and the bigWig
  loader writes each track in blocks of 4 Mi bases:
  - window sums never straddle a block, and a region that touches a block end is carried into the next block, so
    regions, subpeaks and summits are those of the whole-array code;
  - at 2 Mb, the float64 temporaries (16 MB) stay below glibc's 32 MiB mmap threshold and are reused from the heap (8 Mb blocks
    were mmap'd and unmapped each time);
  - freed buffers are returned to the OS with `malloc_trim` (glibc only).
- Element definition reads the coverage dumps with plain file reads instead of resident memory maps. Inputs outside
  the contract of the block-wise code (non-integer, object dtype, mismatched lengths) use the original code.

**Diagnostics**
- `PINTS_PHASE_TIMING=1` makes `peak_calling` write the wall time of each stage to stderr as `[phase]` lines (input
  parsing, element definition, the two testing phases, filtering and pairing). Without it, nothing is printed. It
  never changes the output. The README documents it as a debugging option.

The last commit passes the bugfix PR's bigWig warning and `--peak-width` on to the new block-wise loader and element
definition, which the default path now uses.

The fix commit at the end passes `--stringent-pairs-only` to the pooled per-sample pairing. With several samples
(`--dont-merge-reps`) and `--thread` > 1, that path read `kwargs["stringent_only"]` directly, so the option from the
command line was ignored there. It now uses `_stringent_only(kwargs)`, like the other pairing calls.

### Timing

Measured with the benchmark kit (`bench/`) on cbsubscb19 (AMD EPYC 9554, `performance` governor, 64-CPU allocation,
node-local disk). The measured commit is `f8110c7` (variant `gate1`), the head of the Gate 1 branch before it was split into
B2 and B3, so the numbers cover both PRs. Its `pints/` code equals this branch's except for:
- the bugfix PR's changes, which it does not have: the `get_elbow` arithmetic, option passing, `Poisson.lrt`,
  one-strand pairing and the bigWig warning. Of these, only the bigWig warning runs on these datasets with default
  options, and it is one comparison per interval value;
- where `_concat_result_frames` sits and its docstring;
- setup.py formatting.

The reference is upstream `main` (`81b2dd9`), without the empty-frame dtype fix (A1). That fix alone accounts for
part of the H3 gain (717 → 625 s at 16 threads, variant `pr4`); it was timed on H3 only.

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

Median wall time of three runs (one run for the reference at 1 thread). The reference ran on 2026-10-16/17 and this
code on 2026-10-19, not interleaved. Other users' jobs were on the node: the load recorded before each run was about
150 for the 16- and 64-thread H6 runs of this code, and about 186 for the 16- and 64-thread T3 runs of the reference.
Peak PSS is summed over all PINTS processes. Environment: Python 3.12.13, NumPy 2.1.0, pandas 2.3.1, SciPy 1.15.0,
Numba 0.61.2. The `T` datasets were used while tuning; the `H` datasets were held out.

Peak memory rises on `H1` at 16 and 64 threads (about twice the release's), and on `H2` at 64 threads. At one
thread it is lower everywhere except `H2` and `T3`, where it is about the same (+2% and +4%). On `H1` it grows with
the number of workers (3.5, 9.1, 15.6 and 25.6 GB at 1, 4, 8 and 16 threads). The next PR (B4, I/O and memory)
works on this.

**This PR on its own**, against `pr/03a-gate1-kernels`, on cbsubscb19, three alternating runs each, median wall
time:

| Dataset | 16 threads | 64 threads |
|---|---|---|
| `H3_bam_ct_gm12878_grocap` | 455.5 → 86.9 s (0.191) | 456.7 → 69.4 s (0.152) |
| `T4_bam_reps_hct116` | 446.1 → 35.6 s (0.080) | 449.1 → 26.8 s (0.060) |
| `T5_bam_ct_k562_grocap` | 291.3 → 16.1 s (0.055) | 291.1 → 13.6 s (0.047) |

Outputs are identical. The kernels PR (B2) alone gives 0.56–0.74 at 16 threads on the same datasets (see its
description); the table above with the reference covers both.

### Testing

- New differential tests (Hypothesis where inputs are random) compare the working tree with `pints_ref` bit for bit:
  - `tests/diff/test_differential_gate1.py`: vectorised window sums and regions; contig-parallel BAM parsing and
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
- `tests/test_trackb_fixes.py` gains two tests: the parallel bigWig loader warns about non-integer values, and
  `--peak-width` reaches the block-wise element definition.
- Full suite (`pytest tests`, same environment as above): 211 passed.
- chr22 test data, with bigWig, BAM, and `--annotation-gtf --focused-chrom chr22` input, at 1 and 4 threads: the
  bidirectional, divergent and unidirectional peaks are byte-identical to `main`.
- Benchmark equivalence suite, measured commit `f8110c7` (see Timing): 626 runs over 19 datasets, three configurations
  (default, `--disable-small`, `--sensitive`), 1 to 64 threads and `PYTHONHASHSEED` 0, 99 and 1234. The
  bidirectional, divergent and unidirectional peaks of every run are byte-identical to upstream `main`.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
