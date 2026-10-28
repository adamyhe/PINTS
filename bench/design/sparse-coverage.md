# Design: in-memory sparse coverage (no per-base `.npy` dumps)

Status: proposed, 2026-10-28. Nothing here is implemented yet. This note is the handoff for a session working on the
benchmark node (cbsubscb19, work dir `/workdir/ayh8/PINTS`). Line numbers refer to `perf/io-memory` at `0b91493`.

## 1. Context

| branch | head | what it is |
|---|---|---|
| `main` | `81b2dd9` | upstream hyulab/PINTS, the reference |
| `perf/gate1-bit-identical` | `f8110c7` | gate1: faster, byte-identical to `main` |
| `perf/io-memory` | `0b91493` | gate1 + six I/O and memory commits (`0812507`, `e3a1505`, `7801b61`, `15b8c09`, `70c8c1d`, `0b91493`) |
| `bench/external` | this branch | the benchmark kit (`bench/`); fork only, never for upstream |

Ground rules that still apply:
- **Byte identity.** Every change must reproduce `main`'s `*_bidirectional_peaks.bed`, `*_divergent_peaks.bed` and
  `*_unidirectional_peaks.bed` byte for byte, at any thread count and hash seed. The kit's equivalence suite checks
  this.
- **No upstream PRs.** Don't open anything against hyulab/PINTS until the user says so.
- **Python environments:** use `uv`, not conda.
- **Commits:** author is the user. Commits written by Claude end with
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. Only commits re-authored from Biomni's run carry
  `Co-authored-by: Biomni <biomni@phylo.bio>`.
- **Paired rounds.** Each new `io` commit is timed against fresh gate1 runs in a new paired round:
  `paired_round` in `bench/variants.json`; round k uses repetitions 100k, 100k+1, …. Round 3 is `0b91493`.

## 2. Measurements behind this note

Source: round 3 (`bench_results_cbsubscb19.biohpc.cornell.edu_20261027_1436.tar.gz`) and the profile run
(`profile_results_cbsubscb19_20261028_0410.tar.gz`), both io `0b91493`.

**Round 3, io vs gate1 in the same session:**
- **Equivalence:** all 622 compared variant runs identical.
- **Wall time:** −11 to −19% at 16 threads, −15 to −22% at 64, ±1.5% at 1.
- **H1 at 64 threads:** 61.9 → 49.7 s, peak PSS 20.5 → 6.1 GB.

**Temporary output (profile traces, 64 threads):**

| run | written in total | dense coverage `.npy` | everything else |
|---|---|---|---|
| H1 (2 BAM replicates) | 37.3 GB | 36.8 GB (rep 0, rep 1, merged; 6.1 GB per strand each) | < 0.5 GB |
| H3 (BAM + control BAM) | 37.7 GB | 36.8 GB (sample, control, normalised) | < 0.5 GB |

**H1 input parsing at 64 threads (10.8 s):**

| step | time |
|---|---|
| replicate 1 BAM | 1.7 s |
| replicate 2 BAM | 3.5 s, of which counting reads is 0.4 s (268 region workers) and assembling + writing the arrays is up to 2.9 s per chromosome |
| merge | 5.6 s |

After the first ~12 GB, all writers together sustain about 2.5 GB/s (from `wchar` in `procs.tsv`). That looks like
the kernel throttling writes once unflushed data builds up, but it isn't confirmed. `write_bytes` is not reliable
here, because reaped workers are double-counted. Splitting contigs into regions helped T4 (input parsing
8.0 → 5.8 s) and H3 (4.9 → 3.7 s), but not H1, whose time is now in writing.

**Peak PSS by phase (io, round 3):**

| dataset | 1 thread | 64 threads |
|---|---|---|
| H1 BAM reps | 2.3 GB (filtering) | 6.1 GB (input parsing) |
| H3 BAM + ctrl | **5.1 GB (control normalisation)** | 6.1 GB (control parsing + normalisation) |
| H2 bigWig | **3.1 GB (input parsing)** | 6.6 GB (filtering) |
| T3 bigWig | **3.2 GB (input parsing)** | 5.8 GB (filtering) |
| H6 BAM reps | 2.0 GB | 5.7 GB (input parsing) |
| T1 (chr22 only) | 0.3–0.8 GB | 2.2 GB, which is the worker-pool floor (~35 MB per worker) |

PSS leaves out the page cache. The dumps put up to 37 GB through it, which counts against a SLURM/cgroup memory
limit, so real per-job memory is higher than these numbers.

## 3. Where coverage dumps are written and read today

**Writers** (`pints/io_engine.py`):

| function | line | what it writes |
|---|---|---|
| `get_read_signal` | 891 | per-sample BAM dumps. Allocates `np.zeros(LN, uint32)` for both strands of every contig (≈ 25 GB, only virtual on the fast path; the serial `BamParsers` fallback touches all of it). Also has the `.npy` cache (`load_signal_from_cache`). |
| `parse_bam_by_contig` | 245 | the fast BAM path |
| `_bam_region_worker` | 205 | returns sparse `(positions int32, counts uint32)` per strand |
| `_bam_assemble_worker`, `_save_end_counts` | 223, 155 | densify those and `np.save` |
| `_get_coverage_bw` | 385 | 1-thread bigWig path: dense float64 `bw.values()` per chromosome, 2 GB for chr1 |
| `_bw_chrom_worker`, `_bw_blockwise_sparse` | 531, 484 | parallel bigWig path: block-wise from `bw.intervals()`. Uses one dtype for every chromosome, from the header max. |
| `get_coverage_bw` | 605 | drops chromosomes not present in both strands, and subtracts their sums |
| `merge_replicates_bw` | 1611 | iterates `shared_chromosomes`, a set, so the order depends on the hash seed |
| `_merge_rep_chrom_worker`, `_merge_rep_strand` | 1564, 1504 | int32 sum; dtype from the chromosome max |
| `normalize_using_input` | 1110 | the 1-thread `strand_atom` path loads whole chromosomes as float64 (H3's 5.1 GB) |
| `_normalize_chrom_worker` | 1065 | block-wise, only used with more than one thread. Computes `max(exp − input·sf, 0)`, then a dtype from the max. |

**Readers** (`pints/calling_engine.py` unless noted):

| reader | line | how it reads |
|---|---|---|
| `merge_covs` | 510 | whole chromosome, int32 sum across replicates |
| `cut_peaks_dry_run` (`--annotation-gtf`) | 539 | through `merge_covs` |
| `_load_coverage` | 770 | memory-maps the dump, so workers share it through the page cache |
| `_subpeak_block_stats` (phase 1) | 910 | windows |
| `_check_window_block` (phase 2) | 1027 | windows: peak ± flank, at most 10 kb |
| `_blockwise_element_definition` | ~1944 | block-wise scan of the whole chromosome, with blocks aligned to the 100-bp window |
| `unified_element_definition` | 2039 | — |
| `_coverage_cost` | 1316 | file size, used as the scheduling weight |
| `pints/extension_engine.py` (`pints_boundary_extender`) | 76, 95 | calls `get_read_signal` / `get_coverage_bw` itself, then `np.load` |

**Flow in `peak_calling`:**

| step | line |
|---|---|
| `parse_input_files` | 2631, defined at 2448 |
| replicate merge inside it | 2541 |
| pool forked (`_make_pool`, 1147) | 2654 |
| element definition | — |
| control `parse_input_files` | 2677 |
| `normalize_using_input` | 2693 |
| evaluation, filtering, pairing | — |
| `release_coverage_files` | 1803 |
| `housekeeping` | 1851, called at 2900 |

Because control parsing and normalisation happen after the fork, the workers can only see their results through
files.

**Dumps per input type** (≈ 12 GB per stage genome-wide at int16, both strands):

| input | stages written | dense files | benchmark datasets |
|---|---|---|---|
| BAM, one sample | parse | ~12 GB | T1_bam, T2 |
| bigWig, one sample | parse | ~12 GB | T1_bw, T3, H2, H4, H5, H7_* |
| BAM replicates, merged | per-replicate + merged | ~36 GB | H1, T4, H6 |
| BAM + control | sample + control + normalised | ~36 GB | H3, T5, T1_bam_ct |
| bigWig replicates, or bigWig + control | as above | ~36 GB | **none** |
| replicates + control | n replicates + merged + control + normalised | ~48 GB for 2 replicates | **none** |
| `--dont-merge-reps` | one per replicate (+ normalised if there is a control) | 12 GB × replicates × stages | **none** |

## 4. Design

### 4.1 The store

Each coverage dict maps a chromosome to a `CoverageTrack` instead of a path:

```
CoverageTrack(length, dtype, starts int32[], ends int32[], values dtype[])   # sorted, non-overlapping runs
  .window(a, b) -> dense np.ndarray[dtype] of length b - a   # two searchsorted calls + scatter
  .blocks(block) -> dense blocks aligned like NPY_BLOCK       # for element definition
  .sum(), .max(), .is_zero, .nbytes
```

- **BAM tracks** have runs of length 1. They come straight from `_bam_region_worker` output, concatenated in region
  order. Duplicates at region boundaries are summed exactly as `_save_end_counts` does today.
- **bigWig tracks** keep the intervals as runs. Apply the same element-wise transform as `_bw_blockwise_sparse`:
  `nan_to_num`, flip negatives, `astype(header dtype)`. Drop runs that become 0, so `is_zero` and the support
  stay exact.
- **`dtype`** is the dtype the dense dump would have had, chosen by the same rule at each stage. It matters
  downstream, because readers do arithmetic in it.

Expected size: on the chr22 test BAM, 0.11% of bases carry reads. Genome-wide, the store is probably under 1–2 GB
for any kit dataset. This is an estimate; log `nbytes` per sample in the first commit and check it.

### 4.2 The stages

| stage | today | with the store |
|---|---|---|
| BAM parse | sparse region output → densify → `np.save` | keep the region output; skip assembly. Totals and index-statistics checks unchanged. |
| bigWig parse | intervals → dense file | intervals → runs |
| merge replicates | block-wise file pass, int32 sum | union of positions, int32 sum, then the same dtype rule. Keep the set iteration order and the zero-track handling. |
| control normalisation | block-wise float pass | compute only on the sample's support: with `sf ≥ 0` and finite, `max(exp − input·sf, 0)` is 0 wherever `exp` is 0. Same float64 arithmetic per element, same dtype rule. Fall back to the dense path for NaN, infinite or negative `sf`. |
| element definition | block-wise scan of files | `.blocks()` with the same block alignment |
| phase 1/2, `cut_peaks_dry_run`, `merge_covs` | memory-mapped windows or whole arrays | `.window()` or a whole-chromosome densify |
| scheduling (`_coverage_cost`) | file size | `nbytes`, or the dense size, to keep the current task order |

### 4.3 Process layout

- **Fork order.** Move control parsing and normalisation before `_make_pool`. Workers then inherit every store
  copy-on-write, with no pickling. Only the log order changes.
- **Pools inside parsing** (`parse_bam_by_contig`, the bigWig loader) already return small sparse results.
- **Tracks are read-only after the fork**, so shared pages stay shared; numpy's buffers are separate from the object
  headers that refcounts touch.

### 4.4 What stays file-based

- **Public functions keep returning file paths.** That covers `get_read_signal` and `get_coverage_bw`, which are used
  by `pints_boundary_extender` and possibly by outside users. The store comes from new internal functions.
- **The serial `BamParsers` fallback**, used when the index statistics don't match, builds tracks from its arrays
  instead of writing them. Better: count with `_count_read_ends` per contig, to avoid touching 25 GB.
- **Readers accept either a path or a track** during the migration, so each commit can switch one stage.

### 4.5 Open decision for the user: the `.npy` cache

`get_read_signal` and `_get_coverage_bw` reuse `.npy` files they find in the output directory ("Loading cache from
previous result"). `housekeeping` deletes them after a successful run, so in practice the cache only helps after a
crash. It can also silently reuse stale files from a different input that used the same prefix. With no dumps, the
cache is never written. Choose one:
- (a) drop it, which changes upstream behaviour and needs a note in the PR; or
- (b) keep reading old dumps if present, loading them into tracks.

Ask before implementing.

### 4.6 Expected effect (to be measured)

- **Files:** no coverage dumps, so 12–48 GB less written per run, and no deletion work in `housekeeping` or the
  background deleter.
- **Memory:** no dense per-worker buffers, and no page-cache load. The 64-thread peak should move toward pool floor
  plus store, roughly 3–4 GB.
- **Time:** H1's input phase (10.8 s) and H3's control parsing + normalisation (~9.5 s) should each drop to a few
  seconds.
- **The cost:** filling in windows replaces memory-mapped reads in phases 1 and 2. Watch those phases in the paired
  round.

## 5. Quick wins first (on `perf/io-memory`, before the store)

These are independent of the store and small enough to go into the current PR:

1. **1-thread normalisation.** Have `normalize_using_input` call `_normalize_chrom_worker` inline when
   `threads <= 1`, keeping the dict order, zero tracks and warnings the same as `strand_atom`. Expected for H3 at
   1 thread: 5.1 → ~1.6 GB; the rest of that run stays under 1.6 GB.
2. **1-thread bigWig loading.** Have `_get_coverage_bw` use `_bw_blockwise_sparse`, with the dense fallback, like
   `_bw_chrom_worker`. Keep the read-count accumulation order and the cache branch. Expected for H2/T3 at 1 thread:
   3.1 → ~2 GB (estimate).

Each needs a differential test against the old code, like those in `tests/diff/test_differential_gate1.py`.

## 6. Implementation plan for the store

New branch `perf/sparse-coverage` on top of `perf/io-memory` (after the quick wins). One byte-identical commit per
step:

1. **`CoverageTrack`** plus `track_from_npy` / `track_to_npy`. Readers (`_load_coverage`, `merge_covs`,
   `_blockwise_element_definition`, `_coverage_cost`) accept either form. Property tests: `window`, `blocks`, `sum`
   and `max` against the dense array, including empty tracks, runs crossing block edges, and windows past the
   chromosome ends.
2. **BAM** produces tracks: region output, the contig worker, and the serial fallback. Differential test: tracks
   equal `np.load` of the current dumps (values and dtype) for all twelve layout/end/rc combinations, plus read
   totals.
3. **bigWig** produces tracks, including the chromosome intersection in `get_coverage_bw` and the read-count
   adjustment.
4. **Sparse merge.** Same dicts (including order), dtypes, `pl_rc`/`mn_rc` and zero tracks as
   `merge_replicates_bw`.
5. **Sparse normalisation.** Same result, plus the `sf` edge cases and chromosomes missing from the control.
6. **Process layout.** Control before the pool, dumps no longer written, `release_coverage_files` and the background
   deleter become no-ops for tracks. Decide the cache (4.5).

**Local checks before every push:**
- `pytest -m "not slow" tests/unit tests/diff tests/test_em.py tests/functional_test.py`, in a uv venv built from
  `bench/requirements.txt` plus `pytest` and `hypothesis`.
- The kit's equivalence suite on T1_bw, T1_bam and T1_bam_ct.

## 7. Benchmark plan

1. **Add kit datasets for the untested paths**, built from data the kit already fetches:
   - bigWig replicates (H2 + H5 bigWigs, or T3 + its rep2);
   - bigWig + control;
   - BAM replicates + control;
   - `--dont-merge-reps` on H1 and on T4.

   Run them at least in the equivalence suite, against `ref`.
2. **New variant** `sparse` (`perf/sparse-coverage`) in `bench/variants.json`, with equivalence on every dataset.
   Extend the `paired` suite (or add a second one) to time `sparse` against `io` alternately at 1, 16 and 64
   threads.
3. **Profile** (`bench/run profile --variant sparse --modes trace`) on H1 and H3. Check:
   - written bytes in `files.tsv` (only BED, tabix and pybedtools files should remain);
   - per-phase PSS;
   - store `nbytes`;
   - phase 1/2 times against io.
4. Record `cgroup_peak_gb` (`exec --systemd-scope`), if the node supports it. Unlike PSS, it includes the page cache,
   so it shows what the dumps cost today.

## 8. Not addressed here

- **The worker-pool floor**, about 35 MB PSS per worker (2.2 GB at 64 threads).
- **Filtering, pairing and QC memory** (data frames). This is the 64-thread peak on bigWig and T datasets.
- **`pints_boundary_extender`** keeps its file-based loading.
