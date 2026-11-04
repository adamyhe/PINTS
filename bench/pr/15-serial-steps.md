## Spread the steps after evaluation over the worker pool

Stacked on the opt-in fixes PR (`pr/08-trackb-flags`); review only the last commit. It merges cleanly with PRs
9–14 in any order.

**The output is byte-identical to the base branch.** This PR changes only how the work is split across processes.

### Why

At 64 threads on `H3_bam_ct_gm12878_grocap`, about a third of a run came after the evaluation phase, and during that time only
one or two processes were busy:

| Step | Time (H3, 64 threads) | Why it was slow |
|---|---|---|
| Independent filtering | about 2.6 s per call, two calls per strand | 200 separate `Series.quantile` calls took half the time |
| Re-reading the two strand tables | 2.3 s | read one after the other in the main process |
| Pairing | 4.0 s | one task per strand, each indexing the whole opposite-strand candidate file before its first lookup |

### What changes

- **Pairing is split by chromosome.** With a pool, `pair()` splits each strand's significant peaks into runs of
  consecutive chromosomes. Each pool task gets at least `_PAIR_CHUNK_LINES` (2000) peaks. A task builds an index of
  only its own chromosomes' candidates, read through the candidate file's tabix index
  (`PeakBedIndex.from_tabix`). The tasks' output parts are then joined in order. Each output line depends only on
  its query peak and the candidates on the same chromosome, so the joined files equal those of one call over the
  whole strand.
  - A chromosome missing from the candidate index gets no entry. Looking it up raises the same `ValueError` that
    tabix raises, and `merge_opposite_peaks` already treats that as "no opposite peak".
  - Chromosome names that a tabix region string could misread (for example `HLA-A*01:01`), and records outside
    the format `from_file` accepts, fall back to the whole-file index as before.
  - Without a pool (`--thread 1`), pairing runs as before.
- **Independent filtering makes one quantile call.** `independent_filtering` now calls `quantile` once with all 200
  quantiles. pandas computes a scalar quantile as the one-element list `[q]`, so each value is bit for bit the one a
  separate call returns. In a separate check over 120,600 comparisons, none differed.
- **The strand tables are read on the pool together.** A pickled DataFrame is an exact copy of the one read in
  the worker.

### Timing

`bench/run ab --a tip --b serial --fresh` measures the stack tip (`pr/08-trackb-flags`) against this branch:
- on cbsubscb19 (AMD EPYC 9554);
- with three alternating repetitions per cell, both variants in one session;
- reporting the median wall time.

The output hashes matched in every cell.

| Dataset | Threads | Before | After | After / before | Filtering + pairing |
|---|---|---|---|---|---|
| `H1_bam_reps_karpas422` | 16 | 43.2 s | 41.3 s | 0.956 | 8.1 → 6.3 s |
| `H1_bam_reps_karpas422` | 64 | 27.5 s | 25.6 s | 0.929 | 7.8 → 5.9 s |
| `H3_bam_ct_gm12878_grocap` | 16 | 56.0 s | 52.8 s | 0.943 | 12.6 → 9.6 s |
| `H3_bam_ct_gm12878_grocap` | 64 | 37.2 s | 34.0 s | 0.914 | 12.4 → 9.3 s |
| `T4_bam_reps_hct116` | 16 | 26.7 s | 26.2 s | 0.979 | 3.1 → 2.4 s |
| `T4_bam_reps_hct116` | 64 | 18.2 s | 17.7 s | 0.971 | 3.0 → 2.4 s |
| `T5_bam_ct_k562_grocap` | 16 | 9.4 s | 9.3 s | 0.995 | 1.7 → 1.4 s |
| `T5_bam_ct_k562_grocap` | 64 | 7.4 s | 7.1 s | 0.953 | 1.7 → 1.3 s |

In the H3 64-thread profile after this change:
- each `independent_filtering` call takes 0.5–0.6 s;
- the pairing tasks take about 0.07 s each.

What remains after evaluation is mostly in the main process:
- merging and sorting the three output files: about 1.4 s;
- writing the two significant-peak files: about 0.5 s.

### Testing

- `tests/diff/test_differential_serial_steps.py` (new; Hypothesis property tests) checks:
  - independent filtering against the reference implementation in `pints_ref`, exactly;
  - scalar against list quantiles for integer, float, object, empty and single-value Series;
  - the per-chromosome index against tabix, including a chromosome missing from the file;
  - that region-like chromosome names fall back to the whole-file index;
  - that chunked pairing writes the same files, byte for byte, as per-strand pairing for chunk sizes 1, 3, 40 and
    10⁶, and leaves no part files;
  - that tables read on the pool equal a read in the parent.
- The new tests and `tests/diff/test_differential_pairing.py` pass under pandas 2.3.1 and 3.1.0. The full suite
  passes under pandas 2.3.1 (459 passed).
- chr22 test data, with bigWig, BAM, and `--annotation-gtf --focused-chrom chr22` input, at 1, 4 and 8 threads,
  under pandas 2.3.1 and 3.1.0: the bidirectional, divergent and unidirectional peaks are byte-identical to the
  base branch. No part files are left in the output directory.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
