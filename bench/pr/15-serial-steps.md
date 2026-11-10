## Spread the steps after evaluation over the worker pool, and balance the window tests

Stacked on the dependency PR (`pr/12-fewer-deps`); review only the last three commits.

**The output is byte-identical to the base branch.** This PR changes only how the work is split across processes.

### Why

At 64 threads on `H3_bam_ct_gm12878_grocap`, about a third of a run came after the evaluation phase, and during that time only
one or two processes were busy:

| Step | Time (H3, 64 threads) | Why it was slow |
|---|---|---|
| Independent filtering | about 2.6 s per call, two calls per strand | 200 separate `Series.quantile` calls took half the time |
| Re-reading the two strand tables | 2.3 s | read one after the other in the main process |
| Pairing | 4.0 s | one task per strand, each indexing the whole opposite-strand candidate file before its first lookup |
| Merging the strands | about 1.4 s | the bidirectional, divergent and singleton outputs were merged one after the other |

The window tests (phase 2) also ended with a tail. Their chunks were ranked by number of subpeaks. A chrEBV chunk
costs about twice as much per subpeak as a chunk of a large chromosome, and chrEBV is a single segment, so it
cannot be split. Its two chunks were dispatched 472nd of 590 and ran alone for about 1.3 s after every other chunk
had finished.

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
- **The three strand merges run on the pool.** The bidirectional, divergent and singleton merges are independent:
  each reads its own two files and writes its own output. Without a pool, they run one after the other as before.
- **Window-test chunks are dispatched by total subpeak width.** Width predicts a chunk's run time better than the
  number of subpeaks: on H3 at 16 threads, it underestimates a chrEBV chunk by 1.35× instead of 2.1×, which is
  enough to send it first. Results are still returned in task order.

### Timing

The benchmark kit is [`bench/`](https://github.com/adamyhe/PINTS/tree/bench/external/bench) on the fork's `bench/external` branch: the variant names below are pinned
to commits in [`bench/variants.json`](https://github.com/adamyhe/PINTS/blob/bench/external/bench/variants.json), and the datasets are described in
[`bench/datasets.json`](https://github.com/adamyhe/PINTS/blob/bench/external/bench/datasets.json).

`bench/run ab --a deps --b post --fresh` measures the base branch (`pr/12-fewer-deps`) against this branch:
- on cbsubscb19 (AMD EPYC 9554);
- with three alternating repetitions per cell, both variants in one session;
- reporting the median wall time.

The output hashes matched in every cell.

| Dataset | Threads | Before | After | After / before | Filtering + pairing |
|---|---|---|---|---|---|
| `H1_bam_reps_karpas422` | 16 | 45.8 s | 43.7 s | 0.953 | 8.1 → 5.9 s |
| `H1_bam_reps_karpas422` | 64 | 27.6 s | 25.5 s | 0.922 | 7.8 → 5.6 s |
| `H3_bam_ct_gm12878_grocap` | 16 | 55.9 s | 53.7 s | 0.961 | 12.7 → 8.9 s |
| `H3_bam_ct_gm12878_grocap` | 64 | 37.5 s | 33.3 s | 0.888 | 12.4 → 8.4 s |
| `T4_bam_reps_hct116` | 16 | 26.6 s | 26.1 s | 0.979 | 3.1 → 2.4 s |
| `T4_bam_reps_hct116` | 64 | 18.4 s | 17.8 s | 0.965 | 3.0 → 2.4 s |
| `T5_bam_ct_k562_grocap` | 16 | 9.4 s | 8.9 s | 0.942 | 1.7 → 1.3 s |
| `T5_bam_ct_k562_grocap` | 64 | 7.6 s | 7.1 s | 0.929 | 1.7 → 1.2 s |

Rerun at the branch heads on 2026-11-09 (`bench/run ab --a x12 --b x15 --fresh`: `pr/12-fewer-deps` at `789516e`
against this branch at `7a18b38`), with another user's single-core jobs on the node's other cores, which makes both
sides slower: after / before 0.880–0.983 (H1 0.930 and 0.925, H3 0.934 and 0.880, T4 0.983 and 0.958, T5 0.953 and
0.932 at 16 and 64 threads), filtering + pairing 8.0 → 6.0 s on H1 and 12.7 → 8.8 s on H3 at 64 threads, identical
output in every cell.

Measured separately:
- the strand-merge commit alone gives 0.970–1.010 (H3 at 64 threads: 37.2 → 36.1 s); the rest comes from the
  first commit;
- the chunk-ordering commit gives 0.962 on H3 at 64 threads (32.9 → 31.6 s, window tests 12.98 → 11.80 s) and
  0.978–1.004 elsewhere, with identical output. Only H3 had a chrEBV-sized tail.

In the H3 64-thread profile of the first commit:
- each `independent_filtering` call takes 0.5–0.6 s;
- the pairing tasks take about 0.07 s each;
- writing the two significant-peak files, in the main process, takes about 0.5 s.

### Testing

- `tests/diff/test_differential_serial_steps.py` (new; Hypothesis property tests) checks:
  - independent filtering against the reference implementation in `pints_ref`, exactly;
  - scalar against list quantiles for integer, float, object, empty and single-value Series;
  - the per-chromosome index against tabix, including a chromosome missing from the file;
  - that region-like chromosome names fall back to the whole-file index;
  - that chunked pairing writes the same files, byte for byte, as per-strand pairing for chunk sizes 1, 3, 40 and
    10⁶, and leaves no part files;
  - that tables read on the pool equal a read in the parent;
  - that `pair()` on a pool writes the same files as without one, including when one or both strands have no
    significant peaks;
  - `tests/unit/test_run_tasks.py`: tasks are dispatched costliest first, and results come back in task order.
- The new tests and `tests/diff/test_differential_pairing.py` pass under pandas 2.3.1 and 3.1.0. The full suite
  passes under pandas 2.3.1 (469 passed).
- chr22 test data, with bigWig, BAM, and `--annotation-gtf --focused-chrom chr22` input, at 1, 4 and 8 threads,
  under pandas 2.3.1 and 3.1.0: the bidirectional, divergent and unidirectional peaks are byte-identical to the
  base branch. No part files are left in the output directory. At 1 and 4 threads under pandas 2.3.1, they are
  also byte-identical to `main`.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
