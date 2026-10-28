# Audit: remaining disk, memory and speed work, and how to do it

Status: 2026-10-28. Baseline is `perf/io-memory` at `0b91493` plus the two 1-thread quick wins of
[`sparse-coverage.md`](sparse-coverage.md) §5. Sources: round 3, the io profile runs
(`profile_results_cbsubscb19_20261028_0410.tar.gz`) and the experiments in §3. Byte identity with `main` (§1 of the
sparse-coverage note) applies to everything here.

## 1. Where the cost is now

**Phases, io at 64 threads (round 3 medians, s):**

| phase | H1 (49.7 s) | H2 (46.0 s) | H3 (65.0 s) | runs on |
|---|---|---|---|---|
| input parsing (+ control + normalisation on H3) | 12.2 | 2.4 | 3.7 + 9.6 | main + short pools |
| element definition | 3.0 | 2.4 | 2.6 | pool |
| phase 1 (setup, stats, thresholds) | 2.4 | 2.6 | 3.3 | pool |
| phase 2 setup | 5.4 | 5.9 | 8.0 | **main, serial** |
| phase 2 tests | 12.0 | 13.8 | 15.5 | pool |
| input normalisation/eval | 3.1 | 3.1 | 3.8 | main |
| filtering + pairing + QC | 10.5 | 14.6 | 17.6 | **mostly main, serial** |

**At 1 thread** phase 2 tests are 66–72% of the run (356–434 s), and ZIP fitting alone is 300 s of H1's 647 s and
420 s of H3's 848 s (cProfile). 2.1 M fits on H3, ≈ 17.6 EM iterations each, every iteration a handful of numpy
calls on small arrays: the cost is interpreter and dispatch overhead, not arithmetic.

**Main-process time after evaluation on H3 (py-spy, 64 threads):** waiting for the background subpeak-bed
compression (`_AsyncBed.get`), `SubpeakIndex.__init__` (Python lists from split lines, the bulk of phase 2 setup),
DataFrame construction from nested lists (`_bins_to_df`), `read_csv`, then `pair` (8.6 s, of which
`merge_opposite_peaks` 3.8 s and the rest pybedtools `cat`/`count` and `to_csv`), and `housekeeping`'s `rmtree`
(≈ 1.5 s). About 8 s between the end of evaluation and the start of pairing has no span; add spans before working
on it.

**Pool workers at 64 threads (worker-seconds, H1 / H3):** element definition 129 / 108, bgzip + tabix of the
per-chromosome subpeak beds (`index_bed_file`) 88 / 186, normalisation – / 184, replicate merge 101 / –, BAM
regions + assembly 84 / 97.

**Memory:** after the store, the 64-thread peak is filtering on bigWig and T datasets (5.8–6.6 GB) and the worker
pool floor (≈ 35 MB PSS per worker, 2.2 GB at 64).

**Disk:** after the store, the per-chromosome subpeak beds and their `.gz`/`.csi` (written, compressed, indexed,
then deleted) and pybedtools temp files in pairing.

## 2. Candidates

| # | change | where it pays | exactness risk | method |
|---|---|---|---|---|
| A | sparse coverage store (designed) | input/control 10–13 s at 64 t; 12–48 GB disk; page cache | moderate, fully specified | manual |
| B | keep subpeak beds in memory: no bgzip/tabix of per-chromosome beds; pairing builds its overlap index from the lines | 88–186 worker-s, the wait before filtering, ≈ 30 s at 1 t, disk | low: same records, tabix fallback stays for unsupported inputs | manual |
| C | phase 2 setup off the main process, or numpy-backed `SubpeakIndex` | 5.4–8.0 s serial at 64 t | low | manual |
| D | pairing: run the two strands as two pool tasks; replace pybedtools `cat`/`count`/merge with in-memory code | 4.3–8.6 s serial at 64 t | moderate: bedtools merge/collapse semantics must be replicated byte for byte | manual, with differential tests on BED bytes |
| E | `_bins_to_df` and filtering I/O: build columns as arrays, skip the CSV round trip where the frame is in hand | part of the 8 s gap | moderate: dtypes and float parsing must not change | manual, after spans |
| F | `housekeeping` deletion in the background, after outputs are written | ≈ 1.5 s tail | none | manual |
| G | ZIP EM and local-environment (LER) kernel | 66–72% at 1 t, 15–25% at 64 t | high: the reason Gate 2 exists | **search loop, narrowly**, see §4 |
| H | BAM counting (pysam per-read loop) | 37–53 s at 1 t; ≈ 0.4 s at 64 t | moderate (end coordinates, flags) | defer |
| I | pool floor: `gc.freeze()` before fork, defer heavy imports in workers | up to ≈ 2 GB at 64 t | none | manual, quick experiment |
| J | filtering memory at 64 t (bigWig/T) | the 64-thread peak on bigWig and T datasets | low | manual, profile first |

## 3. Experiments behind the choices (2026-10-28, cbsubscb19, numpy 2.1.0, AVX-512)

- **Memoising fits doesn't pay.** On T1_bam (chr22), 42 357 ZIP fits see 39 856 distinct windows: 5.9% repeats.
  Fits are also not a function of the count histogram alone (masked sums follow numpy's pairwise tree over zero
  positions), so a coarser key would not be exact.
- **A Numba EM can't be byte-identical.** Numba's `exp` differs from numpy's in the last bit on 4.6% of inputs and
  `log` on 0.1% (1 ulp). That rules out compiling the whole EM loop for Gate 1; it is Gate 2 territory.
- **Batching across windows can be.** numpy's `exp`, `log` and `xlogy` give the same bits for an element whether
  called on a scalar, a length-1 array, short chunks or a long array. So an EM that iterates many windows at once
  (struct of arrays, per-window convergence mask, per-window sums in numpy's pairwise order as `_masked_sum`
  already does) is exact on this node. **Check on an AVX2-only node** before relying on it: numpy dispatches
  different SIMD kernels there.

## 4. Decision: manual work, dependency swaps, or search loops

**Manual for A–F, I, J.** These are structural changes (where data lives, which process does what, which tool
writes which file). Each one is a single design with an exactness argument. The equivalence suite can verify them
but gives no gradient to search along, and a wrong step shows up only as a failed equivalence run hours later.
They also need to stay readable for an upstream PR. Order: A (designed), B, C, D, then E/F once spans exist; I and
J are small experiments that can go in between.

**Dependency swaps: no, as a strategy.** Removing tools is worth it (bedtools/pybedtools in pairing and subpeak
handling, bgzip/tabix for internal files: B, D), because it removes subprocesses, temp files and pickling.
Swapping libraries is not:
- pandas CSV I/O → pyarrow/polars changes float parsing and formatting, which reach the output BEDs;
- pysam → another BAM reader helps only 1-thread runs (H) and re-opens the end-coordinate and flag questions;
- Numba for transcendental-heavy code is not exact (§3).

**A search loop only for G, and only once the harness exists.** G is the one place where many small,
interacting implementation choices (batching, array layout, how LER calls sort/unique/var, where to cut Python
overhead) decide the result, and where a cheap, trustworthy oracle can exist. Harness first:
1. a corpus of real `check_window` / `remove_peaks_in_local_env` inputs captured from H1 and H3 (size to be measured);
2. the oracle: bitwise-equal fit outputs on the whole corpus, then `tests/diff`, then the kit's equivalence suite
   on T1_bam, T2, T4 (minutes). Nothing reaches a paired round without passing all three;
3. the objective: corpus wall time at 1 thread, with a paired round to confirm accepted candidates;
4. guardrails: no fast-math, no Numba transcendental functions, no changes outside `stats_engine` and the
   `check_window` path, one candidate per commit, a person reviews each accepted commit, the AVX2 check before
   merging, and a fixed budget (iterations or hours) agreed beforehand.

Without that harness, a loop would optimise against 2.5-hour paired rounds and test suites that don't cover real
window shapes, which is slow and invites overfitting to the kit. If the batched EM turns out to be a clean
single design once prototyped, do it by hand and keep the loop for tuning.

**Expected order of gains.** At 64 threads the serial main-process work (A input, C, D, E: ≈ 25–35 s of 50–65 s)
matters more than the kernel. At 1 thread only G moves the needle. Both are needed; the manual items come first
because they are cheaper and lower risk.

## 5. Gate 2 follow-up (2026-10-28, `perf/gate2-fast-math` at `6ef1c59`, on top of io `4ccecaf`)

- **Kernel speed.** The existing compiled fit is 27x faster than the exact fit on real background windows
  (14 µs vs 386 µs, T2 corpus) and within 1.3e-11 relative; its log-likelihood within 2.1e-14. Speeding it up
  further changes little; the cost had moved to the fits Gate 2 didn't cover.
- **Certified LER (`f515411`).** LER fits use the kernel when every decision is certified to equal the exact one
  (margin above the error bound; bit-exact results carry none), else the exact fit. LER output stays identical
  to Gate 1's, and Gate 2 output is byte-identical to the previous Gate 2 on T1_bw, T1_bam, T1_bam_ct, T2, T4.
  Few-candidate decisions certify 99.5% of the time. Many-candidate calls certify ~40%: masking sub-peaks
  with equal counts leaves identical histograms, the kernel then returns identical values, and exact values
  that differ by ulps decide the tie. Refitting only the candidates near the threshold's order statistics
  refits almost all of them, so it was dropped. Only exact values resolve these ties, i.e. the batched exact
  EM (§3, item G) would remove the rest of this cost, in Gate 1 and Gate 2 alike.
- **Kernel LRT (`6ef1c59`).** Sensitive mode's three LRT fits and log-likelihoods come from the kernel, and
  chi2.sf(x, 2) is exp(-x/2). 20x faster per test; drift budget passes on T1_bam, T2, T4 (default and
  sensitive): max relative dp 6.4e-12, no flips, Jaccard 1.0.
- **Indicative timing** (1 thread, shared node): T2 34.4 → 28.1 s, T4 137.2 → 102.5 s against the previous
  Gate 2. The kit's gate2 runs will give the real numbers.
- **Warm-up.** Numba's cache misses on specialisations that use default arguments (0.2–0.45 s per dtype per
  run); kernels take every argument explicitly.
