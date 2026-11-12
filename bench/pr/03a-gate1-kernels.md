## Exact fast paths for the ZIP fit, the small-peak p-values and independent filtering

Stacked on B1 (`pr/03-test-harness`); review only the last five commits.

**The output is byte-identical to the base branch.** Each fast path performs the same floating-point operations as
the code it replaces, in the same order, so the functions return the same values bit for bit, not only the same peak
files. For inputs where that cannot be guaranteed, the original code runs.

### What changes for users

- No new command-line options, and the output does not change.
- A new optional extra, `pip install 'pyPINTS[fast]'`, installs Numba (`numba>=0.56`). Without it, `ZIP.fit`
  runs the original EM loop and the output is the same; the other changes below do not need Numba. The run time
  without Numba was not measured.
- With Numba installed, the kernels are compiled the first time a window is fitted and cached on disk
  (`@njit(cache=True)`), so later runs load them instead of compiling again. This PR adds no warm-up, so each worker
  process compiles or loads them on its first fit; B3 (parallel engine) compiles them once before its pool forks.

### Why

At one thread, most of a PINTS run is the window tests, and most of their time goes to `ZIP.fit`: millions of
small EM fits, each iteration a handful of NumPy and SciPy calls on short arrays. The cost is call overhead rather
than arithmetic. `scipy.stats.poisson.pmf` alone spends about 50 µs per call in argument parsing and support checks
before it computes anything.

Stratified filtering has two similar costs: the small-peak p-values are computed row by row with
`DataFrame.apply`, and independent filtering runs the BH correction 200 times and builds a DataFrame for each of the
200 quantile thresholds, although only one of them is kept.

### What changes

- **ZIP log-likelihood** (`_poisson_pmf`). `ZIP._log_likelihood_core` calls a copy of the arithmetic inside
  `poisson.pmf` (`exp(xlogy(k, mu) - gammaln(k + 1) - mu)`, clipped to [0, 1], with `k` promoted to int64 as SciPy
  does), without SciPy's generic wrapper. A `k` that is not a 1-D integer array, a non-scalar `mu`, an empty `k`
  and a negative count go to `poisson.pmf` as before. A scalar `mu` that is negative or NaN returns an array of NaN
  of the length of `k` directly, without calling SciPy.
- **ZIP EM** (`ZIP._fit_fast`). It runs the same EM iterations as `ZIP.fit` without building the `hat_z` and
  indicator arrays:
  - `hat_z` is one value on the zero positions and 0 elsewhere, so its sum is a masked sum. `_masked_sum`
    computes that sum with NumPy's own summation order (blocks of 8192 elements, each summed pairwise with eight
    accumulators), so the result equals `np.where(mask, a, b).sum()` bit for bit. It is compiled with Numba
    (no fast-math).
  - The log-factorials of the window's distinct counts are computed once per fit (`_LLCache`).
  - The fast path returns `None`, and `fit` then runs the original loop from the start, if the window is not a
    1-D integer array, has a negative count or a non-finite sum, or if an iteration leaves the range in which the
    two forms are equal (`hat_z` outside [0, 1], `lamb` < 0). It is not used with `debug=True`.
- **Numba is optional.** `setup.py` gains an extra, `pip install 'pyPINTS[fast]'` (`numba>=0.56`). Without Numba,
  `import numba` fails, `_HAS_NUMBA` is false and `ZIP.fit` runs the original loop. The output is the same either
  way; only the speed differs. Numba is kept optional because its wheels often lag new Python releases and it adds
  llvmlite to the install. The kernel is compiled on first use and cached (`cache=True`).
- **Small-peak p-values** (`stratified_filtering`). The row-wise `apply(cached_pois)` becomes one `poisson.sf`
  call over the small peaks. `poisson.sf` is element-wise, so each value is the one the cached scalar call
  returned. Empty tables and a non-numeric `reads` column keep the row-wise path.
- **Independent filtering** (`independent_filtering`). Read counts are integers, so many of the 200 quantile
  thresholds are equal, and equal thresholds select the same windows. BH is now computed once per distinct
  threshold, and only the selected subset is built as a DataFrame, in the same way as before. A NaN threshold is
  never reused. The quantiles tested, windows remaining and windows rejected are unchanged, and so is the log line.

### Timing

Against B1 (`pr/03-test-harness`; no `pints/` changes relative to A1), on an AMD EPYC 9554 server at 16 threads,
three alternating runs each, median wall time, Numba installed:

| Dataset | Wall time |
|---|---|
| `H3_bam_ct_gm12878_grocap` | 616.9 → 454.7 s (0.737) |
| `T4_bam_reps_hct116` | 796.0 → 446.7 s (0.561) |
| `T5_bam_ct_k562_grocap` | 394.6 → 292.2 s (0.741) |

Outputs are identical. 64 threads were not measured: the base's wall time hardly changes between 16 and 64 threads
(see B3 (parallel engine)). Most of the remaining time goes to work that B3 (parallel engine) parallelises; B3 on
its own takes these runs to 0.05–0.19 of this PR's time.

### Testing

- `tests/diff/test_differential.py` (new; Hypothesis property tests) compares the working tree with the reference
  copy `pints_ref`. Results must be equal bit for bit, with the same dtypes and Python types. It covers:
  - interval merging, sliding windows, peak cutting, quasi-max-score segments, IQR thresholds and `get_elbow`;
  - `ZIP.fit` (including very short windows), the ZIP log-likelihood and the ZIP likelihood-ratio test;
  - local environment refinement, the per-chromosome window tests and element definition with replicates;
  - stratified filtering (with and without `--disable-small` and QC) and independent filtering;
  - pairing, BAM parsing and bigWig loading.

  `PINTS_DIFF_EXAMPLES` (default 25) scales the number of examples per test.
- `tests/diff/test_differential_kernels.py` (new) checks:
  - that the fast EM path equals `pints_ref`'s `ZIP.fit` for window lengths around the 8-, 128- and 8192-element
    boundaries of NumPy's summation, for signed, unsigned and float dtypes, for negated windows, and with
    non-default starting values;
  - that integer windows take the fast path;
  - that with `_HAS_NUMBA` false, the result equals the fast path's and `pints_ref`'s.
- Full suite (`pytest tests`, Python 3.12.13, NumPy 2.1.0, pandas 2.3.1, SciPy 1.15.0, Numba 0.61.2): 137 passed.
- chr22 test data, with bigWig, BAM, and `--annotation-gtf --focused-chrom chr22` input, at 1 and 4 threads: the
  bidirectional, divergent and unidirectional peaks are byte-identical to `main`.

The bit-for-bit equality relies on SciPy's `poisson.pmf` arithmetic and NumPy's summation order, which neither
library documents. The differential tests compare against the reference copy, which calls the installed SciPy and
NumPy, so a change in either library would make them fail rather than go unnoticed.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
