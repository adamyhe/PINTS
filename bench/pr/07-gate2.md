## Opt-in `PINTS_FAST_MATH=1`: compiled ZIP fit for the window tests

Stacked on the serial-steps PR (`pr/15-serial-steps`); review only the last five commits.

**With `PINTS_FAST_MATH` unset (the default), the output is byte-identical to the base branch and to `main`.**
With `PINTS_FAST_MATH=1`, p-values and padj change by a small relative amount, within a drift budget fixed before
testing (below). Nothing else in the stack depends on this PR, so it can be left out or merged later.

### Why

The window tests fit a zero-inflated Poisson (ZIP) model to every tested window's background with `ZIP.fit`, an EM
loop in numpy. The compiled path from PR 3a (`ZIP._fit_fast`) is exact: it replays numpy's summation order so the
results stay bit-identical, and that limits how much it can do. A fit allowed to differ in the last bits can be
much faster. On real background windows from T2, this PR's kernel is 27× faster than the exact fit (14 µs vs
386 µs) and agrees with it to 1.3e-11 relative; its log-likelihood agrees to 2.1e-14 (measured on `6ef1c59`, whose
kernel is the same code; see the drift section).

Bit-identical output is the right default, so this is opt-in. Users who accept results that differ at about the
1e-11 level can switch it on.

### What changes

The kernel (`_zip_em_nb` and `_zip_fit_window_nb` in `stats_engine`) runs the whole of `ZIP.fit` in one numba
function: histogram, moment initialisation, EM and log-likelihood. It uses the same iteration, initialisation and
stopping rule as `ZIP.fit`. It differs in two places: the `hat_z` sums are closed-form (`n_zero * v`), and `exp`
and `log` come from LLVM's libm instead of numpy's SIMD kernels. So results agree to about 1e-15 relative per
iteration rather than bit for bit. Where `PINTS_FAST_MATH=1` uses it:

- **The final background fit of each tested window** (`fit_background` in `check_window`). Its λ and π give the
  window's p-value, which is continuous, so the drift is small and stays small.
- **The likelihood-ratio test in `--sensitive` mode** (`lrt_background`). Its three fits and log-likelihoods come
  from the kernel, and `chi2.sf(x, 2)` is computed as `exp(-x / 2)` (`_chi2_sf_df2`, within 1e-12 relative of
  scipy). This is 20× faster than `lrt` on T2 windows, with p within 2e-12.
- **Local environment refinement (LER), only where certified.** LER makes discrete decisions (mask a sub-peak or
  not). A first version used the kernel here too and failed the drift budget: near-ties flipped decisions, and the
  flips changed peak calls. Now the kernel's result is used only when every decision is certified to equal the
  exact fit's, and the exact fit is used otherwise. LER output is therefore identical to the exact path's.
  - The kernel reports near-ties of its own branches: the stopping test `|ll - prev| < stop_diff` (within 1e-12
    relative), and the moment initialisation `m >= s2` (within 1e-9). An initialisation near-tie is resolved by
    passing in numpy's moments, which equal `ZIP.zip_moment_estimators` bit for bit.
  - Results with no zeros, all zeros, or π = 0 throughout are bit-exact (λ = sum / n, as in `ZIP.fit`) and
    carry no error.
  - Many candidates (`bgIQR`, more than three): one compiled call (`_ler_mus_nb`) fits all of them from a single
    histogram. The decisions are certified when every λ is further from the threshold q1 − IQR than its error
    bound (1e-8 relative) plus three times the largest bound. Identical histograms give identical kernel values,
    and only exact values break those ties, so about 60% of these calls fall back on T2.
  - Few candidates: certified when `floor(mu_pk)`, the p-value cutoff and the eLER threshold are each further than
    the error bound from deciding the other way (99.5% of decisions on T2).

What stays exact even with the flag set:
- the fit of each candidate peak itself (`mu_peak`, `pi_peak` in the evaluation loop);
- every LER decision (above);
- the thresholds computed before evaluation;
- the Poisson model (`--model Poisson`), and any ZIP fit with user-supplied initial values or in debug mode.

The fitted background λ and π do feed two comparisons in `--sensitive` mode: whether a peak uses the LRT p-value
(`pi_0 > 0` and the fold-change cutoff). These use the kernel's values, so a value within the drift of a cutoff
could go either way. The drift runs below include `--sensitive`.

**Fallbacks.** Each kernel entry point checks its inputs and calls the exact code when it can't be used: numba not
installed, debug mode, user-supplied initial values, a window that is not a 1-D integer array, negative counts,
or counts above 2^20 (the histogram size limit). `PINTS_FAST_MATH=1` without numba therefore runs the exact path,
and the output is byte-identical to the default. numba is in the existing `fast` extra
(`pip install 'pyPINTS[fast]'`). `numba_warmup` compiles the kernels once before the pool starts, so workers
inherit them.

`PINTS_FAST_MATH` is read once, when `pints.stats_engine` is imported.

### Drift budget

Fixed before testing (Biomni's pre-registered budget; `compare_pair` in `bench/report.py`). Fast and exact runs of
the same dataset and config must have:
- identical candidate coordinates in the per-strand candidate tables;
- p and padj within 1e-6 relative;
- significance flips (padj against 0.1) only within 1e-4 relative of the cutoff, and at most 0.01% of the
  significant candidates;
- a Jaccard index of at least 0.999 for each final set (bidirectional, divergent, unidirectional).

Measured on the old Gate 2 branch, `perf/gate2-fast-math` at `6ef1c59`, against Gate 1 (`f8110c7`), on
cbsubscb19 at 32 threads. Its Gate 2 changes are the same as this PR's: the five commits add and remove the same
lines, and `git diff --stat 4ccecaf 6ef1c59 -- pints` (Gate 2 over its exact base) shows the same two files and
the same 417 insertions and 28 deletions as this PR. Only the base differs: `6ef1c59` sits on the older I/O branch
(`4ccecaf`), without PRs 5 to 15 (`git diff --stat 6ef1c59 pr/07-gate2 -- pints`: 6 files, 1845 insertions, 422
deletions). Results, 19 datasets × 3 configs (default, `--disable-small`, `--sensitive`):

- **57 of 57 pass.**
- Candidate coordinates identical in every run.
- Largest relative difference in p: 1.6e-11. In padj: 1.6e-11, except H4 (`H4_bw_liver`, `--sensitive`): 4.2e-7.
- No significance flips.
- Jaccard of the final sets: unidirectional 1.0 everywhere; bidirectional and divergent at least 0.99947
  (T5, `--disable-small`), and 1.0 on 10 of the 19 datasets in the default config.

TODO: explain why some final bidirectional and divergent sets differ (Jaccard below 1) when no candidate's padj
crosses 0.1. Ties in a later step (pairing) are a guess, not checked.

TODO: rerun the drift report on this branch (`4fc7517`), on the current base.

### Timing

With `PINTS_FAST_MATH=1`, against `pr/15-serial-steps` (`0357468`, this branch's base), on cbsubscb19 (AMD
EPYC 9554), three alternating runs each, median wall time:

| Dataset | 16 threads | 64 threads |
|---|---|---|
| `H1_bam_reps_karpas422` | 41.4 → 23.5 s (0.567) | 25.0 → 17.5 s (0.699) |
| `H3_bam_ct_gm12878_grocap` | 51.5 → 32.0 s (0.621) | 31.7 → 23.7 s (0.747) |
| `T4_bam_reps_hct116` | 22.4 → 11.4 s (0.510) | 14.2 → 9.3 s (0.658) |
| `T5_bam_ct_k562_grocap` | 9.4 → 6.6 s (0.702) | 6.7 → 5.5 s (0.819) |

All of the gain is in the phase-2 window tests (H1 at 16 threads: 22.7 → 5.3 s; H3: 28.0 → 8.7 s); the other
steps are unchanged. Outputs differ from the base, as expected with the flag on (see the drift report above).
Without the flag, this branch's output is identical to `main`.

### Testing

- `tests/unit/test_gate2_fastmath.py` (new, 26 tests; skipped without numba):
  - the kernel fit matches `ZIP.fit` within 1e-9 relative (λ, π, log-likelihood) with the same convergence flag, on
    two seeds of about 400 random windows each, int8 to int64; `fit_background` uses the kernel only when the flag is on;
  - the flag is off by default;
  - fallbacks: user-supplied initial values, debug mode, float windows and 2-D windows never call the kernel and
    return `fit`'s result;
  - certified LER: `bgIQR` and `pkIQR` `remove_peaks_in_local_env` give the same environment, count and cache with
    the flag on and off, on 6 seeds of sparse coverage with many tied sub-peaks, and both certified paths are used;
  - `_ler_mus_nb` matches the exact fits (bit for bit where it reports exact), and rejects negative values;
  - `_chi2_sf_df2` matches `chi2.sf(x, 2)` within 1e-12 relative from 1e-300 to 1400, and handles NaN and large x;
  - the kernel LRT matches `lrt` within 1e-9 relative, and falls back in debug mode and for float input.
- Full suite: 495 passed (base branch: 469), with numba 0.61.2 and pandas 2.3.1.
- chr22 test data, with bigWig, BAM, and `--annotation-gtf --focused-chrom chr22` input, at 1 and 4 threads, flag
  unset: the bidirectional, divergent and unidirectional peaks are byte-identical to `main`.
- On the old branch, the certified-LER commit left the Gate 2 output byte-identical to the commit before it on
  T1_bw, T1_bam, T1_bam_ct, T2 and T4, at 1 and 8 threads.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
