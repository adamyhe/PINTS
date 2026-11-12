## Opt-in `--fast-math`: compiled ZIP fit for the window tests

Stacked on B9 (`pr/15-serial-steps`); review only the last eight commits. `gate2` in the branch name is our label
for this opt-in, tolerance-bounded change; `gate1` (B2, B3) was the bit-identical work.

**Without `--fast-math` (the default), the output is byte-identical to the base branch and to `main`.** With it,
p-values and padj change by small relative amounts. Over 57 runs (19 datasets × 3 configurations): p within 1.6e-11
relative; padj within 1.6e-11 except one run (H4, `--sensitive`) at 4.2e-7, against a budget of 1e-6 fixed before
testing; no significance flips. Tied p-values can change which opposite-strand peak is paired, so the bidirectional
and divergent sets can differ slightly (Jaccard index at least 0.99947). Nothing else in the stack depends on this
PR, so it can be left out or merged later.

### What changes for users

- A new `pints_caller` option, `--fast-math`, off by default. It needs numba, from the existing `fast` extra
  (`pip install 'pyPINTS[fast]'`); without numba it logs a warning and the output is the default one.
- With `--model Poisson` it does nothing: every kernel entry point checks `isinstance(stat_tester, ZIP)`.
- The README describes the option (see below).

### Why

The window tests fit a zero-inflated Poisson (ZIP) model to every tested window's background with `ZIP.fit`, an EM
loop in numpy. The compiled path from B2 (exact kernels), `ZIP._fit_fast`, is exact: it replays numpy's summation
order so the results stay bit-identical, and that limits how much it can do. A fit allowed to differ in the last
bits can be much faster. On real background windows from T2, this PR's kernel is 27× faster than the exact fit
(14 µs vs 386 µs) and agrees with it to 1.3e-11 relative; its log-likelihood agrees to 2.1e-14.

Bit-identical output is the right default, so this is opt-in, for users who accept the differences above.

### What changes

The kernel (`_zip_em_nb` and `_zip_fit_window_nb` in `stats_engine`) runs the whole of `ZIP.fit` in one numba
function: histogram, moment initialisation, EM and log-likelihood. It uses the same iteration, initialisation and
stopping rule as `ZIP.fit`. It differs in two places: the `hat_z` sums are closed-form (`n_zero * v`), and `exp`
and `log` come from LLVM's libm instead of numpy's SIMD kernels. So results agree to about 1e-15 relative per
iteration rather than bit for bit. Where `--fast-math` uses it:

- **The final background fit of each tested window** (`fit_background` in `check_window`). Its λ and π give the
  window's p-value, which is continuous, so the difference is small and stays small.
- **The likelihood-ratio test in `--sensitive` mode** (`lrt_background`). Its three fits and log-likelihoods come
  from the kernel, and `chi2.sf(x, 2)` is computed as `exp(-x / 2)` (`_chi2_sf_df2`, within 1e-12 relative of
  scipy). This is 20× faster than `lrt` on T2 windows, with p within 2e-12.
- **Local environment refinement (LER), only where certified.** LER makes discrete decisions (mask a sub-peak or
  not). A first version used the kernel here too and failed the budget: near-ties flipped decisions, and the flips
  changed peak calls. Now the kernel's result is used only when every decision is certified to equal the exact
  fit's, and the exact fit is used otherwise. LER output is therefore identical to the exact path's.
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

What stays exact even with the option set:
- the fit of each candidate peak itself (`mu_peak`, `pi_peak` in the evaluation loop);
- every LER decision (above);
- the thresholds computed before evaluation;
- the Poisson model (`--model Poisson`), and any ZIP fit with user-supplied initial values or in debug mode.

The fitted background λ and π do feed two comparisons in `--sensitive` mode: whether a peak uses the LRT p-value
(`pi_0 > 0` and the fold-change cutoff). These use the kernel's values, so a value within the kernel's difference
of a cutoff could go either way. The drift runs below include `--sensitive`.

**Fallbacks.** Each kernel entry point checks its inputs and calls the exact code when it can't be used: numba not
installed, debug mode, user-supplied initial values, a window that is not a 1-D integer array, negative counts,
or counts above 2^20 (the histogram size limit). `numba_warmup` compiles the kernels once before the pool starts,
so workers inherit them.

**Option.** `--fast-math` is a `pints_caller` option in the Testing group, next to the opt-in fixes
(`--iqr-scale`, `--fix-ler-flanks`). `peak_calling` passes it to `stats_engine.configure_fast_math` on every run,
before the pool starts, and the pool initializer passes it to the workers as it does those fixes. The first five
commits switch the kernels on with an environment variable, `PINTS_FAST_MATH=1`, read at import; the sixth
replaces it with the option, so the variable does nothing in the final branch.

**README.** The last commit adds a paragraph on `--fast-math` after the option list. It says that a genome-wide run
is about 1.7–2× faster at 16 threads and 1.3–1.5× at 64 with the `fast` extra; in the branch-head table below, the
genome-wide datasets (H1, H2, H3, H6, T3) take 0.492–0.584 of the time at 16 threads and 0.658–0.762 at 64. Partial
datasets vary: T4 gains more (0.467 and 0.631), T2 and T5 less. It also says that the output is no longer
byte-identical, that p and padj change by up to about 1e-6 (the budget; the largest measured change is 4.2e-7),
which can change a few calls near the FDR cutoff (none changed in the drift runs), that tied peaks can be paired
differently, and that the option has no effect without Numba.

**Review order** (`git log pr/15-serial-steps..pr/07-gate2`):

1. `06ae7b1` the kernel (`_zip_em_nb`, `_zip_fit_window_nb`), behind `PINTS_FAST_MATH=1`;
2. `ff45a22` `fit_background`: the kernel only for each tested window's final background fit;
3. `672affb` unit tests;
4. `89252c3` certified kernel fits in LER (`_ler_mus_nb` and the certification helpers);
5. `52f1469` the likelihood-ratio test from the kernel in `--sensitive` mode;
6. `717ec61` the `--fast-math` option replaces the environment variable;
7. `1ea6776` names the code after the option (comments, docstrings and test names only);
8. `480f6e5` README.

### Risks

- About 200 lines of numba code (four compiled functions) and about 200 lines of Python around them (entry points,
  fallbacks and the LER certification) to maintain.
- The output with `--fast-math` depends on the last digits of LLVM's `exp` and `log`. The drift runs were all on one
  AMD EPYC 9554 server; other CPUs or numba versions were not compared.

### Drift budget

The benchmark kit is [`bench/`](https://github.com/adamyhe/PINTS/tree/bench/external/bench) on the fork's
`bench/external` branch: `bench/variants.json` pins each measured variant to a commit, and `bench/datasets.json`
describes the datasets (the IDs below, such as `H4_bw_liver`).

Fixed before testing (`compare_pair` in
[`bench/report.py`](https://github.com/adamyhe/PINTS/blob/bench/external/bench/report.py)). Fast and exact runs of
the same dataset and configuration must have:
- identical candidate coordinates in the per-strand candidate tables;
- p and padj within 1e-6 relative;
- significance flips (padj against 0.1) only within 1e-4 relative of the cutoff, and at most 0.01% of the
  significant candidates;
- a Jaccard index of at least 0.999 for each final set (bidirectional, divergent, unidirectional).

Results, 19 datasets × 3 configurations (default, `--disable-small`, `--sensitive`), on an AMD EPYC 9554 server at
32 threads (PINTS runs on one thread unless `--thread` is given):

- **57 of 57 pass.**
- Candidate coordinates identical in every run.
- Largest relative difference in p: 1.6e-11. In padj: 1.6e-11, except H4 (`H4_bw_liver`, `--sensitive`): 4.2e-7.
- No significance flips.
- Jaccard of the final sets: unidirectional 1.0 everywhere; bidirectional and divergent at least 0.99947
  (T5, `--disable-small`), and 1.0 on 29 of the 57 runs (10 of the 19 datasets in the default configuration).
- Without the option, the output is byte-identical to `main`.

**Why some bidirectional and divergent sets differ.** No candidate crosses the cutoff, but 143 of the 147 changed
pairs come from the choice of partner. When at most one opposite-strand candidate under a peak is significant,
pairing takes the one with the smallest raw p-value (`np.argmin(opposite_pvals)` in `merge_opposite_peaks`, as on
`main`). Peaks of the same height in the same local background have mathematically equal p-values. The exact path
leaves them apart in the last digits (at most 1.9e-14 relative), and the compiled fit can order them the other way
or make them equal, in which case `argmin` takes the first. Either way, the pair's extent changes. Example (T5,
default, chr21): two minus-strand candidates of 3 reads have p = 0.5092152243427122 and 0.5092152243427108 on
the exact path and 0.5092152243427138 for both with the option on, so the pair starts at 44025833 in one run and
at 44025740 in the other. The other 4 changes (H5, `--sensitive` and `--disable-small`) keep the coordinates, so
the Jaccard index is unaffected, but report a different forward summit: two plus-strand peaks in the region have
the same summit height. The order that decides between them was not traced. Every unidirectional difference is in
the printed padj (column 5) only.

<details><summary>Commits measured</summary>

The drift runs and the first timing table were measured at `fd10aab`, an earlier version of this branch with the
same kernel code, from before the option replaced the environment variable (`PINTS_FAST_MATH` off and on, the same
code paths as `--fast-math` now). The kernel micro-benchmark under "Why" was measured at `6ef1c59`, an earlier
version of the likelihood-ratio commit with the same kernel code. The branch-head timings were measured at
`4d17c9d`, the head before the commit that renamed the code after the option (comments, docstrings and test names
only) and the README commit.

</details>

### Timing

This branch without and with `--fast-math`, on 8 datasets, on an AMD EPYC 9554 server on 2026-11-09 (at `4d17c9d`,
see above): three alternating repetitions per cell, median wall time. Another user's single-core jobs ran on the
node's other cores throughout, so both sides are slower than in the earlier measurement below (H1 at 16 threads:
53.2 s here, 41.4 s there without the option), and the window tests slow most:

| Dataset | 16 threads | 64 threads | Window tests, 16 threads |
|---|---|---|---|
| `H1_bam_reps_karpas422` | 53.2 → 28.3 s (0.531) | 28.6 → 19.5 s (0.683) | 31.8 → 6.6 s |
| `H2_bw_caco2` | 57.0 → 28.0 s (0.492) | 30.4 → 20.0 s (0.658) | 36.4 → 7.5 s |
| `H3_bam_ct_gm12878_grocap` | 65.7 → 38.4 s (0.584) | 36.3 → 26.4 s (0.728) | 38.8 → 11.7 s |
| `H6_bam_reps_bw_a549r2_karpas` | 41.1 → 22.9 s (0.556) | 22.7 → 16.6 s (0.731) | 23.6 → 5.4 s |
| `T2_bam_k562` | 6.6 → 5.1 s (0.765) | 5.2 → 4.8 s (0.923) | 2.4 → 0.8 s |
| `T3_bw_gw_a549` | 27.9 → 15.4 s (0.553) | 16.0 → 12.2 s (0.762) | 15.0 → 3.2 s |
| `T4_bam_reps_hct116` | 28.5 → 13.3 s (0.467) | 16.0 → 10.1 s (0.631) | 18.6 → 3.1 s |
| `T5_bam_ct_k562_grocap` | 11.1 → 7.1 s (0.640) | 7.5 → 6.1 s (0.808) | 5.7 → 1.8 s |

H1, H2, H3, H6 and T3 are genome-wide; T2, T4 and T5 cover a few chromosomes. Peak PSS changes by −0.43 to
+0.14 GB, and CPU time (summed over processes) falls by 45–70%.

<details><summary>Earlier measurement, with the environment variable, against the base branch</summary>

With the kernels on (`PINTS_FAST_MATH=1` at `fd10aab`), against this branch's base at the time, on the same server,
three alternating runs each, median wall time:

| Dataset | 16 threads | 64 threads |
|---|---|---|
| `H1_bam_reps_karpas422` | 41.4 → 23.5 s (0.567) | 25.0 → 17.5 s (0.699) |
| `H3_bam_ct_gm12878_grocap` | 51.5 → 32.0 s (0.621) | 31.7 → 23.7 s (0.747) |
| `T4_bam_reps_hct116` | 22.4 → 11.4 s (0.510) | 14.2 → 9.3 s (0.658) |
| `T5_bam_ct_k562_grocap` | 9.4 → 6.6 s (0.702) | 6.7 → 5.5 s (0.819) |

All of the gain is in the phase-2 window tests (H1 at 16 threads: 22.7 → 5.3 s; H3: 28.0 → 8.7 s); the other
steps are unchanged.

</details>

### Testing

- `tests/unit/test_fast_math.py` (new; 10 test functions, 27 cases with the parametrised ones; skipped without
  numba):
  - the kernel fit matches `ZIP.fit` within 1e-9 relative (λ, π, log-likelihood) with the same convergence flag,
    on two seeds of about 400 random windows each, int8 to int64; `fit_background` uses the kernel only when the
    option is on;
  - the option is off by default; `configure_fast_math` turns it on and off, and leaves it off without numba;
  - fallbacks: user-supplied initial values, debug mode, float windows and 2-D windows never call the kernel and
    return `fit`'s result;
  - certified LER: `bgIQR` and `pkIQR` `remove_peaks_in_local_env` give the same environment, count and cache
    with the option on and off, on 6 seeds of sparse coverage with many tied sub-peaks, and both certified paths
    are used;
  - `_ler_mus_nb` matches the exact fits (bit for bit where it reports exact), and rejects negative values;
  - `_chi2_sf_df2` matches `chi2.sf(x, 2)` within 1e-12 relative from 1e-300 to 1400, and handles NaN and large x;
  - the kernel LRT matches `lrt` within 1e-9 relative, and falls back in debug mode and for float input.
- Full suite at `480f6e5`: 496 passed (base branch: 469), with numba 0.61.2 and pandas 2.3.1.
- chr22 test data, with bigWig, BAM, and `--annotation-gtf --focused-chrom chr22` input, at 1 and 4 threads,
  without `--fast-math`: the bidirectional, divergent and unidirectional peaks are byte-identical to `main`.
- In an earlier version of this branch, with the kernels on, the certified-LER commit left the output
  byte-identical to the commit before it (kernel only for the final background fit) on T1_bw, T1_bam, T1_bam_ct,
  T2 and T4, at 1 and 8 threads.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
