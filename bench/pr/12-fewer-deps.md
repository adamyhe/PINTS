## Fewer dependencies: drop statsmodels and requests, make matplotlib optional

Stacked on the opt-in fixes PR (`pr/08-trackb-flags`); review only the last three commits.

The branch also contains a merge of `pr/09-drop-biopython` (group A), because this PR changes the same
`install_requires` lines. That merge drops out once PR 9 is merged upstream and the stack is rebased onto `main`.

**The output is byte-identical to the base branch and to `main`.** The p-value adjustment repeats statsmodels'
operations, so the adjusted p-values are bit for bit the same.

### Why

`pip install pyPINTS` pulled in three large packages for small jobs:

| Package | What PINTS used it for | What it brings with it |
|---|---|---|
| statsmodels | `multipletests` in independent filtering | patsy, packaging; version 0.14.6 requires numpy>=1.22.3,<3, scipy>=1.8 and pandas>=1.4, not 2.1.0 |
| requests | reading PyPI's version of pyPINTS in `check_update` | urllib3, idna, charset-normalizer, certifi |
| matplotlib | the `--debug` diagnostic plots only | contourpy, cycler, fonttools, kiwisolver, pillow, pyparsing |

statsmodels' version bounds on numpy and pandas are the ones most likely to conflict with other packages in a user's
environment. After this PR, all the remaining dependencies are packages PINTS uses for its main work.

`install_requires` in `setup.py`:

| `main` | Base branch (`pr/08-trackb-flags`) | This PR |
|---|---|---|
| numpy>=1.19.2 | numpy>=1.19.2 | numpy>=1.19.2 |
| pandas>=1.1.5 | pandas>=1.1.5 | pandas>=1.1.5 |
| scipy>=1.5.2 | scipy>=1.5.2 | scipy>=1.5.2 |
| pysam>=0.16.0.1 | pysam>=0.16.0.1 | pysam>=0.16.0.1 |
| pyBigWig | pyBigWig | pyBigWig |
| requests | requests | |
| statsmodels>=0.12.1 | statsmodels>=0.12.1 | |
| biopython | biopython | (dropped by PR 9) |
| matplotlib | matplotlib | extra `plots` |
| pybedtools>=0.8.1 | (dropped by PR 6) | |

Extras: `main` has none, the base branch has `fast` (numba), and this PR adds `plots` (matplotlib).

### What changes

- **statsmodels.** `stats_engine.adjust_pvalues` replaces `multipletests(...)[1]` in `independent_filtering`. It
  repeats statsmodels' operations for the four methods `pints_caller` offers: `fdr_bh`, `bonferroni`, `fdr_tsbh`
  and `fdr_tsbky` (the two-stage methods as `fdrcorrection_twostage` with one iteration). For no p-values it raises
  `ZeroDivisionError`, as statsmodels does. `independent_filtering` already catches that error and falls back to
  Bonferroni. An unknown method raises `ValueError`. statsmodels is now needed only by the test that compares
  against it.
- **requests.** `check_update` reads `https://pypi.org/pypi/pyPINTS/json` with `urllib.request`, with the same
  10 s timeout and up to three retries, waiting 0, 6 and 12 s. Any failure is still logged, not raised. The
  epigenomic annotation download already used urllib.
- **matplotlib.** It moves to the `plots` extra: `pip install 'pyPINTS[plots]'`. It was already imported only
  inside the `--debug` branches. Without it, `pints_caller --debug` now stops before any work, with the message
  "--debug saves diagnostic plots, which need matplotlib: pip install 'pyPINTS[plots]'". Before, a missing matplotlib
  gave an `ImportError` after the computation. Runs without `--debug` never import matplotlib. The `--debug`
  help text names the extra.
- **Tests.** The run instructions in `tests/conftest.py` add statsmodels.

### Timing

Not measured. None of the changed code is on a hot path: the adjustment runs the same numpy operations as before,
and statsmodels and matplotlib were already imported lazily.

### Testing

- `tests/unit/test_adjust_pvalues.py` (new) compares `adjust_pvalues` with statsmodels' `multipletests` exactly
  (`assert_array_equal`, same dtype):
  - for each of the four methods and alpha 0.05, 0.1 and 0.5;
  - on uniform p-values, a signal and noise mixture, many ties, exact zeros and ones, all significant, none
    significant, a single p-value, and NaN;
  - on a pandas Series with a non-default index;
  - and checks that both raise `ZeroDivisionError` for no p-values.
- `tests/unit/test_check_update.py` (new): the URL and timeout, the waits between retries, the "latest version"
  and "out-dated" messages, and that being offline only logs an error.
- `tests/unit/test_debug_needs_matplotlib.py` (new): with matplotlib hidden, `peak_calling(..., output_diagnostics=True)`
  exits before doing any work, and the error names `pyPINTS[plots]`.
- Full suite: 460 passed (base branch: 438), with statsmodels 0.14.6, numba 0.61.2 and pandas 2.3.1.
- chr22 test data, with bigWig, BAM, and `--annotation-gtf --focused-chrom chr22` input, at 1 and 4 threads: the
  bidirectional, divergent and unidirectional peaks are byte-identical to `main`.
- Clean install: `pip install` of this branch into a fresh virtual environment (Python 3.12) installs only numpy,
  pandas, pyBigWig, pysam and SciPy (plus pandas' python-dateutil and tzdata); matplotlib is not installed. On the
  chr22 bigWigs, `pints_caller` writes the three peak files. With `--debug` it stops at once with
  "--debug saves diagnostic plots, which need matplotlib: pip install 'pyPINTS[plots]'".

🤖 Generated with [Claude Code](https://claude.com/claude-code)
