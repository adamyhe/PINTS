## Add a test harness: a frozen reference copy of PINTS and characterisation tests

This PR adds tests only: a frozen copy of upstream PINTS (`tests/reference/pints_ref`) and characterisation tests of
the current engine. Nothing in `pints/` changes, so the output is unchanged and nothing changes for users. Stacked
on #30 (`pr/02-trackb-fixes`); review only the last commit (`cf1ec8b`). A PR from a fork can only target `main`, so
GitHub also lists #30's commits, #31's setup.py commit and its merge (see the stacking note below). This PR is a
draft until #30 merges; the branch will then be rebased onto `main`, and only the harness commit will remain.

### Why

The performance PRs that follow (B2–B10 in the overview on #30) change how PINTS computes its results. With default
options, each claims output byte-identical to the current release. That has to be testable on more inputs than a
few real datasets, and the code being tested against has to stay fixed while the working tree changes. This PR adds
the fixed copy and the shared test helpers on their own, so they can be reviewed before any optimisation depends on
them.

### What changes

- **`tests/reference/pints_ref`** is a copy of upstream PINTS `81b2dd9` (`calling_engine`, `extension_engine`,
  `io_engine`, `qc_engine`, `stats_engine`). The only edits are:
  - imports renamed from `pints` to `pints_ref` (6 lines);
  - the one fix of #30 that changes output with default options: when only one strand has pairs, they are merged on
    their own (`pair()`).

  Nothing in `pints/` imports it. The reference is not meant to be edited. Later PRs compare against it, so it
  should change only when an intentional output change is accepted.
- **`tests/conftest.py`** puts `tests/reference` and the repository root on `sys.path`, sets `OMP_NUM_THREADS` and
  `OPENBLAS_NUM_THREADS` to 1 unless already set, gives each test its own pybedtools temp directory, and registers
  two markers (`spec` and `slow`).
- **`tests/helpers.py`** builds synthetic inputs: zero-inflated coverage with peaks, bgzipped and tabix-indexed
  subpeak and candidate beds, BAM files (with the expected 5′/3′ positions for each library layout), and bigWigs.
- **`tests/unit`** characterises the current engine:
  - `test_core_algorithms.py`: interval merging, sliding windows, quasi-max-score segments, peak cutting, the ZIP and
    Poisson models, IQR thresholds, `get_elbow`, independent filtering;
  - `test_io.py`: BAM parsing for each library layout, strand option and `--exp-type` preset, coverage caching and
    dtypes, bigWig loading, input normalisation, replicate merging;
  - `test_ler_calling.py`: local environment refinement, window tests, element definition, stratified filtering
    and pairing.
- **`tests/diff/test_reference.py`** checks that `pints_ref` is loaded from its own directory and not from the
  working tree, and that `independent_filtering` and `ZIP.fit` give identical results in both. These two functions are
  still the same in both, so the checks pass trivially here. They start to matter in B2 (exact kernels).

The tests in this PR need pytest on top of PINTS's own dependencies. The install line in the `tests/conftest.py`
docstring (`pip install pytest hypothesis numba`) covers the whole series: the differential tests added in B2 (exact
kernels) use `hypothesis`, and `numba` is the optional extra B2 adds.

### Stacking note

The branch also carries two commits that are not part of this PR. They come from #31 (`pr/09-drop-biopython`) and
drop out when the branch is rebased onto `main` after #30 and #31 merge:

| Commit | What it is |
|---|---|
| `80187c8` | setup.py with one dependency per line, the first commit of #31 |
| `4f317bc` | The merge commit that brings `80187c8` in. It changes nothing else |

The later PRs need it: the one-per-line setup.py keeps those that add or remove dependencies from conflicting with
#31.

### Testing

- Full suite (`pytest tests`, Python 3.12.13, NumPy 2.1.0, pandas 2.3.1, SciPy 1.15.0): 114 passed.
- `pints/` is unchanged from the base branch (#30), so the peak files are #30's.
- chr22 test data, with bigWig, BAM, and `--annotation-gtf --focused-chrom chr22` input, at 1 and 4 threads: the
  bidirectional, divergent and unidirectional peaks are byte-identical to `main`.
- `pints_ref` was compared file by file against `git show 81b2dd9:pints/<module>.py`: the only differences are the
  import renames and the pairing fix listed above.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
