## Add a test harness: a frozen reference copy of PINTS and characterisation tests

Stacked on the bugfix PR (`pr/02-trackb-fixes`); review only the last commit (`45639d7`). `git rev-list --count
pr/02-trackb-fixes..pr/03-test-harness` gives 3. The other two bring in a commit of another PR and disappear when the
stack is rebased onto `main` after group A has merged:

| Commit | What it is |
|---|---|
| `80187c8` | setup.py with one dependency per line, the first commit of the Biopython PR (`pr/09-drop-biopython`) |
| `eb5b1c6` | The merge commit that brings `80187c8` in. It changes nothing else |

The stack needs it: the one-per-line setup.py keeps later PRs that add or remove dependencies from conflicting with
the Biopython PR.

**The commit under review adds only files under `tests/`. It changes nothing in `pints/`, so the output is that of the
base branch.**

### Why

The performance PRs that follow (B2–B10) change how PINTS computes its results. With default options, each claims
output byte-identical to the current release. That has to be testable on more inputs than a few
real datasets, and the code being tested against has to stay fixed while the working tree changes. This PR adds the
fixed copy and the shared test helpers on their own, so they can be reviewed before any optimisation depends on them.

### What changes

- **`tests/reference/pints_ref`** is a copy of upstream PINTS `81b2dd9` (`calling_engine`, `extension_engine`,
  `io_engine`, `qc_engine`, `stats_engine`). The only edits are:
  - imports renamed from `pints` to `pints_ref` (6 lines);
  - the one fix of the bugfix PR that changes output with default options: when only one strand has pairs, they
    are merged on their own (`pair()`).

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
  still the same in both, so the checks pass trivially here. They start to matter in the next PR.

The harness needs only pytest. The differential tests added in B2 also need `hypothesis`.

### Testing

- Full suite (`pytest tests`, Python 3.12.13, NumPy 2.1.0, pandas 2.3.1, SciPy 1.15.0): 114 passed.
- chr22 test data, with bigWig, BAM, and `--annotation-gtf --focused-chrom chr22` input, at 1 and 4 threads: the
  bidirectional, divergent and unidirectional peaks are byte-identical to `main`. (`main` ignores `--focused-chrom`,
  so for the GTF case it was given the chromosome through `highlight_chromosome`, the name its engine reads.)
- `pints_ref` was compared file by file against `git show 81b2dd9:pints/<module>.py`: the only differences are the
  import renames and the pairing fix listed above.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
