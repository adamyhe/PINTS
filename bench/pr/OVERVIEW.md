## Overview: fixes, speed-ups and fewer dependencies for PINTS

Hi Li, thanks for agreeing to review these. This comment is a map of everything we're proposing, so each PR can be
read in context. There are 14 PRs in two groups. They are labelled A1–A4 and B1–B10 here, because their GitHub numbers aren't known yet:

- **Group A: four independent fixes against `main`**, which can merge in any order;
- **Group B: a stack of ten performance PRs**, to be opened one at a time after group A has merged.

**With default options, every PR leaves `pints_caller`'s output byte-identical to the current release.** The
exceptions are in A1: three options that had no effect before now do what their help text says, and runs where
pairs are found on only one strand, which used to write malformed bidirectional and divergent files, now write
correct ones. A3 fixes the other scripts. Everything else new is opt-in.

### What it adds up to

`main` against the top of the stack (default options), on an AMD EPYC 9554 server. The output files are identical
in every run.

| Dataset (genome-wide) | `main`, 16 / 64 threads | Stack, 16 / 64 threads | Peak memory, 16 threads |
|---|---|---|---|
| PRO-cap BAM, 2 replicates (Karpas-422) | 635 / 572 s | 52 / 28 s | 12.4 → 3.5 GB |
| PRO-cap bigWig (Caco-2) | 585 / 666 s | 56 / 28 s | 6.0 → 3.4 GB |
| GRO-cap BAM with input control (GM12878) | 723 / 727 s | 61 / 34 s | 12.6 → 4.0 GB |
| PRO-cap bigWig (A549) | 383 / 386 s | 25 / 15 s | 5.7 → 2.3 GB |

PINTS runs single-threaded unless `--thread` is given. That is 10–15× faster at 16 threads and 20–27× at 64 on whole-genome data (`main` gains almost nothing from more
threads), with 38–85% less peak memory across 8 datasets. Another user's jobs shared the machine during these runs,
so all absolute times are somewhat high. The opt-in `--fast-math` (last PR) is a further 1.7–2× at 16 threads, at
the cost of exact reproducibility.

A run also no longer writes about 37 GB of per-base `.npy` coverage files on a genome-wide BAM input, and
`pip install pyPINTS` gives a working install without bedtools, bgzip or tabix on the `PATH`.

### Group A: fixes against `main`

| | What | Effect on users |
|---|---|---|
| [A1](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/02-trackb-fixes.md) | NumPy 2.5 crash in `get_elbow` (`np.cross` on 2-D vectors); `--focused-chrom`, `--stringent-pairs-only` and `--peak-width` were ignored; `--model Poisson` LRT; one-strand pairing; integer columns of result frames turning into `object` | A fresh `pip install` currently gets NumPy 2.5 and then fails on every run. The three options now take effect. |
| [A2](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/09-drop-biopython.md) | `pints_normalizer` reads the spike-in FASTA with pysam instead of Biopython; `setup.py` lists one dependency per line | One dependency fewer. Please merge without squashing (the stack contains these commits). |
| [A3](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/10-script-fixes.md) | `pints_normalizer` (did not run), `pints_visualizer`, `pints_counter`, `pints_boundary_extender`, `parse_gtf` | Listed per script in the PR |
| [A4](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/13-pandas4.md) | pandas 4 deprecations; functional tests accept pandas 3's `str` dtype | No output change under pandas 2.3 or 3.1 |

### Group B: the performance stack

Each PR builds on the one before it and passes the full test suite on its own. The branches are on our fork; the
links go to each PR's description.

| # | PR | What it does | Effect |
|---|---|---|---|
| B1 | [Test harness](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/03-test-harness.md) | Frozen copy of the current release as `tests/reference/pints_ref`, plus characterisation tests. No change to `pints/` | The reference every later PR is tested against |
| B2 | [Exact kernels](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/03a-gate1-kernels.md) | `ZIP.fit`, the small-peak p-values and independent filtering without per-call NumPy/SciPy overhead; optional Numba (`pip install 'pyPINTS[fast]'`) | 0.56–0.74× wall time at 16 threads, bit-identical values |
| B3 | [Parallel engine](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/03b-gate1-engine.md) | BAM/bigWig parsing, element definition and the window tests spread over one persistent worker pool, in chunks | Most of the speed-up (8–34× over the release at 16 and 64 threads, together with B2). Peak memory rises on large replicate BAMs; B4 brings it down |
| B4 | [I/O and memory](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/04-io.md) | Follow-ups: faster parsing, less copying, smaller intermediates, pool start-up | 0.61–0.75×; H1 peak memory 21 → 5 GB |
| B5 | [Sparse coverage](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/05-sparse-coverage.md) | Coverage held in an in-memory sparse store instead of per-base `.npy` files | 0.85–0.94×; no `.npy` files; old dumps in the output directory are no longer reused |
| B6 | [No bedtools](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/06-bedtools-free.md) | BED sort, merge, intersect and concatenation in Python (`pints/bed_ops.py`) | Drops pybedtools and the bedtools, bgzip and tabix executables |
| B7 | [Opt-in fixes](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/08-trackb-flags.md) | `--iqr-scale`, `--fix-ler-flanks`, `--strict-independent-filtering`, `--stringent-pairs-padj`: places where the code differs from the paper, the help text or the cited method | Off by default; they change calls only when set |
| B8 | [Fewer dependencies](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/12-fewer-deps.md) | statsmodels replaced by a numpy BH adjustment (bit-identical), requests by `urllib`, matplotlib moved to the `plots` extra (only `--debug` needs it) | `install_requires` becomes numpy, pandas, scipy, pysam, pyBigWig |
| B9 | [Steps after evaluation](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/15-serial-steps.md) | Filtering, pairing and strand merging on the pool; window-test chunks dispatched by size | 0.88–0.98× |
| B10 | [`--fast-math`](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/07-gate2.md) | Opt-in compiled ZIP fit for the window tests; results differ in the last digits (p and padj within about 1e-6 relative) | 0.47–0.77× at 16 threads with the option; the default is unchanged. Optional, and nothing depends on it |

Ratios are new wall time ÷ old wall time against the previous PR, at 16 and 64 threads (1 is unchanged). Across
the stack, `pints/` grows by about 3,600 lines, two of them new modules (`coverage_store.py`, `bed_ops.py`), and
`tests/` by about 8,500, of which about 4,150 are the frozen reference copy.

### How we checked that results don't change

- **Frozen reference.** B1 adds an unmodified copy of the release's `pints/` package. Property-based tests
  (Hypothesis) compare each rewritten function against it bit for bit, including dtypes.
- **Test data.** At every PR, the chr22 test data (bigWig, BAM, and `--annotation-gtf`) at 1 and 4 threads gives
  the same peak files as `main`.
- **Real data.** 19 datasets: PRO-cap and GRO-cap, bigWig and BAM, replicates, input controls, Ensembl-style
  chromosome names, a tissue batch; 13 genome-wide and 6 partial. At the top of the stack, each with default,
  `--disable-small` and `--sensitive` options, every output matches `main` byte for byte (57 runs at 32 threads;
  an earlier sweep at B6 did the same at 1, 32 and 64 threads). Every timing run in the descriptions also compares
  output hashes.
- **Benchmark kit.** Timings use a small harness ([`bench/`](https://github.com/adamyhe/PINTS/tree/bench/external/bench)
  on our fork, not part of any PR) that alternates the two versions being compared and hashes the outputs.

### Changes users will notice

- **Install:** removes pybedtools, bedtools, bgzip, tabix, statsmodels, requests and Biopython. matplotlib becomes
  optional (`[plots]`) and Numba is offered as `[fast]`. The bioconda recipe would need the same changes at the next
  release.
- **New options:** the four flags in B7 and `--fast-math` in B10, all off by default. `PINTS_PHASE_TIMING=1` prints
  how long each stage takes.
- **Behaviour:** `.npy` coverage dumps are no longer written or reused (B5).

### Suggested order, and questions for you

1. Group A in any order (A2 without squashing). Then we rebase the stack onto the new `main` and open B1, B2, …
   one at a time, each after the previous one merges. Does that work for you, or would you prefer the stack
   grouped into fewer, larger PRs?
2. B1 keeps a 4,150-line frozen copy of the release in `tests/`. It is what makes the bit-for-bit claims testable.
   Are you happy to keep it in the repository, or should it be removed once the stack has merged?
3. The B7 flags change results when set. Would you want any of them to become the default in a future release?
4. B10 (`--fast-math`) is the only change that alters numbers. Do you want it at all, or would you rather leave
   it out?

🤖 Generated with [Claude Code](https://claude.com/claude-code)
