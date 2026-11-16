## Faster peak calling with less memory and fewer dependencies (same output)

**With default options, the peak files are byte-identical to `main`'s.** This PR makes `pints_caller` 10–15×
faster at 16 threads and 20–27× at 64 threads on whole-genome data, uses 38–85% less peak memory, stops writing
the per-base `.npy` coverage files (about 37 GB on a genome-wide BAM run), and drops pybedtools, the bedtools,
bgzip and tabix executables, statsmodels and requests from the install. It is one PR because its parts rewrite
the same functions in turn; the review order below follows those parts, each with its own detailed description.

Its first commit is the test harness from #34. Two follow-ups will be opened after it merges, because they change
code this PR rewrites: the opt-in fixes (`--iqr-scale`, `--fix-ler-flanks`, `--strict-independent-filtering`,
`--stringent-pairs-padj`) and `--fast-math`, which do not depend on each other. The `--epig-annotation` fix (#35) is independent: it merges
cleanly with this PR in either order.

### What changes for users

- **Speed.** `--thread` now also parallelises BAM and bigWig parsing, element definition, normalisation, replicate
  merging, filtering and pairing; with `--thread 1` everything runs in the main process. PINTS runs on one thread
  unless `--thread` is given.
- **Memory and disk.** Coverage is held in an in-memory sparse store instead of per-base `.npy` files.
- **`.npy` dumps from an earlier run are no longer reused.** `main` reuses dumps it finds in the output directory
  under the same prefix ("Loading cache from previous result"), which can silently pick up a different input's
  coverage. `PINTS_COVERAGE_STORE=0` keeps the dump-based path with that reuse; it gives the same results.
- **Install.** `install_requires` becomes numpy, pandas, scipy, pysam and pyBigWig. Two optional extras:
  `pip install 'pyPINTS[fast]'` (Numba, for the compiled exact kernels; the output is the same without it) and
  `pip install 'pyPINTS[plots]'` (matplotlib, needed only by `--debug`, which now says so before any work instead
  of failing at the end). The update check uses `urllib`. The bioconda recipe would need the same changes.
- **Files.** No `pybdt_*` temporary directory. The final peak BEDs are unchanged; two intermediates that are
  deleted at the end of the run change form.
- **API.** `pints.calling_engine.housekeeping()` loses its `pybedtools_prefix` parameter.
- **Diagnostics.** `PINTS_PHASE_TIMING=1` prints how long each stage takes to stderr (README).
- The non-integer bigWig warning is logged once per chromosome.

### Why

On whole-genome data, `main` spends most of its time in per-window NumPy/SciPy calls on one process and gains almost
nothing from `--thread` (635 s at 16 threads and 572 s at 64 on a two-replicate BAM). It also writes and rereads
per-base coverage dumps, which dominate peak memory and disk use, and needs external bedtools binaries that make
`pip install pyPINTS` incomplete on its own.

### How to review

Every commit keeps the suite passing and the output unchanged, so the PR can be read commit by commit. The two
groundwork commits change nothing on their own: they let the two follow-ups each add their options without touching the same
lines, so either can merge without the other. The parts, in
order (each link is that part's full description, with its timing, risks and tests):

| Part | Commits | What it does | Wall time vs the part before (16 / 64 threads) |
|---|---|---|---|
| [Test harness](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/03-test-harness.md) | `9c08c15` | Frozen copy of the release as `tests/reference/pints_ref` (imports renamed only, plus #30's one-strand pairing fix), characterisation tests. No change to `pints/` | – |
| [Exact kernels](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/03a-gate1-kernels.md) | `0f48867`..`d9b05a0` (5) | `ZIP.fit`, the small-peak p-values and independent filtering without per-call overhead; bit-identical values | 0.56–0.74 at 16 |
| [Parallel engine](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/03b-gate1-engine.md) | `2c796cf`..`fba112c` (17) | Parsing, element definition and the window tests on one persistent worker pool, in chunks | most of the speed-up (with the kernels, 8–34× over the release) |
| [I/O and memory](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/04-io.md) | `d4933c0`..`a9f858b` (13) | Faster parsing, less copying, smaller intermediates, pool start-up | 0.61–0.75 |
| [Sparse coverage](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/05-sparse-coverage.md) | `1033423`..`47c75fe` (10) | In-memory sparse coverage store instead of `.npy` files | 0.85–0.94 |
| [No bedtools](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/06-bedtools-free.md) | `254cf79`, `6dffa5e` | BED sort, merge, intersect and concatenation in Python (`pints/bed_ops.py`) | 0.94–0.99 |
| [Fewer dependencies](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/12-fewer-deps.md) | `14570dd`..`fad2b8d` (4) | statsmodels' BH adjustment repeated in numpy (bit-identical), requests → `urllib`, matplotlib → `plots` extra | – |
| [Steps after evaluation](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/15-serial-steps.md) | `9f239e3`..`a86de61` (3) | Filtering, pairing and strand merging on the pool; window-test chunks dispatched by size | 0.88–0.98 |
| Groundwork for the follow-ups | `3f2c70a`, `f45ef27` | The LER candidate background in one helper (the code `bgIQR` and `pkIQR` each had inline), and `stats_engine.register_option`: options of the statistical tests registered next to their code, applied by `peak_calling` and passed to pool workers. No option is registered here | – |
| Fixes from review | `5a3544f`..`7d7f3b5` (5) | The tabix fallback tests subpeaks in the sorted file's order, as before; phase-1 results collected in one pass (was units × tasks); the control is normalised before the pool forks on the `PINTS_COVERAGE_STORE=0` path too; `pints_counter` reads gzip region files and skips header lines again; the update check does not retry HTTP errors; tidy-up | – |

The linked descriptions were written when these parts were separate, stacked PRs, so their "stacked on" lines and
commit hashes refer to that layout; the commits themselves are the same.

Size: `pints/`, `scripts/`, `setup.py` and the README +3,608 / −473 lines (two new modules, `coverage_store.py`
and `bed_ops.py`); `tests/` +8,121, of which 4,151 are the frozen reference copy.

### Risks and limitations

- **Platform.** The parallel bigWig loader, input normalisation and replicate merging use
  `multiprocessing.get_context("fork")`, so with `--thread` > 1 they assume Linux or macOS. All testing was on Linux.
- **Memory grows with the number of distinct coverage runs**, not with genome length. The benchmark inputs carry
  signal on 0.02–0.07% of bases; a dense input whose value changes at most bases approaches per-base storage and
  is held in memory rather than on disk.
- **Tie order in BED sorting** replays bedtools 2.31 built with libstdc++. A bedtools built against another C++
  library could order ties differently from `main` as tested.
- **Test dependencies.** The tests need pytest, hypothesis and numba (the `tests/conftest.py` docstring), plus
  pybedtools and a bedtools executable for the reference implementation and the functional tests, and statsmodels
  for the p-value adjustment comparison.

### Timing and memory

`main` against this code at the top of the earlier stack (with the follow-ups' options off), default options, on an
AMD EPYC 9554 server; output files identical in every run. Another user's jobs shared the machine, so absolute times
are somewhat high. Benchmark kit: https://github.com/adamyhe/PINTS/tree/bench/external/bench (`bench/variants.json`
pins each measured variant to a commit, `bench/datasets.json` describes the datasets).

| Dataset (genome-wide) | `main`, 16 / 64 threads | This PR, 16 / 64 threads | Peak memory, 16 threads |
|---|---|---|---|
| PRO-cap BAM, 2 replicates (Karpas-422) | 635 / 572 s | 52 / 28 s | 12.4 → 3.5 GB |
| PRO-cap bigWig (Caco-2) | 585 / 666 s | 56 / 28 s | 6.0 → 3.4 GB |
| GRO-cap BAM with input control (GM12878) | 723 / 727 s | 61 / 34 s | 12.6 → 4.0 GB |
| PRO-cap bigWig (A549) | 383 / 386 s | 25 / 15 s | 5.7 → 2.3 GB |

### Testing

- Output equivalence: 19 real datasets (PRO-cap and GRO-cap, bigWig and BAM, replicates, input controls,
  Ensembl-style chromosome names, a tissue batch; 13 genome-wide), each with default, `--disable-small` and
  `--sensitive` options: every output matched `main` byte for byte (57 runs at 32 threads, on the earlier stack's
  top, which contains this PR's code).
- Property-based tests (Hypothesis) compare each rewritten function with the frozen `pints_ref` bit for bit,
  including dtypes.
- At `3e4c610`, on `main` (`2b72d87`), 2026-11-16: full suite 458 passed (Python 3.12.13, NumPy 2.1.0, pandas 2.3.1,
  SciPy 1.15.0, Numba 0.61.2). `pints.yulab.org` was returning 502, so the three tests in `tests/functional_test.py`
  ran with the `k562` annotation lookup answered locally (same bigBed file); #35 removes that dependency.
- chr22 test data, with bigWig, BAM, BAM with an input control (`--ct-bam tests/datasets/input.chr22.bam`), and
  `--annotation-gtf --focused-chrom chr22` input, at 1 and 4 threads: the bidirectional, divergent and
  unidirectional peaks are byte-identical to `main`'s. With `PINTS_COVERAGE_STORE=0` (bigWig, and BAM with the
  control, at 1 and 4 threads) they are too, with no fork warnings in the log.
- The fallback when the in-memory subpeak index is not supported: the new test shuffles the element-definition
  files so that it is taken, and compares with `pints_ref`. It failed before `5a3544f` and passes now.
- `pints_visualizer`, `pints_counter`, `pints_boundary_extender` and `pints_sample_qc` on chr22: all 24 output files
  are byte-identical to `main`'s. The one run that fails on `main` (`pints_counter` with only `-b`) succeeds here.
- Both follow-ups on top of this branch (opt-in fixes, `--fast-math`) give byte-identical peaks to the earlier
  stack's top for each opt-in option, `--fast-math`, and their combinations.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
