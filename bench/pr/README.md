# Upstream PR plan

The fixes (group A) are independent of each other and of the performance stack, so they can merge in any order.
Open the stack once group A has merged, after rebasing it onto the new `main` (the rebase drops the stack's
copies of PRs 1 and 0). Each stack PR builds on the one before it. With default options, every branch produces
output byte-identical to `main`.

Branch names keep the numbers they had while being built. The "Upstream" column is the order to open them in.

## Group A: fixes on `main`

| Upstream | Branch | What |
|---|---|---|
| A1 | `pr/01-empty-frame-dtype` | Filtering skips empty per-chromosome frames when concatenating results, so the result columns keep their dtypes |
| A2 | `pr/02-trackb-fixes` | NumPy 2.5 crash; `--focused-chrom`, `--stringent-pairs-only` and `--peak-width` were ignored; `--model Poisson`; non-integer bigWig warning; one-strand pairing |
| A3 | `pr/09-drop-biopython` | setup.py lists one dependency per line; Biopython replaced with `pysam.FastxFile` |
| A4 | `pr/10-script-fixes` | `pints_normalizer`, `pints_visualizer`, `pints_counter`, `pints_boundary_extender` and `parse_gtf` fixes |
| A5 | `pr/13-pandas4` | pandas 4 deprecations on `main`'s code |

## Group B: performance stack

| Upstream | Branch | Base | What |
|---|---|---|---|
| B1 | `pr/03-test-harness` | A2 (+ A1, setup.py commit) | Frozen reference copy `pints_ref`, fixtures, characterisation tests; no change to `pints/` |
| B2 | `pr/03a-gate1-kernels` | B1 | Exact statistical kernels (ZIP EM, Poisson p-values, independent filtering), optional numba |
| B3 | `pr/03b-gate1-engine` | B2 | Parallel engine: BAM/bigWig parsing, chunked evaluation, persistent pool, bounded memory |
| B4 | `pr/04-io` | B3 | I/O and memory |
| B5 | `pr/05-sparse-coverage` | B4 | In-memory sparse coverage store instead of per-base `.npy` dumps |
| B6 | `pr/06-bedtools-free` | B5 | BED sort/merge in Python; drops pybedtools and the bedtools, bgzip and tabix binaries |
| B7 | `pr/08-trackb-flags` | B6 | Opt-in flags: `--iqr-scale`, `--fix-ler-flanks`, `--strict-independent-filtering`, `--stringent-pairs-padj` |
| B8 | `pr/12-fewer-deps` | B7 (+ A3) | Drops statsmodels and requests; matplotlib becomes the optional `plots` extra |
| B9 | `pr/15-serial-steps` | B8 | Steps after evaluation on the pool: pairing by chromosome, one quantile call, parallel table reads and strand merges; window-test chunks dispatched by subpeak width |
| B10 | `pr/07-gate2` | B9 | Opt-in `PINTS_FAST_MATH=1`: compiled ZIP fit (results change within a stated tolerance when enabled) |

Gate 2 comes last because it is the most likely to be debated, and nothing else depends on it.

## History

The restructure on 2026-11-04 made these changes:
- moved Gate 2 from after PR 6 to the end;
- split the test harness out of PR 3a;
- folded the two stack-only pandas 4 lines (previously `pr/14-pandas4-stack`) into PRs 4 and 6;
- merged `pr/16-parallel-merge` into PR 15 (and, on 2026-11-08, `pr/17-eval-balance`);
- merged `pr/10-normalizer-cli` and `pr/11-script-fixes` into `pr/10-script-fixes`;
- dropped `pr/00-setup-one-per-line` as a separate PR. Its commit is in A3 and in the stack.

Old heads are kept as local tags `archive/2026-11-04/<branch>`.
