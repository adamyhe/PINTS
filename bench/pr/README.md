# Upstream PR plan

The fixes (group A) are independent of each other and of the performance stack, so they can merge in any order.
Open the stack once group A has merged, after rebasing it onto the new `main` (the rebase drops the stack's
copy of the setup.py commit). Each stack PR builds on the one before it. With default options, every branch produces
output byte-identical to `main`.

Branch names keep the numbers they had while being built. The "Upstream" column is the order to open them in.

## Group A: fixes on `main`

| Upstream | Branch | What |
|---|---|---|
| A1 | `pr/02-trackb-fixes` | Result frames keep their dtypes; NumPy 2.5 crash; `--focused-chrom`, `--stringent-pairs-only` and `--peak-width` were ignored; `--model Poisson`; non-integer bigWig warning; one-strand pairing |
| A2 | `pr/09-drop-biopython` | setup.py lists one dependency per line; Biopython replaced with `pysam.FastxFile` |
| A3 | `pr/10-script-fixes` | `pints_normalizer`, `pints_visualizer`, `pints_counter`, `pints_boundary_extender` and `parse_gtf` fixes |
| A4 | `pr/13-pandas4` | pandas 4 deprecations on `main`'s code; functional tests accept pandas 3's `str` dtype |

## Group B: performance stack

| Upstream | Branch | Base | What |
|---|---|---|---|
| B1 | `pr/03-test-harness` | A1 (+ setup.py commit) | Frozen reference copy `pints_ref`, fixtures, characterisation tests; no change to `pints/` |
| B2 | `pr/03a-gate1-kernels` | B1 | Exact statistical kernels (ZIP EM, Poisson p-values, independent filtering), optional numba |
| B3 | `pr/03b-gate1-engine` | B2 | Parallel engine: BAM/bigWig parsing, chunked evaluation, persistent pool, bounded memory |
| B4 | `pr/04-io` | B3 | I/O, memory and pool follow-ups |
| B5 | `pr/05-sparse-coverage` | B4 | In-memory sparse coverage store instead of per-base `.npy` dumps; old dumps in the output directory are no longer reused |
| B6 | `pr/06-bedtools-free` | B5 | BED sort/merge in Python; drops pybedtools and the bedtools, bgzip and tabix binaries |
| B7 | `pr/08-trackb-flags` | B6 | Opt-in flags: `--iqr-scale`, `--fix-ler-flanks`, `--strict-independent-filtering`, `--stringent-pairs-padj` |
| B8 | `pr/12-fewer-deps` | B7 (+ A2) | Drops statsmodels and requests; matplotlib becomes the optional `plots` extra |
| B9 | `pr/15-serial-steps` | B8 | Steps after evaluation on the pool: pairing by chromosome, one quantile call, parallel table reads and strand merges; window-test chunks dispatched by subpeak width |
| B10 | `pr/07-gate2` | B9 | Opt-in `--fast-math`: compiled ZIP fit (results change within a stated tolerance when enabled) |

Gate 2 comes last because it is the most likely to be debated, and nothing else depends on it.

## Stack tip against `main`

The benchmark kit is [`bench/`](https://github.com/adamyhe/PINTS/tree/bench/external/bench) on the fork's `bench/external` branch: the variant names below are pinned
to commits in [`bench/variants.json`](https://github.com/adamyhe/PINTS/blob/bench/external/bench/variants.json), and the datasets are described in
[`bench/datasets.json`](https://github.com/adamyhe/PINTS/blob/bench/external/bench/datasets.json).

`bench/run ab --a ref --b x07 --reps 1`: `main` (`81b2dd9`) against the tip of the stack (`pr/07-gate2` at
`4d17c9d`, default options, so without `--fast-math`), on cbsubscb19 (AMD EPYC 9554) on 2026-11-09, one run of
each per cell. The three output files are byte-identical in all 16 cells. Another user's single-core jobs ran on
the node's other cores until 20:10, which covers the H rows and the start of the T rows. That load slows both sides,
so these times are higher than on an idle node.

| Dataset | Threads | `main` | Tip | Tip / `main` | Peak PSS |
|---|---|---|---|---|---|
| `H1_bam_reps_karpas422` | 16 | 634.9 s | 52.4 s | 0.083 | 12.42 → 3.46 GB |
| `H1_bam_reps_karpas422` | 64 | 571.9 s | 28.2 s | 0.049 | 12.41 → 4.39 GB |
| `H2_bw_caco2` | 16 | 585.3 s | 55.8 s | 0.095 | 6.03 → 3.38 GB |
| `H2_bw_caco2` | 64 | 666.3 s | 28.1 s | 0.042 | 7.09 → 4.40 GB |
| `H3_bam_ct_gm12878_grocap` | 16 | 723.3 s | 60.5 s | 0.084 | 12.57 → 4.00 GB |
| `H3_bam_ct_gm12878_grocap` | 64 | 726.6 s | 33.9 s | 0.047 | 12.68 → 5.12 GB |
| `H6_bam_reps_bw_a549r2_karpas` | 16 | 557.3 s | 40.0 s | 0.072 | 12.22 → 2.93 GB |
| `H6_bam_reps_bw_a549r2_karpas` | 64 | 451.2 s | 21.1 s | 0.047 | 12.21 → 4.34 GB |
| `T2_bam_k562` (chr1, 21, X) | 16 | 285.9 s | 6.4 s | 0.022 | 6.83 → 1.04 GB |
| `T2_bam_k562` (chr1, 21, X) | 64 | 290.0 s | 5.0 s | 0.017 | 6.84 → 1.43 GB |
| `T3_bw_gw_a549` | 16 | 382.7 s | 25.4 s | 0.066 | 5.68 → 2.34 GB |
| `T3_bw_gw_a549` | 64 | 385.6 s | 14.5 s | 0.038 | 6.96 → 2.98 GB |
| `T4_bam_reps_hct116` (chr17–22) | 16 | 793.7 s | 22.1 s | 0.028 | 7.48 → 2.33 GB |
| `T4_bam_reps_hct116` (chr17–22) | 64 | 788.3 s | 14.3 s | 0.018 | 7.48 → 2.82 GB |
| `T5_bam_ct_k562_grocap` (chr19–22) | 16 | 390.5 s | 9.2 s | 0.024 | 6.91 → 1.23 GB |
| `T5_bam_ct_k562_grocap` (chr19–22) | 64 | 389.0 s | 6.7 s | 0.017 | 6.91 → 1.99 GB |

Unmarked datasets are whole-genome. `main` gains little from more threads, so the speed-up is 10.5–15× at 16 threads
and 20–27× at 64 threads on the whole-genome sets (36–58× on the partial T sets). Peak PSS falls by 38–85%. With `--fast-math`, the tip is a further 0.47–0.77× at 16 threads and
0.63–0.92× at 64 threads ([`bench/pr/07-gate2.md`](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/07-gate2.md)).

## History

The restructure on 2026-11-04 made these changes:
- moved Gate 2 from after PR 6 to the end;
- split the test harness out of PR 3a;
- folded the two stack-only pandas 4 lines (previously `pr/14-pandas4-stack`) into PRs 4 and 6;
- merged `pr/16-parallel-merge` into PR 15 (and, on 2026-11-08, `pr/17-eval-balance`);
- merged `pr/10-normalizer-cli` and `pr/11-script-fixes` into `pr/10-script-fixes`;
- dropped `pr/00-setup-one-per-line` as a separate PR. Its commit is in A3 and in the stack.

Old heads are kept as refs `refs/archive/2026-11-04/<branch>` (not tags: versioneer would read them). All
`refs/archive/*` refs are pushed to the fork (`git fetch origin 'refs/archive/*:refs/archive/*'`), so every commit
a description cites can be opened there; `refs/archive/cited/<sha>` holds benchmarked commits that no branch
contains.

On 2026-11-08 the `.npy` cache was dropped (last commit of PR 5; the stack above it was rebased), and `pr/13` gained the pandas 3 test fix.
Also on 2026-11-08: PR 1 (`pr/01-empty-frame-dtype`) was folded into PR 2, the duplicated pandas 4 commit in PR 4 was
squashed, and PR 3b gained the `--stringent-pairs-only` fix for pooled per-sample pairing. Old heads are kept as
`refs/archive/2026-11-08/<branch>`.

On 2026-11-09: PR 5's element-definition fix (stack above rebased); Gate 2's environment variable became
`--fast-math`; commit messages on PR 6 and Gate 2 reworded; the `Co-authored-by: Biomni` trailers removed from every
branch (trees unchanged; this rewrites PR 2 and the whole stack, and PR 12 again merges `pr/09`'s own commit).
Old heads are kept as `refs/archive/2026-11-09/<branch>` and `refs/archive/2026-11-09-msg/<branch>`.
