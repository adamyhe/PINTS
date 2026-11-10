## Opt-in fixes: `--iqr-scale`, `--fix-ler-flanks`, `--strict-independent-filtering`, `--stringent-pairs-padj`

Stacked on the bedtools-free PR (`pr/06-bedtools-free`); review only the last commit.

This adds four flags to `pints_caller`. Each corrects a place where the code differs from the PINTS paper, the help
text, or the method it cites. **All four are off by default, and with none of them set, the output is byte-identical
to the current release.** Each flag changes only what it names.

| Flag | What it changes | Default (current behaviour) |
|---|---|---|
| `--iqr-scale k` | Fence multiplier of the IQR outlier rule in local environment refinement (Q1 − k·IQR, Q3 + k·IQR) | `1.0`; the paper and Tukey's rule use 1.5 |
| `--fix-ler-flanks` | The 2.5 kb background of a refinement candidate is clipped to the chromosome, and the tested peak is masked wherever it lies | The background is empty within 2.5 kb of a chromosome start, and the mask offsets are right only when the tested peak lies inside the candidate's flank |
| `--strict-independent-filtering` | Windows removed by independent filtering get padj = 1 (Bourgon et al. 2010) | Removed windows keep their padj from BH over all windows, so some of them can still be called significant |
| `--stringent-pairs-padj` | Like `--stringent-pairs-only`, but the opposite peak must pass by padj instead of its raw p-value | Not available; `--stringent-pairs-only` uses the raw p-value |

### Why opt-in rather than default

Turning any of these on changes which peaks PINTS reports, and published analyses and pipelines depend on the current
results. Whether a correction becomes the default is a decision for the maintainers, ideally made with a version bump and a
changelog note. This PR keeps that decision separate: it makes the corrected behaviour available and tested without
changing anyone's results.

Effect on the chr22 test data (`tests/datasets/{pl,mn}.chr22.bw`; bidirectional / divergent / unidirectional peaks):

| Setting | Peaks |
|---|---|
| none (= current release) | 1047 / 986 / 130 |
| `--iqr-scale 1.5` | 1035 / 976 / 129 |
| `--fix-ler-flanks` | 1047 / 986 / 130 (unchanged here; covered by unit tests) |
| `--strict-independent-filtering` | 1047 / 986 / 130 (unchanged here; covered by unit tests) |
| `--stringent-pairs-padj` | 415 / 398 / 907 |
| all four | 416 / 399 / 892 |

Suggested defaults, for the maintainers to decide:

- **`--fix-ler-flanks`**: a plain bug (a wrong or empty background window). A good candidate for the default in the next
  minor release.
- **`--strict-independent-filtering`**: matches the method PINTS cites. Also a good candidate.
- **`--iqr-scale`**: depends on whether the paper (1.5) or the released code (1.0) is meant to be the reference. Worth
  deciding explicitly, since it changes calls.
- **`--stringent-pairs-padj`**: a stricter definition rather than a fix. Best kept opt-in.

### Implementation

- `stats_engine.FIXES` holds the two settings used inside the statistical tests (`iqr_scale`, `fix_ler_flanks`).
  `configure_fixes` sets them in `peak_calling`, and the pool initializer passes them to worker processes, so they
  don't rely on inheritance at fork.
- `IQR.get_outlier_threshold(scale=None)` reads `FIXES["iqr_scale"]` unless given a scale.
- `ler_candidate_background` replaces the two inline copies of the flank code. With the flag off it is the previous
  code.
- `independent_filtering(..., strict_independent_filtering=...)` and `merge_opposite_peaks(..., stringent_pairs_padj=...)`
  take the other two as keyword arguments. `--stringent-pairs-padj` implies `--stringent-pairs-only`.

### Timing

The benchmark kit is [`bench/`](https://github.com/adamyhe/PINTS/tree/bench/external/bench) on the fork's `bench/external` branch: the variant names below are pinned
to commits in [`bench/variants.json`](https://github.com/adamyhe/PINTS/blob/bench/external/bench/variants.json), and the datasets are described in
[`bench/datasets.json`](https://github.com/adamyhe/PINTS/blob/bench/external/bench/datasets.json).

With no flags set, the code paths are the previous ones. `bench/run ab --a x06 --b x08 --fresh`
(`pr/06-bedtools-free` at `efe3a3f` against this branch at `7da8a80`) on cbsubscb19 (AMD EPYC 9554) on 2026-11-09,
three alternating repetitions per cell, median wall time, with another user's single-core jobs on the node's other
cores. Identical output in every cell:

| Dataset | 16 threads | 64 threads |
|---|---|---|
| `H1_bam_reps_karpas422` | 53.7 → 53.6 s (0.998) | 30.9 → 30.9 s (1.000) |
| `H3_bam_ct_gm12878_grocap` | 72.4 → 68.9 s (0.952) | 41.2 → 41.3 s (1.003) |
| `T4_bam_reps_hct116` | 29.6 → 29.5 s (0.999) | 16.8 → 16.9 s (1.007) |
| `T5_bam_ct_k562_grocap` | 11.7 → 11.8 s (1.009) | 7.9 → 8.2 s (1.030) |

### Testing

- `tests/test_trackb_flags.py` (11 tests): CLI defaults and flags; workers receive the settings; IQR fences; LER
  flanks near the chromosome start and end and with the tested peak inside or outside the flank, plus a check that
  the flag-off path equals the previous code; strict independent filtering; padj-based pairing.
- Full suite: 453 passed.
- chr22, bigWig and BAM input: with no flags, output is byte-identical to `main`. With each flag, and with all of
  them, output is byte-identical to the reference implementation of these fixes.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
