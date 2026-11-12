## Fix a crash under NumPy 2.5, three ignored options, `--model Poisson`, one-strand pairing and result dtypes

**With default options, the output is byte-identical to the current release, except in one case:** a run in which
pairs are found on only one strand. There, `main` wrote malformed bidirectional and divergent files (peak ends in
the chromosome column); this branch writes well-formed ones. Apart from that, results change only for runs that set
`--focused-chrom`, `--stringent-pairs-only` or `--peak-width`, which had no effect before and now do what their
help text says. `--model Poisson` runs, which used to crash, now run.

Six independent fixes, one commit each, plus a regression test for the first. The rows below are in commit order;
please review the commits one by one.

| Commit | Problem | Fix |
|---|---|---|
| Filtering: skip empty per-chromosome frames when concatenating results | A chromosome with candidate peaks but no testable peak returns an empty result frame, and all its columns have object dtype. `peaks_single_strand` concatenated it with the other chromosomes' frames, so every integer column of the strand table (`start`, `end`, `reads`, `summit`, `ler_*`, `non_zeros`) became object. `independent_filtering` then ran its 200 quantiles per call on arrays of Python ints | `_concat_result_frames` leaves out the empty frames before `pd.concat`. They contribute no rows, so the values and their order are unchanged. If every frame is empty, it concatenates them as before |
| Add a regression test for dtype-preserving concatenation of result frames | | `tests/test_concat_result_frames.py` |
| Pass --focused-chrom, --stringent-pairs-only and --peak-width on to the engine | The CLI stores them as `focused_chrom`, `stringent_pairs_only` and `peak_width`. The engine read `highlight_chromosome` and `stringent_only`, and never used `peak_width`. So the `--annotation-gtf` dry run always used chr1, the opposite peak was never required to be significant, and `find_peaks` got no width | The engine accepts both spellings, so API callers can keep using `highlight_chromosome` and `stringent_only`. `peak_width` is passed to `find_peaks` only when it is set |
| Implement the likelihood-ratio test for the Poisson model | `Poisson` did not implement `StatModel`'s abstract `lrt()`, so every `--model Poisson` run raised `TypeError` | `lrt()` compares separate Poisson rates with one shared rate (chi-square, 1 df), the same way ZIP's does for its two parameters |
| Warn when a bigWig holds non-integer values | Normalised (CPM/RPM) bigWigs were cast to integers without any message | A warning names the file and chromosome. The values are still truncated, so results do not change |
| Merge pairs when only one strand has any | When pairs were found on one strand only, `pair()` used that strand's raw 16-column file as the merged result. It then parsed it with the 8 column names of a merged file, so the bidirectional and divergent files had peak ends in the chromosome column | That strand's pairs are now merged on their own, as in the two-strand case |
| Compute get_elbow's cross product without np.cross on 2-D vectors | NumPy 2.5 removed `np.cross` for 2-D vectors. A fresh `pip install` gets NumPy 2.5, and then every peak calling run fails with `ValueError: Both input arrays must be (arrays of) 3-dimensional vectors` | The z component is written out with the same arithmetic NumPy used, so elbows are bit-for-bit unchanged |

### What changes for users

- **NumPy 2.5**: peak calling no longer fails with the `np.cross` error above.
- **`--focused-chrom X` with `--annotation-gtf`**: the dry run that learns `--donor-tolerance` now uses chromosome X
  instead of chr1. Previously, if the data had no chr1, the dry run was skipped without any message and the default
  tolerance was used. With the default (`chr1`), nothing changes.
- **`--stringent-pairs-only`**: peaks now pair only with a significant opposite peak, as the help text says.
  Previously the option had no effect.
- **`--peak-width`**: peaks now have to meet the width requirement in `find_peaks`. Previously the option had no
  effect.
- **Pairs on one strand only**: the bidirectional and divergent files are now well-formed BED (see the lead).
- **`--model Poisson`** runs instead of raising `TypeError`.
- **Normalised bigWig input**: a warning is logged once per chromosome that holds non-integer values, so a
  normalised track gives one warning line per chromosome. Results are unchanged.

### Limitations

- `--model Poisson` is covered only by a unit test of `lrt()`. It was not compared end to end with `main`, which
  cannot run it.

### Testing

- `tests/test_concat_result_frames.py` (3 tests) checks:
  - integer columns stay integer when empty frames are mixed in (the empty frame is built the same way as in
    `check_window_chromosome`);
  - the values and row order equal those of a plain `pd.concat`;
  - an all-empty input still returns an empty frame with the original columns.
- `tests/test_bugfixes.py` (6 tests) checks:
  - each of the three options reaches the engine (or `find_peaks`);
  - Poisson's `lrt()`;
  - the non-integer bigWig warning;
  - pairs found on one strand only.
- `tests/test_elbow.py` (2 tests) checks the elbow of a known curve, and that the result equals `np.cross` on 2-D
  vectors wherever NumPy still provides it.
- Full suite at `558bece`: 15 passed.
- chr22 test data (`tests/datasets`), default options, bigWig and BAM input, at 1 and 4 threads: bidirectional,
  divergent and unidirectional peaks are byte-identical to `main`.
- `--annotation-gtf --focused-chrom chr22` on the chr22 data, with `main` given the chromosome under the name its
  engine reads:
  - `main` and this branch both load 2917 annotated TSSs;
  - both learn a donor tolerance of 0.54;
  - all `get_elbow` calls have identical inputs and outputs;
  - the peak files are byte-identical.

  Unmodified `main` skips the dry run, because chr22 data has no chr1.
- Fresh environments with NumPy 2.5.4, SciPy 1.18.1, pysam 0.24.1 and pyBigWig 0.3.29, once with pandas 3.1.1 and
  once with pandas 2.3.3: `main` fails on the chr22 bigWig and BAM data with `ValueError: Both input arrays must be
  (arrays of) 3-dimensional vectors, but they are 2 and 2 dimensional instead.`; this branch (`558bece`) completes
  in all four runs.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
