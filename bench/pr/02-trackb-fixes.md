## Fix a NumPy 2.5 crash, three ignored options, `--model Poisson` and one-strand pairing

Against `main`. Five independent fixes, one commit each. Review them commit by commit.

**With default options, the output is byte-identical to the current release.** Results change only for runs that
set `--focused-chrom`, `--stringent-pairs-only` or `--peak-width`. Those options had no effect before, and now they
do what their help text says.

| Commit | Problem | Fix |
|---|---|---|
| Compute `get_elbow`'s cross product without `np.cross` on 2-D vectors | NumPy 2.5 removed `np.cross` for 2-D vectors. A fresh `pip install` gets NumPy 2.5, and then every peak calling run fails with `ValueError: Both input arrays must be (arrays of) 3-dimensional vectors` | The z component is written out with the same arithmetic NumPy used, so elbows are bit-for-bit unchanged |
| Pass `--focused-chrom`, `--stringent-pairs-only` and `--peak-width` on to the engine | The CLI stores them as `focused_chrom`, `stringent_pairs_only` and `peak_width`. The engine read `highlight_chromosome` and `stringent_only`, and never used `peak_width`. So the `--annotation-gtf` dry run always used chr1, the opposite peak was never required to be significant, and `find_peaks` got no width | The engine accepts both spellings, so API callers can keep using `highlight_chromosome` and `stringent_only`. `peak_width` is passed to `find_peaks` only when it is set |
| Implement the likelihood-ratio test for the Poisson model | `Poisson` did not implement `StatModel`'s abstract `lrt()`, so every `--model Poisson` run raised `TypeError` | `lrt()` compares separate Poisson rates with one shared rate (chi-square, 1 df), the same way ZIP's does for its two parameters |
| Warn when a bigWig holds non-integer values | Normalised (CPM/RPM) bigWigs were cast to integers without any message | A warning names the file and chromosome. The values are still truncated, so results do not change |
| Merge pairs when only one strand has any | When pairs were found on one strand only, `pair()` used that strand's raw 16-column file as the merged result. It then parsed it with the 8 column names of a merged file, so the bidirectional and divergent files had peak ends in the chromosome column | That strand's pairs are now merged on their own, as in the two-strand case |

### Results that change, and why

Only runs that set one of these options are affected:

- **`--focused-chrom X` with `--annotation-gtf`**: the dry run that learns `--donor-tolerance` now uses chromosome X
  instead of chr1. Previously, if the data had no chr1, the dry run was skipped without any message and the default
  tolerance was used. With the default (`chr1`), nothing changes.
- **`--stringent-pairs-only`**: peaks now pair only with a significant opposite peak, as the help text says.
  Previously the option had no effect.
- **`--peak-width`**: peaks now have to meet the width requirement in `find_peaks`. Previously the option had no
  effect.

The `--model Poisson` and one-strand pairing fixes affect only runs that used to crash or produce malformed output.

### Testing

- `tests/test_trackb_fixes.py` (6 tests) checks:
  - each of the three options reaches the engine (or `find_peaks`);
  - Poisson's `lrt()`;
  - the non-integer bigWig warning;
  - pairs found on one strand only.
- `tests/test_elbow.py` (2 tests) checks the elbow of a known curve, and that the result equals `np.cross` on 2-D
  vectors wherever NumPy still provides it.
- Full suite: 12 passed.
- chr22 test data (`tests/datasets`), default options, bigWig and BAM input: bidirectional, divergent and
  unidirectional peaks are byte-identical to `main`.
- `--annotation-gtf --focused-chrom chr22` on the chr22 data, with `main` given the chromosome under the name its
  engine reads:
  - `main` and this branch both load 2917 annotated TSSs;
  - both learn a donor tolerance of 0.54;
  - all `get_elbow` calls have identical inputs and outputs;
  - the peak files are byte-identical.

  Unmodified `main` skips the dry run, because chr22 data has no chr1.
- NumPy 2.5.4 (with pandas 3.1, SciPy 1.18): `main` fails on the chr22 data with the error above. A later branch
  that includes this commit runs to completion. `tests/test_elbow.py` covers this change directly.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
