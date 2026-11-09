## Avoid pandas 4 deprecations; accept pandas 3's `str` dtype in the functional tests

Against `main`. Two commits. Review them commit by commit.

**The output is byte-identical to the current release, under pandas 2.3 and 3.1.**

| Commit | Problem | Fix |
|---|---|---|
| Avoid pandas 4 deprecations: fill split TSS columns with a string, drop `inplace=True` in `peak_bed_to_gtf` | Under pandas 3.1, `pair()` raises a `Pandas4Warning` for each `fillna(value=-1)` on the `str` columns from `str.split(",", expand=True)` (four calls), and `peak_bed_to_gtf` for `sort_values(..., inplace=True)` (PDEP-8). pandas 4 will make both errors | Fill with `"-1"`, which `astype(int)` turns into the same `-1`. Assign the sorted frame instead of sorting in place; the frame is local, so nothing relied on the mutation |
| Functional tests: accept pandas 3's string dtype for text columns | pandas 3 reads text columns with the new `str` dtype instead of `object`, so `test_peak_calling` and `test_boundary_extension` failed their `is_object_dtype` checks, although the output files were correct | `is_string_dtype`, which accepts both. Under pandas 2 it is true for an `object` column whose values are strings |

The other two `inplace=True` calls on `main` (the `sort_values` in `stratified_filtering` and the TSS
`drop_duplicates` in `cut_peaks_dry_run`) sit next to lines that the performance stack rewrites. They are fixed there
instead, so this PR and the stack merge in either order.

### Testing

- chr22 bigWig data (`tests/datasets`), 4 threads, default options, `--annotation-gtf`, and `--debug`: the
  bidirectional, divergent and unidirectional peaks are byte-identical to `main`, under pandas 2.3.1 and 3.1.0.
- `Pandas4Warning`s per run under pandas 3.1.0 go from 50 (51 with `--debug`) on `main` to 2 on this branch. The 2
  left come from `stratified_filtering`'s `sort_values`. Under pandas 2.3.1 there are none on either branch.
- `tests/functional_test.py` (3 tests) under pandas 3.1.0: 2 fail on `main` with "dtype for the Chromosome column
  should be object"; all 3 pass on this branch. They also pass under pandas 2.3.1, and on the performance stack's
  tip merged with this branch under both versions.
- `main`'s test suite passes on the first commit under pandas 2.3.1 (4 passed).
- `--debug` also writes `t_1_peaks.gtf.gz`, the one file the `peak_bed_to_gtf` change touches. Decompressed and
  without the version tag, it is identical to `main`'s (58,338 lines) under pandas 2.3.1 and 3.1.0.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
