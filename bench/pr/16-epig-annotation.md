## `--epig-annotation`: use local bigBed files directly, skip the annotation when the server lookup fails, sort the marks

**Runs without `--epig-annotation` are unchanged.** With it, a run no longer crashes when the PINTS web server cannot
be reached, a local bigBed file works without the server, and the annotated output is the same from run to run.
Two commits, independent of every other open PR.

### What changes for users

- **A path to a local bigBed file** given to `--epig-annotation` is used as is, without contacting
  `pints.yulab.org`. The docstring of `get_epig_annotation` already said a path was accepted, but every value was
  sent to the server first, so a local file failed whenever the server did.
- **When the lookup or the download fails**, or the server has no annotation for the biosample, PINTS logs a
  warning ending in "skipping epigenomic annotation" and writes the peak calls without the annotation. Before,
  `get_epig_annotation` returned `""` and `annotate_tre_with_epig_info` then raised `FileNotFoundError`, so the run
  stopped after pairing. `pints.yulab.org` has returned `502 Bad Gateway` since at least 2026-11-15, so at the
  moment every run with `--epig-annotation k562` (or another biosample name) stops this way.
- **Each peak's marks are written in sorted order** (for example always `CTCF,DNase`, where a run
  could write either order). They came from a Python `set`, whose order depends on the process's hash seed: two identical runs on
  the chr22 test data wrote different text for about a quarter of the annotated peaks. Which marks a peak gets is
  unchanged.
- `--epig-annotation`'s help text mentions both forms (a biosample name or a path) and the fallback.

### What changes

- `c50bd9c` `get_epig_annotation`:
  - returns an existing file path at once;
  - otherwise asks the server as before, and turns every failure into a warning and `""`;
  - restores the previous socket timeout afterwards (the 30 s default stayed set for the whole process after a
    failure);
  - downloads to `<target>.part` and renames it when complete, so an interrupted download is not reused by the
    next run.

  In `peak_calling`, the annotation block runs only when `get_epig_annotation` returned a file. Tests:
  `tests/test_epig_annotation.py`.
- `8dc24d1` `annotate_tre_with_epig_info` joins the marks with `",".join(sorted(cache))`, with a test that runs it
  under three hash seeds.

### Testing

- `tests/test_epig_annotation.py` at `8dc24d1`, 7 passed:
  - a local file is used without network access;
  - a 502 skips the annotation with a warning and restores the timeout;
  - an unknown biosample skips with a warning;
  - an existing download is reused, and a failed one leaves no file;
  - end to end on chr22 with the server down, the peak files equal a run without `--epig-annotation`;
  - the marks are identical and sorted under three hash seeds.
- Full suite at `c50bd9c` with the server returning 502: 21 passed. Before this PR, the three tests in `tests/functional_test.py`
  fail on `main` while the server is down, because `setUp` passes `epig_annotation="k562"`. With the lookup
  answered locally (a reply naming `ENCFF948HPU.bigBed`, which the test links into place), they pass with the
  annotation applied.
- chr22 bigWig data at 4 threads:
  - with `--epig-annotation` and the server down, the peak files equal a run without the option;
  - without the option, they equal `main`'s;
  - with the annotation applied (`tests/datasets/ENCFF948HPU.bigBed`, given as a path, or through a lookup answered
    locally), two runs are byte-identical, and equal to `main`'s output once the marks in `main`'s output are
    sorted. Two runs of `main` differed in 297 of 1,129 bidirectional lines.
- It merges cleanly with the engine PR and with both follow-ups (opt-in fixes, `--fast-math`).

🤖 Generated with [Claude Code](https://claude.com/claude-code)
