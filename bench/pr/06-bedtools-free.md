## Drop pybedtools and the bedtools, bgzip and tabix binaries

Replaces PINTS's few bedtools operations with a small pure-Python module. **The output is byte-identical to the
base branch.** PINTS no longer needs pybedtools, or the bedtools, bgzip and tabix executables, at run time.
Stacked on B5 (`pr/05-sparse-coverage`); review the last two commits.

### What changes for users

- **Dependencies.** `setup.py` no longer lists `pybedtools`, and `pints_caller` no longer exits at startup when
  `bgzip`, `tabix` or `bedtools` is not on the `PATH`.
- **Files.** PINTS no longer creates the `pybdt_*` temporary directory in the output directory. The `.gz`
  intermediates are still BGZF files written by pysam (`pysam.tabix_compress`, which `BedTool.bgzip` also called)
  with CSI indexes (`pysam.tabix_index(..., csi=True)`), as before. The final peak BEDs are plain text and
  unchanged.
- **API.** `pints.calling_engine.housekeeping` loses its `pybedtools_prefix` parameter: it is now
  `housekeeping()`. Code that passed the prefix must drop the argument.

### Why

pybedtools is a wrapper around the bedtools executable, which users have to install separately, and `pints_caller`
exited at startup if `bgzip`, `tabix` or `bedtools` was not on the `PATH`. That made a `pip install` of PINTS
incomplete on its own, and pybedtools itself can be awkward to build. PINTS used it only for simple operations on
small BED files: sort, merge, concatenate and intersect. Compression and indexing already went through pysam.

pybedtools also writes temporary files into a `pybdt_*` directory in the output directory: an H1 (two BAM
replicates) or H3 (BAM with control) run at 64 threads created 7–8 such files (8–24 MB; profile details under
Timing).

### What changes

- **`pints/bed_ops.py` (new)** reproduces what PINTS got from pybedtools and bedtools 2.31, on lists of lines:
  - `bedtools sort`: chromosomes in byte order, then starts. For tied starts, bedtools' order comes from
    libstdc++'s `std::sort`, which is not stable; `sort_lines` replays that algorithm (introsort with its
    heapsort fallback and insertion-sort threshold), so ties come out in the same order. GTF records are sorted on
    their start column; the format is told from the records, as bedtools does, not from the file name.
  - `bedtools merge`, with `-s`, `-d` and the `collapse` and `distinct` operations;
  - `bedtools intersect` with `-v`, `-u` or `-c`, optionally `-s`;
  - `BedTool.cat`, with and without `postmerge`.
- **The default path uses it**: `pair()`, `index_bed_file` (sort, then `pysam.tabix_compress` and
  `pysam.tabix_index`) and the emptiness check of the pairing outputs.
- **So do the other paths that used pybedtools:**
  - borrowing elements from other replicates (`inferring_elements_from_other_reps`, used with `--dont-merge-reps`
    unless `--dont-borrow-info-reps` is given);
  - merging the relaxed calls of `--epig-annotation`;
  - the TSS windows of the alpha dry run (`--annotation-gtf`);
  - the promoter filter of `pints_boundary_extender --promoter-bed`;
  - `pints_counter`'s union of regions.
- **The startup check for `bgzip`, `tabix` and `bedtools` is removed**, and PINTS no longer creates the `pybdt_*`
  directory.
- **`setup.py` no longer lists `pybedtools`.** pybedtools and bedtools stay test dependencies: the reference
  implementation `pints_ref` and the functional tests use them.

One difference from bedtools: `bedtools merge -s` reorders its output through a heap. `bed_ops` returns the same
set of merged records, in sorted order. Only the dry run uses strand-aware merge, and its statistics do not depend
on the order of the windows.

Next to the rewritten line in `cut_peaks_dry_run`, `drop_duplicates` assigns its result instead of working in
place (`inplace=True` is deprecated in pandas 3.1 and removed in pandas 4; the frame is local).

### Risks

- **Order of tied starts.** `sort_lines` replays the order that bedtools 2.31, built with libstdc++, gives records
  with the same chromosome and start. On `main`, that order comes from whichever bedtools build is installed, so
  a bedtools built against another C++ library could order ties differently from both `main` as tested and this
  branch. The comparison tests ran against bedtools 2.31 (named in `tests/unit/test_bed_ops.py`'s docstring); the
  exact build was not recorded in the repository.
- **Test dependencies are not declared.** The reference implementation, the functional tests and an autouse
  fixture in `tests/conftest.py` import pybedtools, but nothing installs it, and the conftest docstring lists only
  `pip install pytest hypothesis numba`. To run the tests, also install pybedtools and a bedtools executable, for
  example `pip install pybedtools` plus bedtools 2.31 from bioconda (`conda install -c bioconda bedtools=2.31`).

### Timing

The benchmark kit is [`bench/`](https://github.com/adamyhe/PINTS/tree/bench/external/bench) on the fork's
`bench/external` branch: [`bench/variants.json`](https://github.com/adamyhe/PINTS/blob/bench/external/bench/variants.json)
pins each measured variant to a commit, and
[`bench/datasets.json`](https://github.com/adamyhe/PINTS/blob/bench/external/bench/datasets.json) describes the
datasets. PINTS runs on one thread unless `--thread` is given; both columns use `--thread`.

`bench/run ab --fresh` on 2026-11-08, B5 against this branch (earlier heads with the same code, see provenance),
on an AMD EPYC 9554 server (64 CPUs), three alternating repetitions per cell, median wall time, identical output in
every cell:

| Dataset | 16 threads | 64 threads |
|---|---|---|
| `H1_bam_reps_karpas422` | 43.7 → 42.9 s (0.981) | 28.0 → 27.4 s (0.976) |
| `H3_bam_ct_gm12878_grocap` | 56.7 → 55.4 s (0.977) | 37.8 → 36.6 s (0.968) |
| `T4_bam_reps_hct116` | 23.6 → 23.5 s (0.994) | 15.4 → 15.0 s (0.975) |
| `T5_bam_ct_k562_grocap` | 10.2 → 9.7 s (0.951) | 7.8 → 7.3 s (0.941) |

The gain is in filtering, pairing and QC (H3 at 16 threads: 13.2 → 12.6 s), where the BED operations no longer
go through temporary files and a bedtools process, and in the input-parsing phase, which is 0.2–0.3 s shorter in
every cell. The benchmark times that phase from process launch, so it includes module imports and the startup
check. On B5, `pints/calling_engine.py` imports pybedtools when it is loaded, and the startup check looks up the
three executables; both are gone here. Whether they account for the 0.2–0.3 s was not traced.

<details><summary>Provenance of the measured commits</summary>

- Timing run: B5 at `80a791b` against this branch at `6171301`. They have the same `pints/`, `scripts/`,
  `setup.py` and `README.md` as `769beba` and `efe3a3f`, the previous heads of B5 and of this branch, and those
  differ from the current heads (`0e60047`, `0a21e1f`) only in comments, test names and B3's phase-timer fix,
  which affects only stderr when `PINTS_PHASE_TIMING` is set.
- The `pybdt_*` counts come from a profile of `6d2056c`, an earlier version of B4 (B4's code before the rebase onto
  A1). In the profile of `adada92` (an earlier version of this branch's first commit, on the store), the H1 and H3
  runs at 64 threads created no `pybdt_*` files.
- Benchmark equivalence of `ca581dd`: all 78 runs (19 datasets; default, `--disable-small` and sensitive
  configurations; 1, 32 and 64 threads) are byte-identical to `main`. Its `pints/` code differs from this branch
  only by the pandas line above, the same line in stratified filtering, and B5's commit that stops reading old
  `.npy` dumps, which does nothing in a fresh output directory.
- The full-suite count below was taken at `efe3a3f`.
- These commits live under `refs/archive/*` on the fork (github.com/adamyhe/PINTS).

</details>

### Testing

- `tests/unit/test_bed_ops.py`: sort, merge (with `collapse` and `distinct`) and `merge -d` against bedtools on
  60 random sets of records with many tied starts; strand-aware merge and `intersect -v/-u/-c` (with and without
  `-s`) against bedtools on 60 more, including zero-length records, `.` strands and an unsorted `-a`; `cat`
  against pybedtools on 20; the heapsort fallback and insertion-sort thresholds; GTF sorting; reading of headers
  and gzipped files.
- `tests/diff/test_differential_bedtools_free.py`: borrowing elements from other replicates, boundary extension
  with a promoter BED (plain and gzipped) and `pints_counter`'s regions (`--max-dist` 0, 1 and 250) give the same
  files as the reference implementation, which uses pybedtools.
- Full suite: 423 passed (265 on B5). Most of the new cases are parametrised: `test_bed_ops.py` runs 60 + 60 + 20
  seeded record sets plus five single tests, and `test_differential_bedtools_free.py` 8 seeds, 2 promoter-BED
  forms and 3 `--max-dist` values.
- chr22 test data, with bigWig, BAM, and `--annotation-gtf --focused-chrom chr22` input, at 1 and 4 threads: the
  bidirectional, divergent and unidirectional peaks are byte-identical to `main`.
- Clean install: `pip install` of this branch into a fresh virtual environment (Python 3.12), then `pints_caller` on
  the chr22 bigWigs at 4 threads with a `PATH` holding only the environment and three coreutils, so no bedtools,
  bgzip or tabix can be found. The run completes and writes the three peak files.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
