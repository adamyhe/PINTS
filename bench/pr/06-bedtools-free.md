## Drop pybedtools and the bedtools, bgzip and tabix binaries

Stacked on the sparse coverage PR (`pr/05-sparse-coverage`); review only the last two commits.

**The output is byte-identical to the base branch.** PINTS no longer needs pybedtools, or the bedtools, bgzip and
tabix executables, at run time.

### Why

pybedtools is a wrapper around the bedtools executable, which users have to install separately, and `pints_caller`
exited at startup if `bgzip`, `tabix` or `bedtools` was not on the `PATH`. That made a `pip install` of PINTS
incomplete on its own, and pybedtools itself can be awkward to build. PINTS used it only for simple operations on
small BED files: sort, merge, concatenate and intersect. Compression and indexing already went through pysam.

pybedtools also writes temporary files into a `pybdt_*` directory in the output directory. In the profile of
`6d2056c` (the code of `pr/04-io`, see the sparse coverage PR), an H1 or H3 run at 64 threads created 7–8 such
files (8–24 MB).

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

### Timing and memory

In the profile of `adada92` (an earlier version of this branch's first commit, on the store), the H1 and H3 runs
at 64 threads created no `pybdt_*` files.

Against `pr/05-sparse-coverage` (`769beba`), on cbsubscb19 (AMD EPYC 9554), three alternating runs each, median
wall time:

| Dataset | 16 threads | 64 threads |
|---|---|---|
| `H1_bam_reps_karpas422` | 43.7 → 42.9 s (0.981) | 28.0 → 27.4 s (0.976) |
| `H3_bam_ct_gm12878_grocap` | 56.7 → 55.4 s (0.977) | 37.8 → 36.6 s (0.968) |
| `T4_bam_reps_hct116` | 23.6 → 23.5 s (0.994) | 15.4 → 15.0 s (0.975) |
| `T5_bam_ct_k562_grocap` | 10.2 → 9.7 s (0.951) | 7.8 → 7.3 s (0.941) |

Outputs are identical. Most of the gain is in input parsing and in filtering, pairing and QC (H3 at 16 threads:
13.2 → 12.6 s), where the BED operations no longer go through temporary files and a bedtools process.

### Testing

- `tests/unit/test_bed_ops.py`: sort, merge (with `collapse` and `distinct`) and `merge -d` against bedtools on
  60 random sets of records with many tied starts; strand-aware merge and `intersect -v/-u/-c` (with and without
  `-s`) against bedtools on 60 more, including zero-length records, `.` strands and an unsorted `-a`; `cat`
  against pybedtools on 20; the heapsort fallback and insertion-sort thresholds; GTF sorting; reading of headers
  and gzipped files.
- `tests/diff/test_differential_bedtools_free.py`: borrowing elements from other replicates, boundary extension
  with a promoter BED (plain and gzipped) and `pints_counter`'s regions (`--max-dist` 0, 1 and 250) give the same
  files as the reference implementation, which uses pybedtools.
- Full suite at this branch's head (`efe3a3f`): 423 passed.
- chr22 test data, with bigWig, BAM, and `--annotation-gtf --focused-chrom chr22` input, at 1 and 4 threads: the
  bidirectional, divergent and unidirectional peaks are byte-identical to `main`.
- Benchmark equivalence of `ca581dd`: all 78 runs (19 datasets; default, `--disable-small` and sensitive
  configurations; 1, 32 and 64 threads) are byte-identical to `main`. Its `pints/` code differs from this branch
  only by the pandas line above, the same line in stratified filtering, and the sparse coverage PR's last commit
  (old `.npy` dumps are not read), which does nothing in a fresh output directory.
- Clean install: `pip install` of this branch into a fresh virtual environment (Python 3.12), then `pints_caller` on
  the chr22 bigWigs at 4 threads with a `PATH` holding only the environment and three coreutils, so no bedtools,
  bgzip or tabix can be found. The run completes and writes the three peak files.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
