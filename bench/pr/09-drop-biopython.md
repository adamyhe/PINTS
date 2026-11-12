## Drop Biopython; list `setup.py` dependencies one per line

**The output is unchanged.** `pints_normalizer` finds the same spike-in records, counts and scales. Two commits;
please review them one by one.

| Commit | Problem | Fix |
|---|---|---|
| setup.py: one dependency per line | `install_requires` was two long lines, so any two changes that add or remove different dependencies conflict | One dependency per line. The list is otherwise unchanged |
| Drop Biopython: read the spike-in FASTA with pysam | `pints_normalizer` used Biopython only to collect the record IDs of the spike-in FASTA, and every install had to pull in Biopython for it | `pysam.FastxFile` gives the same names (the header up to the first whitespace), and pysam is already a dependency. As a side effect, the spike-in file may now also be FASTQ or gzipped: nothing restricts the input to plain FASTA. `biopython` is removed from `install_requires` |

Please merge without squashing: later PRs contain these commits (B1 in the overview contains the `setup.py` commit, B8
both), and a squashed copy would conflict with them in `setup.py`. Merged as they are, both commits merge cleanly with
every other open branch.

### Testing

- `pints_normalizer`'s command line does not run on `main` (fixed in A3), so the check calls
  `spikein_normalization` directly, with bigWig writing stubbed out.
- Spike-in FASTA with two records whose headers carry a description after a tab and after a space; BAMs
  `tests/datasets/chr22.bam` and `input.chr22.bam`: total reads, spike-in reads, spike-in percentages and spike-in
  scales logged by this branch equal those logged by `main` with Biopython 1.88.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
