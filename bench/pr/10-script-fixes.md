## Fix `pints_normalizer`, `pints_visualizer`, `pints_counter`, `pints_boundary_extender` and `parse_gtf`

Against `main`. Two commits: the first makes `pints_normalizer` run, the second fixes the other scripts. Review
them commit by commit.

**`pints_caller`'s output is byte-identical to the current release.** The other scripts write the same files as
before, except where a fix makes an option take effect or moves an output (listed below).

### `pints_normalizer` (first commit)

The script could not run at all: it read `args.output_prefix` (the option is `--output-dir`), and it called
`bam_to_bigwig`, which it never imported. Behind that:

| Problem | Fix |
|---|---|
| RPM passed `output_dir=` instead of `output_pref=`, and converted each BAM twice | Each BAM is converted once |
| Spike-in joined the output directory with the full BAM path, so an absolute path wrote the bigWigs next to the BAM | Both methods now write `<output-dir>/<bam name>_{pl,mn}.bw` |
| Spike-in dropped `--reverse-complement` and `--mapq-threshold` | Both methods call `pints_visualizer`'s `bam_to_bigwig` with every option forwarded |
| No default for `--output-dir` | It defaults to the current directory. A missing directory stops the script with a message |
| The log file name used minutes (`%M`) where the month (`%m`) was meant | `%m` |

### The other scripts (second commit)

| Script | Problem | Fix |
|---|---|---|
| `pints_visualizer` | Without `-o` it crashed with `TypeError` (`None + "_pl.bw"`) | The output defaults to the BAM file name without `.bam`, in the current directory |
| `pints_visualizer` | The log file name used minutes where the month was meant | `%m` |
| `pints_counter` | `--max-dist` was parsed but never passed on, so it had no effect | It is passed to `main` |
| `pints_boundary_extender` | `--save-to` was ignored: the extended elements and the temporary coverage files were written next to the peak calls | They go to `--save-to` when it is given. Without it, they stay next to the peak calls as before |
| `pints_boundary_extender` | `--bw-pl` and `--bw-mn` had to be paths of the same length (the check compared the lengths of two strings) | The check is removed. The numbers of files are already checked |
| `pints_boundary_extender` | The file name tag of an asymmetric extension was built from a set, so `--div-ext-left 60 --div-ext-right 200` gave `60_200` or `200_60` depending on the hash seed | The tag is always `left_right` (`60_200`). A symmetric extension still gives one value (`60`) |
| `parse_gtf` (`pints_sample_qc`, `pints_caller --annotation-gtf`) | It writes a cache, `<gtf>_pd_cache.csv`, next to the GTF. In a read-only folder this failed with `PermissionError`. A cache older than the GTF was still used | The cache is skipped when it cannot be written, and rebuilt when it is older than the GTF |

### Behaviour that changes

- **`pints_boundary_extender --save-to DIR`**: the extended elements (`*_element_*bp.bed`, `*_e.bed`) are now
  written to `DIR` instead of next to the peak calls. `--save-to` now defaults to unset; the log still goes to the
  current directory when it is unset. A `DIR` that does not exist stops the script with a message.
- **Asymmetric extensions**: output names are always `<left>_<right>`, for example
  `t_1_divergent_peaks_element_60_200bp.bed`. Before, some runs named them `200_60`, depending on the hash seed.
- **`pints_counter --max-dist N`** now merges features up to `N` bp apart, as its help text says. With the default
  (0), nothing changes.
- **`parse_gtf`**: a `_pd_cache.csv` older than its GTF is now ignored and overwritten.
- **`pints_visualizer` and `pints_normalizer` logs** are named with the month. A run at 00:02 on 3 November 2026
  used to write `visualizing_2026_02_03_00_02_06.log`; it now writes `visualizing_2026_11_03_00_02_06.log`.

### Testing

- `pints_normalizer`, chr22 BAMs (`tests/datasets/chr22.bam`, `input.chr22.bam`), `-e PROcap`:
  - `main` fails with `AttributeError: 'Namespace' object has no attribute 'output_prefix'`;
  - `-m RPM`: each bigWig has the same chromosomes and intervals as `pints_visualizer -c` on the same BAM;
  - `-m Spikein`, with and without `-r`: each bigWig has the same chromosomes and intervals as `pints_visualizer`
    with the logged spike-in scale as `-n` and the same `-r`;
  - nothing is written next to the BAMs.
- Other scripts on chr22, `main` against this branch, each run on its own copy of `main`'s peak calls:
  - `pints_visualizer` without `-o`: `main` crashes; this branch writes `chr22_{pl,mn}.bw`, byte-identical to
    `main`'s output with `-o`. With `-o`, the bigWigs are byte-identical to `main`'s;
  - `pints_counter`: byte-identical to `main` by default (1178 lines). `--max-dist 1000` gives 734 lines on this
    branch and 1178 on `main`;
  - `pints_boundary_extender` with bigWig input, with and without `--promoter-bed`: elements byte-identical to
    `main`'s. With `--save-to` (bigWig with `--promoter-bed`, and BAM), they are written to that folder only; the
    bigWig ones are byte-identical to `main`'s. No `.npy` files are left behind;
  - a `--bw-mn` path of a different length: `main` stops with "Must provide the same amount of bigwig files for
    both strands"; this branch gives the same elements as with the other path;
  - an asymmetric extension under `PYTHONHASHSEED` 1 to 4: `main` names the files `60_200` once and `200_60` three
    times; this branch always gives `60_200`, byte-identical to `main`'s files;
  - `pints_sample_qc` with the GTF in a read-only folder: `main` fails with `PermissionError`; this branch
    completes;
  - `parse_gtf` ignores a cache older than the GTF and uses a fresh one.
- `main`'s test suite passes on this branch (4 passed).
- `pints_caller` on chr22, bigWig and BAM input at 1 and 4 threads: peaks byte-identical to `main`. With
  `--annotation-gtf --focused-chrom chr22` (both sides given the chromosome under the name `main`'s engine reads),
  also byte-identical.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
