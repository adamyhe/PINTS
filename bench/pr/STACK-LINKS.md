## Reading ahead: the B1–B10 stack, step by step

B1 is open now as draft #34. The other nine will be opened one at a time, because a PR from a fork can
only target `main`, and GitHub would show each one with everything below it. Until then, each link below shows
only that step's own diff against the one before it, on our fork. The descriptions are the ones they will be
opened with.

| | Step | Diff | Commits | Lines + / − | Description | Note |
|---|---|---|---|---|---|---|
| B1 | Test harness | [`02-trackb-fixes...03-test-harness`](https://github.com/adamyhe/PINTS/compare/pr/02-trackb-fixes...pr/03-test-harness) | 1 | 5,014 / 2 | [03-test-harness.md](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/03-test-harness.md) | Also shows #31's setup.py commit (2 lines). Open as draft #34 |
| B2 | Exact kernels | [`03-test-harness...03a-gate1-kernels`](https://github.com/adamyhe/PINTS/compare/pr/03-test-harness...pr/03a-gate1-kernels) | 5 | 688 / 27 | [03a-gate1-kernels.md](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/03a-gate1-kernels.md) |  |
| B3 | Parallel engine | [`03a-gate1-kernels...03b-gate1-engine`](https://github.com/adamyhe/PINTS/compare/pr/03a-gate1-kernels...pr/03b-gate1-engine) | 17 | 2,570 / 126 | [03b-gate1-engine.md](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/03b-gate1-engine.md) |  |
| B4 | I/O and memory | [`03b-gate1-engine...04-io`](https://github.com/adamyhe/PINTS/compare/pr/03b-gate1-engine...pr/04-io) | 13 | 1,088 / 192 | [04-io.md](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/04-io.md) |  |
| B5 | Sparse coverage | [`04-io...05-sparse-coverage`](https://github.com/adamyhe/PINTS/compare/pr/04-io...pr/05-sparse-coverage) | 10 | 1,230 / 127 | [05-sparse-coverage.md](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/05-sparse-coverage.md) |  |
| B6 | No bedtools | [`05-sparse-coverage...06-bedtools-free`](https://github.com/adamyhe/PINTS/compare/pr/05-sparse-coverage...pr/06-bedtools-free) | 2 | 725 / 111 | [06-bedtools-free.md](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/06-bedtools-free.md) |  |
| B7 | Opt-in fixes | [`06-bedtools-free...08-trackb-flags`](https://github.com/adamyhe/PINTS/compare/pr/06-bedtools-free...pr/08-trackb-flags) | 1 | 304 / 34 | [08-trackb-flags.md](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/08-trackb-flags.md) |  |
| B8 | Fewer dependencies | [`08-trackb-flags...12-fewer-deps`](https://github.com/adamyhe/PINTS/compare/pr/08-trackb-flags...pr/12-fewer-deps) | 5 | 233 / 39 | [12-fewer-deps.md](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/12-fewer-deps.md) | Includes the merge of #31, so #31's Biopython change shows too |
| B9 | Steps after evaluation | [`12-fewer-deps...15-serial-steps`](https://github.com/adamyhe/PINTS/compare/pr/12-fewer-deps...pr/15-serial-steps) | 3 | 324 / 19 | [15-serial-steps.md](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/15-serial-steps.md) |  |
| B10 | `--fast-math` | [`15-serial-steps...07-gate2`](https://github.com/adamyhe/PINTS/compare/pr/15-serial-steps...pr/07-gate2) | 8 | 626 / 32 | [07-gate2.md](https://github.com/adamyhe/PINTS/blob/bench/external/bench/pr/07-gate2.md) |  |

Commits are counted without merges. Line counts are `git diff --shortstat` between consecutive branch heads.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
