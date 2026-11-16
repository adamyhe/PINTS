## Update: the performance work is now one PR, plus independent ones

Thanks for merging #30–#33. We have changed how the rest is split, so that every open PR applies to `main` on its
own and none of them conflicts with another. Before, the ten performance PRs (B1–B10 above) formed a stack that had
to be merged in order, one at a time.

Most of that stack rewrites the same functions in turn (`calling_engine.py`, `stats_engine.py`, `io_engine.py`), so
those parts cannot apply to `main` separately. They are now one PR. The two optional changes become follow-ups,
independent of each other:

| PR | What | Applies to | Output with default options |
|---|---|---|---|
| #36 | **Engine**: test harness (was #34), exact kernels, parallel engine, I/O and memory, sparse coverage store, no bedtools, fewer dependencies, steps after evaluation (B1–B6, B8, B9) | `main` | byte-identical to `main` |
| #35 | **`--epig-annotation`**: a local bigBed path is used without contacting the server; a failed server lookup logs a warning and skips the annotation instead of crashing (`pints.yulab.org` has returned 502 since at least 2026-11-15) | `main` | unchanged |
| later | **Opt-in fixes** (B7): `--iqr-scale`, `--fix-ler-flanks`, `--strict-independent-filtering`, `--stringent-pairs-padj` | the engine PR | unchanged (all off by default) |
| later | **`--fast-math`** (B10) | the engine PR | unchanged (off by default) |

- #36 and #35 can merge in either order.
- The two follow-ups will be opened once the engine PR has merged. They change code that the engine PR rewrites, but
  not each other's: each merges on its own, in either order. Applying both gives the same code as the top of the
  old stack, apart from how the options reach the worker processes.
- #34 is closed, since the test harness is the engine PR's first commit.

The engine PR's description has a review order by part. Each part links to its detailed description (timings,
risks, tests), which still describes the old layout. The questions in the overview above still apply: the frozen
reference copy in `tests/`, whether any opt-in fix should become the default, and whether you want `--fast-math`.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
