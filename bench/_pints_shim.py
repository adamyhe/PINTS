"""Run one implementation's unmodified scripts/pints_caller.

    python _pints_shim.py <impl_checkout> <keep_intermediate 0|1> <pints_caller args...>

keep_intermediate=1 disables PINTS's final housekeeping so the per-strand candidate tables
(bench_*_pl.bed.gz / *_mn.bed.gz, with p-values and padj) survive for the Gate 2 drift comparison.

The multiprocessing start method is pinned to fork, the default on Linux up to Python 3.13: the reference
passes its model objects to Pool workers as fork-inherited module globals and fails under spawn/forkserver
(probe B10; macOS, Linux from Python 3.14). The optimised branches create fork contexts explicitly anyway.
"""
import multiprocessing
import os
import runpy
import sys

multiprocessing.set_start_method("fork", force=True)

impl, keep = sys.argv[1], sys.argv[2] == "1"
sys.path[0] = impl          # replace this script's directory with the checkout
if keep:
    import pints.calling_engine as _ce
    _ce.housekeeping = lambda *a, **k: None
script = os.path.join(impl, "scripts", "pints_caller")
sys.argv = [script] + sys.argv[3:]
runpy.run_path(script, run_name="__main__")
