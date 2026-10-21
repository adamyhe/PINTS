"""Run one implementation's unmodified scripts/pints_caller with profiling hooks (used by bench/profiling.py).

    python _pints_profile_shim.py <impl_checkout> <spans 0|1> <cprofile_out|-> <pints_caller args...>

spans=1 wraps the functions in SPANS so that every call, in the parent and in forked pool workers, writes
"[span] <pid> B|E <perf_counter> <name> <context>" to stderr. perf_counter is CLOCK_MONOTONIC on Linux, so the
runner can line the spans up with its memory samples. The wrappers only add timing; arguments and return values
pass through unchanged. Names missing from a checkout are skipped.

cprofile_out profiles the parent process with cProfile (meant for 1-thread runs, where all work is in the parent).
"""
import functools
import multiprocessing
import os
import runpy
import sys
import time

multiprocessing.set_start_method("fork", force=True)

impl, spans, cprof = sys.argv[1], sys.argv[2] == "1", sys.argv[3]
sys.path[0] = impl

# (module, attribute, label). Pool workers resolve functions by module attribute, so patching the module is enough.
SPANS = [
    ("pints.io_engine", "get_read_signal", "bam: get_read_signal"),
    ("pints.io_engine", "parse_bam_by_contig", "bam: parse_bam_by_contig"),
    ("pints.io_engine", "_bam_contig_worker", "bam: contig worker"),
    ("pints.io_engine", "get_coverage", "bam: get_coverage"),
    ("pints.io_engine", "merge_replicates_bw", "reps: merge_replicates_bw"),
    ("pints.io_engine", "_merge_rep_chrom_worker", "reps: chrom worker"),
    ("pints.io_engine", "normalize_using_input", "ctrl: normalize_using_input"),
    ("pints.io_engine", "_normalize_chrom_worker", "ctrl: chrom worker"),
    ("pints.io_engine", "index_bed_file", "io: index_bed_file"),
    ("pints.calling_engine", "parse_input_files", "input: parse_input_files"),
    ("pints.calling_engine", "unified_element_definition", "elements: unified_element_definition"),
    ("pints.calling_engine", "_element_definition_task", "elements: task"),
    ("pints.calling_engine", "evaluate_strand_jobs", "eval: evaluate_strand_jobs"),
    ("pints.calling_engine", "_subpeak_phase1", "eval: phase1 chrom"),
    ("pints.calling_engine", "_segment_chunks", "eval: segment_chunks"),
    ("pints.calling_engine", "_index_bed_task", "io: index_bed_task"),
    ("pints.calling_engine", "stratified_filtering", "filter: stratified_filtering"),
    ("pints.stats_engine", "independent_filtering", "filter: independent_filtering"),
    ("pints.calling_engine", "pair", "pair: pair"),
    ("pints.calling_engine", "merge_opposite_peaks", "pair: merge_opposite_peaks"),
    ("pints.calling_engine", "open_peak_index", "pair: open_peak_index"),
    ("pints.calling_engine", "inferring_elements_from_other_reps", "pair: other reps"),
    ("pints.calling_engine", "get_epig_annotation", "annot: get_epig_annotation"),
    ("pints.calling_engine", "annotate_tre_with_epig_info", "annot: annotate_tre_with_epig_info"),
    ("pints.calling_engine", "on_the_fly_qc", "qc: on_the_fly_qc"),
    ("pints.calling_engine", "housekeeping", "cleanup: housekeeping"),
    ("pints.calling_engine", "_release_free_heap", "cleanup: release_free_heap"),
    ("pybedtools.helpers", "call_bedtools", "ext: bedtools"),
    ("pysam", "tabix_index", "ext: tabix_index"),
    ("pysam", "tabix_compress", "ext: tabix_compress"),
]


def _emit(kind, name, ctx):
    os.write(2, ("[span] %d %s %.4f %s\t%s\n" % (os.getpid(), kind, time.perf_counter(), name, ctx)).encode())


def _ctx(args):
    vals = []
    for a in args:
        items = a if isinstance(a, (tuple, list)) and len(a) <= 12 else (a,)
        vals += [x for x in items if isinstance(x, (str, int)) and not isinstance(x, bool)]
    out = []
    for x in vals[:4]:
        if isinstance(x, str):
            x = os.path.basename(x) if os.sep in x else x
            x = x.split()[0] if x.strip() else x
        out.append(str(x)[:60])
    return ",".join(out) or "-"


def _wrap(name, f):
    @functools.wraps(f)
    def w(*a, **k):
        ctx = _ctx((list(a[0])[:4],) if name == "ext: bedtools" and a else a)   # e.g. bedtools,sort,-i,<file>
        _emit("B", name, ctx)
        try:
            return f(*a, **k)
        finally:
            _emit("E", name, ctx)
    return w


def _install():
    import importlib
    for modname in sorted({m for m, _, _ in SPANS}):
        try:
            importlib.import_module(modname)
        except ImportError:
            pass
    import pybedtools.bedtool  # noqa: F401  (imports call_bedtools by name)
    done = []
    for modname, attr, label in SPANS:
        mod = sys.modules.get(modname)
        orig = getattr(mod, attr, None) if mod else None
        if orig is None or not callable(orig) or hasattr(orig, "py_func"):   # missing, or a Numba dispatcher
            continue
        w = _wrap(label, orig)
        for m in list(sys.modules.values()):
            n = getattr(m, "__name__", "") or ""
            if n.split(".")[0] in ("pints", "pybedtools", "pysam"):
                for k, v in list(vars(m).items()):
                    if v is orig:
                        setattr(m, k, w)
        done.append(label)
    sys.stderr.write("[profile] spans installed: %s\n" % "; ".join(done))


if spans:
    _install()
script = os.path.join(impl, "scripts", "pints_caller")
sys.argv = [script] + sys.argv[4:]
if cprof != "-":
    import cProfile
    prof = cProfile.Profile()
    os.register_at_fork(after_in_child=prof.disable)
    prof.enable()
    try:
        runpy.run_path(script, run_name="__main__")
    finally:
        prof.disable()
        prof.dump_stats(cprof)
else:
    runpy.run_path(script, run_name="__main__")
