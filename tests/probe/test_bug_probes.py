"""Bug / deviation probes against the CURRENT code (pints == pints_ref at the start of the project).

Convention
----------
* Each probe asserts the *correct / documented* behaviour.
* Confirmed bugs are marked ``xfail(strict=True)``: the suite stays green while the bug exists and turns
  red (XPASS) the moment behaviour changes -- so Track A (bit-identical) can never fix a bug by accident,
  and Track B fixes must flip exactly the intended probe(s) (and only behind their flag).
* Probes whose "correct" behaviour is a spec question (paper is ambiguous) are written as
  *characterisation* tests (plain asserts of current behaviour) and tagged ``spec``.
"""
import os
import runpy
import sys
import tempfile
import multiprocessing as mp

import numpy as np
import pandas as pd
import pysam
import pytest

import pints.calling_engine as ce
import pints.stats_engine as se
import pints.io_engine as io
from helpers import write_subpeak_tabix, make_bigwig

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CALLER = os.path.join(ROOT, "scripts", "pints_caller")
CHR = "chrP"


def strict(bug_id, why, raises=AssertionError):
    # `raises` guards against a probe "xfailing" for an unrelated reason (e.g. a harness error)
    return pytest.mark.xfail(strict=True, raises=raises, reason="%s: %s" % (bug_id, why))


# ------------------------------------------------------------------ CLI wiring (B1, B2, B12)
def _cli_kwargs(monkeypatch, tmp_path, extra):
    """Run scripts/pints_caller's argparse + dispatch with peak_calling stubbed out; return kwargs."""
    captured = {}

    def fake_peak_calling(*a, **kw):
        captured.update(kw)

    monkeypatch.setattr(ce, "peak_calling", fake_peak_calling)
    bw = tmp_path / "x.bw"
    bw.write_bytes(b"")
    argv = ["pints_caller", "--bw-pl", str(bw), "--bw-mn", str(bw), "--save-to", str(tmp_path),
            "--file-prefix", "t", "--dont-check-updates"] + extra
    monkeypatch.setattr(sys, "argv", argv)
    old_hook = sys.excepthook
    try:
        runpy.run_path(CALLER, run_name="__main__")
    finally:
        sys.excepthook = old_hook
    return captured


# B1 fixed (Track B, plain CLI wiring): --focused-chrom was stored as `focused_chrom`; engine reads
# `highlight_chromosome`
def test_B1_focused_chrom_reaches_engine(monkeypatch, tmp_path):
    kw = _cli_kwargs(monkeypatch, tmp_path, ["--focused-chrom", "chr7"])
    assert ce._highlight_chromosome(kw) == "chr7"
    assert ce._highlight_chromosome(_cli_kwargs(monkeypatch, tmp_path, [])) == "chr1"


# B2 fixed (Track B, plain CLI wiring): --stringent-pairs-only was stored as `stringent_pairs_only`;
# pair() reads `stringent_only`
def test_B2_stringent_pairs_only_reaches_engine(monkeypatch, tmp_path):
    kw = _cli_kwargs(monkeypatch, tmp_path, ["--stringent-pairs-only"])
    assert ce._stringent_only(kw) is True
    assert ce._stringent_only(_cli_kwargs(monkeypatch, tmp_path, ["--stringent-pairs-padj"])) is True
    assert ce._stringent_only(_cli_kwargs(monkeypatch, tmp_path, [])) is False


_ALL_DESTS = ("bw_pl bw_mn ct_bw_pl ct_bw_mn ct_bam bam_parser seq_rc merge_replicates adjust_method fdr_target "
              "close_threshold stringent_pairs_only min_len_opposite_peaks mapq_threshold small_peak_threshold "
              "window_size_threshold keep_sticks annotation_gtf tss_extension focused_chrom donor_tolerance "
              "ce_trigger top_peak_threshold min_mu_percent peak_distance peak_width div_size_min summit_dist_min "
              "iqr_strategy enable_eler eler_lower_bound sensitive fc_cutoff init_dens_cutoff init_height_cutoff "
              "epig_annotation relaxed_fdr_target chromosome_startswith output_chrom_size disable_qc strict_qc "
              "output_diagnostics borrow_info_reps thread_n").split()
_UNREAD = {}  # B1, B2, B12 fixed in Track B (all CLI dests now reach the engine)


@pytest.mark.parametrize("dest", [pytest.param(d, marks=strict(_UNREAD[d], "CLI dest never read by the engine"))
                                  if d in _UNREAD else d for d in _ALL_DESTS])
def test_cli_dest_is_consumed(dest):
    src = ""
    for f in ("calling_engine.py", "io_engine.py", "stats_engine.py", "qc_engine.py", "extension_engine.py"):
        with open(os.path.join(os.path.dirname(ce.__file__), f)) as fh:
            src += fh.read()
    import re
    pat = r"""((get|pop)\(\s*["']%s["'])|(\[\s*["']%s["']\s*\])|(\b%s\s*=)""" % (dest, dest, dest)
    assert re.search(pat, src), dest


# ------------------------------------------------------------------ models (B13, B5)
# B13 fixed (Track B): Poisson did not implement the abstract lrt(); `--model Poisson` raised TypeError
def test_B13_poisson_model_instantiable():
    se.Poisson()


def test_B13_poisson_lrt():
    from scipy.stats import poisson, chi2
    rng = np.random.default_rng(3)
    a, b = rng.poisson(2.0, 200), rng.poisson(3.0, 300)
    ll = lambda x: poisson.logpmf(x, x.mean()).sum()
    w, p = se.Poisson().lrt(a, b)
    w_ref = 2 * (ll(a) + ll(b) - ll(np.concatenate([a, b])))
    assert np.isclose(w, w_ref) and np.isclose(p, chi2.sf(w_ref, 1))
    assert se.Poisson().lrt(np.zeros(5), np.zeros(5)) == (0., 1.)
    ps = [se.Poisson().lrt(rng.poisson(2, 100), rng.poisson(2, 300))[1] for _ in range(400)]
    assert 0.02 < np.mean(np.array(ps) < 0.05) < 0.09


@pytest.mark.spec
def test_B5_sf_characterisation():
    """ZIP.sf == poisson.sf(floor(mu_peak); mu_bg): ignores pi_peak/pi_bg and var; constant on [k, k+1).

    The paper only says 'the estimated densities are then used to determine the statistically significant
    peaks', so this is recorded as current behaviour (spec ambiguity), not a bug.
    """
    z = se.ZIP()
    a = z.sf(3.01, 1.0, 0.1, 1.0, 1.0, 0.0)
    b = z.sf(3.99, 9.0, 0.9, 1.0, 5.0, 0.95)
    assert a == b
    assert z.sf(0.2, 1, 0, 0.01, 1, 0) == z.sf(0.99, 1, 0, 0.01, 1, 0)


@pytest.mark.spec
def test_D1_iqr_fence_scale():
    """Paper (Methods): outliers are > Q3 + 1.5 x IQR. Code: scale = 1.0 on both fences."""
    data = np.array([1., 2., 3., 4., 5.])
    q1, q3 = np.quantile(data, [.25, .75])
    lo, hi = se.IQR.get_outlier_threshold(data, direction=0)
    assert hi == q3 + 1.0 * (q3 - q1)      # current behaviour
    assert hi != q3 + 1.5 * (q3 - q1)      # deviation from described algorithm


@strict("D1", "IQR fence uses 1.0 x IQR, paper states 1.5 x IQR (fixed only with --iqr-scale 1.5)")
def test_D1_iqr_fence_matches_paper():
    data = np.array([1., 2., 3., 4., 5.])
    q1, q3 = np.quantile(data, [.25, .75])
    lo, hi = se.IQR.get_outlier_threshold(data, direction=0)
    assert hi == q3 + 1.5 * (q3 - q1)


@pytest.fixture
def fixes():
    saved = dict(se.FIXES)
    yield se.configure_fixes
    se.FIXES.clear()
    se.FIXES.update(saved)


def test_D1_iqr_scale_flag(fixes):
    data = np.array([1., 2., 3., 4., 5.])
    q1, q3 = np.quantile(data, [.25, .75])
    fixes(iqr_scale=1.5)
    lo, hi = se.IQR.get_outlier_threshold(data, direction=0)
    assert lo == q1 - 1.5 * (q3 - q1) and hi == q3 + 1.5 * (q3 - q1)
    assert se.IQR.get_outlier_threshold(data, direction=0, scale=1.0)[1] == q3 + (q3 - q1)


# ------------------------------------------------------------------ LER background window (B4)
class RecordingZIP(se.ZIP):
    def __init__(self):
        super().__init__()
        self.calls = []

    def fit(self, window):
        self.calls.append(np.asarray(window).copy())
        return super().fit(window)


def _ler_single_candidate(cov, peak, ler, mu=5.0):
    """Run bgIQR.remove_peaks_in_local_env for a tested peak with exactly one LER candidate.

    Returns the array handed to fit() as the LER candidate's background window.
    """
    L = len(cov)
    cs, cend = peak
    f = 5000
    q = (max(cs - f, 0), max(cs - 1, 0), min(cend + 1, L), min(cend + f, L))
    with tempfile.TemporaryDirectory() as d:
        gz = write_subpeak_tabix(os.path.join(d, "sp.bed"), CHR,
                                 [(ler[0], ler[1], "ler", mu, 0.1, "+", 0.2, 1, ler[0], 1)])
        tbx = pysam.TabixFile(gz, index=gz + ".csi")
        z = RecordingZIP()
        se.bgIQR.remove_peaks_in_local_env(z, tbx, CHR, *q, 5, 1.0, cov, 0.1, {},
                                           disable_ler=False, enable_eler=False, peak_threshold=1e9)
    # calls: [fit(peak region of LER candidate), fit(background of LER candidate)]
    assert len(z.calls) == 2
    return z.calls[1], q


def _expected_ler_bg(cov, ler, q):
    """Documented intent: 2.5 kb flanks of the LER candidate (clipped to the chromosome), with the
    region of the peak under test masked."""
    L = len(cov)
    b, s = ler
    idx = np.r_[max(b - 2500, 0):b, s:min(s + 2500, L)]
    keep = ~((idx >= q[1]) & (idx < q[2]))
    return cov[idx[keep]]


@strict("B4a", "(fixed only with --fix-ler-flanks) mask offsets assume the tested peak lies inside the LER's 2.5 kb flank; when it does not, "
               "negative offsets wrap and mask an unrelated part of the right flank")
def test_B4a_ler_mask_far_candidate():
    cov = np.ones(40_000, dtype=np.int32)
    peak, ler = (15_000, 15_050), (19_000, 19_050)
    got, q = _ler_single_candidate(cov, peak, ler)
    assert len(got) == len(_expected_ler_bg(cov, ler, q))


@strict("B4b", "(fixed only with --fix-ler-flanks) coverage[b-2500:b] with b < 2500 is a negative-start slice -> empty left flank")
def test_B4b_ler_flank_near_chromosome_start():
    cov = np.ones(40_000, dtype=np.int32)
    peak, ler = (3_000, 3_050), (1_000, 1_050)
    got, q = _ler_single_candidate(cov, peak, ler)
    assert len(got) == len(_expected_ler_bg(cov, ler, q))


@pytest.mark.parametrize("peak,ler", [((15_000, 15_050), (19_000, 19_050)),   # B4a
                                      ((3_000, 3_050), (1_000, 1_050)),       # B4b
                                      ((15_000, 15_050), (16_000, 16_050)),   # control
                                      ((39_000, 39_050), (38_000, 38_050)),   # near chromosome end
                                      ((100, 150), (1_000, 1_050))])          # tested peak at chromosome start
def test_B4_fixed_with_flag(fixes, peak, ler):
    fixes(fix_ler_flanks=True)
    rng = np.random.default_rng(0)
    cov = rng.integers(0, 5, 40_000).astype(np.int32)
    got, q = _ler_single_candidate(cov, peak, ler)
    assert np.array_equal(got, _expected_ler_bg(cov, ler, q))


def test_B4_control_case_inside_flank():
    """Control: LER candidate within 2.5 kb of the tested peak, far from the chromosome start -> correct."""
    cov = np.ones(40_000, dtype=np.int32)
    peak, ler = (15_000, 15_050), (16_000, 16_050)
    got, q = _ler_single_candidate(cov, peak, ler)
    assert len(got) == len(_expected_ler_bg(cov, ler, q))


# ------------------------------------------------------------------ LER cache order dependence (B7)
def _find_cache_dependence(max_seed=150):
    """Search for a query whose result changes when the shared cache was populated by an earlier query."""
    from helpers import synth_coverage
    for seed in range(max_seed):
        rng = np.random.default_rng(seed)
        cov = synth_coverage(rng, length=30_000, n_peaks=60)
        L = len(cov)
        starts = np.sort(rng.choice(L - 100, 40, replace=False))
        rows = [(int(s), int(s + rng.integers(5, 60)), "p%d" % i, float(np.round(rng.gamma(1.5, 2.0), 6)),
                 0.1, "+", 0.2, 1, int(s), 1) for i, s in enumerate(starts)]
        with tempfile.TemporaryDirectory() as d:
            gz = write_subpeak_tabix(os.path.join(d, "sp.bed"), CHR, rows)
            tbx = pysam.TabixFile(gz, index=gz + ".csi")
            qs = []
            for a in rng.integers(2_000, L - 2_000, 6):
                a = int(a)
                qs.append((max(a - 5000, 0), a - 1, min(a + 51, L), min(a + 5050, L)))
            for i in range(1, len(qs)):
                cache = {}
                for j in range(i):
                    se.bgIQR.remove_peaks_in_local_env(se.ZIP(), tbx, CHR, *qs[j], 5, 1.0, cov, 0.1, cache,
                                                       enable_eler=True, peak_threshold=5.0)
                warm, _ = se.bgIQR.remove_peaks_in_local_env(se.ZIP(), tbx, CHR, *qs[i], 5, 1.0, cov, 0.1,
                                                             cache, enable_eler=True, peak_threshold=5.0)
                cold, _ = se.bgIQR.remove_peaks_in_local_env(se.ZIP(), tbx, CHR, *qs[i], 5, 1.0, cov, 0.1,
                                                             {}, enable_eler=True, peak_threshold=5.0)
                if warm.shape != cold.shape or not np.array_equal(warm, cold):
                    return seed, i
    return None


@pytest.mark.spec
def test_B7_ler_cache_makes_results_order_dependent():
    """The LER decision for a subpeak is cached after the first evaluation (made relative to the FIRST tested
    peak's neighbourhood) and reused for later peaks.  Results therefore depend on evaluation order;
    deterministic serially, but any parallel split must preserve per-chromosome order and cache sharing."""
    assert _find_cache_dependence() is not None


# ------------------------------------------------------------------ input parsing (B3, B14)
# B3 fixed (Track B): non-integer bigWig values were silently truncated; now a warning is logged
def test_B3_bigwig_float_values_not_silently_truncated(tmp_path, caplog):
    v = np.zeros(1000)
    v[100:110] = 0.5
    bw = make_bigwig(tmp_path / "f.bw", [("chr1", 1000)], {"chr1": v})
    cov, rc = io._get_coverage_bw(bw, "", str(tmp_path), "f")
    arr = np.load(cov["chr1"])
    warned = any("integer" in r.getMessage().lower() or "float" in r.getMessage().lower()
                 for r in caplog.records)
    assert warned or arr.sum() == 5


def test_B14_normalize_drops_chromosomes_without_control(tmp_path, logger):
    """Characterisation: treatment chromosomes absent from the control are dropped from the analysis
    (only a warning).  Reported in the Track B impact report; may be intended."""
    a = tmp_path / "a.npy"; np.save(a, np.ones(10, dtype=np.int8))
    b = tmp_path / "b.npy"; np.save(b, np.ones(10, dtype=np.int8))
    pl, mn = io.normalize_using_input({"chr1": str(a), "chrX": str(a)}, {"chr1": str(b), "chrX": str(b)},
                                      {"chr1": str(b)}, {"chr1": str(b)}, 0.5, str(tmp_path), "n", logger)
    assert set(pl) == {"chr1"} and set(mn) == {"chr1"}


# ------------------------------------------------------------------ element definition edge cases
@strict("O1", "sliding_window drops the final partial (<100 bp) window, so signal in the last "
              "len % 100 bp of a chromosome can never become a peak candidate")
def test_O1_signal_at_chromosome_end_is_scanned():
    cov = np.zeros(1_050, dtype=np.int32)
    cov[1_020:1_030] = 10
    total = sum(w for w, _ in ce.sliding_window(cov))
    assert total == cov.sum()


@strict("O2", "cut_peaks sums events over [m0, m1] inclusive but divides by (m1 - m0); with rel_height=1 both "
              "ends are usually zero-valued bases, so this only bites when signal touches the window end")
def test_O2_cut_peaks_density_uses_consistent_interval(monkeypatch):
    captured = {}
    real = ce.quasi_max_score_segment

    def spy(candidates, **kw):
        captured["c"] = list(candidates)
        return real(candidates=candidates, **kw)

    monkeypatch.setattr(ce, "quasi_max_score_segment", spy)
    w = np.zeros(60, dtype=np.int32)
    w[10:14] = [2, 8, 8, 2]
    w[57:60] = [3, 9, 3]   # signal touching the end of the merged window
    ce.cut_peaks(w, donor_tolerance=0.3, ce_trigger=3, join_distance=1)
    for s, e, ev, dens in captured["c"]:
        assert ev == w[s:e].sum() and dens == ev / (e - s)


# ------------------------------------------------------------------ multiprocessing (B10)
def _worker_has_globals(_):
    import pints.calling_engine as m
    return m.stat_tester is not None and m.iqr_obj is not None


@strict("B10", "stat_tester / iqr_obj are module globals set in the parent; with the spawn/forkserver start "
               "method (forkserver is the Linux default from Python 3.14) workers see None")
def test_B10_workers_see_model_objects_under_spawn():
    ce.stat_tester, ce.iqr_obj = se.ZIP(), se.bgIQR()
    try:
        with mp.get_context("spawn").Pool(1) as p:
            assert all(p.map(_worker_has_globals, [0]))
    finally:
        ce.stat_tester = ce.iqr_obj = None


# ------------------------------------------------------------------ replicate borrowing (B9)
@strict("B9", "inferring_elements_from_other_reps indexes sigs[i] / output names by position in the list of "
              "EXISTING files, so a missing file in one replicate misaligns (or crashes with IndexError)",
        raises=(AssertionError, IndexError))
def test_B9_replicate_alignment_with_missing_file(tmp_path):
    prefix = str(tmp_path / "s")
    line = "chr1\t100\t400\t.\t.\t.\n"
    # rep1: bidirectional only (no unidirectional file); rep2: everything
    open(prefix + "_1_bidirectional_peaks.bed", "w").write(line)
    open(prefix + "_1_divergent_peaks.bed", "w").write(line)
    open(prefix + "_2_bidirectional_peaks.bed", "w").write("chr1\t5000\t5300\t.\t.\t.\n")
    open(prefix + "_2_divergent_peaks.bed", "w").write("chr1\t5000\t5300\t.\t.\t.\n")
    open(prefix + "_2_unidirectional_peaks.bed", "w").write("chr1\t150\t200\t.\t.\t+\n")
    ce.inferring_elements_from_other_reps(prefix, 2)
    # rep1 has no singletons of its own, so nothing should be borrowed into rep1
    rep1 = open(prefix + "_1_bidirectional_peaks.bed").read().strip().splitlines()
    assert len(rep1) == 1


# ------------------------------------------------------------------ numpy forward-compat
def test_F1_get_elbow_uses_deprecated_2d_cross():
    """np.cross on 2-D vectors is deprecated (NumPy 2.0) and will be removed -> future breakage."""
    with pytest.warns(DeprecationWarning):
        se.get_elbow(np.arange(10.), np.arange(10.) ** 2)


def test_B12_peak_width_reaches_find_peaks():
    """B12 fixed (Track B): --peak-width is passed to scipy.signal.find_peaks(width=); None keeps legacy calls."""
    w = np.zeros(200, dtype=np.int32)
    w[20] = 9                                   # 1-bp spike
    w[100:130] = np.r_[np.arange(1, 16), np.arange(15, 0, -1)]   # broad peak
    assert ce.cut_peaks(w, 0.3, 3) == ce.cut_peaks(w, 0.3, 3, peak_width=None)
    legacy = ce.cut_peaks(w, 0.3, 3)
    wide = ce.cut_peaks(w, 0.3, 3, peak_width=5)
    assert any(s <= 20 < e for s, e in legacy) and not any(s <= 20 < e for s, e in wide)
    assert any(s <= 115 < e for s, e in wide)


def test_cli_kwargs_do_not_collide_with_explicit_engine_arguments(monkeypatch, tmp_path):
    """Regression: CLI-only keys must not duplicate keyword arguments that the engine passes explicitly
    together with **kwargs (e.g. pair(..., stringent_only, **kwargs))."""
    kw = _cli_kwargs(monkeypatch, tmp_path, ["--stringent-pairs-only", "--focused-chrom", "chr7"])
    assert "stringent_only" not in kw and "highlight_chromosome" not in kw
