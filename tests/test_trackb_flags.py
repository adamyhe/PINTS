#  PINTS: Peak Identifier for Nascent Transcripts Sequencing
#  Copyright (c) 2019-2023 Yu Lab.
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
"""Opt-in behaviour fixes: --iqr-scale, --fix-ler-flanks, --strict-independent-filtering and
--stringent-pairs-padj. Each default reproduces the previous results; each flag changes only what it names."""
import os
import runpy
import sys
import tempfile

import numpy as np
import pandas as pd
import pysam
import pytest

import pints.calling_engine as ce
import pints.io_engine as io
import pints.stats_engine as se
from helpers import write_subpeak_tabix

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CALLER = os.path.join(ROOT, "scripts", "pints_caller")
CHR = "chrP"


@pytest.fixture
def fixes():
    saved = dict(se.FIXES)
    yield se.configure_fixes
    se.FIXES.clear()
    se.FIXES.update(saved)


def _cli_kwargs(monkeypatch, tmp_path, extra):
    """Run scripts/pints_caller's argument parsing with peak_calling stubbed out; return its keyword arguments."""
    captured = {}
    monkeypatch.setattr(ce, "peak_calling", lambda *a, **kw: captured.update(kw))
    bw = tmp_path / "x.bw"
    bw.write_bytes(b"")
    monkeypatch.setattr(sys, "argv", ["pints_caller", "--bw-pl", str(bw), "--bw-mn", str(bw), "--save-to",
                                      str(tmp_path), "--file-prefix", "t", "--dont-check-updates"] + extra)
    old_hook = sys.excepthook
    try:
        runpy.run_path(CALLER, run_name="__main__")
    finally:
        sys.excepthook = old_hook
    return captured


def test_cli_defaults_and_flags(monkeypatch, tmp_path):
    kw = _cli_kwargs(monkeypatch, tmp_path, [])
    assert (kw["iqr_scale"], kw["fix_ler_flanks"], kw["strict_independent_filtering"],
            kw["stringent_pairs_padj"]) == (1.0, False, False, False)
    kw = _cli_kwargs(monkeypatch, tmp_path, ["--iqr-scale", "1.5", "--fix-ler-flanks",
                                             "--strict-independent-filtering", "--stringent-pairs-padj"])
    assert (kw["iqr_scale"], kw["fix_ler_flanks"], kw["strict_independent_filtering"],
            kw["stringent_pairs_padj"]) == (1.5, True, True, True)
    assert ce._stringent_only(kw) is True  # --stringent-pairs-padj implies the stringent filter


def _worker_fixes(_):
    return dict(se.FIXES)


def test_workers_receive_fixes(fixes):
    """The settings reach pool workers through the initializer, not only by inheritance at the fork."""
    fixes(iqr_scale=1.5, fix_ler_flanks=True)
    pool = ce._make_pool(2)
    try:
        se.FIXES.update(iqr_scale=1.0, fix_ler_flanks=False)  # the parent changing later does not matter
        assert pool.map(_worker_fixes, range(2)) == [{"iqr_scale": 1.5, "fix_ler_flanks": True}] * 2
    finally:
        pool.terminate()


# ------------------------------------------------------------------ --iqr-scale
def test_iqr_scale(fixes):
    data = np.array([1., 2., 3., 4., 5.])
    q1, q3 = np.quantile(data, [.25, .75])
    assert se.IQR.get_outlier_threshold(data, direction=0) == (q1 - (q3 - q1), q3 + (q3 - q1))  # legacy 1.0
    fixes(iqr_scale=1.5)
    assert se.IQR.get_outlier_threshold(data, direction=0) == (q1 - 1.5 * (q3 - q1), q3 + 1.5 * (q3 - q1))
    assert se.IQR.get_outlier_threshold(data, direction=1, scale=1.0) == (None, q3 + (q3 - q1))


# ------------------------------------------------------------------ --fix-ler-flanks
class RecordingZIP(se.ZIP):
    def __init__(self):
        super().__init__()
        self.calls = []

    def fit(self, window):
        self.calls.append(np.asarray(window).copy())
        return super().fit(window)


def _ler_background(cov, peak, ler, mu=5.0):
    """The background window that bgIQR.remove_peaks_in_local_env fits for a tested peak's one LER candidate."""
    L = len(cov)
    cs, cend = peak
    q = (max(cs - 5000, 0), max(cs - 1, 0), min(cend + 1, L), min(cend + 5000, L))
    with tempfile.TemporaryDirectory() as d:
        gz = write_subpeak_tabix(os.path.join(d, "sp.bed"), CHR,
                                 [(ler[0], ler[1], "ler", mu, 0.1, "+", 0.2, 1, ler[0], 1)])
        tbx = pysam.TabixFile(gz, index=gz + ".csi")
        z = RecordingZIP()
        se.bgIQR.remove_peaks_in_local_env(z, tbx, CHR, *q, 5, 1.0, cov, 0.1, {},
                                           disable_ler=False, enable_eler=False, peak_threshold=1e9)
    assert len(z.calls) == 2  # fit(LER candidate), fit(its background)
    return z.calls[1], q


def _intended_background(cov, ler, q):
    """2.5 kb flanks of the LER candidate, clipped to the chromosome, without the region of the tested peak."""
    b, s = ler
    idx = np.r_[max(b - 2500, 0):b, s:min(s + 2500, len(cov))]
    return cov[idx[~((idx >= q[1]) & (idx < q[2]))]]


LER_CASES = [((15_000, 15_050), (19_000, 19_050)),  # tested peak outside the candidate's flank
             ((3_000, 3_050), (1_000, 1_050)),  # candidate within 2.5 kb of the chromosome start
             ((15_000, 15_050), (16_000, 16_050)),  # tested peak inside the flank: legacy is already right
             ((39_000, 39_050), (38_000, 38_050)),  # near the chromosome end
             ((100, 150), (1_000, 1_050))]  # tested peak at the chromosome start


@pytest.mark.parametrize("peak,ler", LER_CASES)
def test_fix_ler_flanks(fixes, peak, ler):
    cov = np.random.default_rng(0).integers(0, 5, 40_000).astype(np.int32)
    legacy, q = _ler_background(cov, peak, ler)
    fixes(fix_ler_flanks=True)
    fixed, _ = _ler_background(cov, peak, ler)
    assert np.array_equal(fixed, _intended_background(cov, ler, q))
    if ler == (16_000, 16_050):
        assert np.array_equal(legacy, fixed)
    if ler in ((19_000, 19_050), (1_000, 1_050)) and peak[0] > 1000:
        assert len(legacy) != len(fixed)


@pytest.mark.parametrize("peak,ler", LER_CASES)
def test_ler_candidate_background_legacy_is_unchanged(peak, ler):
    """With the flag off, the helper is the previous inline code."""
    cov = np.random.default_rng(1).integers(0, 5, 40_000).astype(np.int64)
    b, s = ler
    qel, qsr = peak[0] - 1, peak[1] + 1
    w = np.concatenate((cov[b - 2500:b], cov[s:s + 2500]), axis=None)
    if qsr < b:
        w[qel - (b - 2500):qsr - (b - 2500)] = -1
    elif s < qel:
        w[qel - (s - 2500):qsr - (s - 2500)] = -1
    assert np.array_equal(se.ler_candidate_background(cov, b, s, qel, qsr, fixed=False), w[w >= 0])


# ------------------------------------------------------------------ --strict-independent-filtering
def _ifdf(rng, n=3000):
    reads = rng.negative_binomial(2, 0.2, n) + 1
    signal = rng.random(n) < 0.15
    pval = np.where(signal, rng.beta(0.05, 5, n), rng.random(n))
    pval = np.where(reads < 5, rng.random(n), pval)
    df = pd.DataFrame({"chromosome": "c", "start": np.arange(n) * 10, "end": np.arange(n) * 10 + 5,
                       "name": ".", "padj": 1.0, "strand": "+", "reads": reads, "pval": pval})
    df.loc[df.index[:40], ["reads", "pval"]] = [1, 1e-12]  # very low counts, tiny p-values
    return df


def test_strict_independent_filtering(logger):
    """Bourgon et al. 2010: hypotheses removed by the filter are not tested. By default they keep their padj
    from BH over all windows, so some of them are significant."""
    df = _ifdf(np.random.default_rng(5))
    legacy = se.independent_filtering(df.copy(), fdr_target=0.1, logger=logger)
    strict = se.independent_filtering(df.copy(), fdr_target=0.1, logger=logger, strict_independent_filtering=True)
    assert (legacy.padj.iloc[:40] < 0.1).all() and not (strict.padj.iloc[:40] < 0.1).any()
    changed = legacy.padj.to_numpy() != strict.padj.to_numpy()
    assert changed.any() and (strict.padj[changed] == 1.).all()


# ------------------------------------------------------------------ --stringent-pairs-padj
def _peaks(rows, strand):
    out = [dict(chromosome="chr1", start=s, end=e, name=".", padj=padj, strand=strand, reads=reads,
                pval=pval, mu_0=.1, pi_0=.1, mu_1=1., pi_1=.1, var_1=.1, ler_1=0, ler_2=0, ler_3=0,
                non_zeros=3, summit=(s + e) // 2, summit_val=5) for s, e, padj, pval, reads in rows]
    return pd.DataFrame(out, columns=ce.COMMON_HEADER)


def _pair_starts(stringent, **kw):
    # the first pair's minus peak is significant by p-value (0.05) but not by padj (0.5)
    pl = [(1000, 1050, 1e-5, 1e-6, 50), (5000, 5050, 1e-5, 1e-6, 9)]
    mn = [(800, 850, 0.5, 0.05, 40), (4800, 4850, 1e-5, 1e-6, 4)]
    with tempfile.TemporaryDirectory() as d:
        pre = os.path.join(d, "s")
        dfs = {"pl": _peaks(pl, "+"), "mn": _peaks(mn, "-")}
        for lab, df in dfs.items():
            fn = pre + "_%s.bed" % lab
            df.to_csv(fn, sep="\t", index=False, header=False)
            io.index_bed_file(fn)
        (bid, _, _), _ = ce.pair(pre, dfs, 0.1, stringent, **kw)
        return list(pd.read_csv(bid, sep="\t", header=None).iloc[:, 1])


def test_stringent_pairs_padj():
    assert _pair_starts(False) == [800, 4800]
    assert _pair_starts(True) == [800, 4800]  # --stringent-pairs-only filters on the raw p-value
    assert _pair_starts(True, stringent_pairs_padj=True) == [4800]
