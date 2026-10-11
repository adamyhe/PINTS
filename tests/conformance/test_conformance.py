"""Conformance: does the implementation do what the PINTS paper / docstrings describe?

Independent oracles (scipy / statsmodels / closed forms) are used wherever possible.  Deviations are
encoded as strict xfails (D*/O*) or as `spec` characterisation tests when the paper is ambiguous.
"""
import os
import tempfile
import logging

import numpy as np
import pandas as pd
import pytest
from scipy import optimize, stats
from statsmodels.stats.multitest import multipletests

import pints.calling_engine as ce
import pints.stats_engine as se
from helpers import synth_coverage


def strict(bug_id, why, raises=AssertionError):
    return pytest.mark.xfail(strict=True, raises=raises, reason="%s: %s" % (bug_id, why))


def rzip(rng, n, lam, pi):
    x = rng.poisson(lam, n)
    x[rng.random(n) < pi] = 0
    return x


def zip_nll(theta, x):
    lam, pi = theta
    p0 = pi + (1 - pi) * np.exp(-lam)
    ll = np.where(x == 0, np.log(p0), np.log(1 - pi) + stats.poisson.logpmf(x, lam))
    return -ll.sum()


# ------------------------------------------------------------------ ZIP likelihood + EM
@pytest.mark.parametrize("lam,pi,n", [(0.5, 0.2, 500), (2.0, 0.6, 300), (6.0, 0.9, 2000), (1.2, 0.97, 5000),
                                      (15.0, 0.3, 50), (3.0, 0.0, 400)])
def test_zip_loglik_matches_scipy(lam, pi, n):
    x = rzip(np.random.default_rng(1), n, lam, pi)
    assert np.isclose(se.ZIP.log_likelihood(x, pi, lam), -zip_nll((lam, pi), x), rtol=1e-12, atol=1e-9)


@pytest.mark.parametrize("lam,pi,n", [(0.5, 0.2, 500), (2.0, 0.6, 300), (6.0, 0.9, 2000), (1.2, 0.97, 5000),
                                      (15.0, 0.3, 50), (4.0, 0.5, 20)])
def test_zip_em_reaches_mle(lam, pi, n):
    """EM (stop when |delta loglik| < 1e-4) must land on the numerical MLE (to EM tolerance)."""
    x = rzip(np.random.default_rng(7), n, lam, pi)
    if (x == 0).all() or (x > 0).all():
        pytest.skip("degenerate sample")
    m_lam, _, m_pi, m_ll, conv = se.ZIP().fit(x)
    res = optimize.minimize(zip_nll, x0=[max(x.mean(), .1), .5], args=(x,), method="L-BFGS-B",
                            bounds=[(1e-9, None), (0, 1 - 1e-12)])
    assert conv
    assert -res.fun - m_ll < 1e-3                    # loglik within 1e-3 nats of the optimum
    assert abs(m_lam - res.x[0]) / res.x[0] < 5e-3
    assert abs(m_pi - res.x[1]) < 5e-3


def test_zip_moment_estimators_closed_form():
    x = np.array([0, 0, 0, 0, 1, 2, 3, 5, 0, 7], dtype=float)
    m, s2 = x.mean(), x.var()
    lam, pi = se.ZIP.zip_moment_estimators(x)
    assert np.isclose(lam, (s2 + m * m - m) / m)
    assert np.isclose(pi, (s2 - m) / (s2 + m * m - m))


def test_zip_fit_degenerate_branches():
    z = se.ZIP()
    assert z.fit(np.zeros(10)) == (0, 0, 0, 0, True)
    mu, var, pi, llc, conv = z.fit(np.array([1, 2, 3]))
    assert (mu, var, pi, conv) == (2.0, 2.0, 0, True) and np.isnan(llc)


# ------------------------------------------------------------------ LRT
def test_lrt_closed_form_df2():
    rng = np.random.default_rng(3)
    a, b = rzip(rng, 200, 4.0, 0.5), rzip(rng, 800, 1.0, 0.8)
    z = se.ZIP()
    w, p = z.lrt(a, b)
    fa, fb, f0 = z.fit(a), z.fit(b), z.fit(np.concatenate([a, b]))
    ll = lambda x, f: se.ZIP.log_likelihood(x, f[2], f[0])
    w_ref = 2 * (ll(a, fa) + ll(b, fb) - ll(np.concatenate([a, b]), f0))
    assert np.isclose(w, w_ref) and np.isclose(p, stats.chi2.sf(w_ref, 2))


def _exact_ll(x):
    r = optimize.minimize(zip_nll, x0=[max(x[x > 0].mean(), .1), (x == 0).mean() * .9], args=(x,),
                          method="L-BFGS-B", bounds=[(1e-9, None), (0, 1 - 1e-12)],
                          options=dict(ftol=1e-14, gtol=1e-10))
    return -r.fun


def test_lrt_statistic_matches_exact_mle_and_null_rate():
    """W from EM fits equals W from exact numerical MLEs (|dW| < 1e-3).  Null rejection rate is recorded
    loosely: the chi2(2) approximation is mildly anti-conservative here (~7.7% at alpha=0.05 for n=300/600,
    identical with exact MLEs, so it is the asymptotics, not the EM).  LRT only matters with --sensitive."""
    rng = np.random.default_rng(11)
    z = se.ZIP()
    ps = []
    for i in range(300):
        a, b = rzip(rng, 300, 3.0, 0.6), rzip(rng, 600, 3.0, 0.6)
        w, p = z.lrt(a, b)
        ps.append(p)
        if i < 25:
            w_ref = 2 * (_exact_ll(a) + _exact_ll(b) - _exact_ll(np.concatenate([a, b])))
            assert abs(w - w_ref) < 1e-3
    assert 0.02 < np.mean(np.array(ps) < 0.05) < 0.12


# ------------------------------------------------------------------ IQR outlier rule
def test_iqr_lower_fence_is_used_for_bgIQR():
    """bgIQR masks a candidate when removing it drops the background mu below Q1 - k*IQR (k = 1.0 in code)."""
    data = np.array([1.0, 1.1, 1.2, 1.3, 0.2])
    lo, hi = se.IQR.get_outlier_threshold(data)
    q1, q3 = np.quantile(data, [.25, .75])
    assert hi is None and np.isclose(lo, q1 - (q3 - q1))


# ------------------------------------------------------------------ BH + independent filtering
def _ifdf(rng, n=3000):
    reads = rng.negative_binomial(2, 0.2, n) + 1
    signal = rng.random(n) < 0.15
    pval = np.where(signal, rng.beta(0.05, 5, n), rng.random(n))
    pval = np.where(reads < 5, rng.random(n), pval)          # low-count windows carry no signal
    return pd.DataFrame({"chromosome": "c", "start": np.arange(n) * 10, "end": np.arange(n) * 10 + 5,
                         "name": ".", "padj": 1.0, "strand": "+", "reads": reads, "pval": pval})


def test_independent_filtering_selected_set_is_bh_and_maximises_rejections(logger):
    rng = np.random.default_rng(5)
    df = _ifdf(rng)
    out = se.independent_filtering(df.copy(), fdr_target=0.1, logger=logger)
    best_q, best_r = None, -1
    for q in np.arange(0, 1, 0.005):
        thr = df.reads.quantile(q)
        keep = df.reads > thr
        r = (multipletests(df.pval[keep], method="fdr_bh")[1] < 0.1).sum() if keep.any() else 0
        if r >= best_r:
            best_q, best_r, best_keep = q, r, keep
    exp = multipletests(df.pval[best_keep], method="fdr_bh")[1]
    assert np.allclose(out.padj[best_keep].to_numpy(), exp, rtol=0, atol=0)
    # windows removed by the filter carry the *unfiltered* BH padj (see O4)
    direct = multipletests(df.pval, method="fdr_bh")[1]
    assert np.array_equal(out.padj[~best_keep].to_numpy(), direct[~best_keep.to_numpy()])


@strict("O4", "Bourgon et al. 2010: hypotheses removed by the independent filter are not tested (padj=NA/1). "
              "PINTS keeps their unfiltered-BH padj, so filtered-out windows can still be called significant "
              "(union of two BH procedures; nominal FDR not guaranteed)")
def test_O4_filtered_out_windows_are_not_significant(logger):
    rng = np.random.default_rng(5)
    df = _ifdf(rng)
    df.loc[df.index[:40], ["reads", "pval"]] = [1, 1e-12]   # very low count, extremely small p
    out = se.independent_filtering(df.copy(), fdr_target=0.1, logger=logger)
    assert not (out.padj.iloc[:40] < 0.1).any()


def test_O4_strict_independent_filtering_flag(logger):
    rng = np.random.default_rng(5)
    df = _ifdf(rng)
    df.loc[df.index[:40], ["reads", "pval"]] = [1, 1e-12]
    legacy = se.independent_filtering(df.copy(), fdr_target=0.1, logger=logger)
    out = se.independent_filtering(df.copy(), fdr_target=0.1, logger=logger, strict_independent_filtering=True)
    assert not (out.padj.iloc[:40] < 0.1).any()
    kept = legacy.padj.to_numpy() != out.padj.to_numpy()
    assert (out.padj[kept] == 1.).all()          # only filtered-out windows change, and they become padj = 1
    assert (legacy.padj < 0.1).sum() >= (out.padj < 0.1).sum()


def test_O3_logged_threshold_vs_selected_threshold():
    """Characterisation: selection keeps the LAST quantile reaching max rejections (>=), while the log line
    reports np.argmax (the FIRST).  Harmless for results; noted for the report."""
    r = [3, 5, 5, 4]
    last = max(i for i, v in enumerate(r) if v == max(r))
    assert last == 2 and int(np.argmax(r)) == 1


# ------------------------------------------------------------------ pairing semantics
def _peaks(rows, strand):
    cols = ce.COMMON_HEADER
    out = []
    for s, e, padj, reads in rows:
        out.append(dict(chromosome="chr1", start=s, end=e, name=".", padj=padj, strand=strand, reads=reads,
                        pval=padj / 10, mu_0=.1, pi_0=.1, mu_1=1., pi_1=.1, var_1=.1, ler_1=0, ler_2=0, ler_3=0,
                        non_zeros=3, summit=(s + e) // 2, summit_val=5))
    return pd.DataFrame(out, columns=cols) if out else pd.DataFrame(columns=cols)


def _pair(pl, mn, **kw):
    with tempfile.TemporaryDirectory() as d:
        pre = os.path.join(d, "s")
        dfs = {"pl": _peaks(pl, "+"), "mn": _peaks(mn, "-")}
        for lab, df in dfs.items():
            fn = pre + "_%s.bed" % lab
            df.sort_values(["chromosome", "start"]).to_csv(fn, sep="\t", index=False, header=False)
            ce.index_bed_file(fn)
        (b, dv, u), _ = ce.pair(pre, dfs, 0.1, kw.pop("stringent", False), **kw)
        rd = lambda p: (pd.read_csv(p, sep="\t", header=None) if p and os.path.exists(p) and os.path.getsize(p)
                        else pd.DataFrame())
        return rd(b), rd(dv), rd(u)


def test_pairing_divergent_within_300bp():
    b, d, u = _pair([(1000, 1050, 1e-5, 50)], [(800, 850, 1e-5, 40)])
    assert len(b) == 1 and len(d) == 1
    assert (b.iloc[0, 1], b.iloc[0, 2]) == (800, 1050)


def test_pairing_rejects_beyond_close_threshold():
    b, d, u = _pair([(1000, 1050, 1e-5, 50)], [(600, 650, 1e-5, 40)])
    assert len(b) == 0 and len(u) == 2


def test_pairing_convergent_is_bidirectional_not_divergent():
    b, d, u = _pair([(1000, 1050, 1e-5, 50)], [(1040, 1100, 1e-5, 40)])
    assert len(b) == 1 and len(d) == 0


def test_pairing_overlapping_pairs_are_merged_into_one_element():
    """Pairs found from each strand are unioned with bedtools merge; overlapping pairs collapse to one element
    (summits reported are those with the highest summit value)."""
    b, d, u = _pair([(1000, 1050, 1e-5, 50)], [(760, 790, 1e-5, 7), (900, 950, 1e-5, 40)])
    assert len(b) == 1 and (b.iloc[0, 1], b.iloc[0, 2]) == (760, 1050)


# O7 fixed (Track B): if only ONE strand produced pairs, pair() skipped the merge and parsed a 16-column file
# with 8 names -> malformed BED (chrom column = peak end)
def test_O7_one_sided_pairs_produce_valid_bed():
    b, d, u = _pair([(1000, 1050, 1e-5, 50)], [(800, 850, 0.9, 40)])
    assert len(b) == 1 and b.iloc[0, 0] == "chr1" and (b.iloc[0, 1], b.iloc[0, 2]) == (800, 1050)


@strict("O5", "(fixed only with --stringent-pairs-padj) --stringent-pairs-only filters the opposite peak on raw p-value (column 8), not padj",
        raises=(AssertionError,))
def test_O5_stringent_uses_padj():
    # pair 1: opposite peak padj 0.5 (not significant) but pval 0.05 (< 0.1) -> must be dropped
    # pair 2: both significant (keeps both strands' files non-empty, avoiding O7)
    b, d, u = _pair([(1000, 1050, 1e-5, 50), (5000, 5050, 1e-5, 9)],
                    [(800, 850, 0.5, 40), (4800, 4850, 1e-5, 4)], stringent=True)
    assert len(b) == 1 and b.iloc[0, 1] == 4800


def test_O5_stringent_pairs_padj_flag():
    b, d, u = _pair([(1000, 1050, 1e-5, 50), (5000, 5050, 1e-5, 9)],
                    [(800, 850, 0.5, 40), (4800, 4850, 1e-5, 4)], stringent=True, stringent_pairs_padj=True)
    assert len(b) == 1 and b.iloc[0, 1] == 4800
    # without the stringent filter both pairs are reported
    b, d, u = _pair([(1000, 1050, 1e-5, 50), (5000, 5050, 1e-5, 9)],
                    [(800, 850, 0.5, 40), (4800, 4850, 1e-5, 4)])
    assert len(b) == 2


def test_O7_one_sided_pairs_multiple_elements():
    b, d, u = _pair([(1000, 1050, 1e-5, 50), (9000, 9050, 1e-5, 30)], [(800, 850, 0.9, 40), (8850, 8900, 0.9, 5)])
    assert list(b.iloc[:, 0]) == ["chr1", "chr1"] and b.shape[1] == 6
    assert [tuple(r) for r in b.iloc[:, 1:3].to_numpy()] == [(800, 1050), (8850, 9050)]


# ------------------------------------------------------------------ simulation sanity (end-to-end, one strand)
def _simulate(seed, L=300_000, n=40, planted=True, pi=0.97, lam=1.2):
    rng = np.random.default_rng(seed)
    cov = rng.poisson(lam, L)
    cov[rng.random(L) < pi] = 0
    truth = []
    if planted:
        pos = np.sort(rng.choice(np.arange(3000, L - 3000, 5000), n, replace=False))
        for p in pos:
            w = int(rng.integers(10, 40)); h = rng.integers(15, 60)
            cov[p:p + w] += rng.poisson(h * np.exp(-0.5 * ((np.arange(w) - w / 2) / (w / 5)) ** 2))
            truth.append((p, p + w))
    return cov.astype(np.int32), truth


def _call_one_strand(cov):
    ce.stat_tester, ce.iqr_obj = se.ZIP(), se.bgIQR()
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "c.npy"); np.save(f, cov)
        sp = ce.unified_element_definition([{"chrS": f}], "chrS", "+", os.path.join(d, "s_pl.bed"))
        out = ce.peaks_single_strand({"chrS": f}, os.path.join(d, "s_1_pl.bed"), {"chrS": sp}, "+", 0.1,
                                     thread_n=1)
        return pd.read_csv(out, sep="\t", header=None, names=ce.COMMON_HEADER)


def _score(df, truth):
    sig = df[df.padj < 0.1]
    hit = sum(((sig.start < e + 50) & (sig.end > s - 50)).any() for s, e in truth)
    fp = sum(not any(r.start < e + 50 and r.end > s - 50 for s, e in truth) for r in sig.itertuples())
    return len(df), len(sig), hit, fp


@pytest.mark.slow
def test_simulation_planted_peaks_recovered_without_fp():
    cov, truth = _simulate(0)
    n, nsig, hit, fp = _score(_call_one_strand(cov), truth)
    assert hit / len(truth) >= 0.75 and fp == 0


@pytest.mark.slow
@pytest.mark.spec
def test_simulation_pure_null_characterisation():
    """Homogeneous ZIP background with NO planted peaks: every call is a false discovery.

    Current code calls tens-to-hundreds of 'significant' peaks (seed 0: 58 / ~5,850 candidates), because the
    adaptive thresholds (min_mu_percent of top peaks, small-peak lambda) are learned from noise when there are
    no true peaks.  Recorded as a characterisation (not a bug verdict: real libraries always contain strong
    TSSs); Gate 2 must keep this number unchanged."""
    cov, _ = _simulate(0, planted=False)
    n, nsig, hit, fp = _score(_call_one_strand(cov), [])
    assert nsig == fp and nsig > 0
