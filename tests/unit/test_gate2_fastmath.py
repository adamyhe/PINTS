"""Gate 2 (PINTS_FAST_MATH=1): compiled ZIP EM agrees with the exact (Gate 1) fit within the drift budget,
and LER decisions made from it are certified to equal the exact ones."""
import numpy as np
import pytest

import pints.stats_engine as se

pytestmark = pytest.mark.skipif(not se._HAS_NUMBA, reason="numba not installed")


def _windows(seed, n=400):
    rng = np.random.default_rng(seed)
    out = []
    for i in range(n):
        L = int(rng.integers(1, 6000))
        w = rng.poisson(rng.uniform(0.05, 8), L)
        w[rng.random(L) < rng.uniform(0, 0.995)] = 0
        if i % 5 == 0:
            w[rng.integers(0, L, 3)] += rng.integers(10, 500, 3)
        out.append(w.astype([np.int8, np.int16, np.int32, np.int64][i % 4]) if w.max() < 127 else w.astype(np.int32))
    out += [np.zeros(10, dtype=np.int32), np.ones(10, dtype=np.int32), np.array([0, 3], dtype=np.int16)]
    return out


@pytest.mark.parametrize("seed", [0, 1])
def test_fast_fit_matches_exact(monkeypatch, seed):
    z = se.ZIP()
    exact = [z.fit(w) for w in _windows(seed)]
    fast = [z.fit_fast_math(w) for w in _windows(seed)]
    monkeypatch.setattr(se, "FAST_MATH", True)
    assert [se.fit_background(z, w) for w in _windows(seed)[:50]] == fast[:50]
    monkeypatch.setattr(se, "FAST_MATH", False)
    assert [se.fit_background(z, w) for w in _windows(seed)[:50]] == exact[:50]
    for e, f in zip(exact, fast):
        assert e[4] == f[4]                       # convergence flag
        for a, b in zip((e[0], e[2]), (f[0], f[2])):
            assert abs(a - b) <= 1e-9 * max(abs(a), 1e-300)
        if e[3] is not None and np.isfinite(e[3]):
            assert abs(e[3] - f[3]) <= 1e-9 * max(abs(e[3]), 1.)


def test_flag_off_is_default():
    assert se.FAST_MATH is False or __import__("os").environ.get("PINTS_FAST_MATH") == "1"


@pytest.mark.parametrize("model, window", [
    (se.ZIP(init_mu=2., init_pi=.5), np.array([0, 0, 3, 1, 0, 7], dtype=np.int32)),  # user-supplied init values
    (se.ZIP(debug=True), np.array([0, 0, 3, 1, 0, 7], dtype=np.int32)),              # debug mode
    (se.ZIP(), np.array([0., 0., 3., 1., 0., 7.])),                                  # non-integer dtype
    (se.ZIP(), np.array([[0, 3], [1, 0]], dtype=np.int32)),                          # not 1-D
])
def test_fast_fit_falls_back_to_exact(monkeypatch, model, window):
    called = []
    monkeypatch.setattr(se, "_zip_fit_window_nb", lambda *a: called.append(a))
    assert model.fit_fast_math(window) == model.fit(window)
    assert not called


# ------------------------------------------------------------------ certified LER: decisions equal Gate 1's
def _ler_case(rng):
    """Sparse coverage with short sub-peaks of small integer counts, so tied fits (the hard case) are common."""
    L = int(rng.integers(3000, 30000))
    cov = np.where(rng.random(L) < rng.uniform(0.002, 0.1), rng.poisson(rng.uniform(0.3, 4), L) + 1, 0)
    lines, pos = [], 0
    while True:
        pos += int(rng.integers(5, 400))
        ln = int(rng.integers(1, 30))
        if pos + ln >= L:
            break
        cov[pos:pos + ln] += rng.integers(0, 3, ln)
        lines.append("chrT\t%d\t%d\tp\t%s\t+" % (pos, pos + ln, round(float(cov[pos:pos + ln].mean()), 3)))
        pos += ln
    return cov.astype([np.int8, np.int16, np.int32][int(rng.integers(0, 3))] if cov.max() < 127 else np.int32), lines


@pytest.mark.parametrize("cls", ["bgIQR", "pkIQR"])
@pytest.mark.parametrize("seed", range(6))
def test_certified_ler_equals_exact(monkeypatch, cls, seed):
    from pints.calling_engine import SubpeakIndex
    rng = np.random.default_rng(seed)
    cov, lines = _ler_case(rng)
    index, z, iqr = SubpeakIndex("chrT", lines), se.ZIP(), getattr(se, cls)
    used = {"many": 0, "single": 0}
    for name, key in (("_certified_ler_mus", "many"), ("_certified_ler_single", "single")):
        def wrap(*a, _f=getattr(se, name), _k=key):
            r = _f(*a)
            used[_k] += r is not None
            return r
        monkeypatch.setattr(se, name, wrap)
    for _ in range(40):
        cs = int(rng.integers(0, len(cov) - 50))
        ce = cs + int(rng.integers(5, 50))
        for f in (5000, 2500, 500):
            args = dict(stat_tester=z, bed_handler=index, chromosome="chrT", query_start_left=max(cs - f, 0),
                        query_end_left=max(cs - 1, 0), query_start_right=min(ce + 1, len(cov)),
                        query_end_right=min(ce + f, len(cov)), small_window_threshold=5,
                        peak_in_bg_threshold=float(rng.choice([0.5, 1., 2.])), coverage_info=cov, fdr_target=0.1,
                        disable_ler=False, enable_eler=bool(rng.random() < .5),
                        peak_threshold=float(rng.choice([2., 5., 20.])))
            out = {}
            for fast in (False, True):
                monkeypatch.setattr(se, "FAST_MATH", fast)
                cache = {}
                env, n = iqr.remove_peaks_in_local_env(cache=cache, **args)
                out[fast] = (env.dtype, env.tobytes(), n, sorted(cache.items()))
            assert out[True] == out[False]
    assert used["single"] > 0 and (used["many"] > 0 or cls == "pkIQR")


@pytest.mark.parametrize("seed", [0, 1])
def test_ler_kernel_matches_exact_fits(seed):
    rng = np.random.default_rng(seed)
    z = se.ZIP()
    n_exact = 0
    for it in range(300):
        dt = [np.int8, np.int16, np.int32, np.int64][it % 4]
        la, lb = int(rng.integers(0, 4000)), int(rng.integers(0, 4000))
        dens = rng.uniform(0, 0.3) if it % 7 else 0.
        a = np.where(rng.random(la) < dens, rng.poisson(rng.uniform(0.5, 6), la) + 1, 0).astype(dt)
        b = np.where(rng.random(lb) < dens, rng.poisson(rng.uniform(0.5, 6), lb) + 1, 0).astype(dt)
        nc = int(rng.integers(1, 8))
        st = np.sort(rng.integers(0, la + 1, nc)).astype(np.int64)
        en = np.minimum(st + rng.integers(0, 400, nc), la).astype(np.int64)
        mus, status, cert, exact = se._ler_mus_nb(a, b, st, en, 1e-4, 1000, se._EPS, 1 << 20)
        for c in range(nc):
            le = a.astype(np.int64)
            le[st[c]:en[c]] = -1
            w = np.concatenate((le, b.astype(np.int64)))
            ref = z.fit(w[w >= 0].astype(dt))[0]
            assert status[c] >= 0
            if exact[c]:
                assert mus[c] == ref  # bit for bit
                n_exact += 1
            elif cert[c]:
                assert abs(mus[c] - ref) <= 1e-9 * abs(ref)
    assert n_exact > 0


def test_ler_kernel_rejects_negative_values():
    a = np.array([0, 1, -1, 2], dtype=np.int16)
    _, status, _, _ = se._ler_mus_nb(a, a, np.array([0]), np.array([1]), 1e-4, 1000, se._EPS, 1 << 20)
    assert (status == -1).all()


# ------------------------------------------------------------------ LRT (sensitive mode)
def test_chi2_sf_df2_matches_scipy():
    from scipy.stats import chi2
    for x in np.concatenate([[-5., 0.], np.geomspace(1e-300, 1e-3, 50), np.linspace(1e-3, 1400, 3000)]):
        assert abs(se._chi2_sf_df2(float(x)) - chi2.sf(x, 2)) <= 1e-12 * chi2.sf(x, 2)
    assert se._chi2_sf_df2(1500.) == 0. and np.isnan(se._chi2_sf_df2(float("nan")))


@pytest.mark.parametrize("seed", [0, 1])
def test_fast_lrt_matches_exact(seed):
    rng = np.random.default_rng(seed)
    z = se.ZIP()
    bgs = _windows(seed, n=150)
    for i, bg in enumerate(bgs):
        L = int(rng.integers(1, 300))
        pk = np.where(rng.random(L) < rng.uniform(0., 1.), rng.poisson(rng.uniform(0.5, 20), L), 0).astype(bg.dtype)
        (we, pe), (wf, pf) = z.lrt(pk, bg), z.lrt_fast_math(pk, bg)
        assert abs(we - wf) <= 1e-9 * max(abs(we), 1.)
        if pe > 1e-300 or pf > 1e-300:
            assert abs(pe - pf) <= 1e-9 * max(pe, pf)


def test_fast_lrt_falls_back_to_exact():
    z = se.ZIP(debug=True)
    a, b = np.array([0, 3, 1, 0], dtype=np.int32), np.array([0, 0, 1, 2, 0, 0], dtype=np.int32)
    assert z.lrt_fast_math(a, b) == z.lrt(a, b)
    z = se.ZIP()
    assert z.lrt_fast_math(a.astype(float), b.astype(float)) == z.lrt(a.astype(float), b.astype(float))
