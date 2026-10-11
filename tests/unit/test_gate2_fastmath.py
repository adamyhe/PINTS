"""Gate 2 (PINTS_FAST_MATH=1): compiled ZIP EM agrees with the exact (Gate 1) fit within the drift budget."""
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
