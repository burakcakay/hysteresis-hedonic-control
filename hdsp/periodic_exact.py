"""Semi-analytic long-run average of the trapezoid policy (HDSP-04).

One period of R (ramp speed v, amplitude a, duty q, period T):
  [0, a/v]          R = v t           (slope +v)
  [a/v, qT]         R = a             (slope 0)
  [qT, qT + a/v]    R = a - v(t-qT)   (slope -v)
  [qT + a/v, T]     R = 0             (slope 0)
On a piece with R = R0 + s*tau and fixed beta, everything is closed-form:
  d(tau) = R - E = s/beta + (d0 - s/beta) e^{-beta tau}
  M(tau) = M0 e^{-lam tau} + R0 (1-e^{-lam tau})/lam + s [tau/lam - (1-e^{-lam tau})/lam^2]
The switch R = E (beta_+ <-> beta_-) happens at most once per piece, at
  tau* = ln(1 - d0 beta / s) / beta   (only if d0 and s have opposite signs).
Only J = int eta(M)(R-E) dt needs quadrature (Gauss-Legendre on unit chunks).
"""
import numpy as np

_GL_X, _GL_W = np.polynomial.legendre.leggauss(24)


def _pieces(T, q, a, v, vd=None):
    vd = v if vd is None else vd
    tu, td = a / v, a / vd
    return [(tu, 0.0, v), (q * T - tu, a, 0.0), (td, a, -vd), ((1 - q) * T - td, 0.0, 0.0)]


def _sub(d0, s, dur, th):
    """Split a piece of length dur into sub-pieces of constant beta.
    Returns list of (tau_start, length, beta, d_start) and d at the end."""
    out = []
    tau = 0.0
    d = d0
    while True:
        if d > 0 or (d == 0 and s > 0):
            beta = th.beta_p
        elif d < 0 or (d == 0 and s < 0):
            beta = th.beta_m
        else:  # d == 0 and s == 0: stays on R = E
            beta = th.beta_p
        rem = dur - tau
        cross = None
        if s != 0 and d != 0 and np.sign(d) != np.sign(s):
            ts = np.log(1 - d * beta / s) / beta
            if ts < rem:
                cross = ts
        if cross is None:
            out.append((tau, rem, beta, d))
            d_end = s / beta + (d - s / beta) * np.exp(-beta * rem)
            return out, d_end
        out.append((tau, cross, beta, d))
        tau += cross
        d = 0.0


def _M(M0, R0, s, tau, lam):
    e = np.exp(-lam * tau)
    return M0 * e + R0 * (1 - e) / lam + s * (tau / lam - (1 - e) / lam ** 2)


def _period(M0, E0, T, q, a, v, th, integrate=True, vd=None):
    M, E, J = M0, E0, 0.0
    for dur, R0, s in _pieces(T, q, a, v, vd):
        if dur <= 0:
            continue
        subs, d_end = _sub(R0 - E, s, dur, th)
        if integrate:
            for (t0, L, beta, ds) in subs:
                if L <= 0:
                    continue
                n = max(1, int(np.ceil(L / 1.0)))
                edges = np.linspace(0, L, n + 1)
                mid = 0.5 * (edges[1:] + edges[:-1])[:, None]
                half = 0.5 * (edges[1:] - edges[:-1])[:, None]
                tl = (mid + half * _GL_X[None, :]).ravel()        # local time in sub-piece
                w = (half * _GL_W[None, :]).ravel()
                d = s / beta + (ds - s / beta) * np.exp(-beta * tl)
                Mt = _M(M, R0, s, t0 + tl, th.lam)
                J += np.sum(w * d / (1 + th.alpha * Mt))
        M = _M(M, R0, s, dur, th.lam)
        E = R0 + s * dur - d_end
    return M, E, J


def G_exact(T, q, a, th, v=None, tol=1e-15, max_periods=200, vd=None):
    """dG = G - K0/D of the trapezoid (T, q, a, up-ramp speed v, down-ramp speed vd
    (default vd = v)); -inf if infeasible.  q*T = time from start of ramp-up to start
    of ramp-down."""
    v = th.v_max if v is None else v
    vd = v if vd is None else vd
    if a <= 0 or q * T < a / v - 1e-12 or (1 - q) * T < a / vd - 1e-12:
        return -np.inf
    # mean R and movement cost per period
    meanR = (a * (q * T - a / v) + a * a / (2 * v) + a * a / (2 * vd)) / T
    move = (2 * a * th.mu + th.kappa * a * (v + vd) / 2) / T
    if T > 1e5:  # contrast term is O(1/T): use the asymptote (optimisers only probe this region)
        return -th.c * meanR - move
    MT0, _, _ = _period(0.0, 0.0, T, q, a, v, th, integrate=False, vd=vd)
    M0 = MT0 / (1 - np.exp(-th.lam * T))
    E0 = 0.0
    for _ in range(max_periods):
        _, E1, _ = _period(M0, E0, T, q, a, v, th, integrate=False, vd=vd)
        if abs(E1 - E0) < tol:
            E0 = E1
            break
        E0 = E1
    _, _, J = _period(M0, E0, T, q, a, v, th, vd=vd)
    contrast = J / (th.D * T)
    return contrast - th.c * meanR - move
