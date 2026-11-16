"""Shared definitions for the HDSP problem set.

State (Sigma):  x = (M, E, B),  z = (M, E, B, R),  u = dR/dt
    dM/dt = R - lam*M
    dE/dt = beta_s (R - E),  beta_s = beta_p if R >= E else beta_m
    dB/dt = -a (B-b1)(B-b2)(B-b3) + u_B          (HDSP-01 / HDSP-06)
    eta(M) = 1/(1 + alpha*M)
    K = [B + eta(M)(R - E)] / D,   D = 1 + A + sigma_N^2
Running reward (HDSP-03/04):
    l = [K0 + eta(M)(R-E)]/D - c R - mu|u| - (kappa/2) u^2
"""
from dataclasses import dataclass, replace, asdict
import numpy as np


@dataclass(frozen=True)
class Theta:
    alpha: float = 0.3     # tolerance strength in eta(M)  (a = alpha/lam = 3)
    lam: float = 0.1       # tolerance decay (slow)
    beta_p: float = 1.0    # expectation update when R >= E
    beta_m: float = 0.8    # expectation update when R <  E (losses linger)
    D: float = 1.5         # 1 + A + sigma_N^2
    K0: float = 1.0        # baseline (B frozen) for the control problems
    c: float = 0.005       # running cost of R
    mu: float = 0.01       # |u| cost
    kappa: float = 0.05    # u^2 cost
    v_max: float = 1.0     # rate limit on R
    R_max: float = 1.0
    # bistable B block (HDSP-01 / HDSP-06)
    a: float = 1.0
    b1: float = 0.0
    b2: float = 0.5
    b3: float = 1.0

    def with_(self, **kw):
        return replace(self, **kw)

    def dict(self):
        return asdict(self)


NOMINAL = Theta()


def eta(M, th):
    return 1.0 / (1.0 + th.alpha * M)


def beta_s(R, E, th):
    return np.where(R >= E, th.beta_p, th.beta_m)


def fB(B, th, uB=0.0):
    return -th.a * (B - th.b1) * (B - th.b2) * (B - th.b3) + uB


def K_out(M, E, B, R, th):
    return (B + eta(M, th) * (R - E)) / th.D


# ---------------------------------------------------------------------------
# Periodic trapezoidal input (HDSP-04):  period Tp, duty q, amplitude a,
# ramps at speed v.  R rises 0->a on [0, a/v], holds a until q*Tp,
# falls a->0 on [q*Tp, q*Tp + a/v], stays 0 until Tp.
# ---------------------------------------------------------------------------
def trapezoid(t, Tp, q, amp, v):
    tau = np.mod(t, Tp)
    tr = amp / v
    up = np.clip(tau / tr, 0.0, 1.0) * amp
    down = np.clip(1.0 - (tau - q * Tp) / tr, 0.0, 1.0) * amp
    R = np.where(tau < q * Tp, up, down)
    return R


def periodic_gain(Tp, q, amp, th, v=None, N=3000, settle=3, return_parts=False):
    """Vectorised long-run average reward of the trapezoid policy.

    Returns dG = G - K0/D  (G_constant = K0/D is the best constant policy,
    attained at R = 0, see report).  Arrays Tp,q,amp broadcast together.
    Infeasible combos (ramps do not fit) get -inf.
    """
    v = th.v_max if v is None else v
    Tp, q, amp = np.broadcast_arrays(np.asarray(Tp, float), np.asarray(q, float),
                                     np.asarray(amp, float))
    shape = Tp.shape
    Tp, q, amp = Tp.ravel(), q.ravel(), amp.ravel()
    feas = (q * Tp >= amp / v - 1e-12) & ((1 - q) * Tp >= amp / v - 1e-12) & (amp > 0)
    dt = Tp / N
    lam, al = th.lam, th.alpha

    def rhs(t, M, E):
        R = trapezoid(t, Tp, q, amp, v)
        bs = np.where(R >= E, th.beta_p, th.beta_m)
        return R - lam * M, bs * (R - E), (R - E) / (1 + al * M)

    # pass 1: periodic initial M (linear) ; M(T) from M0=0  ->  M0* = I/(1-e^{-lam T})
    M = np.zeros_like(Tp); E = np.zeros_like(Tp)
    t = np.zeros_like(Tp)
    for k in range(N):
        k1 = rhs(t, M, E)[0]
        k2 = rhs(t + dt / 2, M + dt / 2 * k1, E)[0]
        k3 = rhs(t + dt / 2, M + dt / 2 * k2, E)[0]
        k4 = rhs(t + dt, M + dt * k3, E)[0]
        M = M + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        t = t + dt
    M0 = M / (1 - np.exp(-lam * Tp))
    M = M0.copy(); E = np.full_like(Tp, 0.0)
    acc = np.zeros_like(Tp)
    t = np.zeros_like(Tp)
    for p in range(settle + 1):
        if p == settle:
            acc[:] = 0.0
        for k in range(N):
            a1, b1_, c1 = rhs(t, M, E)
            a2, b2_, c2 = rhs(t + dt / 2, M + dt / 2 * a1, E + dt / 2 * b1_)
            a3, b3_, c3 = rhs(t + dt / 2, M + dt / 2 * a2, E + dt / 2 * b2_)
            a4, b4_, c4 = rhs(t + dt, M + dt * a3, E + dt * b3_)
            M = M + dt / 6 * (a1 + 2 * a2 + 2 * a3 + a4)
            E = E + dt / 6 * (b1_ + 2 * b2_ + 2 * b3_ + b4_)
            acc = acc + dt / 6 * (c1 + 2 * c2 + 2 * c3 + c4)
            t = t + dt
    contrast = acc / (th.D * Tp)
    # exact averages of the cost terms for the trapezoid
    meanR = amp * (q * Tp) / Tp            # area = amp*(q Tp) (ramps symmetric)
    move = (2 * amp * th.mu + th.kappa * amp * v) / Tp
    dG = contrast - th.c * meanR - move
    dG = np.where(feas, dG, -np.inf)
    if return_parts:
        return (dG.reshape(shape), contrast.reshape(shape), (th.c * meanR).reshape(shape),
                move.reshape(shape))
    return dG.reshape(shape)
