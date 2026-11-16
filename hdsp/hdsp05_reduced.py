"""HDSP-05: reduced (slow-tolerance) hysteresis problem.

Singular-perturbation limit lam -> 0 with a = alpha/lam and c_hat = c/lam fixed
(beta, mu, kappa, v_max fixed).  In slow time tau = lam t, S = lam M obeys
dS/dtau = R - S, and E relaxes on the fast scale, so the hedonic term
eta(M)(R-E) only contributes at jumps of R:
    upward jump  dR>0 at S :  + eta(S) dR /(D beta_+)      (exact: int (R-E)dt = dR/beta_+)
    downward jump dR<0 at S:  - eta(S)|dR|/(D beta_-)
    switching cost per jump:  h |dR|,   h = mu + kappa v_max / 2   (ramp at v_max)
Per-cycle net value (relative to the passive baseline K0/D):
    N = -c r T_A + (r/D)[eta(S_L)/beta_+ - eta(S_H)/beta_-] - 2 h r
      = int_{S_L}^{S_H} g(s)/(lam (r-s)) ds + b(S_L) - 2 h r
    g(s) = -c r - lam r (r-s) eta'(s)/(D beta_-),  b(S) = (r/D)(1/beta_+ - 1/beta_-) eta(S)
"""
import json, warnings
warnings.filterwarnings("ignore")
import numpy as np
from scipy.optimize import minimize, brentq
from model import NOMINAL
from periodic_exact import G_exact

th = NOMINAL
A_ = th.alpha / th.lam          # a = alpha/lam
CH = th.c / th.lam              # c_hat
H1 = th.mu + th.kappa * th.v_max / 2


def eta(s): return 1 / (1 + A_ * s)
def deta(s): return -A_ / (1 + A_ * s) ** 2


def J_red(SL, SH, r, lam=th.lam):
    """Reduced average reward (per unit t) of the 2-level hysteresis loop."""
    TA = np.log((r - SL) / (r - SH)) / lam
    TP = np.log(SH / SL) / lam
    N = -CH * lam * r * TA + (r / th.D) * (eta(SL) / th.beta_p - eta(SH) / th.beta_m) - 2 * H1 * r
    return N / (TA + TP)


def g_fun(s, r, lam=th.lam):
    return -CH * lam * r - lam * r * (r - s) * deta(s) / (th.D * th.beta_m)


def b_prime(s, r):
    return (r / th.D) * (1 / th.beta_p - 1 / th.beta_m) * deta(s)


out = {"a=alpha/lam": A_, "c_hat": CH, "h_per_unit_jump": H1}

# --- 1. optimise the 2-level loop (r = R_max and r free) ---------------------
def opt2(r_free=False, lam=th.lam):
    def f(p):
        r = th.R_max / (1 + np.exp(-p[2])) if r_free else th.R_max
        SL = r / (1 + np.exp(-p[0]))
        SH = SL + (r - SL) / (1 + np.exp(-p[1]))
        return -J_red(SL, SH, r, lam)
    best = None
    for p0 in ([-1.0, 0.0, 3.0], [-2.0, 0.5, 2.0], [0.0, -0.5, 4.0], [-3, 1, 1]):
        rr = minimize(f, p0[:3], method="Nelder-Mead", options=dict(xatol=1e-11, fatol=1e-16, maxiter=6000))
        if best is None or rr.fun < best.fun:
            best = rr
    p = best.x
    r = th.R_max / (1 + np.exp(-p[2])) if r_free else th.R_max
    SL = r / (1 + np.exp(-p[0])); SH = SL + (r - SL) / (1 + np.exp(-p[1]))
    return SL, SH, r, -best.fun

SL, SH, r, J = opt2()
phi = lambda s: s * g_fun(s, r)
out["two_level_r=Rmax"] = dict(S_L=SL, S_H=SH, J=J, J_over_lam=J / th.lam,
                               T_A=np.log((r - SL) / (r - SH)) / th.lam, T_P=np.log(SH / SL) / th.lam)
# first-order conditions (spec form and corrected form)
out["FOC_check"] = dict(
    SH_g_SH_minus_rJ=SH * g_fun(SH, r) - r * J,
    SL_g_SL_minus_rJ_spec=SL * g_fun(SL, r) - r * J,
    SL_corrected_minus_rJ=SL * g_fun(SL, r) - th.lam * SL * (r - SL) * b_prime(SL, r) - r * J)
# second-order: Hessian of J at optimum (finite differences)
def Jp(a_, b_): return J_red(a_, b_, r)
h = 1e-4
H = np.array([[(Jp(SL + h, SH) - 2 * J + Jp(SL - h, SH)) / h**2,
               (Jp(SL + h, SH + h) - Jp(SL + h, SH - h) - Jp(SL - h, SH + h) + Jp(SL - h, SH - h)) / (4 * h * h)],
              [0, (Jp(SL, SH + h) - 2 * J + Jp(SL, SH - h)) / h**2]])
H[1, 0] = H[0, 1]
# predicted diagonal from theory: d2J/dSH2 = phi'(SH)/(lam SH (r-SH) T)
T = (np.log((r - SL) / (r - SH)) + np.log(SH / SL)) / th.lam
dphi = lambda s: (phi(s + 1e-6) - phi(s - 1e-6)) / 2e-6
out["hessian"] = dict(H=H.tolist(), det=float(np.linalg.det(H)), eig=np.linalg.eigvalsh(H).tolist(),
                      theory_JHH=dphi(SH) / (th.lam * SH * (r - SH) * T),
                      phi_prime_SH=dphi(SH))
# shape of phi on (0,r): number of sign changes of phi'
ss = np.linspace(1e-4, r - 1e-4, 20001)
dp = np.gradient(phi(ss), ss)
out["phi_shape"] = dict(n_sign_changes_phi_prime=int(np.sum(np.diff(np.sign(dp)) != 0)),
                        argmax_phi=float(ss[np.argmax(phi(ss))]), max_phi=float(phi(ss).max()))

print("two-level", out["two_level_r=Rmax"], out["FOC_check"], flush=True)
SLf, SHf, rf, Jf = opt2(r_free=True)
out["two_level_r_free"] = dict(S_L=SLf, S_H=SHf, r=rf, J=Jf)

# --- 2. staircase policies (m levels) : counterexample search ----------------
def J_stair(levels, durs, lam=th.lam):
    """levels R_1..R_m (R_m may be 0), durations in slow time; periodic."""
    m = len(levels)
    # periodic S0:  S_end = e^{-sum} S0 + sum_k R_k(1-e^{-tau_k}) e^{-sum_{j>k} tau_j}
    S = 0.0; decay = 1.0
    for Rk, tk in zip(levels, durs):
        S = Rk + (S - Rk) * np.exp(-tk); decay *= np.exp(-tk)
    S0 = S / (1 - decay)
    S = S0; C = 0.0; cost = 0.0
    prev = levels[-1]
    for Rk, tk in zip(levels, durs):
        dR = Rk - prev
        C += eta(S) * (dR / th.beta_p if dR > 0 else dR / th.beta_m) / th.D - H1 * abs(dR)
        cost += CH * lam * Rk * tk / lam  # c * R * (tk/lam)
        S = Rk + (S - Rk) * np.exp(-tk)
        prev = Rk
    Ttot = np.sum(durs) / lam
    return (C - cost) / Ttot

# sanity: 2-level staircase reproduces J_red
tA = np.log((r - SL) / (r - SH)); tP = np.log(SH / SL)
out["stair_sanity"] = J_stair([r, 0.0], [tA, tP]) - J

def opt_stair(m, starts=10, seed=0):
    rng = np.random.default_rng(seed)
    best = (-np.inf, None)
    def f(p):
        lv = list(th.R_max / (1 + np.exp(-p[:m - 1]))) + [0.0]
        du = np.exp(p[m - 1:])
        return -J_stair(lv, du)
    for _ in range(starts):
        p0 = np.concatenate([rng.normal(1, 2, m - 1), np.log(rng.uniform(0.1, 3, m))])
        rr = minimize(f, p0, method="Nelder-Mead", options=dict(xatol=1e-9, fatol=1e-14, maxiter=6000))
        if -rr.fun > best[0]:
            best = (-rr.fun, rr.x)
    p = best[1]
    return best[0], (th.R_max / (1 + np.exp(-p[:m - 1]))).tolist(), np.exp(p[m - 1:]).tolist()

stairs = {}
for m in (2, 3, 4):
    Jm, lv, du = opt_stair(m)
    stairs[m] = dict(J=Jm, levels=lv, slow_durations=du, gain_vs_2level=Jm - Jf)
out["staircase"] = stairs
print("staircase", json.dumps(stairs, default=float), flush=True)

# --- 3. singular perturbation: full trapezoid optimum vs reduced, lam -> 0 ----
def full_opt(lam, guess):
    thl = th.with_(lam=lam, alpha=A_ * lam, c=CH * lam)
    def f(p):
        T, q = np.exp(p[0]), 1 / (1 + np.exp(-p[1]))
        return -G_exact(T, q, th.R_max, thl)
    rr = minimize(f, guess, method="Nelder-Mead", options=dict(xatol=1e-8, fatol=1e-14, maxiter=600))
    return -rr.fun, np.exp(rr.x[0]), 1 / (1 + np.exp(-rr.x[1]))

sp = []
for lam in (0.4, 0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125, 0.0015625):
    _, _, _, Jr = opt2(lam=lam)          # reduced J (exactly proportional to lam)
    Ttot = (tA + tP) / lam
    Jfull, Tp, q = full_opt(lam, [np.log(Ttot), np.log(tA / tP)])
    sp.append(dict(lam=lam, J_full=Jfull, J_red=Jr, Jfull_over_lam=Jfull / lam, Jred_over_lam=Jr / lam,
                   rel_err=abs(Jfull - Jr) / abs(Jr), Tp=Tp, q=q, Tp_red=Ttot, q_red=tA / (tA + tP)))
    print(sp[-1], flush=True)
lams = np.array([d["lam"] for d in sp]); errs = np.array([abs(d["J_full"] - d["J_red"]) for d in sp])
pfit = np.polyfit(np.log(lams[-3:]), np.log(errs[-3:]), 1)
out["singular_perturbation"] = dict(table=sp, fitted_exponent_abs_error=pfit[0],
                                    fitted_exponent_rel_error=pfit[0] - 1)

print(json.dumps(out, indent=1, default=float))
json.dump(out, open("results_hdsp05.json", "w"), indent=1, default=float)
