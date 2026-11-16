"""HDSP-01: minimal dimension, equilibria, stability certificates, basins, P1-P4."""
import json
import numpy as np
from scipy.integrate import solve_ivp
from model import NOMINAL, eta, fB, K_out

th = NOMINAL
out = {}

# --- 1. equilibria under constant R = Rbar ---------------------------------
# M* = Rbar/lam, E* = Rbar (=> R-E = 0, eta-term vanishes), B* in roots of fB + uB
def equilibria(Rbar, uB=0.0):
    coeffs = np.array([-th.a, th.a * (th.b1 + th.b2 + th.b3),
                       -th.a * (th.b1 * th.b2 + th.b1 * th.b3 + th.b2 * th.b3),
                       th.a * th.b1 * th.b2 * th.b3 + uB])
    r = np.roots(coeffs)
    Bs = np.sort(r[np.abs(r.imag) < 1e-9].real)
    res = []
    for B in Bs:
        dfB = -th.a * ((B - th.b2) * (B - th.b3) + (B - th.b1) * (B - th.b3) + (B - th.b1) * (B - th.b2))
        # one-sided Jacobians (E-equation is piecewise linear with a kink at R=E)
        Jp = np.diag([-th.lam, -th.beta_p, dfB])
        Jm = np.diag([-th.lam, -th.beta_m, dfB])
        res.append(dict(M=Rbar / th.lam, E=Rbar, B=B, K=B / th.D,
                        spec_plus=np.linalg.eigvals(Jp).real.tolist(),
                        spec_minus=np.linalg.eigvals(Jm).real.tolist(),
                        stable=bool(dfB < 0)))
    return res

Rbar = 0.5
eq = equilibria(Rbar)
out["equilibria_Rbar=0.5"] = eq

# bistability window in u_B for symmetric wells: |u_B| < 2 a w^3 /(3 sqrt 3), w=(b3-b1)/2
w = (th.b3 - th.b1) / 2
uB_crit = 2 * th.a * w ** 3 / (3 * np.sqrt(3))
nroots = []
for uB in np.linspace(-0.1, 0.1, 401):
    nroots.append((uB, len(equilibria(Rbar, uB))))
emp = [u for u, n in nroots if n == 3]
out["uB_crit_formula"] = uB_crit
out["uB_window_numeric"] = [min(emp), max(emp)]

# --- 2. basins: sign(B0 - b2) decides; check by simulation -------------------
rng = np.random.default_rng(0)
def rhs(t, x, Rfun, uBfun=lambda t: 0.0):
    M, E, B = x
    R = Rfun(t)
    bs = th.beta_p if R >= E else th.beta_m
    return [R - th.lam * M, bs * (R - E), fB(B, th, uBfun(t))]

mismatch = 0
N = 400
for _ in range(N):
    x0 = [rng.uniform(0, 20), rng.uniform(0, 1), rng.uniform(-1, 2)]
    sol = solve_ivp(rhs, [0, 300], x0, args=(lambda t: Rbar,), rtol=1e-9, atol=1e-11)
    Bend = sol.y[2, -1]
    predicted = th.b3 if x0[2] > th.b2 else th.b1
    if abs(Bend - predicted) > 1e-4 or abs(sol.y[0, -1] - Rbar / th.lam) > 1e-3:
        mismatch += 1
out["basin_check"] = dict(samples=N, mismatches=mismatch,
                          rule="Basin(x1*)={B<b2}, Basin(x2*)={B>b2}, separatrix {B=b2}")

# Lyapunov certificate on each basin (decoupled): V = (M-M*)^2/lam + (E-E*)^2 + U(B)-U(b_i)
# dV/dt = -2(M-M*)^2 - 2 beta_s (E-E*)^2 - fB(B)^2 <= 0, =0 only at equilibrium points.
Bg = np.linspace(-1, 2, 3001)
U = -np.cumsum(fB(Bg, th)) * (Bg[1] - Bg[0])
out["potential_barriers"] = dict(dU_low=float(np.interp(th.b2, Bg, U) - np.interp(th.b1, Bg, U)),
                                 dU_high=float(np.interp(th.b2, Bg, U) - np.interp(th.b3, Bg, U)))

# --- 3. P4 realisations ------------------------------------------------------
# (a) steady-state paradox: x1* on upper branch at R1=0.2, x2* on lower branch at R2=0.8
out["P4_steady"] = dict(R1=0.2, K1=th.b3 / th.D, R2=0.8, K2=th.b1 / th.D)

# (b) single trajectory from equilibrium, B frozen on upper well, u_B = 0:
#     contrast effect through E.  R: 0.5 -> step 0.9 at t=50 -> ramp to 0.95 by t=150
def Rfun(t):
    if t < 50: return 0.5
    return min(0.9 + 0.05 * (t - 50) / 100, 0.95)
x0 = [0.5 / th.lam, 0.5, th.b3]
sol = solve_ivp(rhs, [0, 300], x0, args=(Rfun,), rtol=1e-10, atol=1e-12, max_step=0.05, dense_output=True)
t1, t2 = 50.5, 300.0
x1, x2 = sol.sol(t1), sol.sol(t2)
K1 = K_out(x1[0], x1[1], x1[2], Rfun(t1), th)
K2 = K_out(x2[0], x2[1], x2[2], Rfun(t2), th)
out["P4_transient_from_equilibrium"] = dict(t1=t1, R1=Rfun(t1), K1=float(K1), t2=t2, R2=Rfun(t2), K2=float(K2),
                                            holds=bool(Rfun(t1) < Rfun(t2) and K1 > K2))

# (c) the same with the E coordinate removed (E == R, instantaneous adaptation): K = B/D constant
out["P4_without_E_from_equilibrium"] = "K(t) = b3/D for all t -> P4 fails (needs E, or u_B coupled to R)"

print(json.dumps(out, indent=1, default=float))
json.dump(out, open("results_hdsp01.json", "w"), indent=1, default=float)
