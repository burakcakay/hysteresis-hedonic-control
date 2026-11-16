"""HDSP-06: stochastic / robust layer.

dx = f(x,R;theta)dt + G dW,  G = diag(sM, sE, sB) constant.
Checks: (1) mean-square bound (explicit Lyapunov), (2) chance constraint P(K<K_min),
(3) robust performance inf_theta J(pi,theta) vs baseline, (4) F-floor (F := R assumed).
"""
import json, itertools, sys, warnings
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize, brentq
from model import NOMINAL, trapezoid
from periodic_exact import G_exact
warnings.filterwarnings("ignore")

th = NOMINAL
SIG = dict(sM=0.05, sE=0.05, sB=0.12)
K_MIN, EPS = 0.5, 0.05
opt04 = json.load(open("results_hdsp04.json"))["optimum"]
Tp_s, q_s, a_s = opt04["Tp"], opt04["q"], opt04["a"]
out = {"policy_from_HDSP04": dict(Tp=Tp_s, q=q_s, a=a_s), "noise": SIG, "K_min": K_MIN, "eps": EPS}

# --- Theta box: +-30% on dynamics params ------------------------------------
BOX = dict(alpha=0.3, lam=0.3, beta_p=0.3, beta_m=0.3, D=0.2)
def vertices():
    keys = list(BOX)
    for signs in itertools.product((-1, 1), repeat=len(keys)):
        yield th.with_(**{k: getattr(th, k) * (1 + s * BOX[k]) for k, s in zip(keys, signs)})

# --- 1. mean-square stability bound ------------------------------------------
# V = M^2 + E^2 + B^2 :  LV <= -k0 V + C   (derivation in report)
def ms_bound(t_):
    Bg = np.linspace(-5, 5, 200001)
    fBv = -t_.a * (Bg - t_.b1) * (Bg - t_.b2) * (Bg - t_.b3)
    CB = np.max(2 * Bg * fBv + Bg ** 2)          # 2 B fB <= -B^2 + CB
    bmin, bmax = min(t_.beta_p, t_.beta_m), max(t_.beta_p, t_.beta_m)
    k0 = min(t_.lam, bmin, 1.0)
    C = t_.R_max ** 2 / t_.lam + bmax * t_.R_max ** 2 + CB + SIG["sM"] ** 2 + SIG["sE"] ** 2 + SIG["sB"] ** 2
    return k0, C, C / k0
bounds = [ms_bound(v) for v in vertices()]
out["mean_square"] = dict(nominal=dict(zip(("k0", "C", "C_over_k0"), ms_bound(th))),
                          worst_C_over_k0_on_box=max(b[2] for b in bounds))

# --- 2. chance constraint ------------------------------------------------------
# (a) B-well occupancy: 1-D stationary density p(B) ~ exp(-2U(B)/sB^2), exact
def p_low(uB, sB=SIG["sB"], t_=th):
    U = lambda B: t_.a * ((B ** 4) / 4 - (t_.b1 + t_.b2 + t_.b3) * B ** 3 / 3
                          + (t_.b1 * t_.b2 + t_.b1 * t_.b3 + t_.b2 * t_.b3) * B ** 2 / 2
                          - t_.b1 * t_.b2 * t_.b3 * B) - uB * B
    w = lambda B: np.exp(-2 * U(B) / sB ** 2)
    lo = quad(w, -3, t_.b2, limit=200)[0]; hi = quad(w, t_.b2, 4, limit=200)[0]
    return lo / (lo + hi)
out["P_lower_well_uB=0"] = p_low(0.0)
u_star = brentq(lambda u: p_low(u) - EPS / 2, 0.0, 0.2)
out["uB_tilt_needed_for_P_low<=eps/2"] = u_star
# (b) policy-induced contrast dips on the upper well (deterministic cycle)
def frac_below(Tp, q, a, t_=th, n_per=6, N=20000):
    tt = np.linspace(0, n_per * Tp, n_per * N)
    dt = tt[1] - tt[0]
    R = trapezoid(tt, Tp, q, a, t_.v_max)
    M = np.zeros_like(tt); E = np.zeros_like(tt)
    for k in range(len(tt) - 1):
        bs = t_.beta_p if R[k] >= E[k] else t_.beta_m
        M[k + 1] = M[k] + dt * (R[k] - t_.lam * M[k])
        E[k + 1] = E[k] + dt * bs * (R[k] - E[k])
    K = (t_.b3 + (R - E) / (1 + t_.alpha * M)) / t_.D
    last = tt > (n_per - 1) * Tp
    return float(np.mean(K[last] < K_MIN)), float(K[last].min())
fb, kmin = frac_below(Tp_s, q_s, a_s)
out["deterministic_cycle_upper_well"] = dict(frac_time_K_below_Kmin=fb, min_K=kmin)

# (c) Monte Carlo of the full SDE under the policy (Euler-Maruyama), uB = u_star tilt
def mc(Tp, q, a, uB, n_paths=400, T=2000.0, dt=0.05, seed=0, t_=th):
    rng = np.random.default_rng(seed)
    n = int(T / dt)
    M = np.full(n_paths, a * q / t_.lam); E = np.zeros(n_paths); B = np.full(n_paths, t_.b3)
    below = 0; cnt = 0; J = np.zeros(n_paths)
    sq = np.sqrt(dt)
    burn = int(0.2 * n)
    for k in range(n):
        tk = k * dt
        R = trapezoid(tk, Tp, q, a, t_.v_max) if a > 0 else 0.0
        bs = np.where(R >= E, t_.beta_p, t_.beta_m)
        K = (B + (R - E) / (1 + t_.alpha * np.maximum(M, 0))) / t_.D
        if k >= burn:
            below += np.sum(K < K_MIN); cnt += n_paths; J += K * dt
        M = M + dt * (R - t_.lam * M) + SIG["sM"] * sq * rng.standard_normal(n_paths)
        E = E + dt * bs * (R - E) + SIG["sE"] * sq * rng.standard_normal(n_paths)
        B = B + dt * (-t_.a * (B - t_.b1) * (B - t_.b2) * (B - t_.b3) + uB) + SIG["sB"] * sq * rng.standard_normal(n_paths)
    Tm = (n - burn) * dt
    meanK = J / Tm
    return below / cnt, float(np.mean(meanK)), float(np.std(meanK) / np.sqrt(n_paths))
for name, (Tp, q, a) in {"baseline_R=0": (Tp_s, q_s, 0.0), "pi_star": (Tp_s, q_s, a_s)}.items():
    for uB in (0.0, u_star):
        pb, mK, se = mc(Tp, q, a, uB)
        out[f"MC_{name}_uB={uB:.4f}"] = dict(P_K_below=pb, mean_K=mK, se=se)
        print(name, uB, pb, mK, se, flush=True)

# --- 3. robust performance over the Theta box ---------------------------------
def robust_value(Tp, q, a):
    vals = []
    for v in vertices():
        vals.append(v.K0 / v.D + G_exact(Tp, q, a, v))
    return min(vals), vals
Jr_star, vals = robust_value(Tp_s, q_s, a_s)
Jr_base = min(v.K0 / v.D for v in vertices())
dg_vertices = [G_exact(Tp_s, q_s, a_s, v) for v in vertices()]
out["robust_nominal_policy"] = dict(J_robust=Jr_star, J_robust_baseline=Jr_base,
                                    improvement=Jr_star - Jr_base,
                                    min_dG_over_vertices=min(dg_vertices),
                                    frac_vertices_dG_positive=float(np.mean(np.array(dg_vertices) > 0)))
# max-min design: maximise min_vertex dG (D only rescales; policy-relevant part)
def neg_minvert(p):
    Tp, q = np.exp(p[0]), 1 / (1 + np.exp(-p[1]))
    return -min(G_exact(Tp, q, th.R_max, v) for v in vertices())
rr = minimize(neg_minvert, [np.log(Tp_s), np.log(q_s / (1 - q_s))], method="Nelder-Mead",
              options=dict(xatol=1e-4, fatol=1e-9, maxiter=200))
Tr, qr = np.exp(rr.x[0]), 1 / (1 + np.exp(-rr.x[1]))
out["maxmin_policy"] = dict(Tp=Tr, q=qr, a=th.R_max, worst_dG=-rr.fun,
                            nominal_dG=G_exact(Tr, qr, th.R_max, th),
                            J_robust=min(v.K0 / v.D + G_exact(Tr, qr, th.R_max, v) for v in vertices()))

# --- 4. F-floor, assuming F := R (Sigma does not define F) ---------------------
ff = {}
for Fmin in (0.2, 0.3, 0.5):
    def negc(p):
        Tp, q = np.exp(p[0]), 1 / (1 + np.exp(-p[1]))
        q = max(q, Fmin / th.R_max)  # <R> = a q >= Fmin with a = R_max
        return -G_exact(Tp, q, th.R_max, th)
    rc = minimize(negc, [np.log(Tp_s), 0.0], method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-10, maxiter=300))
    Tp, q = np.exp(rc.x[0]), max(1 / (1 + np.exp(-rc.x[1])), Fmin / th.R_max)
    # best constant policy with <R> >= Fmin is R = Fmin
    ff[str(Fmin)] = dict(Tp=Tp, q=q, dG_periodic=-rc.fun, dG_constant=-th.c * Fmin)
out["F_floor_assuming_F=R"] = ff

print(json.dumps(out, indent=1, default=float))
json.dump(out, open("results_hdsp06.json", "w"), indent=1, default=float)
