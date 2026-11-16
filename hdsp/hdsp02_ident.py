"""HDSP-02: observability, structural/practical identifiability, input design, CIs, OOS error.

Measurement y_k = K(t_k) + eps_k,  K = [B + eta(M)(R-E)]/D,  D = 1 + A + sigma_N^2.
Sensitivities by central differences on a batched fixed-step RK4 (all perturbed
parameter vectors integrated together on the same grid -> consistent differences).
"""
import json, warnings
import numpy as np
from scipy.optimize import least_squares, minimize
from model import NOMINAL
warnings.filterwarnings("ignore")

th = NOMINAL
DT = 0.02
NAMES = ["alpha", "lam", "beta_p", "beta_m", "A", "sigN2", "a", "b1", "b2", "b3", "M0", "E0", "B0"]
TRUE = dict(alpha=th.alpha, lam=th.lam, beta_p=th.beta_p, beta_m=th.beta_m, A=0.3, sigN2=0.2,
            a=th.a, b1=th.b1, b2=th.b2, b3=th.b3, M0=2.0, E0=0.2, B0=th.b3)


def simulate(P, Rgrid, uBgrid=None):
    """P: dict name -> array (batch). Returns K on the grid, shape (n_t, batch)."""
    al, lam, bp, bm = P["alpha"], P["lam"], P["beta_p"], P["beta_m"]
    D = 1 + P["A"] + P["sigN2"]
    a, b1, b2, b3 = P["a"], P["b1"], P["b2"], P["b3"]
    M, E, B = P["M0"].copy(), P["E0"].copy(), P["B0"].copy()
    n = len(Rgrid)
    uB = np.zeros(n) if uBgrid is None else uBgrid
    out = np.empty((n, M.size))

    def f(M, E, B, R, ub):
        bs = np.where(R >= E, bp, bm)
        return R - lam * M, bs * (R - E), -a * (B - b1) * (B - b2) * (B - b3) + ub

    for k in range(n):
        R = Rgrid[k]
        out[k] = (B + (R - E) / (1 + al * M)) / D
        if k == n - 1:
            break
        Rm_ = 0.5 * (Rgrid[k] + Rgrid[k + 1]); R1 = Rgrid[k + 1]
        um = 0.5 * (uB[k] + uB[k + 1])
        k1 = f(M, E, B, R, uB[k])
        k2 = f(M + DT / 2 * k1[0], E + DT / 2 * k1[1], B + DT / 2 * k1[2], Rm_, um)
        k3 = f(M + DT / 2 * k2[0], E + DT / 2 * k2[1], B + DT / 2 * k2[2], Rm_, um)
        k4 = f(M + DT * k3[0], E + DT * k3[1], B + DT * k3[2], R1, uB[k + 1])
        M = M + DT / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        E = E + DT / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        B = B + DT / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
    return out


def sens(base, names, Rgrid, uBgrid=None, every=25, rel=1e-5):
    """d y(t_k)/d log(p) (or d/dp for params whose true value is 0)."""
    nb = 2 * len(names) + 1
    P = {k: np.full(nb, float(v)) for k, v in base.items()}
    steps = []
    for i, nm in enumerate(names):
        v = base[nm]
        h = rel * (abs(v) if v != 0 else 1.0)
        P[nm][2 * i + 1] += h; P[nm][2 * i + 2] -= h
        steps.append(h / (abs(v) if v != 0 else 1.0))
    Y = simulate(P, Rgrid, uBgrid)[::every]
    S = np.stack([(Y[:, 2 * i + 1] - Y[:, 2 * i + 2]) / (2 * steps[i]) for i in range(len(names))], 1)
    return Y[:, 0], S


def rank_report(S, tol=1e-7):
    sv = np.linalg.svd(S / np.sqrt(S.shape[0]), compute_uv=False)
    r = int(np.sum(sv > tol * sv[0]))
    _, _, Vt = np.linalg.svd(S)
    null = Vt[r:] if r < S.shape[1] else np.zeros((0, S.shape[1]))
    return r, sv, null


T_EXP = 200.0
t = np.arange(0, T_EXP + DT / 2, DT)


def R_square(Tp, q, amp, v=th.v_max):
    tau = np.mod(t, Tp); tr = amp / v
    up = np.clip(tau / tr, 0, 1) * amp
    down = np.clip(1 - (tau - q * Tp) / tr, 0, 1) * amp
    return np.where(tau < q * Tp, up, down)


out = {}
designs = {
    "constant_R=E0(at_equilibrium_E)": np.full_like(t, TRUE["E0"]),
    "single_step_0.2->1": np.where(t < 20, 0.2, 1.0),
    "square_Tp=40_q=0.5": R_square(40, 0.5, 1.0),
}

# --- 1. structural identifiability (generic rank of the sensitivity matrix) ---
rng = np.random.default_rng(3)
cases = {}
for dname, Rg in designs.items():
    for bname, B0, uB in [("B_at_equilibrium", TRUE["b3"], None),
                          ("B_excited_by_uB_pulses", TRUE["b3"], 0.04 * np.sign(np.sin(2 * np.pi * t / 50)))]:
        base = dict(TRUE); base["B0"] = B0
        if dname.startswith("constant"):
            base["E0"] = TRUE["E0"]; base["M0"] = TRUE["E0"] / TRUE["lam"]  # full equilibrium
        _, S = sens(base, NAMES, Rg, uB)
        r, sv, null = rank_report(S)
        nd = []
        for v in null:
            idx = np.argsort(-np.abs(v))[:4]
            nd.append({NAMES[i]: round(float(v[i]), 3) for i in idx if abs(v[i]) > 0.05})
        cases[f"{dname} | {bname}"] = dict(rank=r, n_params=len(NAMES), null_directions=nd,
                                            smallest_sv=[float(x) for x in sv[-4:]])
out["structural"] = cases

# --- 2. observability of the state alone (theta known) -----------------------
obs = {}
for dname, Rg in designs.items():
    base = dict(TRUE)
    if dname.startswith("constant"):
        base["M0"] = TRUE["E0"] / TRUE["lam"] + 1.0  # M off-equilibrium, E0 = R
    _, S = sens(base, ["M0", "E0", "B0"], Rg)
    r, sv, null = rank_report(S)
    obs[dname] = dict(rank=r, sv=[float(x) for x in sv])
out["state_observability"] = obs

# --- 3. Fisher-information input design (identifiable reparametrisation) ------
ID = ["alpha", "lam", "beta_p", "beta_m", "A", "M0", "E0", "B0"]   # D enters via A (sigN2 fixed)
SIG = 0.01      # measurement noise std
EVERY = 25      # sample every 0.5 time units

def fim(Rg):
    _, S = sens(TRUE, ID, Rg, every=EVERY)
    return S.T @ S / SIG ** 2

def score(p):
    Tp, q = np.exp(p[0]), 1 / (1 + np.exp(-p[1]))
    if q * Tp < 1 or (1 - q) * Tp < 1:
        return 1e9
    F = fim(R_square(Tp, q, 1.0))
    return -np.log(max(np.linalg.eigvalsh(F)[0], 1e-300))

des = {}
for dname, Rg in designs.items():
    if dname.startswith("constant"): continue
    ev = np.linalg.eigvalsh(fim(Rg))
    des[dname] = dict(lambda_min=float(ev[0]), energy=float(np.sum(Rg ** 2) * DT))
grid = []
for Tp in (10, 20, 30, 40, 60, 80, 120):
    for q in (0.2, 0.35, 0.5, 0.65, 0.8):
        grid.append((score([np.log(Tp), np.log(q / (1 - q))]), Tp, q))
grid.sort()
rr = minimize(score, [np.log(grid[0][1]), np.log(grid[0][2] / (1 - grid[0][2]))], method="Nelder-Mead",
              options=dict(xatol=1e-3, fatol=1e-4, maxiter=150))
Tp_s, q_s = np.exp(rr.x[0]), 1 / (1 + np.exp(-rr.x[1]))
Rs = R_square(Tp_s, q_s, 1.0)
F_s = fim(Rs)
ev = np.linalg.eigvalsh(F_s)
des["optimal_square"] = dict(Tp=Tp_s, q=q_s, lambda_min=float(ev[0]), energy=float(np.sum(Rs ** 2) * DT),
                             crb_rel_sd={n: float(np.sqrt(np.linalg.inv(F_s)[i, i])) for i, n in enumerate(ID)})
# PRBS-like random multilevel
prbs = np.repeat(rng.choice([0.0, 0.5, 1.0], size=len(t) // 500 + 1), 500)[:len(t)]
prbs = np.convolve(prbs, np.ones(50) / 50, mode="same")  # rate-limit smoothing
ev = np.linalg.eigvalsh(fim(prbs))
des["random_multilevel"] = dict(lambda_min=float(ev[0]), energy=float(np.sum(prbs ** 2) * DT))
# energy efficiency: lambda_min per unit energy for scaled amplitude (amplitude sweep)
eff = []
for amp in (0.25, 0.5, 0.75, 1.0):
    Rg = R_square(Tp_s, q_s, amp)
    eff.append(dict(amp=amp, lambda_min=float(np.linalg.eigvalsh(fim(Rg))[0]), energy=float(np.sum(Rg ** 2) * DT)))
des["amplitude_sweep_optimal_square"] = eff
out["input_design"] = des

# --- 4. minimax (local asymptotic) risk over the Theta box --------------------
risks = []
for _ in range(12):
    pert = {k: TRUE[k] * (1 + rng.uniform(-0.3, 0.3)) if k in ("alpha", "lam", "beta_p", "beta_m", "A") else TRUE[k]
            for k in TRUE}
    _, S = sens(pert, ID, Rs, every=EVERY)
    F = S.T @ S / SIG ** 2
    risks.append(float(np.trace(np.linalg.inv(F))))
out["local_minimax_risk_trace_invFIM_logparams"] = dict(sup_over_samples=max(risks), median=float(np.median(risks)))

# --- 5. estimation on synthetic data + CIs + out-of-sample error --------------
y_true = simulate({k: np.array([float(v)]) for k, v in TRUE.items()}, Rs)[::EVERY, 0]
noise = rng.normal(0, SIG, y_true.size)
y_obs = y_true + noise
n_tr = y_obs.size // 2

def resid(lp, sl):
    P = {k: np.array([float(v)]) for k, v in TRUE.items()}
    for i, nm in enumerate(ID):
        P[nm] = np.array([np.exp(lp[i])])
    y = simulate(P, Rs)[::EVERY, 0]
    return (y - y_obs)[sl] / SIG

lp0 = np.log([TRUE[n] for n in ID]) + rng.normal(0, 0.2, len(ID))
fit = least_squares(resid, lp0, args=(slice(0, n_tr),), x_scale=1.0, xtol=1e-10, ftol=1e-10, max_nfev=200)
Jt = fit.jac
cov = np.linalg.inv(Jt.T @ Jt)
est = {}
for i, nm in enumerate(ID):
    sd = np.sqrt(cov[i, i])
    lo, hi = np.exp(fit.x[i] - 1.96 * sd), np.exp(fit.x[i] + 1.96 * sd)
    est[nm] = dict(true=TRUE[nm], est=float(np.exp(fit.x[i])), ci95=[float(lo), float(hi)],
                   covered=bool(lo <= TRUE[nm] <= hi))
P = {k: np.array([float(v)]) for k, v in TRUE.items()}
for i, nm in enumerate(ID):
    P[nm] = np.array([np.exp(fit.x[i])])
y_pred = simulate(P, Rs)[::EVERY, 0]
te = slice(n_tr, None)
mse_model = float(np.mean((y_pred[te] - y_obs[te]) ** 2))
mse_null = float(np.mean((np.mean(y_obs[:n_tr]) - y_obs[te]) ** 2))
out["estimation"] = dict(params=est, test_mse_model=mse_model, test_mse_null_mean=mse_null,
                         noise_var=SIG ** 2, D_hat=float(1 + np.exp(fit.x[ID.index("A")]) + TRUE["sigN2"]))

print(json.dumps(out, indent=1, default=float))
json.dump(out, open("results_hdsp02.json", "w"), indent=1, default=float)
