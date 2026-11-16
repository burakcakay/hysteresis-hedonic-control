"""HDSP-06: robustness radius -- largest uncertainty box on (alpha, lam, beta_p, beta_m, D)
for which some trapezoid policy still has positive worst-case gain over R = 0."""
import itertools, json, warnings
import numpy as np
from scipy.optimize import minimize
from model import NOMINAL as th
from periodic_exact import G_exact
warnings.filterwarnings("ignore")
KEYS = ["alpha", "lam", "beta_p", "beta_m", "D"]


def verts(rho):
    for s in itertools.product((-1, 1), repeat=len(KEYS)):
        yield th.with_(**{k: getattr(th, k) * (1 + si * rho) for k, si in zip(KEYS, s)})


opt = json.load(open("results_hdsp04.json"))["optimum"]
res = {}
for rho in (0.02, 0.05, 0.08, 0.1, 0.12, 0.15, 0.2, 0.3):
    V = list(verts(rho))

    def f(p):
        T, q = np.exp(p[0]), 1 / (1 + np.exp(-p[1]))
        return -min(G_exact(T, q, th.R_max, v) for v in V)
    best = None
    for p0 in ([np.log(49), np.log(0.25 / 0.75)], [np.log(30), -1.5], [np.log(80), -0.8]):
        r = minimize(f, p0, method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-12, maxiter=800))
        if best is None or r.fun < best.fun:
            best = r
    w = min(V, key=lambda v: G_exact(opt["Tp"], opt["q"], opt["a"], v))
    res[str(rho)] = dict(maxmin_dG=float(-best.fun), Tp=float(np.exp(best.x[0])),
                         q=float(1 / (1 + np.exp(-best.x[1]))),
                         worst_vertex_for_nominal_policy={k: round(getattr(w, k), 4) for k in KEYS})
    print(rho, res[str(rho)], flush=True)
json.dump(res, open("results_hdsp06_radius.json", "w"), indent=1)
