"""HDSP-04: optimal periodic (trapezoid) policy, KKT/Hessian, landscape, ramp-speed extension."""
import json
import numpy as np
from scipy.optimize import minimize
from model import NOMINAL, periodic_gain
from periodic_exact import G_exact

th = NOMINAL
out = {"G_constant": th.K0 / th.D, "note": "values below are dG = G - G_constant"}

# 1. coarse landscape (vectorised RK4) ---------------------------------------
Tp = np.geomspace(4, 400, 60)[:, None, None]
q = np.linspace(0.02, 0.98, 49)[None, :, None]
a = np.linspace(0.1, 1.0, 10)[None, None, :]
dG = periodic_gain(Tp, q, a, th, N=600, settle=2)
i = np.unravel_index(np.argmax(dG), dG.shape)
x0 = (Tp.ravel()[i[0]], q.ravel()[i[1]], a.ravel()[i[2]])
out["coarse_best"] = dict(Tp=x0[0], q=x0[1], a=x0[2], dG=float(dG[i]))
# count local maxima in (Tp,q) at a=R_max to characterise uniqueness
slab = dG[:, :, -1]
locmax = 0
for ii in range(1, slab.shape[0] - 1):
    for jj in range(1, slab.shape[1] - 1):
        w = slab[ii - 1:ii + 2, jj - 1:jj + 2]
        if np.isfinite(slab[ii, jj]) and slab[ii, jj] == np.max(w) and slab[ii, jj] > 0:
            locmax += 1
out["n_local_max_grid_a=Rmax"] = locmax
out["max_dG_vs_a"] = {f"{aa:.1f}": float(np.max(dG[:, :, k])) for k, aa in enumerate(a.ravel())}
out["frac_positive"] = float(np.mean(dG > 0))

# 2. refine (a free in (0, R_max]) ---------------------------------------------
def negG(p):
    T, qq, aa = np.exp(p[0]), 1 / (1 + np.exp(-p[1])), min(th.R_max, np.exp(p[2]))
    return -G_exact(T, qq, aa, th)

p0 = [np.log(x0[0]), np.log(x0[1] / (1 - x0[1])), np.log(x0[2])]
r = minimize(negG, p0, method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-13, maxiter=6000))
T_s, q_s, a_s = np.exp(r.x[0]), 1 / (1 + np.exp(-r.x[1])), min(th.R_max, np.exp(r.x[2]))
G_s = -r.fun
out["optimum"] = dict(Tp=T_s, q=q_s, a=a_s, dG=G_s, G=G_s + th.K0 / th.D)

# 3. first/second-order conditions (central differences on G_exact) ----------
def G3(T, qq, aa):
    return G_exact(T, qq, aa, th)
hT, hq, ha = 1e-3 * T_s, 1e-4, 1e-4
gT = (G3(T_s + hT, q_s, a_s) - G3(T_s - hT, q_s, a_s)) / (2 * hT)
gq = (G3(T_s, q_s + hq, a_s) - G3(T_s, q_s - hq, a_s)) / (2 * hq)
ga = (G3(T_s, q_s, a_s) - G3(T_s, q_s, a_s - ha)) / ha  # one-sided (a at bound?)
HT, Hq = 3e-2 * T_s, 3e-3
f0 = G3(T_s, q_s, a_s)
H = np.zeros((2, 2))
H[0, 0] = (G3(T_s + HT, q_s, a_s) - 2 * f0 + G3(T_s - HT, q_s, a_s)) / HT ** 2
H[1, 1] = (G3(T_s, q_s + Hq, a_s) - 2 * f0 + G3(T_s, q_s - Hq, a_s)) / Hq ** 2
H[0, 1] = H[1, 0] = (G3(T_s + HT, q_s + Hq, a_s) - G3(T_s + HT, q_s - Hq, a_s)
                     - G3(T_s - HT, q_s + Hq, a_s) + G3(T_s - HT, q_s - Hq, a_s)) / (4 * HT * Hq)
ev = np.linalg.eigvalsh(H)
out["first_order"] = dict(dG_dTp=gT, dG_dq=gq, dG_da_left=ga,
                          a_at_bound=bool(abs(a_s - th.R_max) < 1e-9))
out["hessian_Tp_q"] = dict(H=H.tolist(), eig=ev.tolist(), negative_definite=bool(np.all(ev < 0)))

# 4. multistart to probe uniqueness -------------------------------------------
rng = np.random.default_rng(1)
starts = []
for _ in range(16):
    p = [np.log(rng.uniform(8, 200)), rng.uniform(-3, 2), np.log(rng.uniform(0.3, 1.0))]
    rr = minimize(negG, p, method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-13, maxiter=3000))
    starts.append([float(np.exp(rr.x[0])), float(1 / (1 + np.exp(-rr.x[1]))),
                   float(min(th.R_max, np.exp(rr.x[2]))), float(-rr.fun)])
out["multistart"] = starts

# 5. extension: free ramp speed v <= v_max (beyond the (Tp,q,a) family) --------
def negG4(p):
    T, qq = np.exp(p[0]), 1 / (1 + np.exp(-p[1]))
    vv = th.v_max / (1 + np.exp(-p[2]))
    return -G_exact(T, qq, th.R_max, th, v=vv)
r4 = minimize(negG4, [r.x[0], r.x[1], 3.0], method="Nelder-Mead",
              options=dict(xatol=1e-7, fatol=1e-13, maxiter=6000))
out["free_ramp_speed"] = dict(Tp=float(np.exp(r4.x[0])), q=float(1 / (1 + np.exp(-r4.x[1]))),
                              v=float(th.v_max / (1 + np.exp(-r4.x[2]))), dG=float(-r4.fun))

print(json.dumps(out, indent=1, default=float))
json.dump(out, open("results_hdsp04.json", "w"), indent=1, default=float)
