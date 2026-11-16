"""HDSP-04 extension: periodic family with separate up/down ramp speeds (lower bound for HDSP-03)."""
import json, warnings
import numpy as np
from scipy.optimize import minimize
from model import NOMINAL as th
from periodic_exact import G_exact
warnings.filterwarnings("ignore")


def f(p):
    T, q = np.exp(p[0]), 1 / (1 + np.exp(-p[1]))
    v = th.v_max / (1 + np.exp(-p[2])); vd = th.v_max / (1 + np.exp(-p[3]))
    return -G_exact(T, q, th.R_max, th, v=v, vd=vd)


best = None
for p0 in ([np.log(46), -1.07, -0.2, -0.2], [np.log(46), -1.07, 2, -1], [np.log(46), -1.07, -1, 2], [np.log(60), -1.2, 0, 0]):
    r = minimize(f, p0, method="Nelder-Mead", options=dict(xatol=1e-8, fatol=1e-14, maxiter=8000))
    if best is None or r.fun < best.fun:
        best = r
p = best.x
out = dict(Tp=float(np.exp(p[0])), q=float(1 / (1 + np.exp(-p[1]))),
           v_up=float(th.v_max / (1 + np.exp(-p[2]))), v_down=float(th.v_max / (1 + np.exp(-p[3]))),
           dG=float(-best.fun))
print(json.dumps(out, indent=1))
json.dump(out, open("results_hdsp04_asym_ramps.json", "w"), indent=1)
