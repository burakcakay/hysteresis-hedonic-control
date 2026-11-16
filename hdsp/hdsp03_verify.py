"""HDSP-03: re-check saved LP certificates V on several grids and random points."""
import json, sys, warnings
import numpy as np
from hdsp03_lp import Basis, H_grid, l0, bS, psi, th
warnings.filterwarnings("ignore")
out = {}
for d in [int(a) for a in sys.argv[1:]] or [2, 4, 6]:
    w = np.load(f"V_coef_d{d}.npy"); B = Basis(d)
    rec = {}
    for n in (161, 201, 241):
        H, _ = H_grid(B, w, np.linspace(0, th.R_max, n))
        rec[f"grid{n}"] = float(H.max() - th.K0 / th.D)
        del H
    rng = np.random.default_rng(0); best = -np.inf
    for _ in range(20):
        S, E, R = rng.uniform(0, th.R_max, (3, 400000))
        GS, GE, GR = B.grads_points(S, E, R)
        H = l0(S, E, R) + (GS @ w) * th.lam * (R - S) + (GE @ w) * bS(R, E) * (R - E) + psi(GR @ w)
        best = max(best, H.max())
    rec["random_8e6"] = float(best - th.K0 / th.D)
    rec["upper_bound_reported"] = max(rec.values())
    out[d] = rec
    print(d, rec, flush=True)
json.dump(out, open("results_hdsp03_verify.json", "w"), indent=1)
