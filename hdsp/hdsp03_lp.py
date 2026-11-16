"""HDSP-03: upper bound on the optimal long-run average g* via an HJB subsolution.

For any C^1 V and every admissible (z,u) on the compact invariant set
Omega = {0<=S<=Rmax, 0<=E<=Rmax, 0<=R<=Rmax}  (S = lam*M),
    l(z,u) + grad V(z) . F(z,u) <= g     ==>    g* <= g .
(integrate along any admissible trajectory, divide by T, V bounded on Omega).

sup_u is explicit:  sup_{|u|<=v} [p u - mu|u| - kappa u^2/2] = psi(p).
So g_upper(V) = max_z H(z),  H = l0(z) + V_S lam(R-S) + V_E beta_s(R-E) + psi(V_R),
with l0 = [K0 + eta(S)(R-E)]/D - cR.   We minimise max_z H over polynomial V by LP
(cutting planes on u and z), then evaluate max_z H on a fine grid with a Lipschitz margin.
"""
import json, sys, time
import numpy as np
from numpy.polynomial import chebyshev as Ch
from scipy.optimize import linprog
from model import NOMINAL

th = NOMINAL
Rm = th.R_max
aS = th.alpha / th.lam  # eta(S) = 1/(1 + aS*S)


def psi(p):
    ap = np.abs(p)
    v, mu, k = th.v_max, th.mu, th.kappa
    return np.where(ap <= mu, 0.0,
                    np.where(ap <= mu + k * v, (ap - mu) ** 2 / (2 * k), (ap - mu) * v - k * v * v / 2))


def ustar(p):
    ap = np.abs(p)
    v, mu, k = th.v_max, th.mu, th.kappa
    return np.sign(p) * np.clip((ap - mu) / k, 0, v) * (ap > mu)


def l0(S, E, R):
    return (th.K0 + (R - E) / (1 + aS * S)) / th.D - th.c * R


def bS(R, E):
    return np.where(R >= E, th.beta_p, th.beta_m)


class Basis:
    """Tensor Chebyshev basis on [0,Rm]^3, total degree <= d, constant removed."""

    def __init__(self, d):
        self.d = d
        self.idx = [(i, j, k) for i in range(d + 1) for j in range(d + 1) for k in range(d + 1)
                    if 0 < i + j + k <= d]

    def mats(self, x):
        X = 2 * x / Rm - 1
        n = self.d + 1
        T = np.stack([Ch.chebval(X, np.eye(n)[i]) for i in range(n)], 1)
        dT = np.stack([Ch.chebval(X, Ch.chebder(np.eye(n)[i])) * (2 / Rm) for i in range(n)], 1)
        return T, dT

    def grads_points(self, S, E, R):
        """Rows of dV/dS, dV/dE, dV/dR for scattered points (n_pts x n_basis)."""
        TS, dTS = self.mats(S); TE, dTE = self.mats(E); TR, dTR = self.mats(R)
        I = np.array(self.idx)
        GS = dTS[:, I[:, 0]] * TE[:, I[:, 1]] * TR[:, I[:, 2]]
        GE = TS[:, I[:, 0]] * dTE[:, I[:, 1]] * TR[:, I[:, 2]]
        GR = TS[:, I[:, 0]] * TE[:, I[:, 1]] * dTR[:, I[:, 2]]
        return GS, GE, GR

    def coef_tensor(self, w):
        C = np.zeros((self.d + 1,) * 3)
        for c, (i, j, k) in zip(w, self.idx):
            C[i, j, k] = c
        return C

    def grads_grid(self, w, x):
        """Gradient fields on the tensor grid x^3 (fast, via einsum)."""
        C = self.coef_tensor(w)
        T, dT = self.mats(x)
        VS = np.einsum('ijk,ai,bj,ck->abc', C, dT, T, T, optimize=True)
        VE = np.einsum('ijk,ai,bj,ck->abc', C, T, dT, T, optimize=True)
        VR = np.einsum('ijk,ai,bj,ck->abc', C, T, T, dT, optimize=True)
        return VS, VE, VR


def H_grid(basis, w, x):
    VS, VE, VR = basis.grads_grid(w, x)
    S, E, R = np.meshgrid(x, x, x, indexing='ij')
    return l0(S, E, R) + VS * th.lam * (R - S) + VE * bS(R, E) * (R - E) + psi(VR), (S, E, R, VR)


def solve(d, n_samp=22, rounds=40, wbound=2.0, verbose=True):
    B = Basis(d)
    nb = len(B.idx)
    # Chebyshev-Lobatto sample grid
    x = Rm * (1 - np.cos(np.pi * np.arange(n_samp) / (n_samp - 1))) / 2
    S, E, R = [a.ravel() for a in np.meshgrid(x, x, x, indexing='ij')]
    rows_A, rows_b = [], []

    def add(S_, E_, R_, U_):
        GS, GE, GR = B.grads_points(S_, E_, R_)
        A = GS * (th.lam * (R_ - S_))[:, None] + GE * (bS(R_, E_) * (R_ - E_))[:, None] + GR * U_[:, None]
        rhs = -(l0(S_, E_, R_) - th.mu * np.abs(U_) - th.kappa * U_ ** 2 / 2)
        rows_A.append(np.hstack([A, -np.ones((len(S_), 1))]))
        rows_b.append(rhs)

    for u in (-th.v_max, 0.0, th.v_max):
        add(S, E, R, np.full_like(S, u))
    xf = np.linspace(0, Rm, 61)
    hist = []
    w = None
    for it in range(rounds):
        A = np.vstack(rows_A); b = np.concatenate(rows_b)
        cobj = np.zeros(nb + 1); cobj[-1] = 1.0
        bounds = [(-wbound, wbound)] * nb + [(None, None)]
        res = linprog(cobj, A_ub=A, b_ub=b, bounds=bounds, method='highs')
        if res.status != 0:
            print('LP failed', res.message); break
        w, g_lp = res.x[:-1], res.x[-1]
        # exact sup over u on scattered samples + fine grid ; add cuts
        Hf, (Sf, Ef, Rf, VRf) = H_grid(B, w, xf)
        viol = Hf - g_lp
        g_fine = Hf.max()
        hist.append((it, g_lp, g_fine, A.shape[0]))
        if verbose:
            print(f'  d={d} it={it} rows={A.shape[0]} g_lp={g_lp - th.K0/th.D:+.6f} '
                  f'g_fine={g_fine - th.K0/th.D:+.6f}', flush=True)
        if g_fine - g_lp < 2e-6:
            break
        # cuts: worst points on fine grid (with their optimal u), plus optimal-u cuts on sample grid
        flat = np.argsort(viol.ravel())[::-1][:1500]
        sel = flat[viol.ravel()[flat] > 1e-7]
        if len(sel):
            Uc = ustar(VRf.ravel()[sel])
            add(Sf.ravel()[sel], Ef.ravel()[sel], Rf.ravel()[sel], Uc)
        GS, GE, GR = B.grads_points(S, E, R)
        Us = ustar(GR @ w)
        keep = np.abs(Us) > 0
        add(S[keep], E[keep], R[keep], Us[keep])
        # random points
        rng = np.random.default_rng(it)
        Sr, Er, Rr = rng.uniform(0, Rm, (3, 3000))
        GS, GE, GR = B.grads_points(Sr, Er, Rr)
        add(Sr, Er, Rr, ustar(GR @ w))
    return B, w, g_lp, hist


def certify(B, w, n=161):
    """max of H on a fine grid + Lipschitz margin (numerical certificate)."""
    x = np.linspace(0, Rm, n)
    h = x[1] - x[0]
    Hf, _ = H_grid(B, w, x)
    m = Hf.max()
    # Lipschitz constants from finite differences on the grid (x2 safety)
    L = [np.abs(np.diff(Hf, axis=ax)).max() / h * 2 for ax in range(3)]
    margin = sum(Lk * h / 2 for Lk in L)
    return m, margin, L


if __name__ == '__main__':
    degs = [int(a) for a in sys.argv[1:]] or [1, 2, 4, 6, 8, 10]
    out = {"G_constant": th.K0 / th.D, "runs": []}
    for d in degs:
        t0 = time.time()
        B, w, g_lp, hist = solve(d)
        m, margin, L = certify(B, w)
        rec = dict(d=d, n_basis=len(B.idx), dG_lp=g_lp - th.K0 / th.D,
                   dG_upper_grid161=m - th.K0 / th.D,
                   heuristic_lipschitz_margin=margin, seconds=time.time() - t0)
        print(json.dumps(rec), flush=True)
        out["runs"].append(rec)
        np.save(f'V_coef_d{d}.npy', w)
    json.dump(out, open(f'results_hdsp03_d{"_".join(map(str, degs))}.json', 'w'), indent=1, default=float)
