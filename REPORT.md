# Hysteresis in hedonic control — HDSP-01 … HDSP-06 and HDSP-Ω

**Author:** burakcakay · **Computational assistant:** Claude (Anthropic) · **Version:** 0.1.0 (2026-11-16)
Text licensed CC BY 4.0, code MIT. Every number here is reproduced by the scripts in `hdsp/` (see [Reproducibility](#reproducibility)).

Each claim carries one of these tags:

| Tag | Meaning |
|---|---|
| **[Proof]** | a complete mathematical argument is given in this document |
| **[Numerical]** | a reproducible computation; evidence, not a proof |
| **[Ill-posed]** | the problem as stated is degenerate; a corrected formulation is given and solved |
| **[Open]** | not settled here |

---

## 0. Model

State $x=(M,E,B)$, input $R(t)\in[0,R_{\max}]$, control $u=\dot R$, $|u|\le v_{\max}$:

$$
\dot M = R-\lambda M,\qquad
\dot E = \beta_s (R-E),\ \ \beta_s=\begin{cases}\beta_+ & R\ge E\\ \beta_- & R<E\end{cases},\qquad
\dot B = -a(B-b_1)(B-b_2)(B-b_3)+u_B ,
$$

$$
\eta(M)=\frac{1}{1+\alpha M},\qquad K=\frac{B+\eta(M)(R-E)}{D},\qquad D=1+A+\sigma_N^2 .
$$

Running reward (HDSP-03/04): $\ell=\big[K_0+\eta(M)(R-E)\big]/D-cR-\mu|u|-\tfrac{\kappa}{2}u^2$.
It is convenient to use $S=\lambda M\in[0,R_{\max}]$ (so $\dot S=\lambda(R-S)$, $\eta=1/(1+aS)$, $a=\alpha/\lambda$).

**Nominal parameters** (all positive, as Θ requires; chosen so that the slow tolerance has $a=\alpha/\lambda=3$):

| α | λ | β₊ | β₋ | D | K₀ | c | μ | κ | v_max | R_max | a | b₁,b₂,b₃ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.3 | 0.1 | 1.0 | 0.8 | 1.5 | 1 | 0.005 | 0.01 | 0.05 | 1 | 1 | 1 | 0, 0.5, 1 |

β₋ < β₊ means "losses linger": expectations come down more slowly than they go up.

**Undefined symbol.** $F$ (HDSP-06, HDSP-Ω) is not defined in Σ. Where it is needed we take $F:=R$ (a minimum-engagement floor) and say so.

---

## Summary

| Problem | Main result | Status |
|---|---|---|
| HDSP-01 | At every equilibrium $E^*=\bar R$, so the hedonic term vanishes and $K^*=B^*/D$. Hence $n_{\min}=1$, and the one necessary coordinate is $B$. If P4 must be produced from rest by $R$ alone, $n=2$, namely $(E,B)$. Stability holds even though $Df$ does not exist at the equilibrium (kink at $R=E$). Basins are $\{B\lessgtr b_2\}$. | **[Proof]** |
| HDSP-02 | $M$ is unobservable on $\{R=E\}$ for constant inputs. θ is locally identifiable except for $A+\sigma_N^2$ (only $D$ is). $\beta_\pm$ require excursions to both sides of $E$. $(a,b_1,b_2)$ require $B$ to leave its well. $c,\mu,\kappa$ are never identifiable from $(R,y)$. "min energy s.t. rank O = n" has infimum 0 (ill-posed); replaced by a Fisher-information design. | **[Proof]** + **[Numerical]** + **[Ill-posed]** |
| HDSP-03 | The best constant policy is $R=0$, giving $G_{\rm const}=K_0/D$ **[Proof]**. Verification lemma **[Proof]**. Bounds on the optimal gain: $0.002306\le g^*-K_0/D\le 0.00294$. The gap is not closed. | **[Numerical]** / **[Open]** |
| HDSP-04 | Unique optimum in the $(T_p,q,a)$ family: $T_p^*=49.03,\ q^*=0.2492,\ a^*=R_{\max}$ (boundary, KKT). The Hessian is negative definite. $G^*>G_{\rm const}$. The family is **not** globally optimal: freeing the ramp speeds improves the gain by 20%. | **[Numerical]** |
| HDSP-05 | Reduced model derived. The stated condition $S_L g(S_L)=rJ$ is correct **only if β₊ = β₋**; the corrected condition is given. The Hessian is diagonal at critical points. A single hysteresis loop is optimal among all on/off policies. $\lvert J_{\rm full}^*-J_{\rm red}^*\rvert\approx 0.51\,\lambda^2$, i.e. $p=1$ relative. | **[Proof]** + **[Numerical]** |
| HDSP-06 | Mean-square bound **[Proof]**. With symmetric wells, $P(K<K_{\min})\to\tfrac12$ for **every** policy, so the chance constraint is structurally infeasible without a tilt $u_B$ **[Proof]**. The nominal-optimal cycling policy is **not** robust: its robustness radius is about 10%. | **[Proof]** + **[Numerical]** |
| HDSP-Ω | Minimal dimension depends on the target: 1 for the paradox, 2 for the input-driven paradox, 3 for a non-trivial optimal policy. | **[Proof]** (dimension) / **[Open]** (ε→0) |

---

## HDSP-01 — bistability, stability certificates, minimal dimension

**Theorem 1.1 (equilibria) [Proof].** Under constant $R\equiv\bar R$ the equilibria are $M^*=\bar R/\lambda$, $E^*=\bar R$, and $B^*\in Z:=\{B: f_B(B)+u_B=0\}$. Consequently
$$K^*=\frac{B^*}{D}\quad\text{independently of }\bar R,\ \alpha,\ \lambda,\ \beta_\pm .$$
*Proof.* $\dot M=0\Rightarrow M=\bar R/\lambda$. $\dot E=\beta_s(\bar R-E)=0$ with $\beta_s>0$ gives $E=\bar R$, so $R-E=0$ and the $\eta$-term vanishes. ∎

So at steady state neither tolerance ($M$) nor expectation ($E$) changes $K$. Steady-state differences in $K$ can only come from $B$.

**Theorem 1.2 (stability without a Jacobian) [Proof].** The vector field has a kink on $\{R=E\}$, which is exactly where every equilibrium sits, so $Df(x^*)$ does not exist and P2/P3 as written are ill-defined. Two replacements both hold:

1. The Clarke generalized Jacobian is $\partial f(x^*)=\mathrm{conv}\{J_+,J_-\}$ with $J_\pm=\mathrm{diag}(-\lambda,-\beta_\pm,f_B'(B^*))$. Every element is Hurwitz iff $f_B'(B^*)<0$.
2. Lyapunov certificate on each basin:
$$V=\tfrac12(M-M^*)^2+\tfrac12(E-E^*)^2+U(B)-U(B^*),\qquad U'=-(f_B+u_B),$$
$$\dot V=-\lambda(M-M^*)^2-\beta_s(E-E^*)^2-(f_B+u_B)^2\le 0 .$$
Equality holds only at equilibria, so LaSalle gives convergence. The local exponential rate is $\min(\lambda,\beta_-,\beta_+,|f_B'(B^*)|)$.

With $u_B=0$: $f_B'(b_1)=-a(b_1-b_2)(b_1-b_3)<0$, $f_B'(b_3)<0$, and $f_B'(b_2)>0$ (saddle). Nominally the spectra are $\{-0.1,-1.0\text{ or }-0.8,-0.5\}$ at $b_1,b_3$ and $+0.25$ at $b_2$.

**Theorem 1.3 (basins) [Proof].** For $u_B=0$: $\mathrm{Basin}(x_1^*)=\{B<b_2\}$ and $\mathrm{Basin}(x_2^*)=\{B>b_2\}$ (any $M,E$). They are disjoint and separated by $\{B=b_2\}$, the stable manifold of the saddle. *Proof:* $B$ is decoupled and scalar, so it moves monotonically to the root on its side of $b_2$. ∎
Numerical check: 400 random initial conditions, 0 mismatches.
Bistability persists for $|u_B|<u_c$. For symmetric wells $u_c=\frac{2a w^3}{3\sqrt3}$ with $w=(b_3-b_1)/2$; nominally $u_c=0.04811$, which matches the numerical window $\pm0.048$.

**Theorem 1.4 (minimal dimension) [Proof].**
- *Lower bound.* If $\dim x=0$, then $K=H(R)$ and $R_1=R_2\Rightarrow K_1=K_2$, so P1 fails. Hence $n_{\min}\ge1$.
- *The $(M,E)$ part alone can never satisfy P1.* With $B$ frozen, each $\bar R$ has a unique globally stable equilibrium (Theorem 1.2), and in any case $K^*=B/D$ is the same for all $\bar R$.
- *Upper bound.* $x=(B)$ with cubic $f_B$ and $K=B/D$ satisfies P1–P3. P4 holds as well, e.g. upper well at $R_1=0.2$ gives $K=0.667$, lower well at $R_2=0.8$ gives $K=0$.

Therefore **$n_{\min}=1$, and the coordinate must be $B$.**

*Refinement.* Suppose P4 is required along a single trajectory that starts at an equilibrium and is driven by $R$ alone ($u_B$ independent of $R$). With $x=(B)$, $K$ is then constant, so P4 fails. With $x=(E,B)$ the contrast effect works. Example run: $R$ steps 0.5 → 0.9 at $t=50$, then creeps to 0.95. This gives $K(50.5)=0.730$ at $R=0.900$ and $K(300)=0.667$ at $R=0.95$. In that setting $n_{\min}=2$. $M$ is never needed for HDSP-01; it becomes necessary in HDSP-03/05.

---

## HDSP-02 — observability and identifiability

Measurement: $y_k=K(t_k)+\varepsilon_k$. Unknowns: $\theta=(\alpha,\lambda,\beta_+,\beta_-,A,\sigma_N^2,a,b_1,b_2,b_3)$ and $x_0=(M_0,E_0,B_0)$, 13 in total.

**Proposition 2.1 (unobservable set) [Proof].** Let $R$ be constant and write $\varepsilon=R-E$, so $\dot\varepsilon=-\beta_s\varepsilon$. Every Lie derivative has the form $L_f^k h=L_f^k(B)/D+\varepsilon\,P_k(M)$. Hence $\partial_M L_f^k h\propto\varepsilon$, and **the $M$-column of $\mathcal O(x)$ vanishes on $\{R=E\}$**, giving rank ≤ 2. This is exact, not only to first order: if $E_0=R$ then $y\equiv B(t)/D$ for every $M_0$. For a time-varying input $\dot\varepsilon=\dot R-\beta_s\varepsilon$, so any non-constant $R$ restores $\varepsilon\neq0$.
Numerical check (rank of $\partial y/\partial x_0$): constant input at $E=R$ gives rank 2; a step or a square wave gives rank 3.

**Structural identifiability [Numerical — generic rank of the 13-column output-sensitivity matrix].**

| Input | $B$ at its equilibrium | $B$ excited by $u_B$ pulses |
|---|---|---|
| constant $R=E_0$ | rank 3 | rank 7 |
| single step 0.2→1 | rank 7: β₋, a, b₁, b₂, $(A,\sigma_N^2)$, $(M_0,\alpha)$ not separated | rank 11: β₋, $(A,\sigma_N^2)$ |
| square wave $T_p=40$ | rank 9: $(a,b_1,b_2)$, $(A,\sigma_N^2)$ | **rank 12**: only $(A,\sigma_N^2)$ |

**Equivalence classes:**
1. $\{A+\sigma_N^2=D-1\}$ is always unidentifiable; only $D$ is. This one is a proof, since $A$ and $\sigma_N^2$ enter only through $D$.
2. β₋ (resp. β₊) is identifiable only if $R<E$ (resp. $R>E$) occurs.
3. $(a,b_1,b_2)$ are identifiable only if $B$ leaves its well.
4. $(\alpha,M_0)$ need contrast events at two or more distinct tolerance levels.
5. $c,\mu,\kappa$ appear in neither $f$ nor $H$, so they are not identifiable from $(R,y)$ at all. Only inverse optimal control on observed choices could recover them.

**[Ill-posed] "min ∫R² subject to rank O = n".** The rank condition is open and invariant under $R\mapsto\epsilon R$. Worse, for an off-equilibrium $x_0$, $R\equiv0$ already gives full state rank. So the infimum is 0 and is not attained.
**Fix:** maximise $\lambda_{\min}(\mathcal I)$ of the Fisher information in log-parameters. Setting: identifiable set $(\alpha,\lambda,\beta_+,\beta_-,A,M_0,E_0,B_0)$, noise $\sigma=0.01$, horizon 200, sampling every 0.5.

| Design | $\lambda_{\min}(\mathcal I)$ | energy $\int R^2$ |
|---|---|---|
| single step | ≈ 0 (β₋ not excited) | 180.8 |
| square $T_p=40,q=0.5$ | 7.03 | 98.3 |
| random multilevel | 9.39 | 59.0 |
| **optimal trapezoid $T_p=10.8,\ q=0.20$** | **12.00** | 34.9 |
| same, amplitude 0.5 | 11.24 | 9.5 |

Information saturates with amplitude: half the amplitude keeps 94% of $\lambda_{\min}$ for 27% of the energy. Cramér–Rao relative standard deviations at the optimum are α 6.1%, λ 6.6%, β₊ 1.3%, β₋ 1.8%, A 0.7%, B₀ 2.9%, M₀ 11%, E₀ 29%.
Local asymptotic minimax risk (Hájek–Le Cam): $\sup_\theta \mathrm{tr}\,\mathcal I(\theta)^{-1}=0.160$ over a ±30% box (12 samples, log-parameters).

**Estimation on synthetic data [Numerical].** Nonlinear least squares on the first half of the record:
- All 8 true values fall inside their 95% intervals.
- Out-of-sample MSE on the second half is $9.5\times10^{-5}$, about the noise variance $10^{-4}$.
- The constant-mean null model scores $1.09\times10^{-2}$, so $E_{\rm test}[(K_{\rm pred}-K_{\rm obs})^2]<E_{\rm test}[(K_{\rm null}-K_{\rm obs})^2]$ holds by a factor of 115.

Real data would replace the synthetic record in the same pipeline.

---

## HDSP-03 — average-reward optimum and HJB certificate

**Lemma 3.1 (verification) [Proof].** The set $\Omega=[0,R_{\max}]^3\ni(S,E,R)$ is forward invariant and compact. If $V\in C^1(\Omega)$ and $\ell(z,u)+\nabla V(z)\cdot F(z,u)\le g$ on $\Omega\times U$, then $g^*\le g$.
*Proof.* Along any admissible trajectory $\frac1T\int_0^T\ell\,dt\le g-\frac{V(z_T)-V(z_0)}{T}\to g$, since $V$ is bounded on $\Omega$. ∎
The supremum over $u$ is explicit: $\sup_{|u|\le v}[pu-\mu|u|-\kappa u^2/2]=\psi(p)$, a piecewise-quadratic convex function. So $g_{\rm upper}(V)=\max_\Omega\big[\ell_0+V_S\lambda(R-S)+V_E\beta_s(R-E)+\psi(V_R)\big]$.

**Proposition 3.2 (constant policies) [Proof].** Under constant $R$ the long-run average is $K_0/D-cR$, which is maximised at $R=0$. So $G_{\rm const}=K_0/D=0.6667$. Any gain must therefore come from varying $R$.

**Mechanism [Proof, elementary].** Over a period $\langle\dot E\rangle=0$, which gives $\beta_+\langle(R-E)_+\rangle=\beta_-\langle(R-E)_-\rangle$. Positive surprises can be timed for when tolerance is low ($\eta$ large) and negative ones for when it is high ($\eta$ small). That timing is the only source of value, and it requires $M$.

**Bounds on $g^*-K_0/D$** (Chebyshev polynomial $V$ of total degree $d$, LP with cutting planes, then maximum of $H$ on a $161^3$ grid):

| $d$ | LP value | grid-verified upper bound |
|---|---|---|
| 1 | 0.18514 | 0.18517 |
| 2 | 0.008125 | 0.008130 |
| 4 | 0.003733 | 0.003780 |
| 6 | 0.002620 | 0.002944 (max over $161^3$, $201^3$, $241^3$ grids) |

Lower bounds (exact long-run averages of explicit periodic policies):

| Policy | Gain over $K_0/D$ |
|---|---|
| trapezoid with ramps at $v_{\max}$ (HDSP-04) | 0.001924 |
| one common free ramp speed | 0.002271 |
| separate up/down ramp speeds $(0.551, 0.329)$ | **0.002306** |

**Status.** $0.002306\le g^*-K_0/D\le 0.00294$. The best explicit policy therefore achieves at least 78% of the largest possible improvement. The $d=4$ bound (0.00378) is stable across grids and random sampling ($8\times10^6$ points). The $d=6$ polynomial is visibly more oscillatory: its LP value rises from 0.00262 to 0.00294 when checked off the cut points.

- **[Open] $\varepsilon\to0$.** Degree ≥ 8 polynomials oscillate between sample points, and the LP needs regularisation or SOS constraints.
- **Not interval-certified.** The upper bounds are maxima on a grid. Finite-difference Lipschitz margins are far too large to make them proofs; interval branch-and-bound or SOS would be needed.
- **Equality of the sup with the LP value.** On the theory side, $\inf_V g_{\rm upper}(V)=g^*$ is expected for this compact, control-affine problem whose reward is concave in $u$. This follows the linear-programming approach to deterministic long-run average control of Gaitsgory and co-authors. It is not proved here for the non-smooth β-switch.

---

## HDSP-04 — optimal periodic (trapezoid) policy

The family: $R$ rises at $v_{\max}$ to height $a$, holds until $qT_p$, falls at $v_{\max}$, then stays 0.
The long-run average is computed **semi-analytically**. On each linear piece of $R$, $E$ and $M$ are closed-form; the switch $R=E$ happens at $\tau^*=\beta^{-1}\ln(1-d_0\beta/s)$; only $\int\eta(R-E)$ needs Gauss–Legendre quadrature. This agrees with an adaptive DOP853 ODE solution to $10^{-13}$.

**Results [Numerical].**
- **Optimum:** $T_p^*=49.031$, $q^*=0.24918$, $a^*=R_{\max}=1$. $G^*=0.668591$, i.e. $G^*-G_{\rm const}=+0.001924$ (+0.29%).
- **First-order conditions:** $\partial G/\partial T_p=1\times10^{-10}$ and $\partial G/\partial q=1\times10^{-9}$ (≈0).
- **$a^*$ is on the boundary.** $\partial G/\partial a=+0.00344>0$, so "$\partial G/\partial a=0$" does not hold; KKT holds with a positive multiplier.
- **Hessian in $(T_p,q)$:** eigenvalues $-0.0725$ and $-2.6\times10^{-6}$, so $\nabla^2G\prec0$. The $T_p$ direction is very flat.
- **Uniqueness:**
  - One positive local maximum on a 60×49 grid at $a=1$.
  - 13 of 16 random multistarts reach the same point.
  - The other 3 drift to $T_p\to\infty$, where $G\to G_{\rm const}$; this is the trivial "never cycle" limit.
- **Existence [Proof sketch].** $G-G_{\rm const}\to0^-$ as $q\to0$, as $a\to0$, and as $T_p\to\infty$ (contrast is $O(1/T_p)$). Since an interior point has $G>G_{\rm const}$, the maximum is attained on a compact set bounded away from these limits.
- **Optimal amplitude by level:** gain is negative for $a\le0.3$ and grows monotonically to $a=1$.
- **Answer to "$G^*_{\rm periodic}=\sup_{U_{\rm adm}}G$?" — no, not for this family.** Freeing the ramp speeds gives 0.002306, which is 20% better. The remaining gap to the HJB upper bound is given in HDSP-03.

---

## HDSP-05 — reduced hysteresis model, thresholds, singular perturbation

**Lemma 5.1 (exact contrast integral) [Proof].** If $R$ increases monotonically from $R_0$ to $R_1$ starting at $E=R_0$, then $R\ge E$ throughout and
$$\int (R-E)\,dt=\frac{R_1-R_0}{\beta_+}$$
for any ramp shape, once $E$ has relaxed. Symmetrically, $-(R_0-R_1)/\beta_-$ for a decrease.
*Proof.* $\int(R-E)\,dt=\int \dot E/\beta_+\,dt$. ∎

**Reduced model** ($\lambda\to0$ with $a=\alpha/\lambda$ and $\hat c=c/\lambda$ fixed; $E$ fast relative to $S$). Per on/off cycle between $S_L$ and $S_H$ with level $r$:
$$N=\int_{S_L}^{S_H}\frac{g(s)}{\lambda(r-s)}ds+b(S_L)-2hr,\qquad
g(s)=-cr-\frac{\lambda r(r-s)\eta'(s)}{D\beta_-},\qquad b(S)=\frac rD\Big(\frac1{\beta_+}-\frac1{\beta_-}\Big)\eta(S),$$
with $h=\mu+\kappa v_{\max}/2$. Also $T_A=\lambda^{-1}\ln\frac{r-S_L}{r-S_H}$, $T_P=\lambda^{-1}\ln\frac{S_H}{S_L}$, and $J=N/(T_A+T_P)$.
This is the spec's form plus the boundary term $b$, which vanishes iff $\beta_+=\beta_-$.

**Theorem 5.2 (first-order conditions) [Proof].** Let $\varphi(s)=s\,g(s)$ and $\psi(s)=\varphi(s)-\lambda s(r-s)b'(s)$. Then
$$\partial_{S_H}J=\frac{\varphi(S_H)-rJ}{\lambda S_H(r-S_H)T},\qquad
\partial_{S_L}J=-\frac{\psi(S_L)-rJ}{\lambda S_L(r-S_L)T}.$$
Interior optima therefore satisfy $\varphi(S_H)=rJ$ and $\psi(S_L)=rJ$.

**The stated condition $S_Lg(S_L)=rJ$ is correct only when $\beta_+=\beta_-$.** At the nominal β₊=1, β₋=0.8:
- the corrected condition has residual $9.6\times10^{-11}$;
- the stated one has residual $1.3\times10^{-3}$;
- $S_Hg(S_H)=rJ$ holds to $10^{-10}$.

**Theorem 5.3 (second order) [Proof].** At a critical point
$$\partial^2_{S_H}J=\frac{\varphi'(S_H)}{\lambda S_H(r-S_H)T},\qquad \partial^2_{S_L}J=-\frac{\psi'(S_L)}{\lambda S_L(r-S_L)T},\qquad \partial_{S_L}\partial_{S_H}J=0 .$$
So the Hessian is diagonal. It gives a strict local maximum iff $\varphi'(S_H)<0<\psi'(S_L)$: $S_L$ on the rising branch and $S_H$ on the falling branch. When both diagonal entries are negative, $\det\nabla^2J>0$ follows.
Numerical check: $H=\mathrm{diag}(-1.0099,-0.031181)$ with off-diagonal $4\times10^{-8}$. The theory predicts $\partial^2_{S_H}J=-0.0311808$, matching to $10^{-9}$. $\varphi$ is unimodal (one sign change of $\varphi'$, peak at $s=0.185$).

**Theorem 5.4 (existence, characterisation) [Proof].** Assume $\max_S b(S)<2hr$; otherwise infinitely fast cycling looks profitable and the reduction breaks down. Nominally $b<0$. Then $J\to0$ as $S_L\to0$, $J\to-cr$ as $S_H\to r$, and $J\to-\infty$ as $S_H-S_L\to0$. So if some loop has $J>0$, the maximum is attained in the interior. $J^*$ is the unique root of the Dinkelbach function $\Phi(J)=\max_{S_L<S_H}[N-JT]$, which is strictly decreasing.
**Nominal values:** $S_L^*=0.0333$, $S_H^*=0.5820$, $T_A=8.38$, $T_P=28.61$, $J^*=0.005156=0.0516\,\lambda$. With $r$ left free, the optimum is still $r^*=R_{\max}$.

**Theorem 5.5 (single loop is optimal among on/off policies) [Proof].** In the reduced model with $R\in\{0,r\}$, every periodic switching policy has $J\le J^*_{\rm loop}$.
*Proof.*
1. On a closed path in 1-D, the number of up-crossings $n(s)$ of each level $s$ equals the number of down-crossings.
2. Write the up-runs as intervals $[a_i,b_i]$. Up-run $i$ earns $N([a_i,b_i])$ and pays $2hr$. The active time is $\sum_i T_A([a_i,b_i])$. The passive time is $\int n(s)/(\lambda s)\,ds=\sum_i T_P([a_i,b_i])$.
3. Therefore $J=\sum_iN_i/\sum_iT_i\le\max_iN_i/T_i\le J^*_{\rm loop}$.
4. Non-periodic policies are covered by closing a long trajectory segment at bounded cost, which changes the average by $O(1/T)$. ∎

**Multilevel policies [Numerical].** The global-optimality question is "$J_{\rm hysteresis}^*\ge J^*_{\rm all}$?". Optimising 3- and 4-level staircases with levels in $[0,R_{\max}]$ gives no improvement over two levels (≤ $5\times10^{-18}$); the optimisers collapse back to the 2-level loop. No counterexample found; this part is not proved.

**Singular perturbation [Numerical].** "Full" means the best trapezoid in the full model, computed with $\alpha=a\lambda$ and $c=\hat c\lambda$:

| λ | $J_{\rm full}/\lambda$ | $J_{\rm red}/\lambda$ | rel. error | rel. error / λ | $T_p$ / $T_p^{\rm red}$ | $q$ / $q^{\rm red}$ |
|---|---|---|---|---|---|---|
| 0.4 | 0 (no cycling pays) | 0.0516 | 1 | — | ∞ / 9.2 | — |
| 0.1 | 0.0192 | 0.0516 | 0.627 | 6.27 | 49.0 / 37.0 | 0.249 / 0.227 |
| 0.025 | 0.0405 | 0.0516 | 0.215 | 8.60 | 160.0 / 148.0 | 0.237 / 0.227 |
| 0.00625 | 0.0485 | 0.0516 | 0.0599 | 9.59 | 604.2 / 591.9 | 0.230 / 0.227 |
| 0.0015625 | 0.0508 | 0.0516 | 0.0155 | 9.89 | 2379.9 / 2367.5 | 0.2275 / 0.2266 |

- The local convergence exponents rise 0.57 → 0.72 → 0.82 → 0.90 → 0.94 → 0.97 → 0.985. So $\lvert J^*_{\rm full}-J^*_{\rm red}\rvert\approx 9.9\,\lambda\,J^*_{\rm red}\approx0.51\,\lambda^2$. That is **$p=1$ relative ($p=2$ absolute)**.
- The policies converge too. $T_p-T_p^{\rm red}\to12.4$, a constant offset from ramp and expectation transients, so the relative error is $O(\lambda)$, and $q\to q^{\rm red}$.
- Large constant ≈ 10: the asymptotic regime needs $\lambda\ll S_L^*\approx0.03$. For λ ≳ 0.1 the reduced model badly overestimates the value of cycling, and at λ = 0.4 cycling does not pay at all in the full model.
- A rigorous Tikhonov/Fenichel-type error bound is **[Open]**. The mechanism behind the $O(\lambda)$ per-cycle error is the $O(\lambda)$ drift of $S$ during each $E$-transient, which Lemma 5.1 isolates.

---

## HDSP-06 — stochastic and robust layer

Model: $dx=f\,dt+G\,dW$ with $G=\mathrm{diag}(0.05,0.05,0.12)$, $K_{\min}=0.5$, $\varepsilon=0.05$.

**Proposition 6.1 (mean-square boundedness) [Proof].** With $V=M^2+E^2+B^2$,
$$\mathcal LV\le-k_0V+C,\qquad k_0=\min(\lambda,\beta_{\min},1),\qquad C=\tfrac{R_{\max}^2}{\lambda}+\beta_{\max}R_{\max}^2+C_B+\mathrm{tr}\,GG^\top,$$
where $C_B=\max_B[2Bf_B(B)+B^2]<\infty$. The bound uses $2MR\le\lambda M^2+R^2/\lambda$, $2\beta E(R-E)\le-\beta E^2+\beta R^2$, and $2Bf_B\le-B^2+C_B$ (the quartic term dominates). Hence $\sup_tE\|x_t\|^2\le\max(V(x_0),C/k_0)$ for every policy. This is $C/k_0=120.9$ nominally and 238.2 at worst on a ±30% box.

**Proposition 6.2 (the chance constraint is structurally infeasible) [Proof].** $B$ is a 1-D diffusion decoupled from $R$. Its stationary density is $p(B)\propto e^{-2U(B)/\sigma_B^2}$ exactly. For symmetric wells and $u_B=0$, this gives $P(B<b_2)=\tfrac12$. Hence the long-run fraction of time with $K<K_{\min}$ is at least the lower-well mass, $\tfrac12$, for **every** policy π.
Monte Carlo: 0.570 for $R\equiv0$ and 0.570 for π*.

What it takes to fix this [Numerical]:
- A tilt $u_B=0.0306$ brings the lower-well mass to $\varepsilon/2$.
- Even then, in-well noise gives $P(K<K_{\min})=0.061$ for the baseline and 0.064 for π*. A slightly larger tilt or a lower noise level is needed.
- π* adds about 0.3 percentage points through its switch-off contrast dips.
- The noise-free cycle stays above $K_{\min}$ ($\min K=0.515$).

**Robust performance [Numerical].** For the nominal-optimal π* (HDSP-04), take the worst case over the ±30% box on $(\alpha,\lambda,\beta_+,\beta_-,D)$:
- $J_{\rm robust}(\pi^*)=0.5504<J_{\rm robust}(R\equiv0)=0.5556$, so **π* is not robust**.
- The worst vertex has low α, high λ, high β₊, low β₋ and high D.
- The max-min policy over the trapezoid family is $R\equiv0$.

**Robustness radius:**

| box half-width ρ | 2% | 5% | 8% | 10% | ≥ 12% |
|---|---|---|---|---|---|
| best worst-case gain | 0.00149 | 0.00087 | 0.00031 | 0.0000045 | 0 (never cycle) |
| max-min $T_p$, $q$ | 51.5, 0.249 | 56.4, 0.244 | 65.6, 0.224 | 103.7, 0.146 | — |

Robust design lengthens the cycle and shortens the active phase. **$J_{\rm robust}(\pi^*)>J_{\rm robust}(\pi_{\rm baseline})$ holds iff the parameter uncertainty is below about 10–11%.** HDSP-02's input design gives about 6% relative standard deviation on α and λ, i.e. roughly ±12% at 95%. That sits right at the edge, so better identification is what would make cycling worthwhile in practice.

**$F$-floor (assuming $F:=R$) [Numerical].** With $\langle R\rangle\ge F_{\min}$, periodic policies still beat the best constant policy ($R=F_{\min}$):

| $F_{\min}$ | periodic gain | constant gain |
|---|---|---|
| 0.2 | +0.00192 | −0.0010 |
| 0.3 | +0.00184 | −0.0015 |
| 0.5 | +0.00045 | −0.0025 |

---

## HDSP-Ω — synthesis

**Minimal dimension, by target [Proof]:**

| Target | Minimal dimension |
|---|---|
| paradox $R_1<R_2\wedge K_1>K_2$ | $B$ alone (1) |
| paradox produced from rest by the input alone | $(E,B)$ (2) |
| a non-trivial optimal policy (cycling beats constant) | $M$ is also required (3) |

For the last row: without $M$ the gain is $\langle(R-E)\rangle/D$, which is ≤ 0 whenever $\beta_+\ge\beta_-$.

**Constraint status:**

| Constraint | Status |
|---|---|
| structural identifiability | holds modulo $A+\sigma_N^2$, given an exciting input (HDSP-02) |
| observability | holds off $\{R=E\}$ |
| stability | Lyapunov certificates (HDSP-01) |
| mean-square boundedness | proved (HDSP-06) |
| robust feasibility of the $K$ chance constraint | fails without a tilt $u_B$ |
| prediction error | at the noise floor on synthetic data |
| optimality gap | $\le 6.4\times10^{-4}$ absolute (≈ 22% of the gain), not yet → 0 |

**Open items:**
1. Rigorous (interval or SOS) closure of the HJB gap.
2. Global optimality over all continuous $R(t)$ (multilevel policies were tested numerically only).
3. A singular-perturbation theorem matching the observed $p=1$.
4. Fitting to real data, and a definition of $F$.

---

## Reproducibility

```
pip install -r requirements.txt
cd hdsp
python signature.py            # authorship signature check
python hdsp01_bistability.py   # ~5 s
python hdsp02_ident.py         # ~3 min
python hdsp04_periodic.py      # ~2 min
python hdsp05_reduced.py       # ~5 min
python hdsp06_robust.py        # ~5 min (Monte Carlo)
python hdsp03_lp.py 2 4 6      # ~30 min
```

Results are written to `hdsp/results_*.json`; the committed files are the ones quoted in this report. The helper checks `hdsp04_ramps.py`, `hdsp06_radius.py` and `hdsp03_verify.py` reproduce the asymmetric-ramp bound, the robustness radius and the fine-grid LP verification.
