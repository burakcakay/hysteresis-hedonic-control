# hysteresis-hedonic-control

**When does "less" feel better than "more", and when does it pay to cycle?**
A worked study of a small hedonic dynamical system: tolerance $M$, expectation $E$ (with asymmetric updating), a bistable baseline $B$, and output $K$. The study covers minimal models, stability, identifiability, long-run optimal control, and robustness. It answers the problem set **HDSP-01 … HDSP-06, HDSP-Ω** ([`problems/HDSP.json`](problems/HDSP.json)).

Full write-up with proofs: **[REPORT.md](REPORT.md)**. Every claim there is tagged **[Proof]**, **[Numerical]**, **[Ill-posed]** or **[Open]**.

## Key findings

1. **The paradox lives in one coordinate.** At every steady state, expectation catches up with reward ($E^*=\bar R$), so the hedonic term vanishes and $K^*=B^*/D$. "Less reward, more satisfaction" at steady state therefore needs exactly one bistable coordinate, $B$, and $n_{\min}=1$. Producing it from rest with the input alone needs two, $(E,B)$. *[Proof]*
2. **Stability holds even though $Df$ does not exist.** Every equilibrium sits on the kink $R=E$ of the asymmetric expectation update. Stability is certified by the Clarke Jacobian and an explicit Lyapunov function, and the basins are exactly $\{B\lessgtr b_2\}$. *[Proof]*
3. **Identifiability.** Tolerance is invisible whenever $R=E$. All dynamic parameters are identifiable except $A+\sigma_N^2$, but only if the input pushes $R$ both above and below $E$. Cost parameters are never identifiable from (input, output) data. The stated "minimum-energy input with full observability rank" problem has infimum 0; a Fisher-information design replaces it. *[Proof + Numerical]*
4. **Optimal cycling is a hysteresis loop.** The best constant policy is $R=0$. Cycling pays only by timing positive surprises at low tolerance. In the slow-tolerance limit:
   - optimal switching thresholds satisfy $\varphi(S_H)=rJ$;
   - the stated companion condition for $S_L$ is **correct only for symmetric expectation updating**, and the corrected condition is derived;
   - the Hessian is diagonal at critical points;
   - one hysteresis loop beats every other on/off policy;
   - the full model converges with error $\approx0.51\lambda^2$. *[Proof + Numerical]*
5. **Optimality is bracketed, robustness is fragile.** An HJB/LP certificate brackets the optimal gain within $[0.00231,\ 0.00294]$ over the constant policy. The nominal-optimal cycle loses to "do nothing" once parameter uncertainty exceeds about 10%. With symmetric wells, the chance constraint on $K$ is infeasible for **every** policy. *[Numerical + Proof]*

## Reproduce

```bash
pip install -r requirements.txt
cd hdsp
python signature.py
python hdsp01_bistability.py
python hdsp02_ident.py
python hdsp04_periodic.py && python hdsp04_ramps.py
python hdsp05_reduced.py
python hdsp06_robust.py && python hdsp06_radius.py
python hdsp03_lp.py 2 4 6 && python hdsp03_verify.py 2 4 6   # slow (~30 min)
```

Each script writes its results to `hdsp/results_*.json`; the committed JSON files are the ones the report quotes.

| File | Role |
|---|---|
| `hdsp/model.py` | model, nominal parameters, fast vectorised periodic evaluator |
| `hdsp/periodic_exact.py` | semi-analytic periodic evaluator (closed-form E, M; exact switch times) |
| `hdsp/hdsp0*_*.py` | one script per problem; each writes `results_*.json` |
| `hdsp/V_coef_d*.npy` | HJB subsolution coefficients (LP certificates) |

## How this was made

The problem set was formulated by **burakcakay**. Analysis, code and proofs were developed together with Claude (Anthropic) as a computational assistant. Every result is labelled with how much it is worth: a proof, a reproducible computation, or an open question.

## Signature

A visible, documented authorship mark. It is not hidden, and it is not addressed to any program.

$$\sigma(x)=98x^9+117x^8+114x^7+97x^6+107x^5+99x^4+97x^3+107x^2+97x+121,$$

$$\mathfrak S=\sigma(256)=464958426364716893561209 .$$

Write $\mathfrak S$ in base 256 and read the digits as ASCII. `python hdsp/signature.py` does it for you and asserts the result.

## Cite / license

See [`CITATION.cff`](CITATION.cff). Code: MIT ([`LICENSE`](LICENSE)). Text and proofs: CC BY 4.0 ([`LICENSE-docs.md`](LICENSE-docs.md)). Please credit **burakcakay**.

---

## Türkçe özet

Bu repo, hedonik bir dinamik sistem üzerine kurulu HDSP-01…06 ve HDSP-Ω problem setinin çözümüdür. Sistemin parçaları: tolerans $M$, beklenti $E$ (asimetrik güncellenir), iki kararlı durumu olan taban $B$ ve çıktı $K$.

Ana sonuçlar:
- **"Az ödül, çok tatmin" paradoksu.** Dengede beklenti ödüle yetiştiği için hedonik terim sıfırlanır ($K^*=B^*/D$). Paradoks bu yüzden tek bir koordinatla, $B$ ile ortaya çıkar; $n_{\min}=1$ *[İspat]*.
- **Döngüsel politika.** Optimal döngüsel politika tek bir histerezis döngüsüdür *[İspat]*. Problem setindeki $S_L g(S_L)=rJ$ koşulu yalnızca $\beta_+=\beta_-$ olduğunda doğrudur; düzeltilmiş hali türetildi *[İspat]*.
- **Optimumun sınırları.** Optimal kazanç sayısal olarak $[0.00231, 0.00294]$ aralığına sıkıştırıldı *[Sayısal]*.
- **Sağlamlık.** Nominal optimal döngü, parametre belirsizliği %10'u aşınca "hiçbir şey yapmamak" politikasına yenilir *[Sayısal]*.
- **Şans kısıtı.** Simetrik kuyularda $K$ üzerindeki şans kısıtı hiçbir politikayla sağlanamaz *[İspat]*.

Tüm iddialar [REPORT.md](REPORT.md) içinde **[İspat] / [Sayısal] / [Açık]** olarak etiketlidir.
