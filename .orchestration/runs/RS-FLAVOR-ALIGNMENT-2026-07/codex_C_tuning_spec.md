**ANALYSIS/SPEC: F1 Phase-Volume Explanation**

Repo anchors: `C4_LR = -(G_L G_R)/M_KK^2` in [deltaf2.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/deltaf2.py:312), kaon epsilon/mass formulas in [deltaf2.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/deltaf2.py:768), instrumented true phase in [instrument_epsK_phase.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/scripts/instrument_epsK_phase.py:1).

Let

\[
R_K \equiv \frac{2 |M_{12}^{LR}|}{\Delta m_K},\qquad
\theta_K = \arg[(G_L^d)_{12}(G_R^d)_{12}].
\]

Then

\[
|\epsilon_K^{NP}| =
\frac{\kappa_\epsilon}{\sqrt{2}\Delta m_K}|M_{12}^{LR}\sin\theta_K|
=
\frac{\kappa_\epsilon}{2\sqrt{2}}R_K|\sin\theta_K|
\simeq 0.33234\,R_K|\sin\theta_K|.
\]

For anarchic CP, \(\theta_K\) is uniform on \([0,2\pi)\). Therefore

\[
P_{\rm pass}(R_K)
=
P\left(|\sin\theta_K| < \frac{\epsilon_{\max}}{0.33234 R_K}\right)
=
\begin{cases}
1, & R_K \le R_0,\\
\frac{2}{\pi}\arcsin\left(\frac{\epsilon_{\max}}{0.33234 R_K}\right), & R_K > R_0,
\end{cases}
\]

with \(R_0=\epsilon_{\max}/0.33234\). At large \(R_K\),

\[
P_{\rm pass}(R_K) \simeq
\frac{2\epsilon_{\max}}{\pi\,0.33234}\frac{1}{R_K}.
\]

Budgets:

- Observed \(\epsilon_K\): \(\epsilon_{\max}=2.228\times10^{-3}\), \(R_0=6.704\times10^{-3}\), large-\(R\) coefficient \(4.268\times10^{-3}/R_K\).
- Central residual SM-vs-exp budget used by the core deltaf2 path: \(\epsilon_{\max}=|2.228-2.161|\times10^{-3}=6.70\times10^{-5}\), \(R_0=2.016\times10^{-4}\), large-\(R\) coefficient \(1.283\times10^{-4}/R_K\).
- Note: the catalog K001 loose budget is \(3.0371\times10^{-4}\), but `scripts/instrument_epsK_phase.py` uses the core residual unless rerun with an override.

Why F1 follows:

The epsilon cut changes the pre-cut magnitude prior \(p_0(R_K)\) into

\[
p_{\rm surv}(R_K)\propto p_0(R_K)P_{\rm pass}(R_K).
\]

For \(R_K \gg R_0\), this is \(p_0(R_K)/R_K\). So the cut automatically tilts survivors toward smaller \(|C_4|\), even if the original draw was fully anarchic. The 57% magnitude / 7% phase split is therefore not an independent empirical fact; it is the empirical \(p_0(R_K)\) folded with this \(1/R_K\) phase volume.

The positive correlation is also expected. With \(s=|\sin\theta|\),

\[
\epsilon_K \sim R_K s,\qquad \Delta m_K \sim R_K
\]

for the repo’s `ratio_dm_K = |M12|/(\Delta m_K/2)`. Before the cut, \(\mathrm{Cov}(R_Ks,R_K)>0\). After the cut, \(E[\epsilon_K|R_K,\mathrm{pass}]\) rises with \(R_K\) and saturates near \(\epsilon_{\max}/2\), so the correlation remains positive unless the ensemble is dominated only by the saturated tail. A pure phase-cancellation story would instead emphasize large \(R_K\) and tiny \(s\), producing an anti-correlation in \(\epsilon_K/R_K\), not the observed \(\mathrm{corr}(\log\epsilon_K,\log\Delta m_K)\simeq+0.82\).

```python
import numpy as np
import pandas as pd

# Inputs
path = ".orchestration/runs/RS-FLAVOR-ALIGNMENT-2026-07/instrument_S1.parquet"
df = pd.read_parquet(path)

KAPPA_EPS = 0.94
A_EPS = KAPPA_EPS / (2.0 * np.sqrt(2.0))  # 0.332340187...
EPS_EXP = 2.228e-3
EPS_RESID = abs(2.228e-3 - 2.161e-3)

eps_max = EPS_RESID  # use EPS_EXP for the observed-epsilon exercise

# Pre-epsilon population
pre = df[
    (df["passes_pdg"])
    & (df["M_KK_TeV"] <= 3.0)
    & np.isfinite(df["ratio_eps_K"])
    & np.isfinite(df["ratio_dm_K"])
    & np.isfinite(df["Phi_12"])
].copy()

# R_K is exactly the repo Delta-m ratio if LR dominates; validate below.
pre["R_K"] = pre["ratio_dm_K"].astype(float)
pre["sinPhi"] = np.abs(np.sin(pre["Phi_12"].astype(float)))
pre["eps_np"] = pre["ratio_eps_K"] * eps_max
pre["eps_pred"] = A_EPS * pre["R_K"] * pre["sinPhi"]
pre["eps_pass"] = pre["ratio_eps_K"] <= 1.0

# Sanity check: LR phase law should cluster near 1 for C4-dominated kaons.
pre["phase_law_ratio"] = pre["eps_np"] / pre["eps_pred"].replace(0.0, np.nan)
print(pre["phase_law_ratio"].replace([np.inf, -np.inf], np.nan).describe())

# Phase-volume law
x = eps_max / (A_EPS * pre["R_K"])
pre["P_pass_th"] = np.where(
    x >= 1.0,
    1.0,
    (2.0 / np.pi) * np.arcsin(np.clip(x, 0.0, 1.0)),
)

Pbar_emp = pre["eps_pass"].mean()
Pbar_th = pre["P_pass_th"].mean()
Delta_vol = 1.0 / Pbar_th
print({"P_emp": Pbar_emp, "P_th": Pbar_th, "Delta_vol": Delta_vol})

# Bin verification: empirical pass fraction vs arcsin law.
bins = np.geomspace(
    pre.loc[pre["R_K"] > 0, "R_K"].quantile(0.005),
    pre["R_K"].quantile(0.995),
    25,
)
pre["R_bin"] = pd.cut(pre["R_K"], bins=bins, include_lowest=True)

tab = (
    pre.groupby("R_bin", observed=True)
    .agg(
        n=("eps_pass", "size"),
        R_med=("R_K", "median"),
        p_emp=("eps_pass", "mean"),
        p_th=("P_pass_th", "mean"),
    )
    .reset_index()
)
tab["p_emp_err"] = np.sqrt(tab["p_emp"] * (1.0 - tab["p_emp"]) / tab["n"].clip(lower=1))
print(tab)

# Correlation check.
floor = 1e-300
corr_all = np.corrcoef(
    np.log10(np.maximum(pre["ratio_eps_K"], floor)),
    np.log10(np.maximum(pre["ratio_dm_K"], floor)),
)[0, 1]

surv = pre[pre["eps_pass"]].copy()
corr_surv = np.corrcoef(
    np.log10(np.maximum(surv["ratio_eps_K"], floor)),
    np.log10(np.maximum(surv["ratio_dm_K"], floor)),
)[0, 1]
print({"corr_log_all": corr_all, "corr_log_survivors": corr_surv})

# F1 split, matching ledger language.
# "typical" = median R_K of epsilon failers at same M_KK tile.
typ_by_mkk = (
    pre.loc[~pre["eps_pass"]]
    .groupby("M_KK_TeV")["R_K"]
    .median()
)
pre["R_fail_typ"] = pre["M_KK_TeV"].map(typ_by_mkk)
surv = pre[pre["eps_pass"]].copy()

surv["magnitude_survivor"] = surv["R_K"] < 0.1 * surv["R_fail_typ"]
surv["phase_survivor"] = surv["R_K"] >= surv["R_fail_typ"]

print({
    "n_survivors": len(surv),
    "frac_magnitude": surv["magnitude_survivor"].mean(),
    "frac_phase": surv["phase_survivor"].mean(),
    "frac_intermediate": 1.0 - surv["magnitude_survivor"].mean() - surv["phase_survivor"].mean(),
})

# Reweighted-prior proof: these should reproduce observed fractions.
w = pre["P_pass_th"]
Z = w.sum()

mag_region = pre["R_K"] < 0.1 * pre["R_fail_typ"]
phase_region = pre["R_K"] >= pre["R_fail_typ"]

print({
    "pred_frac_magnitude_from_volume": w[mag_region].sum() / Z,
    "pred_frac_phase_from_volume": w[phase_region].sum() / Z,
})
```

Single decisive plot:

- x-axis: \(\log_{10} R_K\), equivalently \(\log_{10}\) `ratio_dm_K`.
- y-axis: \(\log_{10}|\sin\Phi_{12}|\).
- Overlay the epsilon-pass boundary \( |\sin\Phi_{12}|=\epsilon_{\max}/(0.33234R_K)\), plus vertical markers at \(0.1\,R_{\rm fail,typ}\) and \(R_{\rm fail,typ}\).

This plot separates the mechanisms visually: magnitude survivors are on the left, where \(|C_4|\) is small; phase survivors are the bottom-right tail, where \(|C_4|\) is typical/large but \(|\sin\Phi_{12}|\) is tiny. The relative population of those regions is the 57% vs 7% statement.

Three-sentence result:

F1 rigorously establishes that the epsilon_K cut imposes a calculable prior-volume penalty \(P_{\rm pass}(R_K)\sim1/R_K\) on large LR kaon amplitudes. Folding this law with the anarchic pre-cut \(R_K\) distribution predicts that most low-M_KK epsilon_K survivors are magnitude-suppressed, while genuinely phase-aligned, Delta-m_K-loud survivors are a small tail. Therefore flat anarchy does not mainly survive by tuning CP phases; it mostly survives by accidentally drawing small down-sector flavor-violating magnitudes, with phase alignment as the minority path.