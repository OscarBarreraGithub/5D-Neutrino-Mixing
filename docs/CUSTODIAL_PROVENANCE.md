# Custodial RS — implementation provenance ledger

**Purpose.** Establish, equation by equation, *exactly where every piece of our
custodial RS implementation comes from* — reference, equation number, page in the
local source PDF, where it is derived, and whether our code reproduces it
faithfully. This is Phase 0 of the custodial program: provenance **before** the
literature catalog (Phase 1) and the refreshed ours-vs-literature comparison
(Phase 2). Built from three independent source-verification audits against the
local `references/papers/` archive (2026-06-29).

**Status tags used below:**
- ✅ **SOURCED** — exact/faithful reproduction of a published equation (eq# given).
- ◐ **FAITHFUL TRUNCATION** — a deliberate, documented simplification of a
  published equation (conservative and/or flagged in code).
- ⚠️ **PROXY** — repo ansatz with the correct *parametric* form, but the
  coefficients/inputs are not derived from any paper.
- ❌ **UNSOURCED** — no published equation gives this form; it is a repo invention.

---

## 0. Source-paper inventory (all local in `references/papers/`)

| Short | Paper | arXiv / journal | Local file | Used for |
|---|---|---|---|---|
| **CGHNP** | Casagrande, Goertz, Haisch, Neubert, Pfoh | 0807.4937 | `casagrande_goertz_haisch_neubert_pfoh_2008_I.pdf` | oblique S,T,U (§6.2); Z→bb ZMA (§6.4 Eq.170) |
| **ADMS** | Agashe, Delgado, May, Sundrum | hep-ph/0308036 | `agashe_delgado_may_sundrum_2003_custodial.pdf` | custodial SU(2)_R for T (origin) |
| **ACDP** | Agashe, Contino, Da Rold, Pomarol | hep-ph/0605341 | `agashe_contino_darold_pomarol_2006_Zbb.pdf` | P_LR custodial symmetry for Z b_L b̄_L |
| **CPSW-ob** | Carena, Pontón, Santiago, Wagner (oblique) | hep-ph/0701055, PRD76 035006 | `carena_ponton_santiago_wagner_2007_oblique.pdf` | top-partner loop δg_L^b & ΔT (Eqs.13–15) |
| **CPSW-lKK** | Carena, Pontón, Santiago, Wagner (light KK) | hep-ph/0607106 | `carena_ponton_santiago_wagner_2006_lightKK.pdf` | T_top derivation (Eq.35), ΔT (Eq.43) |
| PDG-2025 | Particle Data Group EW review | — | `flavor_catalog/references/EW001/pdg_2025_electroweak_stu.txt` | S coefficient c_S, fit ellipse |

---

## 1. Master provenance table

| # | Code quantity (file) | Formula | Source | Eq# / page | Status | Gap? |
|---|---|---|---|---|---|---|
| 1 | minimal ΔT (`oblique_stu.py`) | $x_1^2\,\pi L/(2\cos^2\theta_W)\,(v/M)^2$ | CGHNP | (147), p.37 | ◐ SOURCED (leading-$L$, drops $-1/2L$, 0.04% conservative) | — |
| 2 | custodial ΔT (`oblique_stu.py`) | $-x_1^2\,\pi/(4\cos^2\theta_W L)\,(v/M)^2$ | CGHNP ← ADMS | (153), p.40 | ✅ SOURCED (exact) | in-code citation missing |
| 3 | ΔS (`oblique_stu.py`) | $c_S\,(v/M)^2,\ c_S=30$ | PDG-2025 (CGHNP S-part consistent ~36.6) | snapshot L35–39 | ✅ SOURCED (prose) | YAML anchor `UNRESOLVED` |
| 4 | ΔU (`oblique_stu.py`) | $0$ | CGHNP / ADMS | with (147),(153) | ✅ SOURCED (tree-level) | — |
| 5 | minimal $g_L^b$ (`rs_ew_couplings.py`) | $+(m_b^2/2M^2)\,B_d$ | CGHNP | (170) 1st line, p.48 | ✅ SOURCED (structure+sign) | convention-dict caveat |
| 6 | minimal $g_R^b$ (`rs_ew_couplings.py`) | $-(m_b^2/2M^2)\,B_Q$ | CGHNP | (170) 2nd line, p.48 | ✅ SOURCED (L↔R symmetric) | convention-dict caveat |
| 7 | custodial $g_L^b$ zeroing | $g_L^b\to0$ on protected mask | ACDP | (3)–(12), p.3–5 | ✅ SOURCED (P_LR, (2,2)$_{2/3}$ embed) | — |
| 8 | custodial $g_L^b$ **residual** | $\kappa_b\,(1/L)\,g_L^{b,\min}$ | — | — | ❌ **UNSOURCED** (ad hoc knob) | **GAP** |
| 9 | custodial $g_R^b$ = 0 (`elementary_zero`) | $g_R^b\to0$ | ACDP | (22)/Table 1, p.8–9 | ◐ SOURCED **only for** $b_R\in(1,1)_{-1/3}$ | **rep not pinned in code** |
| 10 | loop $\delta g_L^b$ singlet | $\frac{\alpha}{16\pi s^2 M_W^2}\frac{m_{q0t}^4}{M_t^2}[1+\dots]$ | CPSW-ob | (14), p.10 | ✅ SOURCED (exact) | — |
| 11 | loop $\delta g_L^b$ bidoublet | $\frac{\alpha}{32\pi s^2 M_W^2}m_t^2[\dots]$ | CPSW-ob | (15), p.10 | ✅ SOURCED (exact) | — |
| 12 | loop ΔT singlet | $T_{\rm top}\,\frac{2m_{q0t}^2}{M_t^2}[\dots]$ | CPSW-ob / -lKK | (13)/(43)+(35) | ◐ SOURCED **but truncated** (omits neg. bidoublet ΔT) | needs override for full |
| 13 | loop $\delta g_R^b$ = 0 | $0$ | CPSW-ob | absent by design, p.10–11 | ✅ FAITHFUL (paper assumes none) | — |
| 14 | loop inputs $\xi,\rho$ knobs | $M_t=\rho_t M,\ m_{q0t}=\xi_s m_t/F$ | — | — | ⚠️ **PROXY** (breaks wavefn correlation) | **GAP** |

---

## 2. Per-sector detail

### 2.1 Oblique S, T, U  (`quarkConstraints/oblique_stu.py`)

**Convention chain (correct).** CGHNP's $M_{KK}$ is the *geometric* scale
$M_{KK}\equiv\Lambda_{IR}=k\epsilon$ (CGHNP Eq.4); the repo's physical
$M_{KK}=x_1\Lambda_{IR}$, $x_1=2.4487$ (`GAUGE_KK_ROOT_NN`). Both T coefficients
are multiplied by $x_1^2$ (`PHYSICAL_MKK_OVER_LAMBDA_IR_SQ`) to convert geo→physical;
$c_S$ is NOT (PDG already calibrates S in physical $M_{KK}$). Numerically verified:
$x_1^2\cdot\pi L/(2\cos^2)=428.80$, $x_1^2\cdot(-\pi/(4\cos^2 L))=-0.17502$ — both
match the code; internal ratio $|T_{\min}/T_{\rm cust}|=2L^2=2450$ holds exactly.

- **T-min** = CGHNP **(147)**: $T=\frac{\pi v^2}{2\cos^2\theta_W M_{KK}^2}(L-\tfrac1{2L})$.
  Repo keeps the leading $L$ (drops $-1/2L$ → 0.04% over-estimate, conservative → tighter floor).
- **T-cust** = CGHNP **(153)**: $T=-\frac{\pi v^2}{4\cos^2\theta_W M_{KK}^2}\frac1L$,
  CGHNP citing **[33]=ADMS**. The $1/L$ suppression *and* the sign flip are both
  genuinely in CGHNP (153). NB: the *closed algebraic form* originates in CGHNP;
  ADMS itself gives it in brane-Higgs-$Z'$ form (ADMS Eq. 4.9/5.6) and establishes
  the no-log-enhancement mechanism, which CGHNP distils into (153).
- **S** = $c_S(v/M)^2$, $c_S=30$ from the PDG-2025 snapshot prose ("$S\approx30v^2/M_{KK}^2$,
  $S\le0.18\Rightarrow M_{KK}\gtrsim3.2$ TeV"). CGHNP S-part (147) gives ~36.6 in
  physical units — same order, consistent.
- **U** = 0 at tree level (both models).

### 2.2 Z → b b̄ tree  (`quarkConstraints/rs_ew_couplings.py`)

- **minimal $g_L^b,g_R^b$** = CGHNP **(170)**, p.48 (§6.4). CGHNP text: "the
  non-universal corrections always reduce the couplings in magnitude," and the fit:
  "$g_L^b$ always larger than SM … $g_R^b$ essentially unaffected." Our code reproduces
  the diagonal bracket $\tfrac1{1-2c}(1/F^2-1+F^2/(3+2c))$ and the flavour-sum
  denominator $1/(1-2c_{di})$ (the B1-fix term) under the internal convention
  dictionary $c_{\rm CGHNP}=-c_{\rm repo},\ F^2_{\rm CGHNP}=2f_{IR,\rm repo}^2$. The
  algebra of that translation is self-consistent; **caveat:** it is asserted in an
  internal "audit slice 3" note, not independently checked against CGHNP App. A.
- **custodial $g_L^b$ zeroing** = ACDP **(3)–(12)**: $P_{LR}\Rightarrow c_1=c_2
  \Rightarrow\delta g_L^b=0$, requiring $q_L\in(2,2)_{2/3}$ (ACDP 7–8). Faithful.

### 2.3 Top-partner loops  (`quarkConstraints/rs_ew_couplings.py`)

- $\delta g_L^b$ **singlet** = CPSW-ob **(14)**, **bidoublet** = **(15)**, ΔT singlet =
  **(13)** with $T_{\rm top}$ = CPSW-lKK **(35)** — all reproduced character-for-character
  (prefactors $16\pi/32\pi$, $s^2M_W^2$, log structure, signs, $N_c=3$).
- $\delta g_R^b$ loop = 0 is **faithful**: CPSW-ob gives loop corrections only to ΔT
  and $Zb_L\bar b_L$; Fig.1 is drawn "assuming no large corrections to $Zb_R\bar b_R$."

---

## 3. The genuinely unsourced / proxy items (the punch-list)

These are the *only* places our custodial implementation is not a faithful copy of a
published equation. Each needs a human to source or to declare as a deliberate proxy.

1. **Custodial $g_L^b$ residual $\kappa_b\,(1/L)\,g_L^{b,\min}$** (item 8) — ❌ **UNSOURCED.**
   Neither ACDP nor CGHNP writes a $1/L$-suppressed P_LR-breaking residual. ACDP's only
   surviving-breaking estimate is *momentum-dependent*, $\delta g/g\sim(\lambda_t/g_{\rm BSM})^2
   \xi_R^{-2}(q^2/\Lambda_{\rm BSM}^2)$ — a different functional form. $\kappa_b$ is a free
   $O(1)$ dial. **To source:** derive the P_LR-breaking operator coefficient in the actual
   production rep ((2,2)$_{2/3}$ for $q_L$) and show it scales as $1/L$ — likely via
   CGHNP-style ZMA with the custodial multiplet, or ACDP Ref.[9] (Contino–Da Rold–Pomarol
   follow-up). Until then: treat as a **phenomenological knob, not a derived residual.**

2. **Top-partner loop inputs $\xi_s,\xi_q,\xi_\chi,\rho_t,\rho_q,\rho_\chi$** (item 14) —
   ⚠️ **PROXY.** The *formulas* (13–15) are exact, but the six KK/mixing masses fed into
   them are re-parametrized as independent $O(1)$ dials ($M_t=\rho_t M_{KK}$,
   $m_{q0t}=\xi_s m_t/F$, …, defaults 1.0). CPSW's central physical claim is that these six
   masses are **not independent** but locked together by the 5D wavefunctions; the repo
   breaks that correlation. **To source:** back-derive the six masses from the actual 5D
   overlap integrals (CPSW App. A) for a benchmark and check the resulting $\delta g_L^b/\Delta T$
   against CPSW Fig.1's correlation band.

3. **Custodial $g_R^b=0$ rep not pinned** (item 9) — ◐ sourced **only** for
   $b_R\in(1,1)_{-1/3}$ (ACDP Eq.22/Table 1), where $T^3_L=T^3_R=0\Rightarrow\delta g_R^b=0$
   by P_C. The code hard-codes the protected branch (`bR_strategy='elementary_zero'`) without
   enforcing/recording the embedding that earns it; ACDP Table 1 lists other embeddings giving
   **nonzero** $\delta g_R^b$ of either sign (incl. the $(2,4)$ used to fit $A_{FB}^b$).
   **To source:** declare $b_R\in(1,1)_{-1/3}$ as the production embedding and cite ACDP
   Eq.(22)/Table 1 row 1.

4. **Singlet-only ΔT truncation** (item 12) — ◐ the singlet ΔT (Eq.13) is exact, but the
   dominant, **negative** bidoublet ΔT is omitted unless `t_sign=-1`/override is set. Honestly
   flagged in code (`TOP_PARTNER_LOOP_FORMULA_SET="...full_Teq_Zbbeq_not_reconstructed"`).
   **To close:** reconstruct the bidoublet ΔT (CPSW notes it is "somewhat complicated,"
   not given in closed form in the main text) or confirm the override path is exercised in
   production rather than left at singlet-only $+$.

---

## 4. Bookkeeping fixes (cosmetic, no physics change)

- **Add in-code citation** to `custodial_rs_plr_t_coefficient`: "CGHNP 0807.4937 Eq.(153),
  custodial SU(2)_R, after Agashe–Delgado–May–Sundrum hep-ph/0308036." (Currently the
  custodial T has *no* in-code reference, though it is fully sourced.)
- **Note the leading-$L$ truncation** of CGHNP (147) on `minimal_rs_t_coefficient`.
- **Fix the `c_S` YAML anchor** in `EW001.yaml` (status `UNRESOLVED`): re-point `anchor_string`
  to a substring that appears in the snapshot (e.g. "30 v"), or document as prose-sourced.
- **Independently verify the CGHNP convention dictionary** ($c\to-c$, $F^2\to2f^2$) against
  CGHNP App. A profile definitions and `warpConfig/wavefuncs.py`, rather than the internal
  audit note. Flag that the repo's UV-localized light singlets ($c>1/2$) sit outside CGHNP's
  stated $-1/2<c<1/2$ expansion window.

---

## 5. Bottom line

**Every closed-form equation in our custodial implementation is faithfully reproduced from a
local, peer-reviewed source PDF** — oblique T (CGHNP 147/153 ← ADMS), oblique S (PDG-2025),
Z→bb tree (CGHNP 170), P_LR zeroing (ACDP 3–12), and the top-partner loops (CPSW 13–15/35).
The custodial T coefficient that looked uncited **is** sourced (CGHNP 153); it only lacks an
in-code comment.

**The implementation's genuine freedom lives in four proxy choices, not in the equations:**
(1) the $\kappa_b/L$ custodial $g_L^b$ residual (unsourced knob), (2) the $\xi/\rho$
top-partner mass inputs (re-parametrization that breaks the wavefunction correlation),
(3) the $g_R^b=0$ embedding (valid only for $(1,1)_{-1/3}$, not pinned), and (4) the
singlet-only ΔT truncation. These four are exactly the targets of the **full custodial
implementation** (roadmap D); see also the Phase-1 literature catalog for how each is
treated in CGHNP / ACDP / CPSW / Blanke / Archer.
