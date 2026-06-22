# Slice 3 audit — RS-EW / Zbb / custodial / oblique

Auditor: Claude (skeptical theoretical-physics code review), 2026-06-10.
Scope: `quarkConstraints/rs_ew_couplings.py`, `quarkConstraints/rs_ew_spectrum.py` (as needed),
`quarkConstraints/oblique_stu.py`, `quarkConstraints/zpole.py`, T010/T011/EW001 (primary),
T014 (secondary — note: prompt said `primary/`, the file actually lives at
`flavor_catalog_constraints/secondary/top_higgs_ew/T014.py`), associated tests, and the
CUSTODIAL-RESEARCH / ZBB orchestration records.

Method: independent re-derivation of the gauge-KK overlap sum via the zero-mode-subtracted
p=0 5D Green's function; direct comparison against the **published CGHNP equations**
(Casagrande–Goertz–Haisch–Neubert–Pfoh, arXiv:0807.4937, full text fetched from ar5iv during
this audit); numerical reproduction of scan-representative points with
`fit_quark_sector(default_quark_targets(), r=0.25, Λ_IR=8168 GeV)` (matches the
ZBB-XCHECK benchmark to ~2%); floor re-computation through the live T010 anchors.

Headline under test: minimal floor 25–30 TeV (physical M_KK) driven by T010/Zbb;
custodial strict floor 2–3 TeV.

**Verdict: the minimal 25–30 TeV floor is NOT defensible as currently computed.** It is
produced by (i) a mistranslated Casagrande fermion-KK admixture whose dominant terms have
the wrong sign and are ~10²–10³ too large at scan points (BLOCKER), interacting with
(ii) a T010 gate that leaves essentially zero headroom on one side because the SM itself
sits at −0.996σ in R_b against a 1.0σ two-sided cut (MAJOR). The custodial-vs-minimal
*qualitative* contrast survives; the custodial branch zeroes the buggy term anyway.

---

## Finding 1 — BLOCKER: Casagrande m_b² fermion-KK admixture is mistranslated from CGHNP (wrong c/F dictionary on the diagonal term; wrong F index in the flavor sum)

**Files/lines:**
- `quarkConstraints/rs_ew_couplings.py:1835-1847` (`_casagrande_zbb_B_profile`):
  ```python
  value = (1.0 / denom_left) * (1.0 / (F * F) - 1.0 + (F * F) / denom_right)
  # denom_left = 1 - 2c ; denom_right = 3 + 2c ; F = repo f_IR
  ```
- `quarkConstraints/rs_ew_couplings.py:957-971` (`build_rs_zbb_fermion_kk_mixing`):
  ```python
  B_d = float(np.dot(row_ratio, profile_b_d))   # row_ratio_i * B(c_d[i], F_d[i])
  B_Q = float(np.dot(column_ratio, profile_b_q))
  delta_g_L_b = +prefactor * B_d ;  delta_g_R_b = -prefactor * B_Q
  ```
- Spec that baked this in: `.orchestration/runs/W2-P6/impl6a_prompt.md:6-9`
  ("g_L uses sum_i |Y_d[2,i]|²/|Y_d[2,2]|² with c_d[i]").

**Published source (fetched this audit, CGHNP 0807.4937, Z→bb̄ ZMA equations):**

```
g_L^b ⊃ +(m_b²/2M_KK²) [ 1/(1−2c_bR) (1/F²(c_bR) − 1 + F²(c_bR)/(3+2c_bR))
                          + Σ_{i=1,2} |Y_d,3i|²/|Y_d,33|² · 1/(1−2c_d_i) · 1/F²(c_bR) ]
g_R^b ⊃ −(m_b²/2M_KK²) [ 1/(1−2c_bL) (1/F²(c_bL) − 1 + F²(c_bL)/(3+2c_bL))
                          + Σ_{i=1,2} |Y_d,i3|²/|Y_d,33|² · 1/(1−2c_Q_i) · 1/F²(c_bL) ]
```

Two independent errors in the repo implementation:

**(1a) Flavor-sum F index.** CGHNP's light-generation sum carries `1/F²(c_bR)` — the
**third-generation** singlet overlap — times the O(1) tower factor `1/(1−2c_d_i)`. The repo
instead applies the full `B(c_d_i, F_d_i)` per light generation, i.e. `1/F²(c_d_i)`.
For UV-localized light down singlets (scan-typical `c_d = [0.649, 0.601, 0.557]`,
`F_d = [0.0021, 0.0094, 0.0333]`) this inflates the i=1,2 terms by `F²(c_d3)/F²(c_d_i)`
≈ 10²–10⁴ per term. Physically: the b_L–KK(d_i) mass mixing is
`∝ v Y^{5D}_{3i} F_Q3 × O(1)` (KK modes are always IR-localized); only the conversion of
`Y²` into `m_b²` introduces a single `1/F²(c_b_R)`.

**(1b) Convention dictionary on the bracket itself.** CGHNP define
`F²(c) = (1+2c)/(1−ε^{1+2c})` with `c_{Q,q} ≡ ±M/k`, so the exact dictionary to repo
variables is `c_CGHNP = −c_repo`, `F²_CGHNP = 2 f_IR,repo²` (proved exactly:
`(1−2c_r)/(1−ε^{1−2c_r}) = 2·(1/2−c_r)/(1−ε^{1−2c_r})`; cross-validated by the gauge piece,
see "verified correct" item V1, where the same dictionary reproduces CGHNP's
`F²(c)/(3+2c)·(L−(5+2c)/(2(3+2c)))` term-by-term, using the identity
`(5+2c)/(2(3+2c)) = 1/2 + 1/(3+2c)`). The correct bracket in repo variables is

```
B_correct(c_r, f_r) = 1/(1+2c_r) · ( 1/(2 f_r²) − 1 + 2 f_r²/(3−2c_r) )
```

The repo coded `1/(1−2c_r)·(1/f_r² − 1 + f_r²/(3+2c_r))` — every factor mistranslated
(c-sign in both denominators, and missing factor 2 in F²). Consequences:
- For `c_r(b_R) > 1/2` (UV-localized b_R, the scan-typical case `c_d3 ≈ 0.56`), the coded
  `1/(1−2c_r)` is **negative**, flipping the sign of δg_L^b relative to CGHNP (whose
  `1/(1−2c_bR) = 1/(1+2c_r) > 0` always in the physical range).
- Magnitudes are off by O(10–30) even on the diagonal term.

**Numerical impact (representative scan point, r=0.25, Λ_IR=8.168 TeV ⇒ M_KK^phys=20 TeV;
fit: c_Q=[0.613,0.564,0.392], c_d=[0.649,0.601,0.557], m_b=2.732 GeV, |Y_d_bulk,33|=0.71,
off-diagonals 0.12–0.35; repo δg matches the recorded ZBB-XCHECK benchmark −2.478e-3 to 2%):**

| quantity | repo | CGHNP-correct | ratio |
|---|---|---|---|
| B_d (drives δg_L^b) | −4.50e4 | +2.36e2 | −191 |
| B_Q (drives δg_R^b) | −5.38e3 | +2.63 | −2042 |
| δg_L^b (fermion piece) | −2.52e-3 | +1.32e-5 | wrong sign, 190× |
| δg_R^b (fermion piece) | +3.01e-4 | −1.47e-7 | wrong sign, 2000× |
| gauge piece δg_L^b (correct, V1) | +8.74e-5 | same | — |

With the correct formula the fermion piece is a ~15% correction to the gauge piece; in the
repo it is ~30× the gauge piece **with opposite sign**. The PHASE2 ledger explicitly states
the 25–30 TeV floor is "Dominated by the m_b² FULL-FLAVOR-SUM Casagrande admixture (NOT the
gauge piece)" (`.orchestration/PHASE2_PROGRAM_LEDGER.md:98`) — i.e., the headline is
dominated by the mistranslated term.

**Floor impact** (re-evaluated through live T010 anchors, scaling couplings as 1/M_KK²):
repo couplings give a T010 crossing at ≈24 TeV at this point (consistent with the quoted
23.5–26.6 TeV); CGHNP-correct couplings give a completely different picture — see Finding 2,
because the corrected floor is then controlled by the gate pathology, not by physics.
Mechanically: the repo's wrong-sign δg_L^b<0 raises R_b; the "floor" at ~24–30 TeV is the
scale where that spurious upward shift re-enters the +1σ window of an R_b measurement the
SM already undershoots by 1σ. The 25–30 TeV headline is an artifact of this bug.

**Process note.** The ZBB-XCHECK reviewer (`.orchestration/runs/ZBB-XCHECK/author_codex_out.md`)
detected the symptom ("the aggressive piece is the m_b² fermion-KK B_d term... code
B_d=−4.43e4 vs compact bracket −9.78e3"; "suspected fix is the B_d/B_Q light-index
bracket/normalization") but the program recorded verdict "CONVENTION-CHOICE". The actual
CGHNP equations (fetched here) settle it: it is a mistake, not a convention. The gate tests
(`tests/test_rs_ew_phase6a_zbb_fermion_mixing.py:135-156`) re-implement the same expression
(`_manual_B` is byte-identical algebra), so they pin the bug rather than check it against
the literature; the pinned values at the test point (c_d3=0.20, near-diagonal Y_d with
off-diagonals ≤4e-5) sit in the one corner where both errors are least visible (sign
correct there; magnitude still ~6× high).

Severity: **BLOCKER**. Confidence: **high** (primary-source equations + exact convention
dictionary cross-validated on the independently verified gauge piece + numerical
reproduction of the scan benchmark).

---

## Finding 2 — MAJOR: T010 veto semantics has ~0.004σ headroom on the low-R_b side (SM sits at −0.996σ of a 1.0σ two-sided gate), so the floor is not robust and the corrected-physics floor becomes a gate artifact

**Files/lines:** `flavor_catalog_constraints/primary/top_higgs_ew/T010.py:653-655`
(`passes = max|pull| ≤ 1.0` with `pull = (predicted − experimental)/σ_combined`), anchors
from `flavor_catalog/processes/top_higgs_ew/T010.yaml` + LEP/SLC snapshot.

Live numbers (this audit): R_b^exp = 0.21629 ± 0.00066; SM-limit prediction (radiator
calibrated to the snapshot SM-fit value) = 0.21562; σ_combined = 6.73e-4 ⇒ **SM pull
= −0.996**. A_b SM pull = +0.627. So the SM passes T010 by 0.004σ, and *any* additional
negative R_b shift of ≳3e-6 absolute (δg_L^b ≳ +4e-6, the standard RS direction) fails the
gate at arbitrarily high M_KK until the shift decays away.

Consequences, quantified through the live anchors at the representative point of Finding 1:
- CGHNP-correct couplings (gauge-dominated, δg_L^b>0 ⇒ R_b decreases): T010 1.0σ "floor"
  = **108 TeV** — obviously not a physics exclusion, just the M_KK at which the NP pull
  decays below the 0.004σ headroom.
- Same couplings with a 2σ gate: floor = **6.8 TeV**.
- Repo (buggy) couplings with a 2σ gate: 19.6 TeV (vs 24.0 at 1σ).

A factor-16 swing in the corrected floor between 1σ and 2σ gating means T010's floor as
currently defined measures the anchor-snapshot tension, not RS physics. Note the snapshot
SM-fit value used (0.21562, old LEP-EWWG-era fit) maximizes the tension; PDG-2024-style SM
fits (R_b ≈ 0.21582) would put the SM pull near −0.7 and partially defuse this, but the
structural problem (one-sided ~0 headroom under a two-sided max-pull-vs-experiment gate)
remains. T011 deliberately uses a "loose edge" budget `|exp−SM| + σ` for exactly this
reason (its docstring says so); T010 does not.

This also makes the custodial branch brittle: custodial T010 passes only because δg is
*exactly* zeroed (pull −0.996 ≤ 1.0 by 0.004); any κ_b/L residual > ~4e-6 in δg_L^b would
flip a custodial point to FAIL at any M_KK.

Recommended fix direction (for the maintainers to decide): gate on the NP shift against an
uncertainty-aware budget (as T011 does), or use a 2σ/95% convention, or include the SM pull
subtraction (Δχ² vs SM). Any of these changes the minimal floor materially once Finding 1
is fixed.

Severity: **MAJOR** (combined with Finding 1 it is what actually sets the headline number).
Confidence: **high** (direct evaluation of the live constraint).

---

## Finding 3 — MAJOR: EW001 minimal ΔT uses the CGHNP geometric-M_KK coefficient but evaluates it with the physical first KK mass ⇒ ΔT underestimated by x₁² ≈ 6.0 (anti-conservative), inconsistent with the S calibration inside the same constraint

**Files/lines:** `quarkConstraints/oblique_stu.py:158-169` (`minimal_rs_t_coefficient` =
`πL/(2c_W²)`), `oblique_stu.py:233-237` (`scale=(v/m_kk)²` with `m_kk` =
`kk_ew_mass_gev`), `EW001.py:338-346` (`_resolve_m_kk_gev` prefers the `kk_ew_mass_gev`
extra, which is the **physical** first gauge KK mass `x₁·Λ_IR`, `rs_ew_spectrum.py:663-667`).

CGHNP (fetched): `T = πv²/(2c_W²M_KK²)·(L − 1/(2L))` and custodial
`T = −πv²/(4c_W²M_KK²L)` with **M_KK ≡ kε = Λ_IR (geometric)**. The repo coefficients match
CGHNP exactly (the dropped −1/(2L) is negligible), but the denominator scale used is the
physical mass 2.45·Λ_IR, so ΔT_minimal is a factor x₁² ≈ 6.005 too small. Meanwhile the
ΔS coefficient c_S = 30 (PDG anchor, `EW001.yaml:156-160`) is calibrated in the *physical*
convention (CORRECTED_PRESCRIPTION.md:28-30 makes this explicit: 2π·x₁² ≈ 36.6 ≈ 30). So
EW001 mixes conventions between its own S and T terms.

Floor impact: minimal-RS EW001 floor moves from ≈6–7 TeV (as coded) to ≈15–17 TeV
(physical) with the convention applied consistently — still below the (buggy) 24–30 TeV
T010 floor, but **above** the corrected gauge-only T010 2σ floor (~7 TeV). In a fully
corrected stack the minimal-model EW floor would plausibly be set by T/Zbb at O(15) TeV
physical, not 25–30. Custodial EW001 is essentially unaffected (T term is −π/(4c_W²L),
tiny in either convention; S dominates and is calibrated in the physical convention), so
the custodial inclusive ~7 TeV statement is not moved by this finding.

Mitigation: EW001 is explicitly tagged a proxy (`OBLIQUE_STU_RS_PROXY_V1`,
NEEDS-HUMAN-PHYSICS) and is excluded from the strict floor lane
(`docs/STATE_OF_PROJECT.md:142`). Documented-proxy status acknowledged — but the
*direction* of the error is anti-conservative and the convention inconsistency inside one
formula is an actual mistake, not an approximation choice.

Severity: **MAJOR** (anti-conservative ×6 in T; does not move current headline because
EW001 is proxy-tagged). Confidence: **high** (primary-source formula + code path).

---

## Finding 4 — MAJOR (claim/process): "custodial strict floor 2–3 TeV" is claimed while the dual-approved custodial prescription forbids that claim without top-partner loop numerics, which are deferred in production

**Files/lines:** `docs/STATE_OF_PROJECT.md:5` and `.orchestration/PHASE2_PROGRAM_LEDGER.md:75`
(claims) vs `.orchestration/runs/CUSTODIAL-RESEARCH/CORRECTED_PRESCRIPTION.md:39-42`:
"Add the FLAG now (required). Can defer numerics for first tree-level validation, but do
NOT claim 2-3 TeV viability without it. At M_KK~2-3 TeV: ΔT_loop ~ O(0.1), δg_L^b|loop ~
1e-3 (relevant for Zbb)".

The production scan path never passes `include_top_partner_loops=True` /
`top_partner_loop_components` (scan harness calls `build_from_rs_ew_inputs` with defaults:
loops deferred; verified in `scripts/run_full_catalog_scan.py` call sites and
`test_t010_t011_custodial_tree_is_active_not_deferred_and_stays_rigorous` asserting
`top_partner_loop_numerics_included is False` while the tag is "rigorous"). At the claimed
2–3 TeV floor, the prescription's own estimates (δg_L^b|loop ~ 1e-3, ΔT_loop ~ 0.1) are
*not* small: δg_L^b ~ 1e-3 is ~10× the gauge piece at 20 TeV scaled to 3 TeV would be
~4e-3 — i.e. loop terms are leading at 2–3 TeV. Tagging custodial T010 results "rigorous"
at 2–3 TeV with the loop sector deferred contradicts the dual-approved spec. (The ΔF=2
origin of the 2–3 TeV strict floor is outside this slice, but the "custodial removes the
Zbb veto down to 2–3 TeV" narrative leans on Zbb being loop-free there.)

Severity: **MAJOR** (claim integrity / tagging, not a formula bug). Confidence: high on
the documented contradiction; the loop-magnitude estimates are the prescription's own.

---

## Finding 5 — MINOR: custodial `bR_strategy="elementary_zero"` zeroes δg_R^b exactly; P_LR does not protect b_R, and "elementary" b_R still has the universal/non-universal gauge shift

`rs_ew_couplings.py:1555-1556` sets `right[2,2]=0` in the custodial branch. P_LR protects
only δg_L^b (correctly documented). For an elementary (UV-localized) b_R the minimal gauge
shift is small but nonzero (a(c)−a_ref ≈ O(1/L) residuals; at 3 TeV custodial points this
is O(1e-5–1e-4) on g_R^b where the A_b budget is 0.02) — numerically negligible, and the
choice matches CORRECTED_PRESCRIPTION.md:35 ("δg_R^b ≈ 0 default"). Documented proxy, not
a mistake; flagged because exact zeroing slightly *understates* custodial effects in A_b
and bakes in "no A_FB^b anomaly explanation" by construction (also documented).

Severity: MINOR. Confidence: high.

---

## Finding 6 — MINOR: T010/T011 metadata copy says custodial protection covers "diagonal down-left" only, while true all-gen bidoublet P_LR would also protect the off-diagonals; `pr1_minimal_offdiag` keeps minimal off-diagonals (conservative for T014/FCNC)

For all-generation `(2,2)_{2/3}` embeddings the P_LR argument zeroes the *entire* left
down-sector gauge-KK coupling matrix (the protection is per-field, basis-independent), so
keeping the minimal off-diagonals (`CUSTODIAL_FCNC_PR1_MINIMAL_OFFDIAG`,
`rs_ew_couplings.py:1558-1600`; test
`test_t014_down_fcnc_offdiagonals_and_result_are_unchanged_in_custodial_branch`) makes the
custodial branch *over-constrained* in Z-FCNC (T014, and Z contributions to ΔF=1) — i.e.
conservative, and apples-to-apples for the Zbb-driven minimal-vs-custodial M_KK comparison
(T014 identical in both branches by construction). The alternative
`all_gen_bidoublet_mass_basis_proxy` mode with κ_fcnc/L residual exists and is the more
faithful representation-aware choice. Documented; consistent with PR1 scoping. No floor
impact (T014 at 2.9e-3 BR limits is far from binding: scan-typical off-diagonal
δg ~ 1e-6–1e-5 gives BR ~ 1e-11–1e-9).

Severity: MINOR (documented, conservative direction). Confidence: high.

---

## Finding 7 — MINOR: assorted small items

- **m_b input to the admixture**: `masses_down[2] = 2.73 GeV` (fit-scale running mass);
  CGHNP use m_b(M_KK) ≈ 2.4 GeV ⇒ ~+29% on the m_b² term. Irrelevant beside Finding 1 but
  worth fixing together. (`rs_ew_couplings.py:944-947`)
- **EW001 volume log**: fixed `DEFAULT_RS_VOLUME_LOG = 35.0` rather than the point's
  L = ln(1/ε) = 34.54 (`oblique_stu.py:47`) — 1.3% on ΔT; cosmetic.
- **a_ref scheme**: the EW-universal subtraction `a_ref = a(c=0.65)` is a documented
  input-scheme choice (`metadata["a_ref_interpretation"]`); fine, but the residual scheme
  dependence for fermions at c ≠ 0.65 in the 0.55–0.8 range is O(few %) of the universal
  piece (a(c) varies by ~1% there, verified numerically) — negligible.
- **T010 prompt path mismatch**: T014 is SECONDARY (`flavor_catalog_constraints/secondary/...`),
  not primary as the audit prompt assumed. Bookkeeping only.
- **Carena top-partner loop proxy** (`rs_ew_couplings.py:1056-1347`): structure resembles
  the leading Carena–Pontón–Santiago–Wagner terms and is honestly labeled
  (`full_Teq_Zbbeq_not_reconstructed`); default-off and never exercised in production
  scans, so it was not audited at depth here. Its (16π s_w² m_W²) normalizations should be
  re-derived against hep-ph/0701055 Eqs. (28)–(30) before it is ever turned on.

---

## Verified correct

- **V1 — Gauge-KK Z-coupling shift (the `_z_delta`/`a(c)` machinery) is exactly right.**
  Independent check: I derived the zero-mode-subtracted p=0 gauge Green's function in
  closed form, G(t,1) = t²log t/2 + (L−½)t²/2 + (1−L)/(4L), giving
  S(c) = F²_CGHNP(c)·[(L−½)/(2(3−2c)) − 1/(2(3−2c)²)] + (1−L)/(4L) in repo variables.
  The spectrum tower sum `a(c) = Σ (x₁²/x_n²)χ_n(1)Ω_n(c)` matches this to ≤1e-4 relative
  across c ∈ [−0.2, 0.8] (numerics run this audit), and the closed form reproduces CGHNP's
  published `F²(c)/(3+2c)(L−(5+2c)/(2(3+2c)))` exactly. So: L-enhancement present; the
  m_Z²/M_KK² prefactor times the (M_KK²/m_n²)-weighted sum gives exactly Σ m_Z²/m_n² —
  **no Λ_IR-vs-physical-M_KK ambiguity in the gauge piece** (the x₁² cancels by
  construction); G(1,1) = (L−1)/2 + 1/(4L) matches the known propagator value.
- Sign/direction: `s_z=−1` ⇒ δg = −g_SM·(m_Z²/m_n²-weighted overlap); for IR-localized b_L
  this gives δg_L^b > 0 (|g_L| reduced, R_b reduced) — the standard RS Zbb problem
  direction.
- EW quantum numbers: `_sm_chiral_z_coupling` gives g_L^b = −½ + s_w²/3, g_R^b = +s_w²/3,
  and the zpole module uses the same convention with sin²θ_eff = 0.2315; b_L T³=−½,
  Q=−1/3 correct everywhere checked.
- Use of geometric Λ_IR in the m_b²/(2Λ_IR²) admixture prefactor is the **correct** CGHNP
  convention (M_KK ≡ kε there); explicitly documented in metadata
  (`M_KK_convention: "geometric Lambda_IR"`). The prefactor and the L↔R B_d/B_Q
  cross-assignment (g_L ← singlet tower, g_R ← doublet tower) match the published
  equations. Only the bracket translation and flavor-sum F-index are wrong (Finding 1).
- m_b factorization dictionary verified: 2·v₁₇₄·|Ȳ_33|·f_Q3·f_d3 = 2.718 ≈ masses_down[2]
  = 2.732 at the representative fit, confirming Ȳ_repo = Y_CGHNP and F_C = √2·f_repo.
- T010/T011 plumbing: reads z_delta_g_L/R_d[2,2] including the admixture when enabled; the
  radiator calibration reproduces the snapshot SM-fit R_b exactly; A_FB = ¾A_eA_b; R_b
  computed as width-weight ratio over u,d,s,c,b with N_c=3. Two-sided max-pull semantics
  as documented (but see Finding 2 for the consequence). T011's loose-edge budget
  `|exp−SM| + σ_comb` on the NP *shift* is sensible and SM-safe.
- T014 width arithmetic: Γ(Z→q_iq̄_j + q_jq̄_i) with charge-state factor 2, N_c=3,
  chirality sum |δg_L|²+|δg_R|², averaged pair radiator, total-width denominator including
  3 ν and 3 charged leptons; massless phase space is fine at the 2.9e-3 BR limit level.
  Pure-NP comparison against direct nonstandard-hadronic-width limits is appropriate.
- Custodial PR1 mechanics match the dual-approved CORRECTED_PRESCRIPTION: all-gen
  down-left diagonal zeroing (gauge *and* fermion-admixture pieces, per ACDRP
  hep-ph/0605341 Eq. (18) T³_R−T³_L argument), κ_b/L residual machinery, minimal_rs branch
  byte-identical when custodial is off (tested), `pr1_minimal_offdiag` keeps T014
  identical across branches (tested).
- Custodial oblique T coefficient −π/(4c_W²L) matches CGHNP's published custodial formula
  exactly (modulo the convention issue of Finding 3, which is subdominant for the
  custodial branch); EW001's correlated-Gaussian (S,T) ellipse algebra
  (inverse covariance, χ²_2dof 95% = 5.99) is correct; PDG 2025 anchors
  (S=0.026±0.075, T=0.047±0.066, ρ=0.90, U fixed 0) load correctly.
- Hermiticity enforcement, mass-basis rotations U†diag(a)U, and the product-minus-SM
  neutral contact tensors are implemented correctly.

---

## Bottom line on the headline claims

1. **Minimal 25–30 TeV floor: REJECTED as currently derived.** It rests on a fermion-KK
   admixture term with the wrong sign and 10²–10³× the correct magnitude (Finding 1,
   primary-source confirmed), filtered through a gate with 0.004σ of SM headroom
   (Finding 2). With CGHNP-correct couplings the T010 floor is either ~108 TeV (current 1σ
   gate — an artifact) or ~7 TeV (2σ gate) at the representative scan configuration;
   a consistent-convention EW001 would sit at ~15–17 TeV but is proxy-tagged. The correct
   minimal Zbb physics (gauge-dominated, δg_L^b ≈ +1.0e-4·(20 TeV/M_KK)²) is ~25× weaker
   than what was scanned.
2. **Custodial strict floor 2–3 TeV:** the Zbb mechanics of the custodial branch are
   internally consistent (and the buggy term is zeroed there), but the 2–3 TeV *claim*
   violates the program's own dual-approved caveat that top-partner loop numerics
   (δg_L^b|loop ~ 1e-3, ΔT ~ 0.1 at 2–3 TeV) must be included first (Finding 4). The
   minimal-vs-custodial qualitative contrast (custodial removes the tree Zbb veto)
   stands.
