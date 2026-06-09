# Custodial branch — CORRECTED prescription (3-agent consensus: 2 codex + Opus, all CHATGPT-CORRECT)

Resolves the contradiction: our earlier Opus panel was WRONG; ChatGPT is right.
Primary source: ACDRP hep-ph/0605341 Eq.(18) — the non-universal gauge-KK vertex
correction is ∝ (T^3_R − T^3_L), which VANISHES for custodial b_L (T^3_L=T^3_R=−1/2).
The P_LR theorem protects the WHOLE non-universal Zb_L vertex correction, INCLUDING
the leading F²(c_Q3)-enhanced gauge-wavefunction piece (which is the RS Zbb problem),
not just the m_b² fermion-mixing admixture. Validated against Carena et al hep-ph/0701055
and Casagrande et al 0807.4937 Eq.(153).

## Implementation spec (custodial branch, ew_model="custodial_rs_plr")

### Representations (metadata)
- Q_L^i ∈ (2,2)_{2/3} for ALL generations i=1,2,3 (not just 3rd) + rotate to mass basis
  — avoids reintroducing unprotected down-left Z couplings via flavor rotations.
- b_L = P_LR eigenstate, T^3_L=T^3_R=−1/2. t_R = (1,1)_{2/3} default. b_R elementary default.

### Z b_L b_L (THE correction): zero BOTH leading tree pieces
- Zero the leading **gauge-profile** shift `z_delta_l_d[2,2]` (the `_z_delta(g_l_d, a_mass["L_d"])`
  entry, the F²(c_Q3)/a(c) piece) — NOT just the m_b² admixture.
- Zero the **m_b² Casagrande fermion-KK admixture** too (protected if bottom Yukawa respects P_LR).
- Apply to all down-type-LEFT entries (all-gen bidoublet), in the mass basis.
- Replace with: δg_L^b ≈ κ_b·(1/L)·δg_L,minimal^b  (κ_b = O(1)), OR zero in the ideal P_LR
  limit + an explicit `custodial_PLR_breaking_residual` flag. PLUS a flagged one-loop
  top-partner term. The κ_b/L term is a model-systematic FLAG, not theorem-mandated.

### Oblique (EW001)
- KEEP S unchanged: repo `c_S` (≈30) is the PHYSICAL first-KK-mass convention; do NOT
  naively replace with the geometric 2π. (2π is the Λ_IR convention; physical m_1≈2.45 Λ_IR
  → ~36 ≈ repo 30.) S becomes the DOMINANT tree oblique constraint after custodial.
- Replace T coefficient: minimal `π L/(2 c_W²)` → custodial `−π/(4 c_W² L)` (the −1/(2L²)
  suppression, sign-flipped). U=0.

### b_R / A_b
- δg_R^b ≈ 0 default (P_LR does NOT protect b_R; minimal b_R is elementary → small).
- A_FB^{0,b} only via an EXPLICIT b_R representation model (positive δg_R^b ~1e-2..2e-2),
  not automatic — do not "solve" the A_FB^b anomaly by hand in the default branch.

### One-loop top-partner (Carena et al Eq.28-30)
- Add the FLAG now (required). Can defer numerics for first tree-level validation, but do
  NOT claim 2-3 TeV viability without it. At M_KK~2-3 TeV: ΔT_loop ~ O(0.1), δg_L^b|loop ~ 1e-3
  (relevant for Zbb); sign model-dependent (bidoublet partners push T negative).

### Conventions
- Keep M_KK convention (physical first-KK m_1≈2.45 Λ_IR vs Λ_IR) EXPLICIT in metadata,
  consistent with the Zbb cross-check finding. All coefficients above must be applied in
  the matching convention.

### Honest omissions (must be flagged in diagnostics)
SU(2)_R gauge KK tower, custodial fermion (custodian) spectrum, exact NC mixing, BKT effects,
and the numerical one-loop top-partner terms (if deferred). The custodial branch is a
representation-aware PROXY, not a full custodial RS-EW sector.
