# Report 11 — Catalog collider_rs / top_higgs_ew / charged_lepton / edm_neutrino

**Structural summary:** CR001=m(g_KK→tt), CR002=m(T_5/3 pair), CR003=m(T pair), CR004=m(B pair), CR005=m(γ/Z_KK→ll), CR006=m(W_KK→ℓν,tb), CR007=m(G_KK→WW/ZZ), CR008=m(T singlet pair), CR009=Λ(llqq contact), CR010=m((T,B) doublet), CR011=σ_fid(VBS W_L W_L), CR012=m(V_KK→VV), CR013=m(G_KK→γγ), CR014=m(top-philic Z′→4t) — all mass-limit proxies. EW001=S,T,U ellipse; EW002/003=CKM SOFT; T001-T008=top FCNC; T010/T011=Z→bb; T012=Z→cc; T014=Z→bs/bd/sd; T015-T017=Z LFV; T018-T020=h LFV. L001=μ→eγ, L002=μ→3e, L003/4/5=μN→eN, L006=muonium, L007/8=τ→ℓγ, L009/10=τ→3ℓ, L023=ν trident. E001=d_e, E002=d_μ, E004=d_n (INFO), E006=d_Hg, E007=d_Ra/Xe, E008=quark cEDM, E009=Weinberg.

### [CRITICAL] EW001 (S,T) anchors do not match any published PDG U=0 fit despite "PDG 2025 Table 10.8" label
- **File:** flavor_catalog/processes/top_higgs_ew/EW001.yaml:82-118
- **Category:** stale-data
- **Claim:** Active HARD ellipse S=+0.026±0.075, T=+0.047±0.066, ρ=0.90 labelled PDG 2025; published PDG U-fixed fits: S=−0.01±0.07, T=+0.04±0.06, ρ=0.92 (2022); S=−0.05±0.07, T=0.00±0.06 (2024) — matches no known edition (stale/fabricated).
- **Evidence:** Numerical impact (re-solving χ²=5.99 with cS=30, cT=428.8): floor = 15.96 TeV (yaml) vs 15.28 (PDG2022) vs 16.44 TeV (PDG2024) — provenance wrong, floor shift ±5%; positive T central mildly anti-conservative vs PDG 2024.
- **Confidence:** high (mismatch), medium (which edition right)
- **Fix:** Re-snapshot the actual PDG U=0 triplet and re-derive the floor.

### [MAJOR] μ→eγ KK-mass convention split: catalog uses physical M_KK (2.45 Λ_IR), legacy uses M_KK=Λ_IR, same prefactor — ×2.45⁴ ≈ 36 in BR
- **File:** flavor_catalog_constraints/rs_ew_builder.py:131; flavorConstraints/muToEGamma.py:31-40; scanParams/scan.py:226-227
- **Category:** convention-inconsistency
- **Claim:** BR = 4e-8·|(ȲȲ†)₁₂|²·(3 TeV/M_KK)⁴ evaluated with physical first gauge-KK mass in the catalog LMFV path but with M_KK=Λ_IR in the legacy scan — same Perez–Randall prefactor, calibrated in only one convention.
- **Confidence:** high (split exists), medium (which side matches PR normalization)
- **Fix:** Pin one M_KK convention for the 4e-8 prefactor and convert the other explicitly.

### [MAJOR] CR005/CR006: SSM Z′/W′ benchmark mass limits applied as HARD vetoes to RS KK EW bosons with √L-suppressed couplings
- **File:** CR005.py:88,392 (CMS 2021 SSM-Z′ 5.15 TeV), CR006.py (W′ SSM 6.0/5.6 TeV)
- **Category:** wrong-formula (benchmark-model mismatch)
- **Claim:** RS KK photon/Z couple to light valence quarks with g≈−g_SM/√L (≈1/5.9), so σ×BR is O(10²–10³) below SSM at fixed mass; HARD-vetoing m<5.15 TeV over-excludes by multiple TeV vs genuine recast (~2–3 TeV reach). Masked in production by the S,T floor, but standalone catalog use misreports direct-search reach.
- **Confidence:** high physics, medium impact
- **Fix:** Downgrade SSM proxies to SOFT/INFO or rescale.

### [MAJOR] EW001 existence floor implemented ≈16.0 TeV, documented as "~18–20 TeV"
- **File:** quarkConstraints/oblique_stu.py:192-260 + docs claims
- **Category:** numerics / convention-inconsistency
- **Claim:** With shipped anchors (c_S=30, c_T=428.8 at L=35, χ²₂,95=5.99) the minimal-RS M_KK floor solves to 15.96 TeV (Λ_IR ≈ 6.5 TeV) — the advertised "existence ~18–20 TeV" is not what the code produces. Custodial lane floor = 5.80 TeV.
- **Confidence:** high (computation), medium (docs may intend different L/CL)
- **Fix:** Recompute/requote the documented floor from shipped anchors or document the parameter set yielding 18–20.

### [MINOR] `build_from_quark_couplings` exports bookkeeping M_KK (=Λ_IR default) as `kk_gluon_mass_gev`
- **File:** point_builder.py:99-103; quarkConstraints/couplings.py:8-11 (xi_KK=1.0)
- **Category:** convention-inconsistency (latent footgun)
- **Claim:** Default xi_KK=1.0 caller gets Λ_IR compared against physical g_KK exclusions (factor 2.45); production avoids only via DEFAULT_XI_KK=2.4487 in the scan script.
- **Confidence:** high
- **Fix:** Require/verify physical-mass flag before exporting.

### [MINOR] EW001 "Gfitter, year: 2026" anchors are the 2018 Gfitter values (diagnostics-only). Fix: tag 2018.

### [MINOR] T010/T011 HARD gate uses an undefined-CL "loose-edge" budget (|exp−SM_fit| + 1σ_combined)
- **File:** T010.py:368-385,693-701; T011.py:569
- **Category:** statistics
- **Claim:** Direction-blind, no stated CL, unlike EW001/CR 95% budgets; improvement shifts >budget are vetoed, worsening up to ~2σ passes. B1 fix itself sound; anchors R_b=0.21629±0.00066, A_FB=0.0992±0.0016, A_b=0.923±0.020 correct; zpole formulas verified.
- **Confidence:** high (description), medium (impact on ~5 TeV floor)
- **Fix:** Signed 95% CL gate or document effective CL.

### [MINOR] EW001 volume log inconsistent: DEFAULT_RS_VOLUME_LOG=35.0 fixed vs actual ln(1.22e19/3000)=35.9, Λ-independent though T∝L (~±1.5%). Fix: compute L from point.

### [NOTE] Collider σ×BR path exists but never populated (dead diagnostic); no Γ/M validity guard, but CR001's active limit is the RS g_KK benchmark itself (broad) so internally consistent there; concern applies to narrow-width benchmarks (CR005 SSM, CR013 k/M̄_Pl=0.1) on broad RS states. Graviton root handling CR007/CR013 (x_G1=3.8317, divide by 2.4487 first) verified correct.

### [NOTE] CR "CMS2026 / PDGLive2026" limits unverifiable and beyond any published Run-2 value
- CR001.yaml: g_KK tt "CMS-B2G-25-009" 5.5 TeV vs Run-2 4.55 TeV; if projection/extrapolation it must not be a HARD budget. Confidence low-medium (offline).

### [NOTE] L- and E-series experimental limits verified current: MEG II 1.5e-13 ✓ (C=1.9365e-3 chain checks out numerically vs PREFAC_BR=4e-8), μ→3e 1.0e-12 ✓, μN→eN Au 7e-13/Ti 6.1e-13 ✓ (Al borrows Au's — placeholder, no direct Al limit exists), τ LFV ✓, d_e=4.1e-30 JILA 2023 ✓, d_n 1.8e-26 ✓, d_Hg 7.4e-30 ✓, ħc conversion ✓. τ→3e 2.5e-8 (vs Hayasaka 2.7e-8) unverified offline.

### [NOTE] Coverage gaps: no Σm_ν/0νββ constraint in edm_neutrino (E003/E005 absent); no (g−2)_μ constraint anywhere; hadronic EDMs (E004 etc.) are INFO bookkeeping with no RS prediction path — no RS EDM floor actually exists in the catalog despite the family name.

### [NOTE] Mixed CLs across HARD constraints (90% LFV/EDM, 95% collider/EW001, ~1σ T010/T011) with uniform ratio≤1 veto — "N passed" has no common statistical meaning.

## Cross-module consistency
- oblique_stu: catalog adapter pure re-export — consistent; x₁² conversion verified numerically; c_S=30 calibration reproduces "M_KK ≥ 3.2 TeV" PDG context exactly.
- zpole: T010/T011/T012 single source of truth, formulas verified — no disagreement.
- collider_resonance: consistent wrapping; fault line is M_KK meaning (xi_KK=1.0 default vs physical) — production papers over via DEFAULT_XI_KK; lepton LFV path genuinely disagrees legacy-vs-catalog ×2.45⁴.
- top_fcnc: partial-width formulas verified against standard results (Γ(t→qZ), Γ(t→qg)=⅔α_s m_t|λ|², Γ(t→qh), Γ_t=1.41 GeV) — no factor-2 found HERE (nb: the sibling rare-decay audit found the ×2 in the *photon/gluon dipole* widths in the same file; this agent verified the Z/h/total widths).
- Trident L023: internally consistent normalization, correct.

## Not reviewed (line-by-line)
CR002–CR004, CR008–CR012, CR014 bodies; T003–T008, T012, T014–T020 bodies; L002–L010, L023 bodies; E002, E006–E009 bodies; adapters mu_e_conversion.py, lfv_three_body*, muonium_conversion.py, zpole_lfv*, higgs_lfv.py, vbs_longitudinal.py, weinberg_operator.py, quark_cedm.py internals; tests. All YAML anchors WERE checked; formula-level bugs inside those bodies may remain.
