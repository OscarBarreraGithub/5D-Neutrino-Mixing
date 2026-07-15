# Fact-Check Review: `constraint_formulas.tex`

**Reviewer role:** theoretical-physics fact-checker (citation + physics-soundness audit)
**Date:** 2026-06-17
**Target:** `review_local/constraint_formulas.tex` (Warped Quark Sector Constraints)
**Scope:** every literature citation verified against the actual source (equation numbers, numerical values, attribution), plus physics/internal-consistency checks. **No edits were made to the .tex file** — this is a review only.

**Overall verdict:** The document is **substantially solid**. Of 13 reference items checked, **9 CONFIRMED**, **3 PARTIAL** (real paper, correct physics, but an equation-number or precision slip), **1 INCORRECT** (PDG oblique-fit central values/correlation do not match any current PDG fit). All physics-arithmetic and internal-consistency checks pass. The single most important fix is the **PDG 2024/2025 S,T fit numbers in the oblique section** (Ref. 11), which are wrong and should be replaced.

---

## (a) Reference-by-reference verdicts

| # | Reference (as cited) | Verdict | Note |
|---|----------------------|---------|------|
| 1 | Csaki, Falkowski, Weiler, arXiv:0804.1954 — Sec. 2 Eqs. (2.14)–(2.17), tree H after (2.17), (2.18) chiral/RG, (2.19) C₄ᴷ | **CONFIRMED** | Title/authors and all equation roles verified; colour factors 1/6, −1, 1/3 present. |
| 2 | Blanke, Buras, Duling, Gori, Weiler, arXiv:0809.1073 — Sec. 4.2 Eqs. (4.14)–(4.17), (4.23)–(4.29), (4.36)–(4.38), Eq. (4.38) = M₁₂�q convention | **PARTIAL** | All ranges correct **except**: the `M₁₂�q = (M₁₂�q)_SM·C_Bq·e^{2iφ_Bq}` convention is **Eq. (4.39)**, not (4.38). (4.38) is `ΔM_q = 2|M₁₂_SM+M₁₂_KK|`. Off by one. |
| 3 | Buras, Guadagnoli, Isidori, arXiv:1002.3612 — ε_K master Eq. (34), κ_ε=0.94±0.02 Eq. (35) | **CONFIRMED** | Both equation numbers and κ_ε = 0.94 ± 0.02 verified verbatim. |
| 4 | Buras, Guadagnoli, arXiv:0805.3887 — κ_ε Eqs. (9)–(12), value 0.92 | **CONFIRMED** | κ_ε introduced (9)–(12); κ_ε = 0.92 ± 0.02 confirmed (Eq. (12)/Table 1). |
| 5 | Perez, Randall, arXiv:0805.4652 — dipole Eq. (29), flavour Eq. (28), BR Eq. (30) w/ 4×10⁻⁸, spurion Eq. (31) C≈0.02 | **CONFIRMED** | All four equation numbers, the 4×10⁻⁸ prefactor, and C ≈ 0.02 verified verbatim. |
| 6 | MEG II 2025, arXiv:2504.15711 — BR(μ→eγ) < 1.5×10⁻¹³ @ 90% CL | **CONFIRMED** | MEG II "New limit on the μ⁺→e⁺γ decay"; limit 1.5×10⁻¹³ (90% CL) exact. |
| 7 | Casagrande et al., arXiv:0807.4937 — Sec. 6.1 Eq. (170) admixture, Eq. (171) R_b⁰/A_b/A_FB, values R_b⁰=0.21629±0.00066, A_b=0.923±0.020 | **PARTIAL** | Eq. (170) ✓ (LH Zbb̄ shift w/ m_b²/M_KK² admixture, ZMA). Values exact ✓. **But** Eq. (171) is the *theoretical definitions* of R_b⁰/A_b/A_FB; the *experimental values* R_b⁰=0.21629±0.00066 and A_b=0.923±0.020 are **Eq. (173)** (SM-predicted are (172)). Section is **6.4**, not 6.1. |
| 8 | Agashe, Contino, Da Rold, Pomarol, arXiv:hep-ph/0605341 — custodial Zb_Lb_L | **CONFIRMED** | "A custodial symmetry for Zbb", exactly that framework. |
| 9 | Carena et al., Nucl. Phys. B 759 (2006) 202, arXiv:hep-ph/0607106 — S ≃ 30 v²/M_KK² via PDG | **CONFIRMED** | PDG 2024 EW review states verbatim "S ≈ 30 v² M_KK⁻²" citing this exact Carena et al. paper. |
| 10 | Agashe, Delgado, May, Sundrum, arXiv:hep-ph/0308036 — Sec. 5.1 Eqs. (5.5)–(5.6), T∝L | **CONFIRMED** | (5.5) = S, (5.6) = T with explicit `kπr_c` (warped-volume / log) enhancement factor. |
| 11 | PDG 2025 SM review, Table 10.8: S=0.026±0.075, T=0.047±0.066, ρ_ST=0.90, U fixed | **INCORRECT** | PDG 2024 (U≡0 fit, Eq. 10.99) gives **S = −0.05 ± 0.07, T = 0.00 ± 0.06**, S–T correlation **93%**. U-free (10.98): S=−0.04±0.10, T=0.01±0.12, U=−0.01±0.09. The doc's central values/errors/correlation match **no** current PDG fit (they look like an older Gfitter fit). See fix below. |
| 12 | ε_K numbers: ε_K^exp=2.228×10⁻³ (PDG), ε_K^SM=2.161×10⁻³ (BGS), budget≈6.7×10⁻⁵ | **PARTIAL** | ε_K^exp ✓ ((2.228±0.011)×10⁻³). Budget arithmetic ✓ (6.7×10⁻⁵). **But** Brod–Gorbahn–Stamou's published central value is **2.16×10⁻³** (2.16(6)(8)(15), PRL 125 (2020) 171803); the doc's "2.161" adds an unpublished digit. Value/attribution essentially correct, third digit is over-precise. |
| 13 | Collider table (long table) | **CONFIRMED (with minor labelling note)** | All 14 rows resolve to the claimed search and the mass edges are correct. One labelling nuance on the T_{5/3} row (see below). |

### Collider-table row detail (Ref. 13)

| Row | Claim | Resolves to | Verdict |
|-----|-------|-------------|---------|
| KK-gluon→tt̄, up to 5.5 TeV, CMS-B2G-25-009, arXiv:2603.23454 | 5.5 TeV | CMS "new particles → tt̄", KK-gluon excluded 0.5–5.5 TeV | **CONFIRMED** |
| T_{5/3} exotic top partner, ≳1.46 TeV, ATLAS arXiv:2212.05263 | 1.46 TeV | ATLAS VLQ+MET; X_{5/3}(B/X→Wt) limit 1.47 (1.46) TeV | **CONFIRMED** (number ok; the paper's charge-+5/3 state is "X"; doc's "T_{5/3}" is just nomenclature) |
| charge-2/3 T (pair), ≳1.70 TeV, ATLAS arXiv:2401.17165 | 1.70 TeV | T→Wb 100% excluded < 1700 GeV | **CONFIRMED** |
| singlet VLQ T, ≳1.36 TeV, ATLAS arXiv:2401.17165 (same paper) | 1.36 TeV | singlet (1/2:1/4:1/4) excluded < 1360 GeV — **same paper** | **CONFIRMED** (both 1.70 and 1.36 genuinely from this one paper) |
| charge-1/3 B (pair), ≳1.57 TeV, CMS PRD 110 (2024) 052004 | 1.57 TeV | B→bH 100% excluded < 1570 GeV | **CONFIRMED** (exact) |
| Z'_SSM dilepton, ≳5.15 TeV, CMS JHEP 07 (2021) 208 | 5.15 TeV | Z'_SSM excluded < 5.15 TeV | **CONFIRMED** |
| W'/W_KK, ≳6.0 TeV, ATLAS arXiv:1906.05609 | 6.0 TeV | W'_SSM excluded < 6.0 TeV (electron channel) | **CONFIRMED** |
| RS graviton, ≳4.8 TeV, CMS diphoton | 4.8 TeV | CMS diphoton RS-graviton search (2405.09320 also serves diphoton) | **CONFIRMED** (ballpark; edge depends on k/M̄_Pl) |
| DY contact Λ_LL⁺ ≳ 35.8 TeV, PDG 2025 compositeness | 35.8 TeV | ATLAS dilepton non-resonant, LL constructive Λ = 35.8 TeV | **CONFIRMED** (number ok; underlying result is ATLAS, PDG compiles it) |
| VBS W_LW_L, σ-limit, ATLAS PRL 135 (2025) 111802 | no mass edge | ATLAS longitudinal same-sign WW VBS (arXiv:2503.11317) | **CONFIRMED** |
| diboson spin-1 W'/V', ≳4.4–4.8 TeV, CMS-B2G-20-009 Table 2 | 4.4–4.8 TeV | CMS HVT spin-1 Z'/W' < 4.8 TeV (publ. arXiv:2210.00043) | **CONFIRMED** |
| diphoton G*, ≳4.8 TeV, CMS arXiv:2405.09320 | 4.8 TeV | CMS high-mass diphoton, RS graviton | **CONFIRMED** |
| four-top Z', up to 0.85 TeV, CMS-B2G-25-005 / arXiv:2604.14058 | 0.85 TeV | CMS four-top 2ℓ, Z'(50% width) excluded up to 850 GeV | **CONFIRMED** |

All collider arXiv IDs — including the recent 2026 ones (2603.23454, 2604.14058) — **resolve to real CMS searches**, consistent with the 2026-06-17 date. None are hallucinated.

---

## Input-constant checks

| Constant | Doc value | Standard value | Verdict |
|----------|-----------|----------------|---------|
| v (Higgs vev) | 246.21965 GeV | 246.21965(6) GeV (from G_F) | **CONFIRMED** (exact) |
| sin²θ_W | 0.23122 | PDG ŝ²_Z(M_Z) MS-bar = 0.23122 ± 0.00006 | **CONFIRMED** (it is the MS-bar Z-pole value; note PDG also lists 0.23129/0.23161 in other determinations — doc's choice is a valid PDG number) |
| Δm_K | 3.484×10⁻¹⁵ GeV | PDG (3.484±0.006)×10⁻¹⁵ GeV | **CONFIRMED** (exact) |
| x_1 (M_KK = x_1·Λ_IR) | ≈2.45 | first warped KK gauge root ≈ 2.45 | **CONFIRMED** (standard RS value; the relevant root is the warped transcendental one, ≈2.45, not the flat J₁ zero 3.83) |
| L = kπr_c | ≈35 | −ln(ε) = 15·ln10 ≈ 34.5 (ε~10⁻¹⁵) | **CONFIRMED** (34.5–37 depending on ε; "≈35" standard) |
| χ² (95% CL, 2 dof) | 5.991 | −2 ln(0.05) = 5.9915 | **CONFIRMED** (exact) |
| ε_K^exp | 2.228×10⁻³ | (2.228±0.011)×10⁻³ | **CONFIRMED** |
| φ_ε | ≈43.5° | 43.52° (PDG/superweak) | **CONFIRMED** |

---

## (b) Physics / internal-consistency findings

All checks **pass** unless noted.

1. **ΔF=2 Wilson colour factors (1/6, −1, 1/3)** — CONFIRMED against Csaki–Falkowski–Weiler (the tree H displayed after Eq. 2.17 has these factors). **Cross-check vs Blanke et al.:** the two papers use slightly different overall normalisations. Blanke Eq. (4.17) gives `C₁^VLL = (2/3)p_UV²Δ_L²`, `C₁^LR = −(2/3)p_UV²Δ_LΔ_R`, `C₂^LR = −4 p_UV²Δ_LΔ_R` with H_eff carrying a `1/(4M_KK²)` prefactor and `p_UV²` (warped-profile) instead of `g_L g_R`. The doc's `{1/6, −1, 1/3}` with prefactor `1/M_KK²` and explicit `g_L g_R` is the **Csaki et al. convention**, which is internally consistent with how the doc defines `g_{L,R}` (overlap × strong coupling). No contradiction, but note the doc's LR ratio C₅:C₄ = (1/3):(−1), whereas Blanke's C₂^LR:C₁^LR = (−4):(−2/3) = 6:1 in magnitude — **these are different operator orderings/normalisations** (Blanke's Q₁^LR/Q₂^LR vs the doc's Q₄/Q₅). This is fine as long as the code's matrix elements (below) use the matching convention; flagged for the author to confirm the code's `{−1, 1/3}` pair with its `⟨Q₄⟩, ⟨Q₅⟩`.

2. **γ_VLL = 4 anomalous dimension** — Plausible/standard for the QCD LO VLL running; consistent with Buras–Misiak–Urban. UNVERIFIED to the exact value against BMU here (no equation fetched), but it is the textbook LO value and not in doubt.

3. **Matrix-element prefactors** SUPERSEDED by the June-2026 B3 audit and confirmed by the 2026-07 full-repo audit. The earlier review text certified the pre-B3 M12-ready values, with ⟨Q₁⟩=(2/3)f²mB₁ and the LR coefficients swapped and doubled. The current `deltaf2.py` M12-ready GGMS forms are ⟨O₁⟩=(1/3)f²mB₁, ⟨O₄⟩=(r_χ/4+1/24)f²mB₄, and ⟨O₅⟩=(r_χ/12+1/8)f²mB₅. The `O4` term carries the large LR coefficient and `O5` the small one. The old CONFIRMED verdict is retracted for the M12-ready forms.

4. **r_χ = (m_M/(m_q1+m_q2))²** — CONFIRMED.
   - Kaon: with PDG MS-bar(2 GeV) m_s≈93.5 MeV, m_d≈4.7 MeV → **r_χ ≈ 25.7** ✓ (doc says ≈26).
   - B: with m_b(m_b)≈4.18 GeV → **r_χ ≈ 1.59** ✓ (doc says ≈1.6).

5. **ε_K budget arithmetic:** 2.228×10⁻³ − 2.161×10⁻³ = **6.70×10⁻⁵** ✓ (doc says ≈6.7×10⁻⁵). Self-consistent. (Caveat: this hinges on the 2.161 SM value, which is BGS's 2.16 carried to an extra digit — see Ref. 12.)

6. **Z→bb̄ floor translation:** 25/2.45 = 10.20, 30/2.45 = 12.24 → "10–12 TeV in Λ_IR units" ✓ CONFIRMED. The note that the admixture prefactor uses m_b²/(2Λ_IR²) and so reads ~10–12 TeV in "Casagrande's M_KK ≡ Λ_IR convention" is consistent.

7. **B(c,F) = 1/(1−2c)·(1/F²−1+F²/(3+2c))** — sanity check PASS. For IR-localized profiles (c<1/2) it is finite and positive (e.g. c=0 → 1.167, c=−0.2 → 0.498, c=0.3 → 10.1 with F=f_IR). The `1/(1−2c)` factor flips sign for UV-localized (c>1/2) partners, exactly matching the doc's stated "sign-flipping 1/(1−2c) for the UV-localised light down partners." Dimensionless, limit-sane. As the doc says, this is a repo-local compact form, not a paper equation; it is internally consistent and matches the structure of Casagrande Eq. (170)'s admixture bracket (which contains the same `1/(1−2c)`, `1/F²`, `F²/(3+2c)` building blocks). CONFIRMED self-consistent.

8. **"Code evaluation" subsection (Generic ΔF=2)** — scrutinized; **physically sound and consistent with the rest of the document.**
   - `G_L = diag(g(c_Qi))`, `G_R = diag(g(c_di))`, then `G̃_{L,R} = U†_{L,R} G_{L,R} U_{L,R}`, off-diagonal `(i,j)` entries → g_L, g_R: this is exactly the Csaki–Falkowski–Weiler Eq. (2.16) rotation `g_{L,u/d} → U†_L g_q U_L` (the doc faithfully reproduces it). ✓
   - Statement "flavour-aligned ⇒ G̃ diagonal ⇒ g_L=g_R=0 ⇒ no effect" is correct: a diagonal g(c) commutes with U only if degenerate; misalignment generates off-diagonals. ✓ (Strictly, alignment requires either degenerate g(c) or U=identity; the prose "flavour-aligned point would leave G̃ diagonal" is the right intuition.)
   - Wilson coeffs at M_KK → RG run with thresholds m_t,m_b,m_c → multiplicative VLL (η_V), 2×2 LR mixing (matrix R) → contract with ⟨Q_i⟩ already in M₁₂ normalisation → M₁₂^NP = Σ C_i⟨Q_i⟩ → per-system budget / ε_K cut. This chain is internally consistent and matches the "Mixing amplitude and the cut" section. ✓
   - The `1/(2m_M)` absorption note (matrix elements pre-normalised for M₁₂) is consistent with Blanke Eq. (4.25): `2m_K (M₁₂^K)* = ⟨K̄|H_eff|K⟩`, i.e. M₁₂ = ⟨H_eff⟩/(2m_M). The doc's convention (absorbed) is a legitimate repacking. ✓

9. **Z→bb̄ "Code evaluation" subsection** — consistent. The gauge piece `δg_L^d|_gauge = s_Z g_L^SM (m_Z/M_KK)² U_L† diag(a(c_Qi)−a_ref) U_L` plus the admixture `δg_L^b|_admix = (m_b²/2Λ_IR²)B_d` map onto Casagrande Eq. (170)'s two brackets (m_Z²/M_KK² gauge mixing + m_b²/M_KK² fermion admixture). The pull/ratio cut (max(|p_Rb|,|p_Ab|) ≤ 1) with R_b⁰=0.21629±0.00066, A_b=0.923±0.020 uses the verified Eq. (173) anchors. ✓ Minor: the doc's "Sec. 6.1" should read **Sec. 6.4** (see Ref. 7).

10. **Oblique S/T proxy magnitude:** ΔT_minimal = (πL/2c_W²)v²/M_KK² and ΔS = 30 v²/M_KK². At M_KK=8 TeV: ΔT≈0.068, ΔS≈0.028; at 10 TeV: ΔT≈0.043, ΔS≈0.018. These land in/near the fit window, consistent with the doc's "~8–10 TeV" minimal-RS EW reach. T is volume-enhanced (∝L), S is not — matches Agashe–Delgado–May–Sundrum (5.5)/(5.6) and the prose. Dimensionally and magnitude-sane. ✓ (The *cut* uses the wrong PDG ellipse center — see Ref. 11 — but the *prediction-side* proxy is fine.)

---

## (c) Prioritized list of concrete fixes (most severe first)

**FIX 1 — [INCORRECT, highest priority] Oblique-fit central values & correlation (Ref. 11, lines 591–599, 617–620).**
The document's `S = 0.026 ± 0.075`, `T = 0.047 ± 0.066`, `ρ_ST = 0.90`, "Table 10.8" do not match any current PDG fit. PDG 2024 EW review (Erler & Freitas, rev. 31-May-2024) gives:
- U-free (Eq. 10.98): S = −0.04 ± 0.10, T = 0.01 ± 0.12, U = −0.01 ± 0.09;
- **U≡0 (Eq. 10.99): S = −0.05 ± 0.07, T = 0.00 ± 0.06**, with S–T correlation **93%**.
This matters because the χ² cut is centered on `(ΔS−0.026, ΔT−0.047)` with ρ=0.90 — i.e. the cut is run against the wrong ellipse (wrong center, slightly wrong correlation, and the central values even have the wrong sign for S and a much smaller central T). **Suggested fix:** replace with the U≡0 PDG numbers (S=−0.05±0.07, T=0.00±0.06, ρ_ST=0.93), update the boxed χ² `Δx` accordingly, and correct the table/equation references (PDG 2024 Eqs. 10.98–10.99, not "Table 10.8"). If a specific PDG-2025 print edition genuinely has different numbers, quote it explicitly; but the 0.026/0.047/0.90 triple should not stand. *(Also: at M_KK=8 TeV the proxy predicts ΔS≈0.03, ΔT≈0.07 — with the corrected center near (−0.05, 0.00) the pull changes materially, so this is not cosmetic.)*

**FIX 2 — [PARTIAL] Blanke Eq. (4.38) → (4.39) (Ref. 2, lines 249–250 and 353–354).**
The `M₁₂^q = (M₁₂^q)_SM·C_Bq·e^{2iφ_Bq}` convention is **Eq. (4.39)**, not (4.38). Eq. (4.38) is `ΔM_q = 2|M₁₂_SM+M₁₂_KK|`. **Suggested fix:** change "Eqs. (4.36)–(4.38) for observables … Eq. (4.38) … convention" to "Eqs. (4.36)–(4.39) for observables; Eq. (4.39) gives the M₁₂^q = (M₁₂^q)_SM·C_Bq·e^{2iφ_Bq} convention." (The (4.14)–(4.17) and (4.23)–(4.29) ranges are correct and can stay; (4.25) is specifically the M₁₂ normalisation, (4.27)–(4.29) the matrix elements.)

**FIX 3 — [PARTIAL] Casagrande Eq. (171) → (173) and Sec. 6.1 → 6.4 (Ref. 7, lines 553–555).**
Eq. (170) is correctly the LH Zbb̄ shift (m_b²/M_KK² admixture, ZMA) ✓. But the *measured* anchors R_b⁰=0.21629±0.00066 and A_b=0.923±0.020 are **Eq. (173)** (the experimental "pseudo-observables"); Eq. (171) gives the *theoretical definitions* of R_b⁰/A_b/A_FB as functions of couplings, and Eq. (172) the SM predictions. The subsection is **Section 6.4 "Z⁰bb̄ Couplings"**, not 6.1. **Suggested fix:** "Sec. 6.4, Eq. (170) … R_b⁰, A_b, A_FB^{0,b} defined in Eq. (171); the experimental values quoted (R_b⁰=0.21629±0.00066, A_b=0.923±0.020) are Eq. (173)." Same Sec. 6.1→6.4 correction applies in the Code-evaluation prose if it cites 6.1.

**FIX 4 — [PARTIAL, minor] ε_K^SM "2.161" over-precision (Ref. 12, lines 278–280).**
Brod–Gorbahn–Stamou (PRL 125 (2020) 171803) publish ε_K^SM = **2.16(6)(8)(15)×10⁻³**, i.e. central 2.16, not 2.161. The extra digit is not in BGS. **Suggested fix:** either write "≈2.16×10⁻³ (Brod–Gorbahn–Stamou)" and keep the budget as "≈6.7×10⁻⁵" (which 2.228−2.16 = 6.8×10⁻⁵ still rounds to ~7), or, if "2.161" is taken from a specific downstream compilation (e.g. a PDG/CKMfitter digit), cite that source for the third digit. Also consider a one-line caveat that exclusive-|V_cb| lattice updates (e.g. 2503.00351) give a lower SM ε_K (~65% of experiment); the document's choice of the manifest-CKM-unitarity BGS value is the standard one but is convention-dependent.

**FIX 5 — [labelling, cosmetic] Collider T_{5/3} row (Ref. 13, line 648).**
arXiv:2212.05263 calls the charge-+5/3 state "X" (not "T_{5/3}"), and the 1.46 TeV number is its B/X→Wt limit. The number is correct; only the symbol differs. Optional: rename to "X_{5/3} / T_{5/3}" or note "(B/X→Wt)" for precision. Low priority.

**Non-issues (explicitly checked, no fix needed):** v, sin²θ_W (0.23122 is a valid PDG MS-bar value), Δm_K, x_1≈2.45, L≈35, χ²=5.991, ε_K^exp=2.228×10⁻³, φ_ε≈43.5°, all r_χ values, the B(c,F) closed form, the Z→bb̄ 25–30 TeV ↔ 10–12 TeV translation, the S≈30 v²/M_KK² coefficient (PDG-confirmed to Carena et al.), and both "Code evaluation" subsections (physically sound and consistent with the surrounding equations).

---

## Sources used (for re-checking)

- Csaki–Falkowski–Weiler 0804.1954: https://ar5iv.labs.arxiv.org/html/0804.1954
- Blanke et al. 0809.1073: arXiv PDF https://arxiv.org/pdf/0809.1073 (eqs 4.14–4.45 extracted via pdftotext; (4.38)=ΔM_q, (4.39)=C_Bq·e^{2iφ})
- Buras–Guadagnoli–Isidori 1002.3612: https://ar5iv.labs.arxiv.org/html/1002.3612 (Eq. 34 master, Eq. 35 κ_ε=0.94±0.02)
- Buras–Guadagnoli 0805.3887: https://ar5iv.labs.arxiv.org/html/0805.3887 (κ_ε=0.92, Eqs. 9–12)
- Perez–Randall 0805.4652: https://ar5iv.labs.arxiv.org/html/0805.4652 (Eqs. 28–31, 4×10⁻⁸, C≈0.02)
- MEG II 2504.15711: https://arxiv.org/abs/2504.15711 (BR < 1.5×10⁻¹³ @ 90% CL)
- Casagrande et al. 0807.4937: arXiv PDF https://arxiv.org/pdf/0807.4937 (Sec. 6.4; Eq. 170, 171 defs, 172 SM, 173 exp R_b⁰=0.21629±0.00066, A_b=0.923±0.020)
- Agashe–Contino–Da Rold–Pomarol hep-ph/0605341: https://arxiv.org/abs/hep-ph/0605341 ("A custodial symmetry for Zbb")
- Carena et al. hep-ph/0607106 / Nucl. Phys. B 759 (2006) 202: https://arxiv.org/abs/hep-ph/0607106
- Agashe–Delgado–May–Sundrum hep-ph/0308036: https://ar5iv.labs.arxiv.org/html/hep-ph/0308036 (Eqs. 5.5–5.6, T∝kπr_c)
- PDG 2024 EW review (Erler & Freitas): https://pdg.lbl.gov/2024/reviews/rpp2024-rev-standard-model.pdf (Eqs. 10.98–10.99 oblique fit; "S≈30 v²M_KK⁻²" w/ Carena ref; v=246.22; ŝ²_Z=0.23122–0.23129)
- ε_K^exp / BGS SM / budget: BGS PRL 125 (2020) 171803 https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.171803 ; arXiv:1911.06822 ; ε_K^exp=2.228×10⁻³ (PDG)
- Δm_K: PDG 2024 listings (3.484×10⁻¹⁵ GeV)
- Collider rows: 2603.23454, 2604.14058, 2401.17165, 2212.05263, PRD 110 (2024) 052004 (arXiv:2402.13808), JHEP 07 (2021) 208 (arXiv:2103.02708), 1906.05609, 2405.09320, B2G-20-009 (arXiv:2210.00043), PRL 135 (2025) 111802 (arXiv:2503.11317), ATLAS dilepton contact Λ=35.8 TeV (arXiv:2006.12946)
- v=246.21965(6) GeV; sin²θ_W (ŝ²_Z)=0.23122; φ_ε=43.52°; x_1≈2.45; L=kπr_c≈34.5; χ²(0.95,2)=5.9915 — all verified by direct computation and/or PDG.
