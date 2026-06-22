# Slice 2 — MFV model construction & KK-gluon couplings (DONE)

Agent id: aa4a9ba5107a5fc73 (~130k tokens, 27 tool uses). Verbatim final report:

---

# Audit report: RS quark-sector model construction & KK-gluon couplings

## FINDINGS

### F1. ε_K depends on arbitrary SVD column phases — no rephasing to a physical CKM phase convention — **BLOCKER**, confidence: high (numerically demonstrated)

**Where:** `quarkConstraints/fit.py:245-258` (`mass_matrix_observables` — raw `numpy.linalg.svd` factors, no phase fixing), consumed by `quarkConstraints/couplings.py:149-152` and `quarkConstraints/deltaf2.py:787-789` (`evaluate_epsilon_k` uses `Im(M12^NP)`).

**Issue:** The SVD factors `U_L_d`, `U_R_d` are defined only up to per-column phases `U_L→U_L P`, `U_R→U_R P` (P diagonal unitary). This freedom leaves masses, |CKM| and the Jarlskog invariant (everything the fitter scores) unchanged, but rotates the mass-basis couplings `D_ij → e^{i(φ_j−φ_i)} D_ij`, hence all four ΔF=2 Wilsons by a common phase `e^{2iΔφ}`. `|M12|` observables (ΔM_K, B_d, B_s, D — all use `abs(m12)`) are invariant, but **ε_K uses Im(M12^NP), which is not**. The code never rephases to the standard (PDG) CKM convention in which `EPSILON_K_SM = 2.161e-3` (Brod–Gorbahn–Stamou) and the NP budget are defined, so the Im/Re split is whatever LAPACK happens to return.

**Demonstration (run in repo):** for the default benchmark fit at r=0.25, g*=3: applying a random physically-equivalent column rephasing to `(U_L_d, U_R_d, ckm)` left masses, |CKM|, J, and B_d/B_s/D ratios bit-identical, but changed the ε_K ratio_to_bound from **17.34 to 105.1** (factor ~6).

**Impact:** Per the methodology note's own appendix, ε_K binds 99.50% of the post-audit RUNA PDG-passing ensemble — i.e., the headline M_KK floors (47.26 TeV p50 at g*=3, etc.) are percentiles of a convention-dependent quantity. M_KK floor per point shifts as √ratio (O(1)–O(few), either direction). Worse conceptually: the 5D-MFV program's central claim (0710.1869: MFV ties CP phases to CKM, protecting ε_K) cannot be exhibited, because any calculable phase alignment is scrambled by the unfixed convention phase. For anarchic random-spurion ensembles the *distribution* may be roughly right (the spurious phase acts like one more random O(1) phase), but point-level pass/fail and any alignment/protection structure in r are not well-defined. Both lanes affected (modern `build_modern_point_couplings` calls the same `compute_quark_kk_gluon_couplings` on the same un-rephased fit result). `tests/test_epsilon_k_physics.py` only tests synthetic hand-built couplings, so this is untested.

**Fix:** after SVD, rephase columns of `U_L_u/U_L_d` (with matched `U_R` columns) so the CKM matrix is in the PDG standard-parameterization phase convention (5 quark rephasings, global phase irrelevant) before building mass-basis couplings; or compute ε_K from a rephasing-invariant combination relative to the SM `(V_ts*V_td)²` phase.

### F2. Default KK-gluon coupling normalization g_s(M_KK)≈1.05 vs literature g_s*≈3–6 — **documented convention, not a hidden bug** (would be MAJOR if undocumented), confidence: high

`quarkConstraints/couplings.py:145-147`: default `g_s = sqrt(4π α_s(M_KK)) ≈ 1.05`; the true first-KK-gluon vertex for IR-localized fermions carries the √(2πkr_c)≈8.5 / γ(c) enhancement (note's own R(c): ≈4 at c=½, ≈8 at c→−∞). The docstring and methodology note §"What g_s* is" document this exhaustively; the paper headline is g*=3, `modern/inputs.py:801` defaults `g_s_star=3.0`. However: legacy `quarkConstraints/scan.py:378` and all three `validation.py` call sites hardcode `g_s_star=None` (perturbative), so any floor quoted off those paths is ~×2.9 weaker than the headline convention. Additional un-pinned O(1): the code's flavor-basis profile is `diag(F²)` with F²=(½−c)/(1−ε^{1−2c}) = f_CFW²/2, and the CFW vertex is ≈ g_s*·f_i f_j·γ(c) with γ(c)~O(1); neither the factor 2 (f² vs F²) nor γ(c) is pinned by any test against a literature benchmark — it is absorbed into the g*∈{1.05, 3, 6} dial. Flag: the "g*=3" headline is a convention choice ~2× below CFW tree-level matching, stated as such in the note.

### F3. Universal (flavor-diagonal) KK-gluon coupling piece omitted entirely — **MINOR now, latent MAJOR for future consumers**, confidence: high

`couplings.py:_mass_basis_overlap` uses only `diag(F²)`; the universal piece ≈ −g_s/√(2πkr_c) ≈ −0.12 g_s is dropped. Because it is ∝ 𝟙, it exactly cancels in mass-basis off-diagonals, and `deltaf2._pair_couplings` reads only i≠j entries → ΔF=2 unaffected (verified). But the *diagonal* entries of `left_*`/`right_*` are wrong in both magnitude and sign for UV-localized quarks (code gives +g_s F² ~ +10⁻⁴ g_s for c≈0.6 vs true ≈ −0.12 g_s, off by ~10³ and sign). Currently no consumer: `collider_resonance.py` is an explicit mass-limit proxy (`NEEDS-HUMAN-PHYSICS` documented) that never touches couplings. Any future σ·BR recast, KK-gluon width, or SM-NP interference computed from these matrices' diagonals would be badly wrong. Recommend a guard comment or adding the −1/(2πkr_c)·𝟙 term.

### F4. ξ_KK bookkeeping (M_KK = Λ_IR) vs physical m_gKK = 2.449·Λ_IR — **documented; one residual approximation**, confidence: medium-high

- `GAUGE_KK_ROOT_NN = 2.448687135269161` **verified numerically**: it is the first root of J₀(x)Y₀(εx) − Y₀(x)J₀(εx) (Neumann–Neumann gauge BCs) at ε = 3000/1.2209e19; J₀(x) = −0.022558 matches Y₀(x)J₀(εx)/Y₀(εx) = −0.022576 at my warp-log precision. `SPIN2_GRAVITON_KK_ROOT = 3.8317` = first zero of J₁ ✓.
- No double conversion found: explicit `M_KK` short-circuits the Λ→M helper and `xi_KK` is recomputed as `M_KK/Λ_IR` (`couplings.py:132-143`; pinned by `tests/test_quark_couplings.py:59-63`). Publication scripts (`export_accepted_quark_scan.py`, `plot_publication_figures.py`, `export_collaborator_*`) apply the root exactly once via `xi_kk=GAUGE_KK_ROOT_NN` or `ratio_rescale = 1/ROOT²`.
- Residual issue: the `1/ROOT²` ratio rescale in `scripts/export_collaborator_5tev_points.py:313` is only exact if α_s(M_KK) and the M_KK→2 GeV RG factors were unchanged; matching at 7.35 TeV vs 3 TeV shifts the LR running enhancement at the several-% level. The script also recomputes exact physical-mass phenomenology alongside, so MINOR.
- Note the default-lane internal inconsistency (documented as "bookkeeping"): with ξ=1 the propagator uses Λ_IR² where the physical KK-gluon mass² is 6× larger — conservative (over-excludes) by ~×2.45 in the floor, partially offsetting F2's perturbative-g_s direction. The two dials are independent and both documented, but their *product* is what fixes the floor; anyone comparing repo floors to literature must set both (the publication scripts do).

### F5. ε_K Im-part aside: arbitrary phases also enter Re/Im of C4 vs C1 interference — subsumed by F1.

## VERIFIED CORRECT (explicitly checked)

1. **MFV chirality assignment** (`model.py:235-237`): `C_u = Y_u†Y_u`, `C_d = Y_d†Y_d` (singlets), `C_Q = r·Y_uY_u† + Y_dY_d†` (doublets) — correct sides; Hermitian and PSD by construction (r≥0 enforced), `eigh` valid.
2. **Eigen-ordering & c-map** (`model.py:214-222`, `BulkMassMap.map_eigenvalues`): eigenvalues sorted ascending, eigenvectors permuted in lockstep (`vecs[:, idx]`); map c = c_uv − (c_uv−c_ir)·λ/(λ+1) is strictly monotone decreasing in λ → gen-3 (largest spurion eigenvalue) gets smallest c (most IR-localized), and f_IR is monotone decreasing in c → largest overlap. Bounded window [0.30, 0.72] keeps f_IR finite, straddles c=½ sensibly.
3. **Bulk-basis rotation sides**: `Y_bulk = R_Q† Y R_u/d` (`model.py:251-252`) — correct for Q̄ Y u with fields `flavor = R · bulk`; pinned by `tests/test_quark_model.py:33-43`.
4. **SVD conventions** (`diagonalization/diag.py:8-18`): wrapper returns `(U, s, Vh.conj().T)` i.e. genuine V, so `U_L† M U_R = diag(s)` — the classic Vh/V conflation is **not** present. Light-to-heavy reorder applied to columns of both factors consistently (`fit.py:227-231`).
5. **Mass-basis coupling rotation** (`couplings.py:45-48`): `D = U† diag(F²) U` with U = the corresponding sector/chirality SVD factor; down-left uses `U_L_d`, up-left uses `U_L_u` (both with F_Q — doublet profile sector-independent), right sectors use `U_R_u/F_u`, `U_R_d/F_d`. Identity `left_up = V_CKM · left_down · V_CKM†` holds and is tested (`tests/test_quark_couplings.py:51-52`).
6. **CKM** = `U_L_u† U_L_d` (`fit.py:249`) — standard convention.
7. **Mass normalization / v convention**: `M = 2·v·F_Q Y F_u` with `V_EWSB = 174` (`fit.py:239`, `warpConfig/baseParams.py:4`). Since f_lit = √2·F_repo, this is identically m = (246/√2)·Y·f_Q f_q — exact match to CFW/literature convention. Same v constant shared with lepton modules. CLAUDE.md's m = 2vk f Y₅ f convention ⇒ code spurion Y = kY₅ (dimensionless), consistent.
8. **f_IR** (`warpConfig/wavefuncs.py`): matches documented f² = (½−c)/(1−ε^{1−2c}); c→½ limit 1/(−2 ln ε) correct; no overflow within the bulk-map window.
9. **ΔF=2 tree matching pattern** (`deltaf2.py:342-345`): C1_VLL = g_L²/(6M²), C4_LR = −g_L g_R/M², C5_LR = +g_L g_R/(3M²) — matches the standard color-octet KK-gluon Fierz results (CFW 0804.1954) including relative signs; ε_K normalization (2/3 f²m B factors, κ_ε/√2ΔM_K, no double factor of 2) separately audited per `docs/STATE_OF_PROJECT.md:29`.
10. α_s evaluated at the resolved M_KK with the 4-loop `qcd` module — right scale for the chosen convention.

## Numerical impact summary on M_KK floors

| Finding | Effect on floors |
|---|---|
| F1 ε_K phases | O(1)–O(6) per point, either direction (demonstrated ×6 in ratio ⇒ ×~2.5 in per-point floor); ensemble percentiles distorted by an unphysical random phase; MFV-protection claims not testable. Dominant constraint (binds 99.5%). |
| F2 g_s default | documented ×2.9 between pert. and g*=3 headline (note quotes both); legacy scan/validation lanes on the weak convention. |
| F3 universal piece | zero on ΔF=2; ~10³ + sign error on diagonal couplings if ever consumed. |
| F4 ξ_KK | documented; ~few-% residual from the 1/ROOT² rescale ignoring α_s/RG scale shift. |

The single actionable defect is F1; everything else checked is either correct or explicitly documented as a convention/approximation in `docs/quark_scan_methodology_note.tex` (gauge-coupling section + scope appendix) and the couplings docstring.
