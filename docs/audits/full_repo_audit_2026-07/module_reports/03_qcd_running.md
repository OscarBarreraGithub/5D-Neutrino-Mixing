# Report 03 — QCD/RG stack (qcd/, quarkConstraints/qcd_running.py, pdg_quark_masses.py, ckm_extraction.py)

**Structural summary:** The `qcd/` package (4-loop β, 4-loop γ_m, ODE running, CKS threshold decoupling) is largely correct: β₀–β₃, γ₀–γ₃, and the α_s decoupling constants c2=11/72, c3 all match the literature, and spot-checks reproduce RunDec-grade values (α_s(2 GeV)=0.3015, α_s(m_b)=0.2246, α_s(m_t)=0.1084, α_s(1 TeV)=0.0885; m_c(3 GeV)=0.989, m_b(M_Z)=2.866, m_t(3 TeV)=134.4). `quarkConstraints/qcd_running.py` (LO ΔF=2 Wilson running) is self-consistent: γ_VLL=+4, the BMU→scalar-basis LR map [[−16,−6],[0,2]] is algebraically correct, and the evolution reproduces the known LO magic-number structure (C4 factor 3.54 from 3 TeV→2 GeV). `ckm_extraction.py` is clean: standard rephasing-invariant β/β_s definitions, exact unitarity, J=3.12×10⁻⁵, sin2β=0.708. Genuine problems: wrong 3-loop mass-decoupling constants d3, a spurious threshold matching triggered by the 163.5 vs 162.5 GeV m_t(m_t) inconsistency, and an understated PDG top-mass uncertainty that a test actively enforces. Numerical impact of all found bugs is sub-percent, but they are real physics/transcription errors.

### [MAJOR] 3-loop MS-bar mass-decoupling constants d3 are wrong (value and n_l-slope sign)
- **File:** qcd/decoupling.py:117-128
- **Category:** wrong-formula
- **Claim:** The tabulated 3-loop ζ_m constants (d3=1.39151/1.36214/1.33277 for n_l=3/4/5, decreasing with n_l) do not match the CKS closed form quoted in the code's own comment, which gives 1.9218/1.9465/1.9713 (increasing with n_l).
- **Evidence:** Code comment gives `d3 = 2951/2916 - 407/864 * zeta(4) + 5/4 * zeta(4) - B4/36 + n_l*(1327/11664 - 2/27*zeta(3))` (first ζ4 should be ζ3 — CKS hep-ph/9708255 Eq. 17). Evaluating the correct form (B4=−1.7628): nl=3: CKS d3 = 1.9218 vs code 1.39151; nl=4: 1.9465 vs 1.36214; nl=5: 1.9713 vs 1.33277 — code slope is −0.0294/n_l vs correct +0.0247/n_l (RunDec: ζ_m = 1 + 0.2060 a² + (1.8476 + 0.0247 n_l) a³). Measured impact: matching-factor error ≈ 1.0–1.1×10⁻³ (relative) at the charm threshold, 2.1×10⁻⁴ at m_b.
- **Confidence:** high
- **Fix:** Replace d3 with the CKS closed form `1.84764 + 0.02473*n_l` (and fix the ζ4→ζ3 typo in the comment).

### [MAJOR] Spurious mass matching when the top's own reference scale (162.5) sits below the threshold table's m_t (163.5)
- **File:** qcd/mass_running.py:187-191, 230-246 (with qcd/decoupling.py:132-180 lacking a no-op guard)
- **Category:** logic-bug
- **Claim:** Running m_t from μ_ref=162.5 (n_f_ref=6) upward past the 163.5 GeV table threshold produces an nf "6→6" crossing that still applies a threshold matching — and because `direction = "up" if nf_next > nf else "down"` evaluates to "down", it applies the *downward* matching factor while running up.
- **Evidence:** `run_msbar_mass(162.5, 162.5, 3000, n_f_ref=6, matching_loops=3)` = 134.4405 vs `matching_loops=0` = 134.4002 → spurious +3.0×10⁻⁴ relative shift (~+40 MeV on m_t(3 TeV)). `match_alpha_s` has an `n_f_from == n_f_to` early return; `match_msbar_mass` has none. Root cause: `qcd/constants.py:19` `M_TOP_MS = 163.5` vs `pdg_quark_masses.py:112-115` t reference 162.5 (PDG 2024).
- **Confidence:** high
- **Fix:** Skip matching when `nf_next == nf`, and set `M_TOP_MS` = 162.5 to agree with the PDG table.

### [MAJOR] PDG top-quark MS-bar mass uncertainty understated ~2.5×, and a test bakes it in
- **File:** quarkConstraints/pdg_quark_masses.py:110-116; tests/test_pdg_quark_masses.py:64-73
- **Category:** stale-data / test-bug
- **Claim:** `sigma_GeV=0.7` for m_t(m_t)=162.5 does not match PDG 2024, which quotes the MS-bar top mass (from cross sections) as 162.5 +2.1/−1.5 GeV; 0.7 looks like a pole-mass/direct-measurement uncertainty pasted onto the MS-bar value, and it tightens the fitter's 2σ window ~2.5×.
- **Evidence:** PDG 2024 Summary Table: m_t(MS-bar from cross-section) = 162.5 +2.1/−1.5 GeV (symmetric avg ≈1.8). Test comment enforces it: `# ... a tight 0.0086 for top` and asserts `rel["t"] < rel["s"]` — with the correct σ≈1.8 the test would fail, i.e. the wrong number is load-bearing in the test.
- **Confidence:** high
- **Fix:** Set σ_t ≈ 1.8 GeV (symmetrized PDG 2024) and update the test's ordering assertion/comment.

### [MINOR] run_msbar_mass anchors alpha_s in a scheme that ignores `n_f_ref`
- **File:** qcd/mass_running.py:171-196
- **Category:** convention-inconsistency
- **Claim:** The α_s anchor is `qcd.running.alpha_s(mu_ref)`, whose flavor scheme is fixed by the scale (via `_n_f_at_scale`), not by the caller's `n_f_ref`, so the γ_m/β integration can start with an α_s in the wrong flavor scheme.
- **Evidence:** For the top: `_alpha_s(162.5)` returns α^(5) (162.5 < 163.5 threshold) but the segment integrates with nf=6; for `run_msbar_mass(m, M_CHARM, 2.0, n_f_ref=3)` the anchor is α^(4)(m_c) but the segment runs with nf=3. Scheme mismatch is O(2×10⁻⁴) at m_t, ~2×10⁻³ at m_c — hidden by the 0.5% test tolerance.
- **Confidence:** high
- **Fix:** After anchoring, apply `match_alpha_s` so the starting α_s scheme equals `n_f_ref`.

### [MINOR] α_s(M_Z) reference disagrees across modules (0.1180 vs 0.1179)
- **File:** qcd/constants.py:8 vs quarkConstraints/qcd_running.py:34
- **Category:** convention-inconsistency / stale-data
- **Claim:** `qcd` uses the PDG 2024 world average 0.1180 while the ΔF=2 Wilson running uses the older 0.1179.
- **Confidence:** high
- **Fix:** Import `qcd.constants.ALPHA_S_MZ` in `qcd_running.py` (LO effect negligible, but unify).

### [NOTE] Stale top pole mass and light-quark σ symmetrization
- **File:** qcd/constants.py:22; quarkConstraints/pdg_quark_masses.py:74-95
- **Category:** stale-data / statistics
- **Claim:** `M_TOP_POLE = 172.69` is the PDG 2022 average (PDG 2024/2025: 172.57±0.29); u-quark σ=0.49 MeV takes only the plus side of PDG's asymmetric +0.49/−0.26 despite the docstring's "averaged" claim (average ≈0.38); d-quark σ=0.43 similar (lower confidence).
- **Confidence:** medium
- **Fix:** Update M_TOP_POLE to 172.57 and recompute symmetrized light-quark σ consistently.

### [NOTE] CKM target δ = 1.196 rad is ~2σ above the PDG global-fit δ
- **File:** quarkConstraints/modern/inputs.py:638-641 (consumed by ckm_extraction.py)
- **Category:** stale-data
- **Claim:** δ=1.196 rad (68.5°) vs PDG fit δ≈1.14-1.15 rad pushes J to 3.118×10⁻⁵ (still inside 1σ) while reproducing sin2β=0.708 exactly; likely a deliberate joint tuning, not a transcription bug.
- **Confidence:** low
- **Fix:** Document the provenance of δ=1.196 or retune to the PDG fit quadruple.

## Cross-module convention inconsistencies
1. **m_t(m_t):** qcd/constants.py (163.5, also the nf=5/6 threshold and the fitter's `mu_common`) vs pdg_quark_masses.py (162.5, PDG 2024) — directly causes the spurious-matching MAJOR; quarkConstraints/qcd_running.py also hardcodes 163.5.
2. **α_s(M_Z):** 0.1180 vs 0.1179.
3. **m_b/m_c thresholds:** qcd/constants.py 4.18/1.27 vs PDG mass table 4.183/1.273 reference scales (negligible; should be one source of truth).
4. Verified consistent: β-function conventions self-consistent in both packages; BMU→scalar LR basis map and γ_VLL sign verified against LO magic-number structure.

## Not reviewed
- qcd/alphaS.ipynb; tests/test_qcd_running.py Section 3 integration fixtures; downstream consumers (deltaf2.py etc.).
