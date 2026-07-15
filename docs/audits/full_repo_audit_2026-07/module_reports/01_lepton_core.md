# Report 01 — Lepton core stack (warpConfig, solvers, neutrinos, diagonalization, yukawa, flavorConstraints, scanParams)

## Structural summary

The lepton-sector stack: `warpConfig/baseParams.py` builds RS geometry (ε = Λ/k, r_c, brane positions); `warpConfig/wavefuncs.py` computes zero-mode overlap factors f_IR/f_UV. `solvers/bessel.py` finds KK tower masses from Bessel cross-product quantization conditions (gauge NN, fermion ++/--). `neutrinos/neutrinoValues.py` holds NuFIT inputs, computes mass spectra (NO/IO) and the PDG-parametrized PMNS matrix; `neutrinos/massConstraints.py` sweeps the lightest mass under the Σm_ν = 0.082 eV bound. `diagonalization/diag.py` provides SVD and Autonne–Takagi factorization. `yukawa/` inverts m_E = 2vk f_L Y f_E (charged leptons) and the universal-limit seesaw m_ν = 2k²v²f_L²f_N²Y²/((f_N^UV)²M_N) (neutrinos), then assembles Y_N_matrix = V_PMNS·diag(Y_N). `flavorConstraints/muToEGamma.py` applies the NDA dipole bound |(Ȳ_NȲ_N†)₁₂| ≤ C·(M_KK/3 TeV)². `scanParams/scan.py` grid-scans all of this with perturbativity/naturalness/LFV filters; `anarchy.py`/`postprocess.py` do prior scoring and CSV reclassification.

Verified numerically: f_IR/f_UV formulas and their c→1/2 limits, the identity f_UV = ε^(1/2−c)f_IR, PMNS unitarity + Jarlskog vs analytic, NO/IO Δm² sign conventions, both Yukawa inversion round-trips (reproduce PDG masses and target ν spectrum to machine precision), Takagi on 200+ random complex-symmetric matrices including degenerate/rank-deficient/near-real cases (max error ~5e-14; mild, bounded error growth near the sqrtm branch cut), LFV C-coefficient arithmetic (0.02 paper / 1.94e-3 MEG II, both consistent with BR = 4e-8·|·|²·(3 TeV/M_KK)⁴), and scan filter logic.

---

### [CRITICAL] KK solver skips and misorders roots for n_roots ≳ 5
- **File:** solvers/bessel.py:153-179 (`_bracket_around`), 264-285 (bracket assembly in `solve_kk`)
- **Category:** numerics / logic-bug
- **Claim:** `solve_kk` silently skips KK roots and can return unsorted or wildly wrong masses once the requested root index is ≳5, because brackets are built with a fixed ±20% relative window around each seed — which exceeds the ~π root spacing for x ≳ 15 — and results are never deduplicated, sorted, or verified against expected root count.
- **Evidence:** Brute-force sign-scan of the solver's own `_F_exact` vs `solve_kk` output (ε for Λ=3 TeV, k=M_Pl):
  ```
  gauge NN  solver: [2.449 5.566 8.700 11.838 14.978 21.259 24.400 27.541]
            true  : [2.449 5.566 8.700 11.838 14.978 18.118 21.259 24.400]   <-- 18.118 skipped
  fermion ++ c=0.58: solver skips 18.198; fermion ++ c=-0.7: skips 16.772;
  fermion -- c=0.58: skips 16.591
  fermion ++ c=0.27 (nu=-0.23): solver: [2.751 5.876 9.011 12.150 52.985 18.431 24.713 21.572]
            true : [2.751 5.876 9.012 12.150 15.290 18.431 21.572 24.713]   <-- 5th "root" is 52.98, output unsorted
  ```
  The 5th returned mass for c=0.27 is off by a factor 3.5 and the array is not monotonic. Cause: for seed x₀ the bracket (0.8x₀, 1.2x₀) has width 0.4x₀ > π when x₀ ≳ 8, so it can contain 2–3 roots with equal end-signs; the expansion loop (`a*=0.8, b*=1.25`) makes it worse, and `brentq` then converges to an arbitrary interior root. No warning is emitted (the "fewer roots" warning fires only when brackets are missing, not wrong). First 3–4 roots (the default `n_roots=3`) are correct in all tested cases.
- **Confidence:** high
- **Fix:** Bracket with a fixed absolute half-width (< π/2) around each seed, or sequentially scan in steps < π; then sort, dedupe (tolerance ~xtol), and assert consecutive roots are ~π apart.

### [MAJOR] Flagship benchmark point violates its own perturbativity bound (docs/examples misleading)
- **File:** yukawa/compute_yukawas.py:220-232 (docstring example), yukawa/__init__.py:6-22, CLAUDE.md "Quick Method"
- **Category:** convention-inconsistency (documentation vs physics)
- **Claim:** The canonical example point (Λ=3 TeV, c_L=0.58, c_E=[0.75,0.60,0.50], c_N=0.27) advertised everywhere as the validation benchmark yields Ȳ_E = [2.94, **4.37, 5.42**], i.e. it fails the |Ȳ| < 4 perturbativity bound the same docstring tells the user to check; the comment "Should be O(1) to O(few)" hides this.
- **Evidence:** Ran `compute_all_yukawas` with the documented inputs: `Y_E_bar = [2.9364, 4.3719, 5.4194]`, `is_perturbative() = False`. The scan test `tests/test_scan.py:37-38` even asserts `not row["perturbative"]` for this exact point — so the code is self-consistent, but every user-facing example presents a non-viable point as the reference workflow without a warning.
- **Confidence:** high (numeric); the "bug" is in the documentation/benchmark choice, not the arithmetic.
- **Fix:** Either change the benchmark c_E to values giving Ȳ < 4, or annotate the example as intentionally non-perturbative.

### [MINOR] τ mass is the pre-2023 world average but labeled "PDG 2024"
- **File:** yukawa/constants.py:6-9
- **Category:** stale-data
- **Claim:** `M_TAU = 1.77686` GeV is the PDG 2022 average; PDG 2023/2024 (after Belle II's 2023 measurement) list m_τ = 1776.93 ± 0.09 MeV, so the value is ~70 keV low while the comment claims "PDG 2024".
- **Evidence:** `M_TAU = 1.77686  # 1.777 GeV` under header `# Charged lepton masses (GeV) - PDG 2024`. m_e and m_μ are correct and current. The shift is 4e-5 relative — negligible for Ȳ, but the label is wrong.
- **Confidence:** medium-high
- **Fix:** Update to `M_TAU = 1.77693` or relabel the comment.

### [MINOR] `_F_exact` returns a fake constant sign at x ≤ 1e-12, and seed-based recovery can't rescue early roots
- **File:** solvers/bessel.py:113-117, 276-285
- **Category:** numerics
- **Claim:** For x ≤ `_MIN_X` the quantization function is replaced by `np.sign(x)*1.0`, an arbitrary positive value that need not match the true small-x sign of F (which is (2/π)ln ε < 0 for ν=0), so any bracket that expands down to `_MIN_X` acquires a spurious sign change; additionally, the fallback linear scan starts at `seeds[-1]`, so a root missed *below* the last seed can never be recovered.
- **Evidence:** `return np.sign(x) * 1.0  # avoid evaluating at 0` vs. numerically F(1e-3) < 0 for ν=0, ε=2.5e-16. `_scan_for_brackets(F, x_start=seeds[-1], ...)` only searches upward of the highest seed while the c=0.27 case above shows the missing root (15.29) lies *below* returned ones.
- **Confidence:** high
- **Fix:** Return the analytic small-x sign (or clamp brackets to a small positive floor above the first root region), and start the fallback scan from the lowest unresolved region, not `seeds[-1]`.

### [NOTE] Seesaw prefactor 2k²v² is in tension with the Dirac-mass convention 2vk by a factor of 2
- **File:** yukawa/neutrino.py:141-147; yukawa/charged_lepton.py:99-105
- **Category:** factor-of-2 (convention, flagged not asserted)
- **Claim:** If the neutrino Dirac mass follows the same IR-brane formula as charged leptons, m_D = 2vk f_L Y_N f_N, then a naive seesaw m_D²/M_eff with M_eff = M_N/(f_N^UV)² gives a prefactor **4**k²v², not the implemented **2**k²v²; the factor 2 is only consistent if the 4D Majorana mass carries an extra 2 (e.g. from ½M_N NN normalization or orbifold-fixed-point brane-term counting).
- **Evidence:** Code: `prefactor = (f_N_UV**2 * M_N) / (2.0 * k**2 * v**2 * f_L**2 * f_N**2)`; charged leptons: `Y_E = m_E / (2.0 * v * k * f_L * f_E_arr)`. The implemented formula matches CLAUDE.md and `scripts/audit_perez_randall_consistency.py` exactly (round-trip verified to machine precision), so the repo is internally consistent — but the origin of 2 vs 4 relative to the Dirac convention is nowhere derived, and the repo's own audit script shows the Perez–Randall Eq. (10) benchmark is off by O(3–4)× in kY_N, which an unresolved factor of 2 in the seesaw prefactor would partially feed into (√2 ≈ 1.4 of it).
- **Confidence:** low-medium (as a defect); high that the ambiguity is undocumented.
- **Fix:** Add a derivation note fixing the Majorana brane-term normalization that produces 2k²v², or reconcile with m_D = 2vk f_L Y f_N explicitly.

### [NOTE] NuFIT "6.1 (2025)" oscillation inputs could not be independently verified
- **File:** neutrinos/neutrinoValues.py:19-66
- **Category:** stale-data (unverified)
- **Claim:** Δm²₂₁ = 7.537e-5 eV², sin²θ₁₂ = 0.3088, δ_CP = 212°(NO)/274°(IO) are attributed to "NuFIT 6.1, IC24 with SK-atm (Nov 2025)"; these differ noticeably from NuFIT 6.0 (Δm²₂₁ = 7.49e-5, sin²θ₁₂ ≈ 0.308) and could not be confirmed this session. All internal usage (Δm² signs, NO/IO branches, legacy `delta_m32_sq_IH ≡ m₁²−m₃²` alias) is self-consistent and numerically correct.
- **Confidence:** low (possible fabricated/rounded dataset values; possibly genuine JUNO-era update)
- **Fix:** Cross-check the five best-fit values against the actual NuFIT 6.1 table before publication use.

### [NOTE] Scan "naturalness" filter silently rejects every point with a massless lightest neutrino
- **File:** scanParams/scan.py:517-521
- **Category:** logic-bug (edge case)
- **Claim:** The naturalness check requires min|Ȳ| ≥ 0.1 over the *combined* Y_E_bar and Y_N_bar arrays, so `lightest_nu_mass = 0` (Y_N₁ = 0) can never pass, and more generally the smallest neutrino Yukawa — set by the input lightest mass, not by model naturalness — dominates the filter.
- **Evidence:** `all_y = np.concatenate([abs(Y_E_bar), abs(Y_N_bar)]); is_natural = (min_y >= lo)`; at the default benchmark Ȳ_N₁ = 0.204, only ~2× above the 0.1 floor.
- **Confidence:** high (behavior), medium (that it's unintended)
- **Fix:** Apply the lower naturalness bound to charged-lepton Yukawas only, or document the coupling to `lightest_nu_mass`.

### [NOTE] Takagi sqrtm branch cut: mild accuracy loss for nearly-real indefinite matrices
- **File:** diagonalization/diag.py:109-120
- **Category:** numerics
- **Claim:** For complex-symmetric inputs that are a tiny imaginary perturbation of a real indefinite matrix, Z has eigenvalues near −1 (sqrtm branch cut) and reconstruction error grows to ~the perturbation scale (observed up to 5e-7 for 1e-10 imaginary noise in one seed) instead of machine epsilon; never catastrophic in 200+ random trials, and the real-path shortcut handles the exactly-real case correctly.
- **Evidence:** Numerical probe: imag scale 1e-8 → error 1.0e-8; one seed at 1e-10 → 4.9e-7. Random complex-symmetric worst case 5e-14; degenerate and rank-deficient cases exact.
- **Confidence:** high (measured), low severity
- **Fix:** None required for current use; if hardened, use the eigendecomposition-of-A†A Takagi algorithm or symmetrize sqrt(Z).

---

## Cross-module convention inconsistencies

1. **Yukawa dimension conventions differ between sectors:** lepton code uses dimensionful 5D Y with m = 2·v·k·f_L·Y·f_E and Ȳ = 2kY; the quark-sector tests (`tests/test_diagnostics.py:45`) use m_q = 2·v·F_Q·Y·F_q with dimensionless Y (no k). Anyone porting formulas across sectors will pick up a stray 2k.
2. **M_KK conventions coexist:** the LFV module defaults to the internal convention M_KK = Λ_IR (xi_KK = 1) while also shipping `GAUGE_KK_ROOT_NN = 2.4487` for the physical first-KK convention; the bound scales as M_KK², so mixing conventions shifts the μ→eγ floor by ×6. It is documented, but the scan CSV column `M_KK` records the internal (unphysical) value by default.
3. **`k` vs `k̄` (reduced Planck mass):** production code uses k = M_Pl = 1.2209e19, while `scripts/audit_perez_randall_consistency.py` uses k = M̄_Pl = 2.435e18 for the "paper-like" geometry — correct per-context but easy to conflate; both constants live in `baseParams.py` with no guard.
4. **Perturbativity strictness:** `is_perturbative` and `postprocess.classify_row` use strict `<` while the μ→eγ check uses `lhs <= rhs`; cosmetic, but boundary points classify differently across filters.

## Not reviewed
- `derivations/` LaTeX contents (only listed, not read) — the seesaw factor-of-2 question above could be settled there.
- Notebooks (`wavefuncsTest.ipynb`, `besselExample.ipynb`, `PMNS.ipynb`, `allowedMass.ipynb`, `muToEGamma.ipynb`, `diagonalizationTest.ipynb`).
- Tests `test_lmfv_lepton_parameters.py`, `test_finite_stats.py`, and the quark-sector/catalog test suite (out of scope).
- Live verification of NuFIT 6.1 and MEG II 2025 numbers against the web.
