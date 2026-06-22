# Slice 6 audit — QCD running package (`qcd/`) and alpha_s consumers

Auditor: Claude (skeptical theoretical-physics code review), 2026-06-10.
Scope: `qcd/beta_function.py`, `qcd/running.py`, `qcd/decoupling.py`, `qcd/mass_running.py`,
`qcd/constants.py`, consumers `quarkConstraints/{qcd_running,couplings,pdg_quark_masses,benchmarks,scales}.py`,
plus the cross-slice check on the chiral-enhancement quark-mass scales in `quarkConstraints/deltaf2.py`.
Tests inspected/run: `tests/test_alpha_s.py`, `tests/test_mass_running.py`, `tests/test_qcd_running.py`,
`tests/test_pdg_quark_masses.py` (68 passed).

**Reference standard used:** CRunDec3 via the `rundec` Python package (installed for this audit,
`pip install --user rundec`, v0.7). All numerical cross-checks below are against CRunDec with
alpha_s(M_Z)=0.1180, MS-bar thresholds at m_c(m_c)=1.27, m_b(m_b)=4.18, m_t(m_t)=163.5, 4-loop
running / 3-loop decoupling, matching at mu = m_h (log-free), i.e. exactly the conventions the
package documents.

---

## Findings

### F1. MAJOR — wrong 3-loop mass-decoupling constants d3 in `qcd/decoupling.py`
- **Where:** `qcd/decoupling.py:112-129` (`_coeffs_msbar_mass`), used by `match_msbar_mass`.
- **Issue:** The tabulated 3-loop MS-bar mass-decoupling constants
  `d3 = 1.39151 (n_l=3), 1.36214 (n_l=4), 1.33277 (n_l=5)` are wrong, in both intercept and the
  **sign of the n_l slope**. The correct CKS value at mu = m_h (verified two independent ways) is

      d3(n_l) = 2951/2916 − (407/864)ζ3 + (5/4)ζ4 − B4/36 + n_l·(1327/11664 − (2/27)ζ3)
              = 1.84763 + 0.02473·n_l
      → 1.9218 (n_l=3), 1.9466 (n_l=4), 1.9713 (n_l=5)

  Extraction of CRunDec's coefficients (`DecMqDownMS` order-by-order at mu=m_h) gives
  d2 = 0.206019 = 89/432 (matches repo) and d3 = 1.9218 / 1.9465 / 1.9712 for n_l = 3/4/5 —
  i.e. the closed form above, **increasing** with n_l. The repo's numbers (decreasing with n_l,
  slope −0.0294) match neither CRunDec, nor the closed form the code's own comment quotes (the
  comment additionally typos `407/864·ζ(4)` for `407/864·ζ(3)`). The claim in the comment that the
  values were "cross-checked against the full literature" is false.
- **Measured impact:** every mass run that crosses the b threshold inherits a 2.2e-4 relative error
  (measured: m_s(2 GeV)→m_s(M_Z): repo 0.05395419 vs CRunDec-exact-ODE 0.05394236, rel 2.19e-4;
  same 2.19e-4 at the 163.5 GeV fitter scale, affecting the u,d,s,c PDG-2024 target masses).
  Crossing the top threshold: 2.5e-5 (decoupling factor 0.99969982 vs 0.99967443). All effects are
  ~10x below the smallest PDG input uncertainty involved (m_b: 1.7e-3 relative) and have no impact
  on any downstream constraint or scan conclusion at present.
- **Severity:** MAJOR (genuine wrong literature constant in a precision-flagship module with a false
  provenance claim; effective matching accuracy silently degrades from 3-loop to 2-loop).
  Numerical impact today: harmless (≤2.2e-4).
- **Confidence:** high (reproduced against CRunDec and the closed form independently).
- **Fix:** replace the three constants by `d3 = 1.84763 + 0.02473*n_l` (or the exact closed form
  with B4 = 16 Li4(1/2) − 13/2 ζ4 − 4 ζ2 ln²2 + 2/3 ln⁴2 ≈ −1.762800), and fix the comment's ζ3 typo.

### F2. MAJOR — `run_msbar_mass` breaks when `n_f_ref` is inconsistent with the threshold table; `pdg_quark_masses_at_scale` crashes below m_b
- **Where:** `qcd/mass_running.py:97-113` (`_ordered_thresholds_between`) + segment builder
  (`qcd/mass_running.py:181-192`); consumer contract at
  `quarkConstraints/pdg_quark_masses.py:132-168` (`pdg_quark_masses_at_scale`).
- **Issue:** The segment builder takes the caller's `n_f_ref` on faith and only crosses thresholds
  strictly between `mu_ref` and `mu_target`. The PDG top input is (m_ref=162.5, n_f_ref=6), but
  162.5 GeV sits **below** the package's top threshold (163.5), so a downward run from the top
  reference never decouples the top: the state stays n_f=6 until it hits the m_b crossing tuple
  (4.18, 4, 5), where `match_alpha_s(n_f_from=6, n_f_to=4)` raises
  `ValueError: Matching only supports adjacent thresholds`.
  **Reproduced:** `pdg_quark_masses_at_scale(mu)` (default flavors include 't') crashes for every
  `mu < 4.18 GeV` — including 2 and 3 GeV — directly contradicting its documented contract
  ("Must be strictly above the charm mass — no path runs below m_c"). For 4.18 ≤ mu < 163.5 the top
  run silently evolves with n_f=6 γ_m/β while anchored to the α^(5)(162.5) value returned by
  `alpha_s` (scheme mixing; measured effect ~1e-6 on the mass — negligible, but mislabeled).
- **Impact:** the production fit path (`mu = 163.5`, `benchmarks.py`/`scales.py`) is unaffected and
  verified correct. But this is exactly the consistent-scale call one would use to fix the deltaf2
  chiral-mass scales at 2 GeV (see F6), so the advertised remedy is currently a landmine unless
  `flavors=` excludes `'t'` (verified: `pdg_quark_masses_at_scale(2.0, flavors=("u","d","s","c","b"))`
  works and returns sane values: m_c(2)=1.0936, m_b(2)=4.9608).
- **Severity:** MAJOR (latent crash + docstring-contract violation on a public helper).
- **Confidence:** high (reproduced; root cause traced).
- **Fix options:** (a) validate/derive n_f at `mu_ref` from the threshold table and insert the top
  decoupling step when the path implies it (a zero-length matching step at 163.5 when starting at
  162.5 with n_f_ref=6 and running down), or (b) at minimum raise a clear error and fix the
  `pdg_quark_masses_at_scale` docstring.

### F3. MINOR — alpha_s at exactly mu = m_threshold returns the unmatched (lower-scheme) coupling
- **Where:** `qcd/running.py:134-156` (strict inequalities `t_ref < t_thresh < t_target`) vs
  `_n_f_at_scale` (`mu >= mass` → n_f above).
- **Issue:** `alpha_s(163.5)` = 0.108447 = α^(5)(m_t) (no 3-loop matching applied), while
  `alpha_s(163.5 + ε)` = 0.108425 = α^(6); `_n_f_at_scale(163.5)` labels the point n_f=6. Pure
  boundary-convention ambiguity (the two schemes genuinely differ by 2e-4 at the threshold), but the
  returned scheme at the exact threshold is inconsistent with the module's own n_f labeling.
- **Impact:** none downstream (KK-scale callers are at multi-TeV; mass running handles its own
  matching). Severity MINOR, confidence high.

### F4. MINOR — stale doctest value in `qcd/running.py`
- **Where:** `qcd/running.py:95-96`: `>>> f"{alpha_s(1000.0):.4f}"` → `'0.0884'`.
- **Issue:** actual value is 0.088507 → `'0.0885'` (CRunDec agrees: 0.088507; the test suite itself
  asserts 0.0885 in `tests/test_alpha_s.py:73-75`). Doctests are not collected by pytest config, so
  it never fails. Documentation typo only.

### F5. MINOR — alpha_s(M_Z) anchor inconsistency between modules
- **Where:** `quarkConstraints/qcd_running.py:34` (`_ALPHA_S_MZ_DEFAULT = 0.1179`, PDG 2022) vs
  `qcd/constants.py:8` (`ALPHA_S_MZ = 0.1180`, PDG 2024).
- **Impact:** ≪ the LL truncation error of the (deliberately) 1-loop WC-running module; cosmetic
  consistency issue. MINOR, confidence high.

### F6. Chiral-enhancement quark-mass scales in `deltaf2.py` (cross-slice check) — part documented approximation, part undocumented; consistent-scale API exists but see F2
What the qcd/ side actually delivers, and the quantified effect per system
(`(m_P/(m_q1+m_q2))²` enters LR matrix elements as r/6+1/4 and r/2+1/12):

- **Kaon** (`deltaf2.py:639-654, 719-724`): masses m_s, m_d at 2 GeV; B_1 at 2 GeV (consistent);
  **B_4/B_5 at 3 GeV** — a **documented** caveat with TODO (`deltaf2.py:643-650`). Quantified with
  `run_msbar_mass`: m_s(3 GeV)=0.08443, m_d(3 GeV)=0.004221 → r_chi = 25.75 (2 GeV masses) vs 31.51
  (3 GeV masses); at fixed B(3 GeV) the code's hybrid underestimates the LR MEs by ~17-18%
  (O4 ratio 0.825, O5 0.818). Partially compensated by evaluating WCs at 2 GeV instead of 3 GeV
  (LO C4 grows ~8% from 3→2 GeV), net ~10% underestimate of LR kaon terms ⇒ **anti-conservative**
  direction for eps_K/Delta m_K. Classification: documented approximation (acknowledged TODO), not a
  hidden mistake; recommend closing the TODO.
- **B_d** (`deltaf2.py:667-668`): m_b(m_b)=4.18 + m_d(2 GeV)=0.00467 — mixed scales as flagged by
  the other slice. Consistent m_d(m_b)=0.003944 (computed via `run_msbar_mass(0.00467, 2.0, 4.18, 4)`).
  Effect on O4 ME: 0.02% — numerically irrelevant because m_b dominates the denominator. MINOR.
- **B_s** (`deltaf2.py:681`): m_b(m_b) + m_s(2 GeV)=0.0934 vs consistent m_s(m_b)=0.0789:
  O4 ME ratio 0.9965, O5 0.9939 — ≤0.6%. MINOR (well inside the ~5-10% bag uncertainties).
- **D0** (`deltaf2.py:694-695`): m_c(m_c)=1.27 + m_u(2 GeV)=0.00216. This is the largest hybrid:
  r_chi = 2.149 (code) vs 2.915 with m_c(2 GeV)=1.090 (vs 3.568 at 3 GeV) → LR O4 ME underestimated
  by ~17% (vs 2 GeV) / ~28% (vs 3 GeV), anti-conservative. **Not documented** (unlike the kaon
  caveat), but B_4_D=B_5_D=1.0 are themselves order-of-magnitude "estimated" inputs, so the scale
  choice is within the stated input uncertainty. Severity MINOR; recommend a caveat comment matching
  the kaon one, or fixing all r_chi masses to the bag scale.
- **Consistent-scale call availability:** YES — `qcd.mass_running.run_msbar_mass(m, mu_ref,
  mu_target, n_f_ref)` covers every needed conversion (verified values above against CRunDec to
  ≤6e-4, the residual being F1's d3 issue), and `pdg_quark_masses_at_scale(2.0|3.0,
  flavors=(...no 't'...))` returns a consistent set. The **default** full-flavor
  `pdg_quark_masses_at_scale(2.0)` crashes (F2) and must be fixed before recommending it as the
  deltaf2 remedy.

### F7. Note — test-coverage gaps (no code defect)
- No external-reference test pins the *magnitude* of the 3-loop mass-decoupling step (would have
  caught F1); `tests/test_mass_running.py` only checks identity/continuity/round-trip properties.
- No test exercises `pdg_quark_masses_at_scale` below 163.5 (would have caught F2).

---

## Verified correct (numerical spot-checks run)

1. **Beta coefficients & convention** (`qcd/beta_function.py`): β0 = 11 − 2n_f/3, β1 = 102 − 38n_f/3,
   β2 = 2857/2 − 5033n_f/18 + 325n_f²/54, β3 = van Ritbergen-Vermaseren-Larin — all match literature;
   e.g. n_f=5: (7.6667, 38.6667, 180.9074, 4826.1563) ✓. Expansion variable is a = αs/4π with
   dαs/dln μ² (the classic αs/π slip is absent; README even documents the pitfall explicitly).
2. **4-loop alpha_s with 3-loop decoupling vs CRunDec** (thresholds crossed both up and down,
   matching at mu = m_h with MS-bar masses incl. top 163.5):
   | mu (GeV) | repo | CRunDec | rel |
   |---|---|---|---|
   | 1.0 | 0.476548 | 0.476548 | 5e-9 |
   | 1.777 | 0.319570 | 0.319570 | 3e-9 |
   | 2.0 | 0.301496 | 0.301496 | 3e-9 |
   | 4.18 | 0.224616 | 0.224616 | 2e-9 |
   | 91.1876 | 0.118000 | (anchor) | — |
   | 1000 | 0.088507 | 0.088507 | 4e-10 |
   | 3000 | 0.079664 | 0.079664 | 9e-10 |
   | 10000 | 0.071819 | 0.071819 | 2e-9 |
   | 50000 | 0.063482 | 0.063482 | 3e-9 |
   This also validates the alpha_s decoupling constants c2 = 11/72 and
   c3 = 564731/124416 − (82043/27648)ζ3 − (2633/31104)n_l, the n_f direction at every crossing, the
   up/down series inversion, and threshold application in both directions. (The audit prompt's rough
   expectations 0.087 @ 1 TeV / 0.073 @ 10 TeV were slightly off; CRunDec sides with the code.)
3. **gamma_m convention and 1-4-loop coefficients** (`qcd/mass_running.py`): d ln m/d ln μ² = −γ_m(a),
   a = αs/4π, γ0 = 4, γ1 = 202/3 − 20n_f/9, γ2/γ3 per Chetyrkin 1997. Verified against CRunDec's
   exact-ODE mode (`AsmMSrunexact`, same truncated-series ODE formulation):
   m_b(m_b)=4.18 → m_b(M_Z) = 2.86326001 (repo) vs 2.86326001 (CRunDec), rel 5e-10;
   m_c(m_c)=1.27 → m_c(3 GeV) = 0.98531209 vs 0.98531209, rel 1e-9. The prompt's "m_b(M_Z) ≈ 2.88"
   corresponds to older inputs; with αs(M_Z)=0.1180 and m_b=4.18 the CRunDec-verified value is 2.863.
   (Repo vs CRunDec's *series-solution* `mMS2mMS` differs by 7e-5 — pure 5-loop truncation-scheme
   ambiguity, both legitimate; the exact-ODE comparison is the apples-to-apples one.)
4. **Mass-decoupling d2 = 89/432** ✓ (CRunDec-extracted 0.206019). Pole↔MS-bar top handling:
   thresholds and default M_TOP are MS-bar (163.5); M_TOP_POLE kept for reference only and never
   used in matching ✓ (consistent with the log-free mu=m_h prescription).
5. **`couplings.py` g_s\*** : `alpha_s(M_KK, precision='high')` = 4-loop MS-bar with n_f=6 above
   163.5 ✓; g_s = sqrt(4π·0.079664) = 1.0006 at 3 TeV, CRunDec-verified. The g_s_star override
   semantics (RS-enhanced coupling vs perturbative fallback) are documented and validated.
6. **`quarkConstraints/qcd_running.py` (LO WC running)**: 1-loop β0 = (33−2n_f)/3 ✓; threshold/n_f
   segment assignment correct at boundaries (n_f=6 only above m_t etc.); Landau-pole guard works;
   LO continuous matching is the correct LL prescription. Deliberate LO design is documented
   (α_LO(2 GeV)=0.267 vs 4-loop 0.301 — inherent LL truncation, documented approximation).
7. **ΔF=2 Wilson evolution** (`evolve_deltaf2_wilsons`): γ_VLL = 6(N−1)/N = 4 ✓ (BMU);
   the scalar-basis LR coefficient ADM [[−16,−6],[0,2]] was re-derived by hand from the BMU operator
   ADM γ^(0) = [[2,12],[0,−16]] via the stated map C_BMU = (−C5/2, C4): T⁻¹·γ^(0)T·T with
   T = [[0,−1/2],[1,0]] reproduces [[−16,−6],[0,2]] exactly ✓. Segment-product evolution reproduces
   direct ODE integration of the LO RGE to 1e-10 (3 TeV → 2 GeV: VLL factor 0.72913; U_LR =
   [[3.53816, 0.89476],[0, 0.85389]] in (C4,C5)) — correct LO magnitudes (VLL suppression ~0.73,
   LR C4 enhancement ~3.5, C4 does not feed C5 at LO ✓).
8. **Scale orthogonality bookkeeping** (`scales.py`, `pdg_quark_masses.py`, `benchmarks.py`): mass
   targets at mu = m_t(m_t) = 163.5 in n_f=5/6 as documented; WC reference scale 3 TeV kept separate;
   the 163.5 GeV path verified: u 0.00119, d 0.00259, s 0.05149, c 0.60229, b 2.73198, t 162.42 GeV —
   all consistent with CRunDec up to F1's ≤2.2e-4.
9. **Tests:** 68/68 pass (`test_alpha_s.py`, `test_mass_running.py`, `test_qcd_running.py`,
   `test_pdg_quark_masses.py`); the asserted anchors (0.1084 @ m_t, 0.0885 @ 1 TeV, 0.0797 @ 3 TeV,
   0.0718 @ 10 TeV) are all CRunDec-confirmed.

## Verdict
The alpha_s machinery (4-loop running, 3-loop decoupling, thresholds both directions, conventions)
is **exactly right** — agreement with CRunDec at the 1e-9 level across 1 GeV-50 TeV, and the γ_m
running is equally exact. The two MAJOR findings are a wrong 3-loop mass-decoupling constant
(numerically harmless today, ≤2.2e-4, but a real literature error with a false provenance comment)
and a latent crash/contract violation in the consistent-scale mass helper below m_b — which happens
to be the very API needed to fix the deltaf2 chiral-mass scale hybrids (themselves ≤1% for B_d/B_s,
~17% for D0 and ~10-18% for the kaon LR terms, the kaon one already documented as a TODO).
No BLOCKERs; no downstream physics conclusion changes at current precision.
