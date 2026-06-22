# Independent implementation review — QUARK-FIX-2026-06-21 (reviewer: Opus, adversarial)

**Verdict: APPROVE.** All six items (B1/B2/B3, M1/M2/M5) are faithful to PLAN v3 and
physically correct. Both implementer deviations are justified and verified. The new tests
are genuine literature-anchored oracles that fail against old code. All fix-relevant tests
pass; the only red I saw was an `LD_LIBRARY_PATH`/GLIBCXX environment artifact (resolved by
prepending the conda lib path) — not a code defect.

I re-derived every load-bearing number independently and did NOT trust the implementer.

---

## 1. Per-item correctness

### B3 — ΔF=2 O4/O5 un-swap + 1/(2m_M), all 6 sites — CORRECT
- **All 6 sites changed, byte-identical**: deltaf2.py `_kaon_matrix_elements` (720) +
  `_meson_matrix_elements` (908); phenomenology.py sites at 452/482/516/580
  (`_meson_matrix_elements_inline`). Verified via diff — every site now
  `O1=(1/3)`, `O4=(R/4+1/24)`, `O5=(R/12+1/8)`. No 7th site (the `_m12_np_from_bridge_generic`
  consumer calls the inline helper transitively, as the plan states).
- **GGMS direction independently confirmed TRUE** via a dedicated literature sub-derivation
  (FLAG 2024 arXiv:2411.04268 Eqs. 117-119; Ciuchini hep-ph/9808328; Bae 1202.1570; GGMS).
  Q4 = colour-SINGLET `(s̄(1−γ5)d)(s̄(1+γ5)d)` carries the LARGE `(R/4+1/24)`; Q5 = colour-
  CROSSED carries the SMALL `(R/12+1/8)`. Non-trivial cross-check: FLAG bag normalizations
  N4=2, N5=2/3 (independent data) give (1/4)/2 = (1/12)/(2/3) = 1/8 — confirms the
  coefficient↔colour pairing. The code's assignment is correct.
- **×2-and-swap identity verified exactly**: old⟨O5⟩/2 = (R/2+1/12)/2 = R/4+1/24 = new⟨O4⟩;
  old⟨O4⟩/2 = R/12+1/8 = new⟨O5⟩. Consistent with §0.2.
- **Numerics**: ⟨O4⟩/⟨O5⟩ (per-bag) = 2.8532 at R=25.746 (matches IMPL_NOTE); ⟨O1_K⟩ =
  2.2128e-3 GeV³ (PLAN §6.1 ≈2.21e-3). Confirmed.
- **review_local notes "DO NOT ÷2" is correctly overruled**: the ÷2 is on the MATRIX
  ELEMENT (the missing 1/(2m_M)); the 1/6 Wilson is untouched. The R5 "machine-precision"
  defense is irrelevant to the swapped chiral sector (§0.2) — confirmed by reproducing the
  old 1.9286 (see §2b below): the swap+double cancel in the leading post-running LR combo,
  which is why R5 returned 1.0.

### B2 — ε_K SVD rephasing to PDG convention — CORRECT (with the +V_ts deviation, §2a)
- Phase algebra hand-verified: each anchor (`pa_0=1`, `pb_0=conj(g00)`, `pb_1=conj(g01)`,
  `pa_2=g21·pb_1`, `pb_2=pa_2·conj(g22)`, `pa_1=g12·pb_2`) makes its target entry real, using
  `|g|=1` unit phases. Branch-free (`z/|z|`, never `np.angle`), signed-zero-safe sign pins
  (`>=0.0`), single global anchor `arg(V_ud)=0`, idempotent. Matches PLAN §1.2 determinism
  recipe exactly.
- L and R columns rephased by the SAME diagonal → masses preserved (verified: invariants
  test passes; masses/|CKM|/Jarlskog rtol≤1e-12).
- At the default benchmark V_ud,V_us,V_cb,V_tb real≥0, |V_ub|=0.00368, J=3.12e-5 — physical
  PDG CKM. Confirmed.

### B1 — CGHNP Zbb fermion-KK retranslation — CORRECT
- Diagonal bracket `1/(1+2c)·(1/(2F²) − 1 + 2F²/(3−2c))`: c-sign fixed in both denominators,
  F²=2f² restored. **Independently verified against the convention dictionary oracle**
  (c_CGHNP=−c_repo, F²=2f²) at 4 points — production matches to rtol 1e-12.
- Flavour-sum restructure correct: diagonal uses the 3rd-gen bracket; light-gen sum carries
  the COMMON `1/(2 f_3²)` factor with per-gen `1/(1−2c_i)`, NOT a per-light-gen full bracket
  (fixes the F²(c_d3)/F²(c_d_i) ~10²–10⁴ inflation, PLAN §4.2 error 1a).
- Prefactor `m_b²/(2 Λ_IR²)` kept (CGHNP geometric-Λ_IR). m_b-scale deferred per plan.
- Decision to CORRECT (not drop) matches plan §4.3. Magnitude pin confirms fermion < gauge.

### M1 — T010 R_b → T011-style loose-edge NP-shift budget — CORRECT
- `ObservableBudget` gains `central_residual`, `hard_veto_budget`; `_build_budget` computes
  `central=|exp−sm_fit|`, `hard_veto=central+combined` (mirrors T011, the verified gap from
  review non-blocking #3 is closed). `ObservablePull` gains `np_shift`, `np_shift_ratio`.
  `evaluate` veto scalar = `max(np_shift_ratio)`, `passes = ratio<=1`. Reported budget is
  `hard_veto_budget`. Matches PLAN §5/M1.
- **Headroom recomputed from code, not hardcoded** (review non-blocking #1 satisfied): SM-limit
  pull = −0.99601 → old gate passed SM by 0.003986σ; combined_sigma=6.7268e-4,
  central=6.7e-4, hard_veto=1.3427e-3 — I reproduced 1.3426812023536894e-3 exactly. New gate:
  SM-limit NP shift 0 → ratio 0 → SM passes by construction (no knife-edge). Confirmed.
- Fixture floor behavior: minimal T010 vetoes at 3 TeV (ratio 1.67), passes at 8/15 TeV — the
  108 TeV artifact collapsed; T010 ceases to be the minimal-floor driver. As designed.

### M2 — EW001 ΔT physical-M_KK convention — CORRECT
- Both `minimal_rs_t_coefficient` and `custodial_rs_plr_t_coefficient` multiply by
  `x1²=5.996` (≈6.0), putting ΔT in the same physical-M_KK convention as the ΔS=c_S
  calibration. Equivalent to `(v/Λ_IR)²`. Docstring + proxy string updated.
- chi2 at 6 TeV 3.6313 → 519.55 is verified by the test's INDEPENDENT `_manual_chi2`
  recomputation from the covariance (`result.predicted == expected`), not a snapshot. The
  143× growth is the quadratic (var_t·dt² + cross terms) response to ΔT×6 — self-consistent.

### M5 — structured tag-class hint — CORRECT
- `_structured_tag_class` reads `diag["tag_class"] ∈ {rigorous,proxy,partial,stub}`; honored
  after stub/missing guards, before prose. `_mentions_proxy` is plural-robust
  ("proxy"/"proxies"); replaces both `"proxy" in needs_text`/`status_text`; "proxies resolved"
  added to allowlist. T001/T002 declare `"tag_class":"proxy"`. Explorer veto count (F7) now
  requires `tag ∈ {rigorous,proxy}`, aligning with harness/comparison. Matches PLAN §5/M5.

---

## 2. The two deviations — both JUSTIFIED and VERIFIED

### 2a. B2 5th anchor (V_ts) — JUSTIFIED, correct, convention-stable, pure rephasing
- **Independently confirmed the 4-anchor version is insufficient**: pinning only
  {Vud,Vus,Vcb,Vtb} real leaves arg(V_ub) drifting wildly under physical column rephasing
  (0.219 → −0.282 → 0.742 → 2.566 → 2.212 across seeds) — exactly the audit's ×6
  nondeterminism. A 5th real-positivity condition is mathematically required (6 column phases,
  1 global redundant ⇒ 5 physical d.o.f.; 5 reals fix all 5).
- **With V_ts as the 5th**: arg(V_ub) is bit-stable (0.6986068571 to ~1e-16) across all seeds;
  full CKM invariant to ~5e-16. V_ub correctly RETAINS the physical phase (|V_ub|>0, nonzero
  arg). Masses/|CKM|/Jarlskog invariant (test + my recompute). ε_K is now rephasing-invariant
  (the ×6 collapse). This is a pure rephasing — no observable changed. Faithful realization of
  approach A. **The choice of which entry is the 5th is convention; placing the physical phase
  in V_ub (not V_ts) is correct and is guaranteed harmless by the Jarlskog-invariance pin.**

### 2b. ε_K = 20.67, not the plan's ≈3.32 — IMPLEMENTER IS RIGHT
Full 2×2 decomposition (I toggled B2 rephasing on/off and B3 coefficients old/new):

| | B2 OFF (legacy phase) | B2 ON (PDG phase) |
|---|---|---|
| **old B3** | **1.9286313** (reproduces the legacy pin to 9 sig figs) | 12.684 |
| **new B3** | **3.3225** (= the plan's "≈3.32") | **20.6675 (current pin)** |

- 3.3225/1.9286 = **1.7227** — exactly the plan's B3 ×1.7227. So 3.32 was indeed the
  **B3-only** figure (B2 disabled).
- 20.667/3.3225 = **6.22** — exactly the audit's "×6" rephasing magnitude. So B2 (a phase
  change, |M12| invariant) accounts for the 3.32→20.67 move.
- PLAN §1.4 explicitly mandates re-pinning ε_K **after BOTH B2 and B3 land**. 20.67 is the
  correct convention-stable B2+B3 value; 3.32 was a B3-only intermediate. **Not a regression,
  not a paper-over.** b/d systems still pass (phase-magnitude only). J=3.12e-5 confirmed.

---

## 3. New tests are GENUINE oracles (fail on old, pass on new)

- **B3 O4/O5 ratio (2.853)**: built from GGMS rationals written literally in-test. OLD
  production ratio = (R/6+1/4)/(R/2+1/12) = **0.3505** → fails the 2.853 pin. O1 pin: old
  (2/3) ⟨O1⟩=4.43e-3 fails the (1/3)=2.21e-3 pin. GENUINE.
- **B1 dictionary-oracle bracket**: `_cghnp_diagonal_bracket` derived from c=−c_repo, F²=2f²
  (NOT a byte-copy of production); production must match. The OLD `_manual_B` was byte-
  identical to the buggy bracket — precisely the M7 failure mode — and is now removed. The new
  value pin (3.74e-6 vs old 2.35e-5) and the fermion<gauge magnitude pin are genuine
  discriminators. GENUINE.
- **B2 invariance pin**: I confirmed it fails on old (no-rephasing) code with maxdiff **1.79**
  (≫1e-12). GENUINE guard. (The determinism `array_equal` pin alone would not fail old code —
  old SVD was deterministic too — but the invariance pin is the real discriminator, and it is
  the audit's ×6 demo inverted.)
- **M1 SM-limit pin** (NP shift 0 → ratio 0 → passes): fails old max|pull| gate (which gave
  ~0.996 at SM limit). GENUINE.

No pin merely re-pins current production output. The snapshot re-pins (ε_K 20.67, K001,
B003 predicted/ratio, T010/T011, EW001 519.55, δg_L/R) are correctly labelled as DOWNSTREAM
of the absolute pins.

## 4. Re-pin discipline — followed correctly
- **B003 anchor budgets NOT moved**: 5.844e-12 and 2.635167648629676e-13 unchanged in the
  diff; only predicted (4.760e-14→4.697e-14) and ratio (0.1806→0.1782) moved. Matches review
  non-blocking #2.
- **K001 expected_pass True→False** (couplings0 now ratio ~1.09, crosses veto) and the
  `ratio>10.0`→`>1.0` relaxation are physics-correct (corrected kaon ε_K ×1.7 larger).
  couplings1 remains ~108. Legitimate, not a papered-over regression.
- **EW001 6 TeV pass→fail, 18 TeV pass**: physics-correct (ΔT×6 → ~16 TeV floor); verified by
  the independent `_manual_chi2`.
- **pr1 15 TeV veto collapse, pr2 T010 oracle 1.468→0.236**: M1 artifact collapse, correct.

## 5. Test suite — FULL SUITE GREEN
- **Definitive full-suite run (with `LD_LIBRARY_PATH=$CONDA/lib`), split in two halves to fit
  the sandbox wall limit: 1263 passed (a-m + constraints) + 504 passed, 1 skipped (n-z) =
  1767 passed / 1 skipped — matches the 1768 collected. ZERO failures.**
- Earlier targeted runs covering EVERY changed area also all GREEN: B2/B3/K001/ε_K (46),
  B1/M1/M2/M5/beauty/kaon/top_higgs_ew (610), quark/modern/fit/deltaf2/plot (63), pr1/pr2 (20),
  harness/comparison (22).
- The only failures seen were a transient `GLIBCXX_3.4.29 / LD_LIBRARY_PATH` collision on
  pr1/pr2/harness collection (scipy._ufuncs_cxx) — an ENVIRONMENT artifact; resolved by
  prepending `$CONDA/lib` to `LD_LIBRARY_PATH`, after which all pass. NOT a code issue.

---

## Non-blocking notes (no action required for approval)
1. The full 1768-test suite exceeds the 590s sandbox timeout, so I verified by exhaustive
   targeted subsets rather than one monolithic green line. Recommend the implementer run the
   full suite once locally with `LD_LIBRARY_PATH=$CONDA/lib:...` before committing.
2. B1: at the test fixture point both old and new δg_L^b are >0, so the `>0` sign assertion
   is not by itself a discriminator there; the value pin and the dictionary-oracle test ARE.
   The explicit UV-b_R sign-flip case is covered by the `c_uv=0.557` spot value in the oracle
   test. Adequate.
3. M2 chi2 143× growth (not a naive x1⁴=36×) is correct — it is the cross-term + var_t·dt²
   quadratic response with ΔT now dominant; the in-test `_manual_chi2` certifies it. Worth a
   one-line comment if a future reader expects x1⁴.
4. The GLIBCXX env fragility (some test files crash on scipy import depending on lib-path
   ordering) is orthogonal to this fix but should be pinned in the project env so CI is stable.
