# Slice 1 audit report — ΔF=2 matching, RG, ε_K

Auditor: Claude (skeptical theoretical-physics code review), 2026-06-10.
Scope: `quarkConstraints/deltaf2.py`, `quarkConstraints/qcd_running.py`,
`tests/test_epsilon_k_physics.py`, `tests/test_quark_deltaf2.py`; context docs
(`docs/STATE_OF_PROJECT.md`, `docs/quark_scan_methodology_note.tex`,
`docs/audits/wilson_rg_reference_values.md`, `docs/audits/bag_param_inventory.md`).
Mode: read-only on repo code; numerical checks run with in-repo modules plus
independent closed-form re-derivations. Repo HEAD: 5328a22.

Builds on the prior slice's findings (SVD phase-convention ambiguity in
`fit.py:245-258`; tree-matching pattern verified vs CFW 0804.1954). I agree
with both prior conclusions; the C4/C5 relative sign (C4/C5 = −3) was
independently re-derived here via the color-Fierz decomposition of octet
exchange and confirmed.

---

## Executive summary

The Wilson matching (given), the RG evolution machinery, and the ε_K master
formula are all correct. The **hadronic matrix elements are not**. The
functions `_kaon_matrix_elements()` and `_meson_matrix_elements()` implement,
numerically exactly:

```
ME_code(O1) = 2 × ME_correct(O1)
ME_code(O4) = 2 × ME_correct(O5)   ← swapped
ME_code(O5) = 2 × ME_correct(O4)   ← swapped
```

i.e. a global missing 1/(2 m_M) state-normalization factor **and** a swap of
the chirally-enhanced coefficients between the two LR operators. For the
kaon (chiral factor R ≈ 25.7) the swap dominates and the two errors do not
cancel: at the pinned regression benchmark, the corrected ε_K^NP is
**×1.72 larger** than the code value, implying ε_K-driven M_KK floors are
currently **~31% too LOW (anti-conservative)**. This directly moves the
paper's headline numbers (47.26 TeV p50 → ≈62 TeV; CFW-matched 23.37 TeV →
≈30.7 TeV; the "factor 2.2 stronger than CFW" claim → ≈2.9).

Finding count: 1 BLOCKER, 1 MAJOR, 2 MINOR.

---

## Finding F1 — BLOCKER: O4/O5 chiral coefficients swapped in all hadronic matrix elements

**Location:**
- `quarkConstraints/deltaf2.py:719-724` (`_kaon_matrix_elements`)
- `quarkConstraints/deltaf2.py:908-913` (`_meson_matrix_elements`)
- Vendored copies: `quarkConstraints/modern/phenomenology.py:450-454, 481-484,
  515-518, 577-580` (kept in sync by `tests/test_quark_deltaf2.py:120-155`,
  so the bug is enforced-consistent across lanes).

**Code expression:**
```python
o4_lr = (m_ratio_sq * (1.0 / 6.0) + 1.0 / 4.0) * fk2_mk * B_4_K_3GEV
o5_lr = (m_ratio_sq * (1.0 / 2.0) + 1.0 / 12.0) * fk2_mk * B_5_K_3GEV
```
with `m_ratio_sq = R = (m_K/(m_s+m_d))² ≈ 25.75`, `fk2_mk = f_K² m_K`.

**Correct result and why.** The code's own docstrings (deltaf2.py:711-717),
the FLAG B4/B5 bag definitions, the CFW-verified matching
(C4 = −g_L g_R/M², C5 = +g_L g_R/(3M²)), and the BMU-mapped ADM in
`qcd_running.py` all define the LR operators in the conventional SUSY/GGMS
basis:

```
O4 = (s̄_α P_L d_α)(s̄_β P_R d_β)   (color singlet-singlet, strongly chirally enhanced)
O5 = (s̄_α P_L d_β)(s̄_β P_R d_α)   (color crossed, weakly enhanced)
```

The standard VSA/bag matrix elements in this basis, in the M12-ready
normalization ⟨O⟩/(2m_K) consistent with `M12 = Σ C_i ⟨O_i⟩` as used in
`_compute_m12_np` (Gabbiani-Gabrielli-Masiero-Silvestrini hep-ph/9604387;
Ciuchini et al. hep-ph/9808328; CFW 0804.1954 eq. (23); FLAG normalizations
reduce to these), are:

```
⟨O1⟩ = (1/3)              f_K² m_K B1
⟨O4⟩ = (R/4  + 1/24)      f_K² m_K B4     ← large coefficient belongs to O4
⟨O5⟩ = (R/12 + 1/8 )      f_K² m_K B5     ← small coefficient belongs to O5
```

The code instead assigns coefficient (R/6 + 1/4) to O4 and (R/2 + 1/12) to
O5. Verified numerically (exact identities):
`code_O4_coeff = 2 × correct_O5_coeff` and `code_O5_coeff = 2 × correct_O4_coeff`.
Note (R/6 + 1/4) is exactly |⟨Q1_LR^BMU⟩| = |−2⟨O5⟩| — the magnitude of the
Fierzed BMU *vector*-LR matrix element. The most plausible origin is that the
matrix elements were transcribed from a BMU-basis table while the Wilson
coefficients, bags, and ADM live in the scalar O4/O5 basis. No consistent
reinterpretation rescues the code: if c4/c5 were secretly BMU coefficients,
the matching pattern would have to be (±g g/(6M²), ±g g/M²), not
(−g g/M², +g g/(3M²)), and the RG would have to use the unmapped BMU ADM.

**Consequences.** The dominant, RG-enhanced coefficient C4 is paired with
the small matrix element and vice versa, and the C4/C5 relative sign then
produces a spurious partial cancellation: at the fitted benchmark the code's
C5 term cancels 19% of the C4 term post-RG, whereas with correct MEs the C5
term is +2.4%. Pre-RG, code underestimates the kaon LR combination
|⟨O4⟩ − ⟨O5⟩/3| by ×2.4 (1.117 vs 5.33 in units of f²m·bag, before the
F2 factor 2 partially compensates).

**Numerical impact on M_KK floors** (computed with the repo's own RG +
fitter; identical Wilson evolution in both columns):

| Quantity (fitted benchmark r=0.25, M_KK=3 TeV) | code | corrected MEs |
|---|---:|---:|
| ε_K ratio_to_budget | 1.9286 (pinned in tests) | 3.3225 |
| ε_K^NP scaling | 1 | ×1.7227 |
| implied ε_K M_KK floor | 1 | ×1.3125 |

Generic anarchic-phase points give ×1.72–1.75 whenever LR dominates Im M12
(LL-only points instead go *down* ×2, F2 alone). Headline consequences
(ε_K is the binding constraint per the methodology note): RUNA p50
47.26 TeV → ≈62 TeV; p95 127.13 TeV → ≈167 TeV; CFW-matched 23.37 TeV →
≈30.7 TeV; the documented "factor 2.2 stronger than CFW" → ≈2.9.
B_d/B_s: R ≈ 1.6, the two errors nearly cancel (code/correct ≈ 1.02–1.08,
floors ≈ 1–4% too strong). D0: code/correct ≈ 0.89–0.93 (floor ≈ 4–6% too
weak). Direction for the kaon is **anti-conservative** — published floors are
understated.

**Severity:** BLOCKER (moves the headline M_KK floor by ~31% in the
unsafe direction; propagates to catalog constraints K001/K002, B001-B004,
C001/C002 via `flavor_catalog_constraints/physics_adapters/deltaf2.py` and to
the modern lane via the vendored copies).

**Confidence:** High. Anchored by (i) three independent literature
conventions (GGMS, Ciuchini et al., CFW eq. 23; FLAG bag normalization),
(ii) the exact ×2-and-swap numerical identity, (iii) internal consistency:
the code's own RG and matching unambiguously fix the basis the MEs must
live in.

---

## Finding F2 — MAJOR: missing 1/(2 m_M) — all matrix elements are 2× the M12-ready values

**Location:** same functions as F1 (`deltaf2.py:722`, `deltaf2.py:910`); the
2m_M division is applied **zero** times anywhere in the M12 chain
(`_compute_m12_np`, `deltaf2.py:733-741`; `compute_m12_np`,
`deltaf2.py:916-935`).

**Code expression:** `o1_vll = (2.0 / 3.0) * fk2_mk * B_1_K` combined with
`M12 = Σ C_i ⟨O_i⟩` directly.

**Correct result and why.** With relativistic state normalization,
⟨K̄0|(s̄γ_μ P_L d)²|K0⟩ = (2/3) f_K² **m_K²** B_K (textbook; (8/3) f²m²B̂ for
the (V−A)⊗(V−A) normalization), and M12 = ⟨K̄|H|K⟩/(2 m_K), giving the
M12-ready value (1/3) f_K² m_K B_K. Structural cross-check performed: with
the (1/3) convention the SM box reproduces the textbook
M12 = (G_F²/12π²) f_K² B̂_K m_K M_W² λ_t² η S0 identically (verified to
1e-16 numerically); with the code's (2/3) it is 2× too large. The
single-power-of-m_K form the code uses is only correct *after* the 2m_K
division — the code has relativistic-magnitude coefficients with
M12-ready power counting.

**This contradicts a documented audit claim.** `docs/STATE_OF_PROJECT.md:29`
states the ε_K normalization "is audited and tested with no extra factor of
two ... Do not 'fix' epsilon_K by multiplying by 2." The warning's direction
is right but the conclusion is wrong: the fix is to *divide* by 2 (or use the
GGMS coefficients). `docs/audits/bag_param_inventory.md:57-60` explicitly
says the prefactors `2/3, 1/6, 1/4, 1/2, 1/12` were "inspected but not
changed" and defers the basis-normalization audit to hole #6, whose artifact
(`docs/audits/wilson_rg_reference_values.md`) only audited the RG evolution.
The matrix-element normalization audit therefore never actually happened.

**Numerical impact:** in isolation, ×2 on every |M12^NP| and ε_K^NP →
M_KK floors ×√2 too *high* (overconservative). In practice it is entangled
with F1; net per-system impacts are given under F1. Fixing only F2 without
F1 (or vice versa) would make the kaon numbers worse, not better — the two
must be corrected together.

**Severity:** MAJOR. **Confidence:** High.

---

## Finding F3 — MINOR: mixed-scale quark masses inside the chiral enhancement factors

**Location:** `deltaf2.py:667-668, 681, 694-695` feeding
`_meson_matrix_elements` (`deltaf2.py:909`).

**Issue:** R_χ = (m_P/(m_q1+m_q2))² is built from m_b(m_b) = 4.18 GeV
combined with m_d(2 GeV)/m_s(2 GeV) for B_d/B_s, and m_c(m_c) = 1.27 GeV
with m_u(2 GeV) for D0. The chiral factor should use both masses at the same
scale, and that scale should match the bag-parameter scale (m_b for the
FLAG B-meson bags, 3 GeV for ETM-style D bags). Using m_b(2 GeV) ≈ 4.9 GeV
vs m_b(m_b) changes R_χ(B) by ~25%; for D the m_c scale choice shifts R_χ
by ~30%. The kaon case is internally consistent at 2 GeV (modulo the
documented 3-GeV-bags caveat). Distinct from the *documented* bag-endpoint
caveat (which concerns the Wilson-coefficient scale): the quark-mass scale
pairing inside R_χ is not mentioned anywhere.

**Impact:** LR matrix elements for B/D shift by O(10-30%); B floors are
subdominant and D bags are order-one estimates anyway, so M_KK floor impact
is a few percent. **Severity:** MINOR. **Confidence:** High on the
inconsistency, medium on assigning it "undocumented" (the audit-doc bag
table touches scales but only for bags, not masses).

---

## Finding F4 — MINOR: regression tests pin the buggy values; conventions are untested

**Location:** `tests/test_quark_deltaf2.py:164`
(`assert np.isclose(epsilon_k_ratio, 1.9286313761001348)`),
`tests/test_epsilon_k_physics.py:182-198, 358-383`.

**Issue:** the test suite checks dimensional analysis, CP structure, M_KK
scaling, and "LR is chirally enhanced" — all of which pass equally with the
swapped/doubled matrix elements. The only sharp numbers pinned are
regression snapshots of the buggy output and the bag *inputs*. There is no
test of the VSA coefficients against an external convention (e.g.
⟨O4⟩/⟨O5⟩ → (R/4+1/24)/(R/12+1/8) ≈ 2.85, or the SM-box structural anchor
used in this audit). Additionally the ME formulas exist in four copies
(deltaf2.py ×2, modern/phenomenology.py ×3 call sites) and a sync test
enforces the constants only — any fix must touch all copies and the pinned
regression numbers simultaneously.

**Severity:** MINOR (audit-process). **Confidence:** High.

---

## Disagreements with prior conclusions / docs

- `docs/STATE_OF_PROJECT.md:29` ("no extra factor of two", "Do not fix
  epsilon_K by multiplying by 2") — **disagree** as detailed in F2.
- `docs/quark_scan_methodology_note.tex` headline attribution (factor 2.2 vs
  CFW from BGS budget + FLAG bags + LO sign convention) — incomplete: with
  corrected MEs the factor becomes ≈2.9, and the current 2.2 partially
  reflects F1+F2 rather than physics inputs.
- Prior slice's two conclusions (phase-convention BLOCKER; CFW matching
  pattern correct): **agree**. The matching relative sign C4/C5 = −3 was
  independently re-derived here via T^a⊗T^a color decomposition + Dirac
  Fierz.

---

## Verified correct (explicitly checked)

1. **BMU→scalar LR ADM map** (`qcd_running.py:54-57`): with
   C_BMU = (−C5/2, C4) and BMU coefficient ADM [[2,0],[12,−16]], the
   similarity transform T⁻¹γT was recomputed analytically and equals the
   code's [[−16,−6],[0,2]]. Mixing direction (C5 feeds C4; C4 never feeds
   C5 at LO) matches BMU hep-ph/0005183.
2. **γ_VLL = 4 and evolution direction** (`qcd_running.py:48`): exponent
   γ/(2β0) with (α_high/α_low) base reproduces the RGI B̂_K exponent
   α_s^{−2/9} (nf=3); C1 is suppressed running down — factor 0.7291
   (3 TeV → 2 GeV), matching `docs/audits/wilson_rg_reference_values.md`
   to 16 digits when re-run.
3. **LR evolution numerics**: C4 unit → 3.5382, C5 unit → (0.8948, 0.8539)
   for 3 TeV → 2 GeV; closed-form check (1/3)(η_{−16} − η_{+2}) = 0.8948
   exact. The LO kaon C4 enhancement ≈3.5 sits at the low edge of the
   literature O(4-8) band — consistent with LO vs NLO magic numbers, and
   LO-only running is a documented approximation.
4. **Threshold handling**: nf per segment (6/5/4 for 3 TeV→m_t→m_b→2 GeV),
   segment-matrix ordering (lowest-applied-last), continuous α_s matching,
   charm threshold correctly not crossed for μ_low = 2 GeV; α_s(2 GeV) =
   0.2672, α_s(3 TeV) = 0.0804 (one-loop, documented).
5. **ε_K master formula** (`deltaf2.py:789`): ε_K^NP = κ_ε|Im M12|/(√2 Δm_K)
   with κ_ε = 0.94 (Buras-Guadagnoli) and Δm_K = 3.484e-15 GeV — structure
   and prefactor (1.908e14 GeV⁻¹) correct; |·| applied once; budget
   |ε_exp − ε_SM| = 6.7e-5 (BGS 2020) is the documented strict central
   choice with a documented override knob (`deltaf2.py:790-798`,
   STATE_OF_PROJECT.md:126).
6. **Running direction and single application**: matching at μ = M_KK,
   one evolution M_KK → μ_had = 2 GeV, matrix elements evaluated once;
   kaon B1 (2 GeV) is correctly paired with 2 GeV Wilsons.
7. **2m normalization bookkeeping intent**: M12 = ⟨H⟩/(2m) is *intended*
   exactly once (single-power-of-m formulas) — the error is in the
   coefficient values (F2), not double counting.
8. **D0 inputs**: Δm_D = 6.25e-15 GeV ↔ x ≈ 0.4% (HFLAV) consistent;
   budget Δm_exp/2 with long-distance-dominated SM is documented; scalar
   bags 1.0 documented as estimates (ETM 2015 cross-check in
   bag_param_inventory.md gives 0.91/0.97 — code choice conservative).
9. **Budget constructions** (`deltaf2.py:949-961`): B_d/B_s
   max(Δm_exp/2, |Δm_exp − Δm_SM|/2) evaluates to Δm_exp/2 — the
   incoherent "|NP| ≤ observed amplitude" criterion is explicitly documented
   in the methodology note (§ "What the FCNC ratio is", eq. ratio-mass) and
   is internally consistent with using |M12^NP| (no SM phase needed).
   CPV in B/D (sin2β, φ_s, D0 |q/p|) is intentionally handled by separate
   catalog constraints B002/B004/C001, not by this module.
10. **f_K = 155.7 MeV convention** consistently the PDG (√2-free)
    normalization expected by the (8/3)/(1/3) formulas; m_s+m_d at 2 GeV
    for the kaon chiral factor.
11. **Documented approximations confirmed as documented** (not findings):
    LO-only running; kaon B4/B5 FLAG 3 GeV used at 2 GeV
    (`deltaf2.py:643-654`); B_d/B_s bags at m_b vs 2 GeV Wilson endpoint
    (documented in `docs/audits/wilson_rg_reference_values.md` endpoint
    table); legacy operator-weight fallback non-default.

---

## Recommended fix (for the orchestrator; no code was modified)

In `_kaon_matrix_elements` and `_meson_matrix_elements` (and the three
vendored copies in `modern/phenomenology.py`):

```python
o1_vll = (1.0 / 3.0) * fk2_mk * B_1
o4_lr  = (m_ratio_sq / 4.0  + 1.0 / 24.0) * fk2_mk * B_4
o5_lr  = (m_ratio_sq / 12.0 + 1.0 / 8.0 ) * fk2_mk * B_5
```

then re-baseline the pinned regression values (ε_K benchmark 1.9286 → 3.3225)
and re-run the RUNA floor scan (expect ε_K-driven floors ×≈1.31). Update
STATE_OF_PROJECT.md:29 and the CFW-comparison paragraphs. Consider adding a
convention test: ⟨O4⟩/⟨O5⟩ = (R/4+1/24)/(R/12+1/8) and the SM-box structural
anchor.
