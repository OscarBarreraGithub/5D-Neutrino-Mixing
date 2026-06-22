# Independent adversarial review of PLAN.md (QUARK-FIX-2026-06-21)

Reviewer: Opus (independent dual-signoff lane). Read-only on production code.
Repo HEAD context: 5328a22. Method: re-derived the disputed physics from first
principles + primary sources, verified every claimed code location against the
actual files, and reproduced the key numbers end-to-end through the repo's own
pipeline (no production code modified; B3 fix tested via in-process monkeypatch).

**VERDICT: APPROVE** (with non-blocking refinements; no BLOCKING issues found.)

The plan's three fix directions (B1 sign-flip + flavour-index, B2 PDG rephasing,
B3 un-swap + ÷2) are each independently confirmed. The single most contestable
claim — the B3 direction and the "R5 returned 1.0000 but is irrelevant"
resolution — survives adversarial re-derivation and is, if anything,
*understated* by the plan.

---

## 1. B3 direction (highest stakes) — CONFIRMED, independently re-derived

### 1.1 The corrected coefficients are literature-universal

I re-derived the SUSY/GGMS-basis VSA matrix elements three independent ways:

- **First-principles colour combinatorics.** With the code's own operator
  definitions (deltaf2.py:710-711, confirmed): O4 = colour-singlet
  `(s̄_α P_L d_α)(s̄_β P_R d_β)`, O5 = colour-crossed. The chirally-enhanced
  pseudoscalar-density² (R) term saturates the SINGLET O4 with the full N_c
  colour factor; the crossed O5 gets it reduced by exactly 1/N_c. Check:
  R-term ratio O5/O4 = 1/12 ÷ 1/4 = 1/3 = 1/N_c (✓); const-term ratio
  O5/O4 = 1/8 ÷ 1/24 = 3 = N_c (✓). This is the textbook colour pattern and it
  *forces* the large R-coefficient onto the singlet.
- **Ciuchini/Gabbiani B-normalization constants** (hep-ph/9808328, confirmed via
  literature search): the VSA B-normalization constants c4=2, c5=2/3 reproduce
  exactly `⟨O4⟩ = (1/24 + R/4) f²m B4`, `⟨O5⟩ = (1/8 + R/12) f²m B5`.
- **Cross-reference** to FLAG reviews, Bona/UTfit (hep-ph/0606167), Buras-Misiak-
  Urban — these coefficients are universal in the literature.

Result: the corrected `⟨O4⟩ = (R/4 + 1/24)`, `⟨O5⟩ = (R/12 + 1/8)`, `⟨O1⟩=(1/3)`
**match the plan §2.2 exactly**. O4/O5 = 2.853 at R=25.75 (plan says ≈2.85 ✓).
The code's `(R/6+1/4)` on O4 and `(R/2+1/12)` on O5 are exactly **backwards and
×2** — verified to the rational digit: `code O4 = 2·correct O5`,
`code O5 = 2·correct O4`, `code O1 = 2·correct O1`.

If a reviewer could cite a standard reference where the singlet carries the SMALL
R-coefficient, the fix would invert — I actively looked and found the opposite in
every source. The fix direction is solid.

### 1.2 The "R5 = 1.0000" defense is genuinely irrelevant to the O4/O5 split — CONFIRMED EXACTLY

This is the plan's load-bearing adjudication (§0.2). I verified it symbolically.
With the matched Wilsons C4 = −gg/M², C5 = +gg/(3M²) (deltaf2.py:342-343,
confirmed) and gL=gR, the *leading chiral (R) term* of M12 pre-running is:

- code coefficients (O4 R-coeff 1/6, O5 R-coeff 1/2): **exactly 0**.
- correct coefficients (1/4, 1/12): −2/9 · gg/M².

More strikingly, the FULL pre-running M12 (general complex gL,gR):
- **correct**: (gL² + gR² − 4R·gL·gR)/(18M²)  — R-term present and dominant.
- **code**:    (gL − gR)²/(9M²)               — **the entire R-term has cancelled**.

So the swap doesn't just mis-weight the chiral enhancement; pre-running it
*annihilates* it, leaving only a non-enhanced (gL−gR)² remnant. The chiral piece
the code reports is resurrected solely by C5→C4 RG mixing. This is a *stronger*
statement than the plan makes and fully explains why a check that didn't isolate
the O4/O5 chiral split (the VLL-sector R5 comparison) could read 1.0000 while the
LR sector was swapped. The plan's resolution of the audit-vs-note contradiction is
correct.

### 1.3 End-to-end magnitude — REPRODUCED through the real RG pipeline

I ran the repo's own `evaluate_delta_f2_constraints(_fit_result(0.25),
M_KK=3000)` and applied the B3 fix via monkeypatch (no edit):
- baseline ε_K ratio = 1.9286313765 (matches the pinned 1.9286313761…)
- fixed ε_K ratio = **3.3224626535** → plan's ≈3.3225 (✓, exact)
- ratio = **×1.7227048642** → plan's ×1.7227 (✓, exact); floor ×√1.7227 = ×1.3125 (✓)

### 1.4 All SIX matrix-element sites located and verified byte-identical

| # | File:lines | Function | Verified |
|---|---|---|---|
| 1 | deltaf2.py:720-722 | `_kaon_matrix_elements` | ✓ (R/6+1/4),(R/2+1/12),(2/3) |
| 2 | deltaf2.py:908-910 | `_meson_matrix_elements` | ✓ same |
| 3 | modern/phenomenology.py:452-454 | `_kaon_m12_np_from_bridge_match` | ✓ same |
| 4 | modern/phenomenology.py:482-484 | `_evaluate_epsilon_k_from_bridge` | ✓ same |
| 5 | modern/phenomenology.py:516-518 | `_evaluate_delta_mk_from_bridge` | ✓ same |
| 6 | modern/phenomenology.py:578-580 | `_meson_matrix_elements_inline` | ✓ same |

The plan's site table (§2.1) is exact, including line numbers. All six are
byte-identical algebra; the "must change in one commit" requirement is sound.
NOTE the plan's §2.1 labels site 6 as `_meson_matrix_elements` generic — the
actual function name is `_meson_matrix_elements_inline` (site 2 is the generic
`_meson_matrix_elements`). Cosmetic; both are correctly located. There is also a
3rd consumer of the generic helper, `_m12_np_from_bridge_generic` (phenomenology.py
~584), but it *calls* `_meson_matrix_elements_inline` rather than re-typing the
algebra, so it is fixed transitively — no 7th site.

---

## 2. B1 Zbb (correct-vs-drop) — CONFIRMED, "correct it" is the right call

### 2.1 Convention dictionary verified exactly

I checked `F²_CGHNP(c_repo) = 2 f_IR,repo²` numerically at the repo's ε for
c_r ∈ {0.45, 0.557, 0.601, 0.649}: ratio = 1.00000000 in every case (the
identity `(1−2c)/(1−ε^{1−2c}) = 2·(½−c)/(1−ε^{1−2c})` is exact). The
`c_CGHNP = −c_repo` map and the missing factor-2 in F² are real.

### 2.2 The sign flip is real for the scan-typical regime

For UV-localized b_R (c_r > ½, scan-typical c_d3≈0.557), the repo bracket
`1/(1−2c)·(…)` has **1−2c < 0 ⇒ negative bracket**, while the corrected CGHNP
bracket `1/(1+2c)·(…)` is **positive**. Verified at c_r=0.557/0.601/0.649: repo
= −9.10e3/−6.97e4/−1.01e6 vs corrected = +2.45e2/+3.20e3/+6.55e4. The sign flip
of δg_L^b and the ~190× magnitude inflation (audit table) are reproduced in
sign and order of magnitude. The flavour-sum index bug (957-971 applies the full
per-light-gen `B(c_d_i,F_d_i)` instead of the common `1/F²(c_b3)`) is confirmed
by reading `np.dot(row_ratio, profile_b_d)` with `profile_b_d =
_casagrande_zbb_B_profile_triplet(c_d, f_d)`.

### 2.3 "Correct, don't drop" is defensible

The corrected term is a genuine ~15% correction to the independently-verified
gauge piece (slice-3 V1). Dropping would still require computing the corrected
magnitude to state any bound, so correcting is strictly more informative. The
plan's fallback (zero + metadata bound + NEEDS-HUMAN tag) is a reasonable
contingency. APPROVE the recommendation to correct.

---

## 3. B2 ε_K rephasing & determinism — APPROVE approach A; one determinism caveat

- Location verified: fit.py:227-231 (`_ordered_dirac_svd`), 245-258
  (`mass_matrix_observables`), raw `numpy.linalg.svd`, `ckm = U_L_u.conj().T @
  U_L_d`, no phase fixing. Exactly as claimed.
- The convention dependence is real: Im(M12) is not rephasing-invariant; |M12|,
  masses, |CKM|, J are. The ×6 demonstration is structurally sound.
- Approach A (PDG rephasing) over B (rephasing-invariant ε_K): I agree A is
  preferable because it also fixes the phase-sensitive gates B002/B004/C002 in
  one place and matches the convention in which EPSILON_K_SM (BGS) is defined.

**Determinism — the one place the plan is under-specified (non-blocking).** The
plan correctly flags determinism as load-bearing but does not pin the exact
recipe. My analysis: the dangerous nondeterminism source is NOT degenerate SVD
subspaces — quark singular values are strongly hierarchical, so SVD columns are
unique up to phase. The real hazards are (i) `np.angle` on a *signed* zero:
`np.angle(-0.0+0j) = π`, not 0, so a pivot real-part that can flip sign at the
±0 boundary is a branch hazard; and (ii) the global-phase fix. The plan must
require: pin each rephased real entry by `sign = +1 if x.real >= 0` with an
explicit, documented tie-break (never let a raw `np.angle` of a near-zero pivot
choose the branch), and pin the global phase by a single fixed entry (e.g.
arg(V_ud)=0). Recommend the implementation also add the rephasing-invariance
regression pin (plan §6.2) AND a determinism pin (same fit run twice / under a
random pre-rephasing must give bit-identical ε_K), which §6.2 implies but should
state explicitly as a paired-scan guard.

---

## 4. Fix coupling / order / anchored pins — CONFIRMED

- "B3's two sub-bugs must be fixed together" is correct and now has a sharp
  mechanism (§1.2 here): fixing only one leaves the chiral term either still
  cancelled (swap-only) or wrong-signed/doubled. The §0.3 direction table (B1
  DOWN, B3 UP) is consistent with my checks; order matters only for
  interpretation, agreed.
- The new ABSOLUTE pins the plan proposes (O4/O5 ratio = (R/4+1/24)/(R/12+1/8);
  singlet>crossed at R≫1; SM-box structural anchor; CGHNP Zbb value computed
  independently in-test from the c=−c, F²=2f² dictionary; rephasing-invariance)
  are genuinely literature-anchored, not snapshots. This directly remediates M7,
  which is real: I confirmed via subagent sweep that the entire `tests/` tree
  contains NO literature-anchored absolute pin for the O4/O5 ratio, an SM-box
  ε_K check, or a CGHNP Zbb value — only regression snapshots, self-referential
  re-implementations (`_manual_B` is byte-identical to the production bracket),
  order-of-magnitude inequalities, and input-constant pins. The 77/77-pass-with-
  three-blockers-live claim is therefore credible.
- Budget anchoring: slice-4 (I-1) affirms budgets are exp/SM-anchored and were
  NOT calibrated to the buggy ⟨O⟩, so the plan's "no budget re-derivation needed,
  only re-run" is correct.

---

## 5. Scope calls (M6/M1/M2/M5/M8/re-scan) — ACCEPTABLE

- **M6 defer**: acceptable. Its direction for the kaon LR piece is conservative
  (~19% would push ε_K floors UP), so deferring cannot manufacture an
  anti-conservative claim. Doing it half-right (without NLO bag scheme-matching)
  could introduce a new inconsistency — better as its own dual-gated item. The
  in-pass doc caveat is the right minimum.
- **M1 (T010 gate) as a separate sign-off**: correct. It changes a *published
  floor* and is a physics-policy choice (1σ vs T011-style uncertainty-aware
  budget). With corrected B1 the floor swings ×16 between 1σ (108 TeV artifact)
  and 2σ (6.8 TeV), so it must be settled with the user, not unilaterally. Right
  call to flag, not bundle.
- **M2 / M5 in scope**: both are cheap, anti-conservative/integrity fixes;
  including them is fine. M5's structured-tag-hint is the durable fix.
- **M8 defer**: acceptable (latent ≤2.2e-4; off the floor path).
- **Re-scan + restatement deferred to a compute-gated follow-up**: correct, given
  the M1 decision materially moves the minimal floor and the site must not
  auto-deploy. The dependency list (§8.1) is complete and includes correcting
  the three review_local notes and STATE_OF_PROJECT.md:29.

---

## 6. The review_local notes are correctly identified as the DEFENDANT — CONFIRMED

Subagent verification confirms all of the plan §0.1 claims:
- constraint_formulas.tex prints O1=2/3 (L111), O4=(R/6+1/4) (L112),
  O5=(R/2+1/12) (L113), "no extra 1/(2m_M)" (L128-130), and the Zbb bracket
  `1/(1−2c)(1/F²−1+F²/(3+2c))` (L439-443).
- epsilon_k_review.tex prints the same coefficients (L693/695/697) and the
  "do not 'fix' the 1/6 / a ×2 would BREAK literature agreement" callout
  (L804-816).
These notes document and defend the buggy code; they are downstream
documentation to be corrected, not independent authority. The plan's framing is
right. The subtlety the plan handles correctly: the 1/6 *Wilson* factor is
genuinely correct (½ propagator × ⅓ colour-Fierz, deltaf2.py:340-341); only the
*matrix element* is ÷2'd and un-swapped. The notes' "do not touch the 1/6"
warning is true but mis-applied — it does not protect the matrix-element
coefficients, which are separately wrong.

---

## Non-blocking suggestions

1. **B2 determinism recipe.** State the exact branch-free phase-pinning recipe
   (signed-zero-safe sign convention + single global-phase anchor) and add an
   explicit determinism regression pin (same point twice → bit-identical ε_K),
   not only the rephasing-invariance pin. (§3 here.)
2. **Path precision.** The K001 test pins are at
   `tests/constraints/primary/kaon/test_K001.py:197-198`
   (unrun=0.26555011996248595, run=2.844741580402456), not `test_K001.py:197-198`
   as written in §2.3/§6.4. Note §2.3 also lists run/unrun in reverse line order
   (line 197 is unrun). Fix the path/ordering in the snapshot-repin list so the
   implementer hits the right file.
3. **Site-6 function name.** §2.1 site 6 is `_meson_matrix_elements_inline`, not
   `_meson_matrix_elements` (which is site 2). Also note `_m12_np_from_bridge_
   generic` consumes the inline helper transitively (no 7th algebra site) — worth
   a one-line note so the implementer doesn't hunt for a missing site or
   accidentally introduce a divergent 7th copy.
4. **D0 ÷2 direction.** For D0 (R≈2.2) the const terms are not negligible vs R;
   the plan's "D0 slightly DOWN ×0.89-0.93" comes from the audit's post-running
   anarchic-phase evaluation. Suggest the impl re-derive the D0/B_d/B_s post-fix
   ratios at the actual benchmark rather than carrying the audit's ballpark into
   test pins (the plan already says "re-derive exact literals at impl" — keep
   that discipline for C002/B003 too).
5. **`_manual_B` replacement.** Make explicit that the replacement CGHNP/GGMS
   in-test check must be computed from the convention dictionary
   (c=−c_repo, F²=2f²), not from the new production line, or it re-pins against
   itself (same M7 failure mode that hid the bug).

---

## Bottom line

Every disputed physics claim re-derived independently lands on the plan's side,
and the two highest-stakes items (B3 direction and the R5-irrelevance
resolution) are confirmed *more strongly* than the plan argues — the swap
exactly cancels the pre-running chiral enhancement. The ×1.7227 ε_K factor is
reproduced to all printed digits through the repo's own RG pipeline. All six B3
sites, both B1 sites, and the B2 site are verified at the claimed lines. The
scope deferrals are safe (conservative or off-path directions). The only gaps
are an under-specified B2 determinism recipe and minor path/naming imprecisions
in the snapshot-repin bookkeeping — all non-blocking.

**VERDICT: APPROVE.**
