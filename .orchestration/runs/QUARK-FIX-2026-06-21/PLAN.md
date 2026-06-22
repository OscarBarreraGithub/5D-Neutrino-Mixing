# QUARK-FIX-2026-06-21 — implementation plan (PLAN v3)

Plan author: Claude/Opus (plan-author lane). **READ-ONLY on production code; this file
is the only thing written.** Dual-signoff gate applies: this plan must be approved by an
independent Codex (gpt-5.4 xhigh) AND iterated to mutual APPROVE before any code lands.

Repo HEAD at planning time: `5328a22` (working tree had only the untracked
`.orchestration/runs/QUARK-AUDIT-2026-06-10/`).

---

## CHANGELOG (v2 → v3)

PLAN v2 was independently reviewed (`review_plan_opus_v2.md`) and **APPROVED on every axis
(M1 gate decision, refinement fold-in, restatement, sequence) EXCEPT one BLOCKER**: the cheap
floor-extraction recipe in §8.2 used a single blanket `M_KK·sqrt(max ratio)` rescale and
claimed it "valid for the FULL binding set." That is wrong for the gates that set the
*inclusive* floor and biased for the strict floor (verified against production code):
- **Collider-resonance CR0xx** (`collider_resonance.py:166–168`): `ratio = m_limit / mass_tev`
  with `mass_tev ∝ M_KK` ⇒ ratio ∝ **1/M_KK (linear)** ⇒ correct rescale is `M_KK·ratio`, not
  `M_KK·sqrt(ratio)`. CR001/002/003/004/008/010/012/013 bind in the inclusive lane.
- **EW001 oblique S/T** (`oblique_stu.py:264–266`): `ratio = chi2/budget` with `chi2` a quadratic
  form in `(dS,dT)`, each `∝ 1/M_KK²` ⇒ ratio ∝ **1/M_KK⁴** ⇒ correct rescale is `M_KK·ratio^(1/4)`.
- **ΔF=2 / T010**: even the "1/M_KK²" gates carry an `α_s(M_KK)`/QCD-running log drift, an
  `ε=Λ/k` f-factor drift, a nonlinear (R_b, A_b) veto ratio, AND an **M_KK-dependent
  perturbativity skip set** (`run_full_catalog_scan.py:830–842`) — all biasing single-tile
  extraction. The cited `rs_anarchy_mkk_min_hist*.py` only rescales a `deltaf2_ratios` field that
  the full-catalog JSONL (`run_full_catalog_scan.py:891–901`, one opaque scalar `ratio`/process)
  does not even contain.

**v3 resolution (BLOCKER fixed — option (a), the cleaner one):** §8.2 is **DROPPED as a floor
source**. ALL headline/inclusive/strict floors come from the **full §8.1 paired-1M grid** — which
is cheap anyway (~30 min wall, ~115 core-hrs) and is the trustworthy native source the
explorer/comparison builders consume. §8.2 is **retained only as an OPTIONAL, strict-ΔF2-ONLY,
order-of-magnitude cross-check**, rewritten with the **correct per-gate rescale exponents** and a
**per-draw perturbativity skip mask**, and **explicitly labelled NOT the headline/inclusive
source**. The blanket "valid for the FULL binding set" claim is struck. The §8 sequence (cheap
cross-check FIRST, then full 1M) is replaced by: **§8.1 full grid is the floor source; §8.2 is an
optional sanity tripwire only**. §9 (Phase II) and the §7 stress-test #8 are updated to match.

**Two non-blocking review notes folded in v3:**
- **§5/M1 headroom number** — the −0.996σ / 0.004σ figures use T010's `sm_prediction` (the
  SM-LIMIT formula R_b), which differs slightly from the YAML SM-fit anchor; the implementer MUST
  **recompute this headroom from code at implementation**, NOT hardcode the audit's "0.004σ". The
  108-vs-6.8 TeV qualitative story does not depend on the exact decimal. Noted in §5/M1.
- **§6.4 test pins** — re-pin C002/B003 *predicted/ratio* lines ONLY, NOT the ⟨O⟩-invariant
  exp/SM-derived anchor-budget literals (B003:151/298, C002 anchor.budget/value/m12_budget_gev),
  which are not ⟨O⟩-dependent and must NOT move. Tightened in §6.4. Also: T010's `_build_budget`
  (`T010.py:325–354`) computes ONLY `combined_sigma` today — it has NO `central_residual`/
  `hard_veto_budget`; the implementer must ADD them. The plan's existing "mirror T011" /
  "factor out of T011" step already covers this; §5/M1 now states it explicitly.

---

## CHANGELOG (v1 → v2)

PLAN v1 was independently reviewed (`review_plan_opus_v1.md`) and **APPROVED with 5
non-blocking refinements**. Scope then **expanded** (user mandate 2026-06-22): the re-scan
and the corrected minimum M_KK floors are now **IN SCOPE**, not deferred.

**A. Reviewer refinements folded in (all 5):**
1. **B2 determinism recipe** — §1.2 now specifies a branch-free, signed-zero-safe phase-pin
   (explicit `sign = +1 if x.real >= 0.0 else -1` tie-break; never `np.angle` a near-zero
   pivot; note `np.angle(-0.0+0j)=π`), a single global-phase anchor (`arg(V_ud)=0`), and a
   NEW determinism regression pin (same point evaluated twice → bit-identical ε_K) added to
   §6.2 alongside the rephasing-invariance pin.
2. **K001 path/line fix** — §2.3 / §6.4 now cite `tests/constraints/primary/kaon/test_K001.py`
   with the correct line→value mapping (197 = unrun 0.26555…, 198 = run 2.84474…); the v1
   reversed order is corrected.
3. **Site-6 naming** — §2.1 now names site 6 `_meson_matrix_elements_inline` (site 2 is the
   generic `_meson_matrix_elements`); a note records that `_m12_np_from_bridge_generic`
   consumes the inline helper transitively (no 7th algebra site).
4. **Test anchoring discipline** — §6.1/§6.3 now state explicitly that the `_manual_B`
   replacement and the CGHNP/GGMS in-test checks MUST be computed from the convention
   dictionary (`c=−c_repo`, `F²=2 f²`) / paper numbers, NOT copied from the new production
   line (else they re-pin against themselves — the same M7 failure mode that hid the bug).
5. **Re-derive post-fix literals at the benchmark** — §2.3 now mandates re-deriving D0/B_d/B_s
   pinned literals (C002/B003) at the actual benchmark point, not the audit's anarchic-phase
   ballpark — same discipline already applied to ε_K.

**B. Expanded scope (re-scan + new floors now IN SCOPE):**
- **M1 (T010 gate policy) RESOLVED HERE** — no longer a separate sign-off. Decision +
  source justification in §5/M1 (rewritten). It sets the corrected minimal floor.
- **§8 rewritten** from "deferred follow-up" to an **in-scope re-scan + restatement plan**:
  exact reproduction commands + compute cost (§8.1), a cheaper floor-extraction option
  (§8.2), the full restatement target list (§8.3), and the deploy gate (§8.4).
- **§9 commit structure** updated: M1-gate folded into the code-fix sequence; explicit
  ordering (all code fixes + tests green FIRST → run scan → restate), one scan-artifacts+
  restatement commit at the end.

---

## 0. Orientation, scope, and a load-bearing correction to the framing

The audit (`QUARK-AUDIT-2026-06-10/`) found three independent BLOCKERs (B1 Zbb, B2 ε_K
SVD phase, B3 ΔF=2 matrix elements) plus Majors M1–M8. This plan fixes the three BLOCKERs,
selected Majors (incl. **M1, the T010 gate, resolved here — §5/M1**), and the test/documentation
debt, and **(v2 mandate) carries it through to the re-scan and the corrected minimum M_KK floors**:
the sequencing is code fixes + tests GREEN first, THEN run the scan, THEN restate every floor
artifact (§8–§9). The re-scan is in scope, not deferred.

### 0.1 The review_local/*.tex notes are the DEFENDANT, not the source of truth

**Important framing correction for the reviewer.** The prompt lists
`review_local/{deltaF2_framework_review, epsilon_k_review, constraint_formulas}.tex` as
"independent derivations / source translations." They are **not** corrected derivations —
they are internal notes that *document and defend the current (buggy) code*. Concretely:

- `deltaF2_framework_review.tex` eqs (O4ME)/(O5ME) and `epsilon_k_review.tex` eqs (O4ME)/(O5ME)
  and `constraint_formulas.tex` lines 111–113 ALL print the disputed code coefficients
  `⟨O4⟩=(R/6+1/4)f²mB4`, `⟨O5⟩=(R/2+1/12)f²mB5`, `⟨O1⟩=(2/3)f²mB1`, and explicitly say
  "no extra 1/(2m_M)" and "DO NOT 'fix' the 1/6; a ×2 would BREAK agreement with literature."
- `constraint_formulas.tex` lines 439–443 print the disputed Zbb `B(c,F)=1/(1-2c)(1/F²-1+F²/(3+2c))`.

So these `.tex` files are **downstream documentation that must be corrected by this fix**, in
exactly the same way as `STATE_OF_PROJECT.md:29`. We do NOT treat them as the authority. The
authority for the corrected formulas is:
- **B3**: GGMS (Gabbiani–Gabrielli–Masiero–Silvestrini, hep-ph/9604387) Eq. (8); ETM
  1310.5461 Eqs. (2)–(3) for the operator/ξ assignment; CFW 0804.1954.
- **B1**: CGHNP (Casagrande–Goertz–Haisch–Neubert–Pfoh, arXiv:0807.4937) Z→bb̄ ZMA equations,
  fetched primary-source by audit slice 3.
- **B2**: PDG standard CKM phase convention; Brod–Gorbahn–Stamou ε_K^SM=2.161e-3 is defined in
  that convention.

### 0.2 Resolution of the audit-vs-note contradiction (the "R5" defense)

The notes claim a past "R5 cross-check" compared the repo's `2|M12^NP|` to CFW 0804.1954
Eqs. 3.3–3.7 and got ratio `1.0000` to machine precision, "proving" the MEs correct. An
independent re-derivation done for this plan (GGMS Eq. (8) exact rational arithmetic, ETM ξ
assignment) **confirms the audit and overturns the R5 defense**:

- GGMS gives the **colour-singlet** Q4 the LARGE coefficient `(R/4 + 1/24)` and the
  **colour-crossed** Q5 the SMALL coefficient `(R/12 + 1/8)` (M12-ready, single power of m,
  all B=1). The code assigns the singlet O4 the small `R/6` and the crossed O5 the large
  `R/2` — exactly **backwards**, and each is **×2** too large:
  `code⟨O4⟩ = (R/6+1/4) = 2·(R/12+1/8) = 2·⟨Q5⟩_correct` and
  `code⟨O5⟩ = (R/2+1/12) = 2·(R/4+1/24) = 2·⟨Q4⟩_correct` — the audit's "×2-and-swap" to the
  rational digit.
- **Why R5 returned 1.0000:** the swap makes the *leading chiral term cancel* in the
  dominant post-running LR combination. Schematically the leading-R contribution to M12 is
  `C4·(R/6) + C5·(R/2)`; with `C4 = −gg/M²`, `C5 = +gg/(3M²)` this is
  `(−R/6 + R/6)·gg/M²·f²m = 0` at leading R *before running mixes C5 into C4*. The R5 check
  therefore validated only the VLL/colour-split (O1) sector — where the (2/3)↔(1/6) repo
  convention IS internally consistent — and never exercised the O4/O5 chiral split. The
  Blanke 0809.1073 Eq. 4.33 "~1.69" mismatch was the real warning, mis-attributed to
  "basis/running." **The R5 "machine-precision agreement" is real but irrelevant to the
  swapped sector; it is not evidence the MEs are correct.**

This is the single most important point for the reviewer to stress-test (see §7). If the
reviewer can produce a genuine literature ME (GGMS/Ciuchini/CFW Eq. 23) in which the
colour-singlet operator carries the SMALL coefficient, the B3 fix direction must be revisited.

### 0.3 Net direction of the three fixes (they move floors in OPPOSITE directions)

| Fix | Channel | Direction of floor change | Magnitude |
|---|---|---|---|
| B1 | minimal Zbb / T010 | DOWN (over-exclusion removed) | corrected fermion piece ~25× weaker; T010 collapses as the minimal driver |
| B2 | ε_K (Im M12) | either way per-point; ensemble currently ill-defined | up to ×6 in ratio at a point until convention fixed |
| B3 | ε_K (kaon LR) | UP (anti-conservative removed) | ε_K^NP ×1.72 ⇒ ε_K floor ×1.31 |
| B3 | B_d/B_s | ~no change | code/correct ≈1.02–1.08 |
| B3 | D0 | slightly DOWN | code/correct ≈0.89–0.93 |

Because B1 lowers and B3 raises, **the order matters only for the *interpretation* and the
re-scan, not for code correctness** — each fix is independent at the code level except that
**B3's two sub-bugs (swap + 1/(2m_M)) must be changed together** (fixing one alone makes the
kaon worse, per §0.2 cancellation). B2 is a prerequisite for any *defensible ε_K statement*
(without it Im M12 is convention-dependent), so B2 lands before the ε_K floors are restated.

---

## 1. FIX B2 — ε_K SVD column-phase convention (`quarkConstraints/fit.py`)

Ordered first per the audit (prerequisite for any ε_K statement). Lowest-risk, self-contained.

### 1.1 Location
- `quarkConstraints/fit.py`, `_ordered_dirac_svd` (lines 227–231) and `mass_matrix_observables`
  (lines 245–258). The raw `numpy.linalg.svd` factors `U_L_u, U_R_u, U_L_d, U_R_d` are used with
  **no phase fixing**; `ckm = U_L_u.conj().T @ U_L_d` (line 249).
- Consumers that inherit the convention dependence: `quarkConstraints/couplings.py:45–48,149–152`
  (mass-basis couplings `D = U† diag(F²) U`), `quarkConstraints/deltaf2.py` ε_K via
  `Im(M12^NP)`, and every phase-sensitive gate (K001, B002, B004, C002).

### 1.2 Wrong vs corrected expression
**Wrong (as coded):** SVD factors are defined only up to per-column phases
`U_L → U_L P`, `U_R → U_R P` (P diagonal unitary). This leaves masses, `|CKM|`, and the
Jarlskog invariant unchanged but rotates `D_ij → e^{i(φ_j−φ_i)} D_ij`, hence each ΔF=2 Wilson
by a common phase `e^{2iΔφ}`. `|M12|` is invariant; **`Im(M12)` (the ε_K input) is not**.
Numerically demonstrated by the audit: a physically equivalent rephasing moved ε_K ratio
17.34 → 105.1 (×6) with masses/|CKM|/J/B_d/B_s/D bit-identical.

**Corrected — chosen approach (A): rephase the SVD factors to the PDG standard CKM convention.**
After `_ordered_dirac_svd` returns `(U_L_u, masses_up, U_R_u)` and `(U_L_d, masses_down, U_R_d)`,
apply a deterministic rephasing so that `CKM = U_L_u† U_L_d` sits in the PDG standard
parameterization (real, positive `V_ud, V_us, V_cb, V_tb`; the single physical phase in `V_ub`).
Concretely, the standard 5-phase + 1-global rephasing freedom of a 3×3 unitary CKM:

1. Compute `V = U_L_u† U_L_d`.
2. Determine the diagonal up-phases `P_u = diag(e^{iα_i})` and down-phases `P_d = diag(e^{iβ_j})`
   (6 phases, 1 global redundant ⇒ 5 physical) that bring `V → P_u† V P_d` into PDG standard
   form (e.g. make row-1 and column-3 phases canonical: `V_ud, V_us` real≥0, `V_cb, V_tb` real≥0,
   and fix the remaining phase so `V_ud V_ub* V_cb V_cb*`-type Jarlskog sign is preserved).
3. Apply the SAME phases to the rotation matrices that build the mass-basis couplings, with the
   right-handed columns rotated by the **matched** phase so the mass matrix stays diagonal and
   positive: `U_L_u → U_L_u P_u`, `U_R_u → U_R_u P_u`, `U_L_d → U_L_d P_d`, `U_R_d → U_R_d P_d`.
   (Rephasing a Dirac field multiplies BOTH its L and R rotation columns by the same phase, so
   `U_L† M U_R = diag(s)` is preserved and masses are untouched.)
4. Recompute `ckm = U_L_u† U_L_d` (now in PDG convention) and return.

**Why approach A over the rephasing-invariant alternative (B):** alternative B (compute ε_K
from `arg[(V_ts* V_td)²]`-relative combination) is also valid and the audit lists it, but A is
preferred because (i) it fixes the convention for ALL downstream phase-sensitive gates
(B002/B004/C002) in one place rather than per-observable, (ii) it makes the MFV phase-protection
structure in `r` physically meaningful (the central 0710.1869 claim — currently untestable), and
(iii) it matches the convention in which `EPSILON_K_SM = 2.161e-3` (BGS) is defined.

**Determinism requirement (load-bearing) — EXACT branch-free recipe (refinement 1).** The
rephasing must be a deterministic function of the input matrices with NO sign/branch ambiguity,
or it reintroduces nondeterminism into a 1M-row paired scan (the minimal/custodial pairing joins
on `(r, mkk, draw_seed)` and would desync if ε_K flipped phase branch). The reviewer flagged the
real hazards: (i) singular values are strongly hierarchical so SVD columns are unique up to phase
(degenerate-subspace nondeterminism is NOT the threat); (ii) the threat is `np.angle` on a
**signed zero** — `np.angle(-0.0+0j) = π`, not 0 — i.e. a pivot real-part that can cross the ±0
boundary is a branch hazard; and (iii) the global-phase fix.

**Recipe (implement exactly):**
1. **Never feed a near-zero pivot to `np.angle`.** To pin the phase of a target entry `z`, use
   the explicit branch-free form: `phase = z / abs(z)` when `abs(z) > tol` (tol e.g. 1e-300 /
   a documented floor), and an explicit documented fallback (multiply by `+1`, i.e. leave the
   column phase as the canonical SVD output) when `abs(z) ≤ tol`. Do NOT route the pivot through
   `np.angle` at all — work with the complex unit `z/|z|` so signed-zero of the real part is
   irrelevant.
2. **Real-entry sign pin (signed-zero-safe).** Where the convention requires a real entry ≥0,
   pin its sign with `sign = +1.0 if x.real >= 0.0 else -1.0` (the `>= 0.0` makes both `+0.0`
   and `-0.0` map to `+1.0`, removing the signed-zero branch). Document this tie-break in-code.
3. **Single global-phase anchor.** After the 5 relative phases are fixed by making
   `V_ud, V_us, V_cb, V_tb` real ≥0, pin the one remaining global phase by a single fixed
   convention: **`arg(V_ud) = 0`** (V_ud real positive). One anchor, no branch.
4. **No `np.angle(0)` leak.** Guard every phase extraction by the `abs > tol` test above; the
   degenerate/zero-entry case takes the documented `+1` fallback rather than raising in the hot
   loop (raising would abort scan draws; a documented canonical fallback keeps determinism).
5. **Idempotence by construction** (already in §1.3): applying the rephasing to an
   already-canonical matrix is a no-op, so re-running is bit-stable.

This recipe is paired with TWO regression pins in §6.2: **rephasing-invariance** (random physical
column rephasing → identical ε_K) AND **determinism** (same point evaluated twice → bit-identical
ε_K), per refinement 1.

### 1.3 Numerical sanity check
- **Invariants preserved:** at the default benchmark (r=0.25, g*=3), `masses_up`, `masses_down`,
  `|CKM|` element-by-element, and Jarlskog `J` must be bit-identical (rtol≤1e-12) before/after.
- **Convention pinned:** after the fix, `V_ud, V_us, V_cb, V_tb` are real and ≥0; `arg(V_ub)` is
  the physical CKM phase (≈ −γ ≈ −65°-ish from the fit, sign per convention).
- **Idempotence:** applying the rephasing twice equals applying it once.
- **ε_K now convention-stable:** the audit's ×6 demonstration must collapse — feeding a randomly
  re-phased `(U_L_d, U_R_d)` through the (now-rephasing) pipeline gives the SAME ε_K ratio.

### 1.4 Interaction with other fixes
B2 changes `Im(M12)` for ε_K (K001) and the phase of B002/B004/C002. It is **orthogonal to B3's
magnitude fix** (B3 changes ⟨O⟩, B2 changes the phase of the couplings), but the two compound in
the K001 ε_K number, so the K001 regression re-pin (§4) must be done AFTER both B2 and B3 land.
The slice-4 note (I-3) requires B2 to be revalidated against all four phase-sensitive gates
(K001, B002, B004, C002), not just K001.

---

## 2. FIX B3 — ΔF=2 matrix elements: O4/O5 swap + missing 1/(2m_M) (fix TOGETHER)

### 2.1 Locations — ALL SIX computation sites (enumerated, verified by grep)
| # | File | Function / lines | Operators |
|---|---|---|---|
| 1 | `quarkConstraints/deltaf2.py` | `_kaon_matrix_elements` lines **720–722** | O1/O4/O5 kaon |
| 2 | `quarkConstraints/deltaf2.py` | `_meson_matrix_elements` lines **908–910** | O1/O4/O5 generic |
| 3 | `quarkConstraints/modern/phenomenology.py` | lines **452–454** (`_kaon_m12_np_from_bridge_match`) | kaon |
| 4 | `quarkConstraints/modern/phenomenology.py` | lines **482–484** (`_evaluate_epsilon_k_from_bridge`) | kaon |
| 5 | `quarkConstraints/modern/phenomenology.py` | lines **516–518** (`_evaluate_delta_mk_from_bridge`) | kaon |
| 6 | `quarkConstraints/modern/phenomenology.py` | lines **578–580** (`_meson_matrix_elements_inline`) | generic |

**Naming note (refinement 3):** site 6 is `_meson_matrix_elements_inline` (NOT `_meson_matrix_elements`
— that is site 2, the generic `deltaf2.py` helper). The modern-lane consumer
`_m12_np_from_bridge_generic` (phenomenology.py ~584) **calls** `_meson_matrix_elements_inline`
rather than re-typing the algebra, so it is fixed **transitively** — there is **no 7th algebra
site**. Do not hunt for a missing site and do not introduce a divergent 7th copy.

All six are byte-identical algebra and **must be changed in the same commit**, or the sync test
`tests/test_quark_deltaf2.py:120–155` (which pins constants, not formulas — it would not catch a
formula divergence, but the modern lane would silently differ) and the two lanes diverge.

### 2.2 Wrong vs corrected expression (repo variables)

Let `R = (m_M/(m_q1+m_q2))²` (`m_ratio_sq` / `r_chi`), `fp2_mp = f_M² m_M`.

**Wrong (all six sites):**
```
⟨O1⟩ = (2/3)        · fp2_mp · B1
⟨O4⟩ = (R/6 + 1/4)  · fp2_mp · B4     # = 2·⟨O5⟩_correct  (small coeff on singlet — backwards)
⟨O5⟩ = (R/2 + 1/12) · fp2_mp · B5     # = 2·⟨O4⟩_correct  (large coeff on crossed — backwards)
```

**Corrected (GGMS Eq. (8), M12-ready single-power-of-m normalization):**
```
⟨O1⟩ = (1/3)        · fp2_mp · B1
⟨O4⟩ = (R/4  + 1/24) · fp2_mp · B4    # colour-SINGLET gets the LARGE coefficient
⟨O5⟩ = (R/12 + 1/8 ) · fp2_mp · B5    # colour-CROSSED gets the SMALL coefficient
```

Source dictionary: GGMS hep-ph/9604387 Eq. (8) with full `(1∓γ5)` currents (single power of
`m_M`, B=1 normalization), operator labelling Q4=singlet `[P_L][P_R]`, Q5=crossed (ETM 1310.5461
Eqs. 2–3). This is simultaneously the **un-swap** (singlet↔large) and the **÷2** (1/(2m_M)
restored): note `(R/6+1/4)/2 = R/12+1/8 = ⟨O5⟩_correct` and `(R/2+1/12)/2 = R/4+1/24 = ⟨O4⟩_correct`
— so the corrected ⟨O4⟩ is the old ⟨O5⟩/2 and vice versa.

**Self-consistency of the (1/3) vs (2/3) for O1:** the audit's structural SM-box check (verified
to 1e-16) shows that with `C1_VLL = gg/(6M²)` and `⟨O1⟩=(2/3)f²mB1`, the textbook
`M12 = (G_F²/12π²)f²B̂ m M_W² λ_t² η S0` comes out **2× too large**; the M12-ready value is
`(1/3)f²mB1`. So the (1/3) is REQUIRED for the VLL sector too — the "colour-in-the-Wilson"
defense in the notes is half-right (the 1/6 Wilson factor is correct) but the matrix element
must be (1/3), not (2/3). **The fix is: keep the Wilson coefficients as-is (verified correct),
divide every ⟨O⟩ by 2, and un-swap O4/O5.** Equivalently, write the corrected coefficients
directly as above.

### 2.3 Numerical sanity check (before/after, from audit + adjudication)
- **Coefficient identities (exact, unit test these):**
  `⟨O4⟩_new / ⟨O5⟩_new = (R/4+1/24)/(R/12+1/8) ≈ 2.85` at R≈25.7.
  `⟨O1⟩_new = (1/3)f²mB1` exactly half the old.
- **ε_K benchmark:** the pinned `tests/test_quark_deltaf2.py:164` ε_K ratio
  `1.9286313761001348` → expected ≈ `3.3225` (×1.7227). (Re-derive the exact new literal during
  impl; do not hand-transcribe.)
- **K001 ratios (`tests/constraints/primary/kaon/test_K001.py`, refinement 2 — correct
  path AND line→value order):**
  - **line 197 = `unrun` = `0.26555011996248595`** → expected ≈ `0.4575` (×1.7227).
  - **line 198 = `run` = `2.844741580402456`** → expected ≈ `4.90` (×1.7227).
  (v1 reversed run/unrun across these two lines; line 197 is the **unrun** value. The full path
  is `tests/constraints/primary/kaon/test_K001.py`, not a top-level `test_K001.py`.) Re-derive
  exact literals at impl — these are the pre-fix snapshots, not anchors.
- **D0 / B_d / B_s literals — re-derive at the benchmark, do NOT carry the audit ballpark
  (refinement 5).** The per-system ratios below are the audit's **post-running anarchic-phase**
  estimates and are fine for direction/order-of-magnitude, but the C002 (D0) and B003 (B_d/B_s)
  pinned test literals MUST be recomputed at the **actual benchmark point** the test uses (same
  discipline already applied to ε_K), not hand-transcribed from this table. For D0 in particular
  (R≈2.2) the constant terms are not negligible vs the R-term, so the ÷2-and-swap shifts the
  ratio in a point-dependent way — re-derive it.
- **Per-system (audit table, DIRECTION ONLY — confirm exact numbers at re-scan/impl):**
  kaon ε_K ×1.72 (floor ×1.31); B_d/B_s ≈×1.02–1.08; D0 ≈×0.89–0.93.
- **SM-box structural anchor:** with corrected ⟨O1⟩=(1/3), the SM box reproduces the textbook
  ε_K^SM master formula to ~1e-16 (this becomes a literature-anchored test, §6).

### 2.4 Interaction
- B3 must be fixed as ONE change (swap + ÷2 together) — §0.2.
- B3 + B2 compound only in the K001/B002/B004/C002 numbers; re-pin after both land.
- **M6 (bag/Wilson scale pairing) decision — see §3.** M6 is the natural companion to B3 but is
  a separate, smaller (~10–19%) systematic; this plan DEFERS the M6 *numerics* but adds explicit
  scale documentation (see §3 for the justification and the partial action taken).

---

## 3. M6 — bag/Wilson scale pairing: DECISION = DEFER numerics, document explicitly

**What M6 is:** kaon `B4/B5` are FLAG MS̄(3 GeV) used at the 2 GeV Wilson point (~19% on O4);
B-meson bags at `m_b` vs 2 GeV Wilsons (~7%); mixed-scale quark masses inside `R_χ` (D0 ~17%).

**Decision: DEFER the RG conversion of bags to the 2 GeV scheme; do NOT bundle into B3.** Reasons:
1. It is a **separate physics question** (scheme/scale of hadronic inputs) from the B3 algebra
   error, and is already documented in-code (`deltaf2.py:643–650` TODO) — it is a known
   approximation, not a transcription bug.
2. Its direction for the kaon LR piece is **conservative** (~19% would push ε_K floors UP, the
   safe way), so deferring it does not create an anti-conservative claim.
3. RG-converting B4/B5 from 3→2 GeV correctly requires the NLO scheme-matching of the bag
   parameters (same scheme as the LO Wilson running) — doing it half-right could introduce a
   NEW inconsistency. Better to do it as its own dual-gated item with the right scheme treatment.

**Action taken in THIS pass (cheap, no numeric change):** update the in-code caveat and the
`docs/audits/bag_param_inventory.md` / B003 docstring so the B-meson `μ=m_b` pairing is flagged
(currently only "Yellow"), and add an explicit note that the quark-mass scale inside `R_χ` is
mixed-scale (slice-1 F3 / slice-4 M-1). Record M6 as an open follow-up item with its ~10–19%
size and conservative-for-kaon direction. **The reviewer should confirm this defer is acceptable
given the kaon direction is safe.**

---

## 4. FIX B1 — Zbb fermion-KK admixture (`quarkConstraints/rs_ew_couplings.py`)

Ordered after B2/B3 per the audit (corrected term is tiny; its removal collapses the spurious
minimal floor). **Highest physics-judgment content** — the corrected term is a ~15% correction to
the gauge piece, so the design decision is *correct it* vs *drop it with a documented bound*.

### 4.1 Locations
- `_casagrande_zbb_B_profile` (lines **1835–1847**): the per-flavour bracket `B(c,F)`.
- `build_rs_zbb_fermion_kk_mixing` (lines **931–971**): the flavour sums `B_d, B_Q` and the
  `1/F²(c_i)` index, plus `m_b` and prefactor.
- Tests that pin the bug: `tests/test_rs_ew_phase6a_zbb_fermion_mixing.py` (`_manual_B`
  re-implementation ~135–156; pinned `delta_g_L_b=2.3519e-5`, `delta_g_R_b=-5.4919e-4` at 173–176;
  T010/T011 pins `t010.predicted=0.21328…`, `t010.ratio=4.4704…` at 277–280).

### 4.2 Wrong vs corrected expression (CGHNP dictionary)

**Convention dictionary (proved exactly in slice 3, cross-validated on the verified gauge piece):**
`c_CGHNP = −c_repo`, `F²_CGHNP = 2 f_IR,repo²` (i.e. `F²_CGHNP(c) = (1−2c_r)/(1−ε^{1−2c_r}) =
2·(1/2−c_r)/(1−ε^{1−2c_r})`).

**CGHNP 0807.4937 (Z→bb̄ ZMA), in repo variables:**
```
g_L^b ⊃ +(m_b²/2M_KK²)[ (1/(1+2c_bR))·( 1/(2 f_bR²) − 1 + 2 f_bR²/(3−2c_bR) )
                         + Σ_{i=1,2} |Y_d,3i|²/|Y_d,33|² · 1/(1−2c_d_i) · 1/(2 f_bR²) ]
g_R^b ⊃ −(m_b²/2M_KK²)[ (1/(1+2c_bL))·( 1/(2 f_bL²) − 1 + 2 f_bL²/(3−2c_bL) )
                         + Σ_{i=1,2} |Y_d,i3|²/|Y_d,33|² · 1/(1−2c_Q_i) · 1/(2 f_bL²) ]
```
Two repo errors:
- **(1a) Flavor-sum F index.** The light-generation sum carries `1/F²(c_b3)` (the 3rd-gen
  singlet overlap), NOT `1/F²(c_d_i)`. Repo applies the full per-light-gen `B(c_d_i, F_d_i)`,
  inflating i=1,2 terms by `F²(c_d3)/F²(c_d_i) ≈ 10²–10⁴`.
- **(1b) Bracket dictionary.** Correct diagonal bracket in repo variables is
  `B_correct(c_r, f_r) = 1/(1+2c_r) · ( 1/(2 f_r²) − 1 + 2 f_r²/(3−2c_r) )`. Repo coded
  `1/(1−2c_r)·(1/f_r² − 1 + f_r²/(3+2c_r))` — c-sign wrong in BOTH denominators and the F²=2f²
  factor missing. For UV-localized b_R (`c_r>1/2`, scan-typical) the repo `1/(1−2c_r)` is
  NEGATIVE ⇒ wrong SIGN of δg_L^b.

**Corrected code (two coupled changes):**
1. `_casagrande_zbb_B_profile(c, F)` → return
   `(1.0/(1.0 + 2.0*c)) * (1.0/(2.0*F*F) - 1.0 + (2.0*F*F)/(3.0 - 2.0*c))`
   (guard `1+2c≠0`, `3−2c≠0`). This is the **diagonal** (3rd-gen) bracket only.
2. `build_rs_zbb_fermion_kk_mixing`: the diagonal term uses the b_R (resp. b_L) bracket; the
   light-generation flavour sum uses `1/(1−2 c_d_i)` (resp. `1/(1−2 c_Q_i)`) times the COMMON
   `1/(2 f_b3²)` factor — NOT a per-light-gen `B`. Restructure `B_d/B_Q` accordingly:
   `B_d = B_correct(c_d3, f_d3) + (1/(2 f_d3²)) · Σ_{i=0,1} row_ratio_i / (1 − 2 c_d_i)` and
   symmetric for `B_Q` with `c_Q`, `f_q3`, `column_ratio`.
   Keep `m_b² /(2 Λ_IR²)` prefactor (verified correct CGHNP geometric-Λ_IR convention).
   Optionally also fix the m_b scale (slice-3 Finding 7: use m_b(M_KK)≈2.4 not 2.73; ~+29%) —
   **defer m_b-scale to M6/follow-up; not required for sign/magnitude correctness.**

### 4.3 DESIGN DECISION: correct the term, do NOT drop it
The audit suggested "consider dropping with a documented bound." **Decision: implement the
corrected CGHNP term, do not drop.** Reasons: (i) it is a genuine ~15% correction to the
verified gauge piece and the machinery already exists; (ii) dropping would leave a documented
*bound* that still needs the corrected magnitude to state, so we must compute it anyway;
(iii) keeping the term, correctly signed, makes the minimal-vs-custodial contrast quantitatively
honest (custodial zeroes it; minimal keeps a small correct +15%). The reviewer should weigh
correct-vs-drop; if review prefers drop, the fallback is to set the fermion piece to zero with a
metadata bound `|δg_L^b,fermion| ≲ 0.15·|δg_L^b,gauge|` and a NEEDS-HUMAN tag — but the
plan's recommendation is to correct it.

### 4.4 Numerical sanity check (audit table, representative point r=0.25, Λ_IR=8.168 TeV)
| quantity | repo (wrong) | CGHNP-correct | check |
|---|---|---|---|
| B_d | −4.50e4 | +2.36e2 | sign flip, ratio −191 |
| B_Q | −5.38e3 | +2.63 | sign flip, ratio −2042 |
| δg_L^b (fermion) | −2.52e-3 | +1.32e-5 | **sign flip, ~190×** |
| δg_R^b (fermion) | +3.01e-4 | −1.47e-7 | sign flip, ~2000× |
| δg_L^b (gauge, unchanged) | +8.74e-5 | +8.74e-5 | corrected fermion is ~15% of gauge |
- Post-fix, corrected fermion piece must have the SM-Zbb sign (δg_L^b>0 for IR-localized b_L,
  reducing R_b) and be ~15% of the gauge piece, not ~30× it.
- The two existing pinned test literals (`2.3519e-5`, `-5.4919e-4`) and the T010 pins
  (`0.21328`, ratio `4.4704`) are pins-of-the-bug and MUST be re-derived (§6). The `_manual_B`
  re-implementation in the test must be replaced by a **literature-anchored** check, not a
  byte-copy of the new code (else it re-pins-against-itself; §6/M7).

### 4.5 Interaction with M1/M2 (see §5)
B1 alone does not set the corrected minimal floor — the T010 *gate policy* (M1) does. With
corrected couplings the T010 1σ "floor" is ~108 TeV (gate artifact) vs ~6.8 TeV (2σ). M1 is now
**resolved in-scope** (§5/M1: adopt the T011-style loose-edge budget on R_b), so B1's commit lands
the corrected couplings and the M1 commit lands the gate that turns them into a sensible floor; the
two are reasoned about together and both precede the re-scan. The CODE fix for B1 is independent of
the M1 gate change (different files), so they remain separate commits.

---

## 5. Majors M1, M2, M5 — scope decisions

### M1 (T010 gate policy) — DECISION RESOLVED HERE: adopt the T011-style loose-edge NP-shift budget on R_b. IN SCOPE, NOT a separate sign-off.

**Scope change (v2):** the user mandate makes the new minimal floor in-scope, and the floor is
*set by* this gate, so M1 can no longer be a downstream proposal — it is resolved here as the
plan author's design call (§B1 of the mandate) and the reviewer will check it.

**The problem (slice-3 Finding 2, source-confirmed).** T010 currently vetoes at
`passes = max|pull| ≤ 1.0` with `pull = (predicted − experimental)/σ_combined`
(`T010.py:653–655`). The SM-limit R_b prediction already sits at pull **−0.996σ** against the
LEP/SLC SM-fit anchor (R_b^exp=0.21629±0.00066, SM-limit=0.21562, σ_comb=6.73e-4) — i.e. the SM
passes by **0.004σ**.

> **Implementer note (review non-blocking #1 — RECOMPUTE, do not hardcode).** The −0.996σ /
> 0.004σ figures here use T010's `sm_pull = (sm_prediction − experimental)/σ_combined`
> (`T010.py:401`), where `sm_prediction` is the **SM-LIMIT formula** R_b, NOT the YAML SM-fit
> anchor — the two differ slightly. The implementer MUST **recompute this headroom number from
> code at implementation time** (evaluate the SM-limit pull directly) rather than transcribing the
> audit's "0.004σ" literal into any test/pin. The 108-vs-6.8 TeV qualitative story and the gate
> decision do NOT depend on the exact decimal, so this is a precision/pinning matter only — but
> the literal must come from code, not from this paragraph.

This is a **two-sided max-pull-vs-experiment** gate with essentially zero
headroom on the low-R_b side: *any* RS shift in the standard direction (δg_L^b ≳ +4e-6) fails it
at arbitrarily high M_KK until the 1/M_KK² shift decays below 0.004σ. With corrected B1 couplings
that yields a **108 TeV** "floor" at the representative point — a pure gate artifact, not physics
(it swings to 6.8 TeV under a 2σ gate: a ×16 sensitivity to the gate convention proves the number
measures anchor-snapshot tension, not RS).

**The decision: change T010's veto scalar for the R_b/A_b/A_FB pseudo-observables to the
T011-style loose-edge budget on the absolute NP shift:**
```
ratio   = |predicted − SM_limit_prediction| / hard_veto_budget
budget  = |exp − SM_fit| + sqrt(σ_exp² + σ_SM²)        # central residual + combined σ
passes  = ratio ≤ 1.0
```
i.e. veto on how far the **RS shift** moves the observable from the **SM-limit** prediction,
normalized by a budget that already absorbs the standing SM↔data tension — exactly the recipe
T011 already implements (`T011.py:481–513` `_build_budget`, :568–569 `np_shift`/`ratio`; its
docstring lines 34–38 state the rationale verbatim: "the veto scalar is the absolute NP shift
relative to the SM-limit pseudo-observable divided by the … loose-edge budget … Thus the SM point
is not vetoed merely because of the legacy anomaly, while a large RS shift is excluded").

**Source justification (the reviewer will check this against CGHNP 0807.4937 §6.4 + slice-3):**
1. **It is the correct statistical object.** CGHNP §6.4 constrains RS via the *shift* the KK tower
   induces in the Z→bb̄ pseudo-observables relative to the SM, read against the combined R_b/A_b/
   A_FB^b electroweak fit — not by re-testing whether the SM itself agrees with data. A
   `max|pred − exp|/σ` gate conflates the standing SM tension (which RS does not author) with the
   NP exclusion; the loose-edge budget cleanly separates them. This is the same logic CGHNP use
   when they quote an M_KK reach from the EW fit's tolerance to the RS contribution.
2. **It is internally consistent with the validated sibling gate.** T011 (the A_b/A_FB^b half of
   the *same* Z→bb̄ measurement) already uses this budget and is dual-approved and rigorous-tagged.
   Having T010 (R_b) use a *different, stricter* statistic on the *same* underlying LEP/SLC data is
   the actual bug: the two halves of one measurement must use one gate convention. Adopting T011's
   makes the Zbb sector self-consistent.
3. **It is conservative where it must be and not artifactual.** The budget `|exp−SM|+σ_comb` is
   *wider* than 1σ_exp (it adds the central residual and σ_SM in quadrature), so it does not
   manufacture spurious exclusion from the SM's own −0.996σ offset; but it still excludes a genuine
   large RS shift. CGHNP's own R_b reach (their published O(few–10) TeV figures) is recovered by a
   tolerance of this order, NOT by a 0.004σ-headroom cut.

**Resulting floor behavior (to be confirmed numerically at re-scan, §8):** with corrected B1
(gauge-dominated δg_L^b ≈ +1.0e-4·(20 TeV/M_KK)², ~25× weaker than the buggy term) AND this gate,
the spurious 108 TeV collapses; the T010 R_b contribution to the minimal floor drops to the
**O(few TeV)** range (the representative point gave ~6.8 TeV under the comparable 2σ comparison,
and the loose-edge budget is of similar width). T010 therefore **ceases to be the minimal-floor
driver**; the corrected minimal strict floor is then set by the corrected ε_K (B3, ×1.31 up) and
EW001 (proxy, ~15–17 TeV if M2 applied) — to be read off the re-scan, not asserted here.

**Fallback / reviewer off-ramp:** if the reviewer judges the loose-edge budget too permissive for
R_b specifically, the documented alternative is a **2σ two-sided cut** (floor ≈6.8 TeV at the
representative point) — still source-defensible and still collapsing the artifact. The plan's
primary recommendation is the T011-style budget for sibling consistency; the 2σ cut is the
fallback. Keeping the bare 1σ gate is NOT acceptable post-B1 (it yields a meaningless 108 TeV).

**Implementation note:** T010 and T011 share the Z-pole machinery; the change is to T010's
`evaluate` veto scalar (`T010.py:648–655`) to reuse the loose-edge budget construction (factor it
out of T011 or mirror `_build_budget`). **Verified gap (review non-blocking #3):** T010's
`_build_budget` (`T010.py:325–354`) today computes ONLY `combined_sigma = sqrt(σ_exp²+σ_sm²)` — it
does **NOT** carry a `central_residual = |exp − sm_fit|` or a `hard_veto_budget`. The implementer
must **ADD** the `central = |exp − sm_fit|` term and the `hard_veto_budget = central + combined`
field to T010's budget builder (mirroring `T011._build_budget`, :481–515). The plan's "factor out
of / mirror T011" instruction already covers this; it is called out here so the implementer does
not assume the field exists. Add a literature-anchored pin: at the SM-limit prediction
the NP shift is 0 ⇒ ratio 0 ⇒ SM passes by construction (no 0.004σ knife-edge); and a known RS
shift of a fixed magnitude maps to a known ratio. Re-pin the T010 snapshot literals
(`t010.predicted`, `t010.ratio`) — they change with B1 anyway (§6.4).

### M2 (EW001 ΔT Λ_IR-vs-M_KK convention) — DECISION: IN SCOPE (small, clean fix).
`oblique_stu.py:158–169,233–237` + `EW001.py:338–346`: the CGHNP geometric-Λ_IR ΔT coefficient
`πL/(2c_W²)` is evaluated with the PHYSICAL `M_KK = x₁·Λ_IR`, so ΔT is `x₁²≈6.0` too small
(anti-conservative), while the ΔS coefficient `c_S=30` is calibrated in the physical convention —
EW001 mixes conventions inside one formula. Fix: evaluate ΔT in a single, documented convention
consistent with the ΔS calibration (use `(v/Λ_IR)²` for the ΔT term, i.e. divide the scale by
`x₁²`, OR equivalently rescale the coefficient). EW001 is proxy-tagged and excluded from the
strict floor, so this does NOT move the strict headline, but it is a genuine convention bug with
an anti-conservative direction and is cheap to fix. **Include it; it makes the inclusive ~7 TeV
custodial floor and the ~6–7→15–17 TeV minimal EW001 statement honest.** (Cosmetic M2 sub-items —
`DEFAULT_RS_VOLUME_LOG=35.0` vs point L=34.54, 1.3% — defer.)

### M5 (tag_result substring bug "proxies" vs "proxy") — DECISION: IN SCOPE (cheap integrity fix).
`scripts/run_full_catalog_scan.py:984–987` (`"proxy" in needs_text` misses the plural "proxies")
mis-tags T001/T002 (evaluated+passing) as `partial`→`hard_not_evaluated` on 100% of rows ⇒
`coverage_complete=False` everywhere, and a FAILING partial-HARD would silently never veto. Fix:
replace the brittle substring match with a structured boolean tag hint (slice-5 F8 recommends a
`tag_hint`/diagnostic flag rather than prose matching), OR at minimum make the proxy/rigorous
detection robust to plurals and to "proxy pending rigorous" phrasings. **Recommendation: add a
small structured diagnostic on the constraint result (e.g. `tag_class ∈ {rigorous,proxy,partial,
stub}`) and have `tag_result` read it, with the prose match as a deprecated fallback.** This is
the more durable fix and prevents the latent silent-non-veto. It does not change any current
veto count (T001/T002 pass by 7–9 orders of magnitude) but fixes the coverage metric and the
brittleness. Include the matching explorer/comparison consistency fix (slice-5 F7: explorer
counts any active+evaluated HARD failure regardless of tag, unlike the harness — align them).

### M3, M4 — documentation/claim items: handled in the re-scan/restatement follow-up (§8), NOT
code. M4 (grid-honest brackets "(20,30] TeV", missing r=0.05 floor) and M3 (custodial 2–3 TeV
top-partner-loop caveat) are prose corrections to STATE/ledger/methodology note done when floors
are restated.

### M7 — test strategy: IN SCOPE, mandatory (§6).
### M8 — qcd/ package (d3 3-loop constants ≤2.2e-4; pdg_quark_masses_at_scale crash <4.18 GeV):
**DEFER** — latent (≤2.2e-4 effect) and the crash is outside the ΔF=2/Zbb path exercised by the
floors; record as a separate low-priority item. The reviewer may pull M8's crash-guard in if it
is trivially cheap, but it is not required for floor defensibility.

---

## 6. Test strategy — literature-anchored ABSOLUTE pins (M7)

The audit's M7 is the meta-finding: the suite pins **code-against-itself** (regression snapshots
and `_manual_B`/sync re-implementations), so it passed 77/77 with all three blockers live. Every
fix MUST add at least one **literature-anchored absolute pin** (fixed known input → a number
from the paper/derivation, NOT a snapshot of current output).

**Anchoring discipline (refinement 4 — applies to ALL pins in §6.1 and §6.3):** every "absolute"
pin below MUST compute its expected value from the **convention dictionary / paper numbers**
(GGMS Eq.(8) closed-form rationals; the `c=−c_repo`, `F²=2 f²` CGHNP map; the textbook ε_K
master formula), and must NOT read the expected value off the new production line. A pin that
copies the post-fix production output re-pins the code against itself — the exact M7 failure mode
that let 77/77 pass with three blockers live. Where a pin needs a literature constant, hard-code
the closed-form rational or the cited paper number in the test, with a source comment.

### 6.1 B3 absolute pins (NEW tests, `tests/test_epsilon_k_physics.py` / `test_quark_deltaf2.py`)
- **Coefficient ratio pin:** `⟨O4⟩/⟨O5⟩ = (R/4+1/24)/(R/12+1/8)` evaluated at the repo's kaon R
  must equal the **closed-form rational computed in-test** (≈2.847 at R=25.7) — anchors the
  un-swap. Source: GGMS Eq.(8). Do NOT read the ratio from the production ⟨O⟩ functions; build
  it from the rational coefficients written literally in the test.
- **Singlet-carries-large pin:** assert the colour-singlet O4 coefficient `> ` the crossed O5
  coefficient at R≫1 (guards against re-swap).
- **SM-box structural anchor:** with the corrected `⟨O1⟩=(1/3)f²mB1` and the repo's VLL Wilson
  normalization, an independent textbook SM ε_K box (`G_F²/12π² · f²B̂ m M_W² λ_t² η S0`) must
  reproduce to ~1e-12. Source: Buras ε_K master formula. (This is the pin that would have caught
  the missing 1/(2m_M).)
- **Absolute ME pin:** `⟨O1_K⟩(2 GeV) ≈ (1/3)·f_K²·m_K·B1_K` to a fixed literature number with
  FLAG-tolerance (slice-4 M-2 suggests ≈2.21e-3 GeV³ M12-ready — recompute exact at impl).

### 6.2 B2 absolute pins (NEW, `tests/test_quark_fit.py` / `test_epsilon_k_physics.py`)
- **Rephasing invariance pin:** generate a point, apply a RANDOM physical column rephasing to the
  inputs, run the full ε_K pipeline; assert the ε_K ratio is INVARIANT (rtol≤1e-10). This is the
  exact ×6 demonstration from the audit, inverted into a regression guard.
- **Determinism pin (NEW, refinement 1 — paired-scan guard):** evaluate the SAME benchmark point
  TWICE through the full fit→ε_K pipeline and assert the two ε_K values are **bit-identical**
  (`==`, not just `rtol`). Separately, pre-rephase the inputs by a random physical phase and
  assert bit-identical ε_K to the un-rephased run. This catches any residual `np.angle`/signed-zero
  branch nondeterminism that rephasing-invariance-with-rtol could mask, and is the property the
  1M-row minimal/custodial `(r,mkk,draw_seed)` join depends on.
- **PDG-convention pin:** after the fit, assert `V_ud,V_us,V_cb,V_tb` real and ≥0 and that the
  CKM phase reproduces a known input (the fit's target γ within tolerance).
- **Invariants pin:** masses/|CKM|/Jarlskog unchanged by the rephasing (rtol≤1e-12).

### 6.3 B1 absolute pins (REPLACE `_manual_B`, `tests/test_rs_ew_phase6a_zbb_fermion_mixing.py`)
- **CGHNP value pin (refinement 4, load-bearing):** at a fixed `(c_bR, f_bR)` corresponding to a
  documented CGHNP/scan point, assert `B_correct(c,F)` equals the number from the **CGHNP
  equation**, computed independently IN THE TEST from the `c_CGHNP=−c_repo`, `F²=2 f²` dictionary
  — i.e. the test writes out `1/(1+2c)·(1/(2F²) − 1 + 2F²/(3−2c))` from the convention map and
  compares. It MUST NOT call the production `_casagrande_zbb_B_profile` or byte-copy the new
  bracket line. The old `_manual_B` was byte-identical to the (buggy) production bracket — that is
  precisely why the suite passed with the sign-flipped term live. The replacement check is the
  literature dictionary, not the code.
- **Sign pin:** for a UV-localized b_R (`c_r>1/2`), assert `δg_L^b,fermion > 0` (SM-Zbb sign).
- **Magnitude pin:** corrected fermion piece is within a factor ~few of the gauge piece (assert
  ratio `|fermion|/|gauge| < 1`), guarding against the ~30× inflation.
- **Limit pins (keep):** `m_b=0 ⇒ δg=0`; `δg ∝ 1/Λ_IR²` (ratio 0.25 for 2× Λ_IR) — these are
  physics-correct and stay.

### 6.4 Re-pin (NOT re-baseline-blindly) all snapshot literals
Every snapshot literal enumerated below must be RECOMPUTED from the corrected code and replaced,
in the SAME commit as its fix, with a comment pointing at the literature-anchored test that
justifies the new value (so a future reader sees the snapshot is downstream of an absolute pin,
not the source of truth):
- B3: `test_quark_deltaf2.py:164` (1.9286…→≈3.3225), `test_quark_plot_data.py:33` (same),
  `test_K001.py:127,197,198` (budget unchanged; ratios ×1.7227), the chiral-enhancement ballpark
  bounds `test_epsilon_k_physics.py:182–198` (still pass; tighten if desired),
  `test_modern_phenomenology.py`.
  - **B003 / C002 — re-pin the PREDICTED/RATIO lines ONLY (review non-blocking #2).** For B003,
    re-pin `result.predicted` (≈4.76e-14) and `result.ratio` (≈0.18065, ≈×1.05); do **NOT** touch
    the budget literals at `test_B003.py:151` / `:298` — those are exp/SM-derived **anchor
    budgets** (`budget == 2.635167648629676e-13` etc.), which per slice-4 I-1 are **NOT
    ⟨O⟩-dependent and must NOT move**. For C002, re-pin only the *predicted*/*ratio* lines (D0
    ≈×0.92, the relevant lines among `test_C002.py:226,228,230`); do **NOT** re-pin
    `anchor.budget` / `anchor.value` / `m12_budget_gev`, which are ⟨O⟩-invariant anchors. Re-pin
    every moved literal at the **actual benchmark point**, re-derived (not hand-transcribed). If a
    "re-derived" literal lands on an anchor-budget line, that is a signal the wrong line is being
    re-pinned — anchor budgets are fixed.
- B2: `test_K001.py:197–198` (if B2 changes Im split — recompute after B2+B3), `test_quark_fit.py`
  quotient-invariance tests (should still hold, now with the rephasing applied).
- B1: `test_rs_ew_phase6a…:173–176,277–280` (all change: sign + magnitude + T010 ratio).

---

## 7. Key points the reviewer must stress-test (open questions)

1. **B3 direction (the big one).** Is GGMS Eq. (8) / ETM ξ-assignment correct that the
   colour-singlet O4 carries the LARGE `(R/4+1/24)` and the crossed O5 the small `(R/12+1/8)` in
   the M12-ready normalization? If the reviewer can cite a standard ref where the singlet carries
   the SMALL coefficient, the fix inverts. **The whole 31% anti-conservative correction rides on
   this.** Also confirm the §0.2 cancellation explanation for why R5 returned 1.0000 (i.e. that
   the R5 check never exercised the O4/O5 chiral split).
2. **B3 vs the three review_local notes.** The notes explicitly forbid the ÷2 ("DO NOT fix the
   1/6; a ×2 would BREAK literature agreement"). Confirm the resolution: the 1/6 Wilson is NOT
   touched; only the matrix element is ÷2'd and un-swapped, and this is consistent with — not a
   violation of — the CFW comparison the notes cite (because CFW agreement was on the VLL sector).
3. **B2 approach A vs B and determinism.** Is the explicit PDG-rephasing (A) the right call vs the
   rephasing-invariant ε_K (B)? Is the proposed phase-pinning genuinely branch-free for a 1M-row
   paired scan (no `np.angle(0)` / sign-flip nondeterminism that would break pairing)?
4. **B1 correct-vs-drop.** Should the corrected CGHNP fermion term be implemented (~15% of gauge),
   or dropped with a documented bound? Plan recommends implement.
5. **M1 gate policy (now RESOLVED in-scope, §5/M1).** Is adopting T011's loose-edge NP-shift
   budget for T010's R_b the right call, justified from CGHNP §6.4 + slice-3 (the R_b half of a
   measurement whose A_b/A_FB half already uses this budget must use the same convention; a
   `max|pred−exp|/σ` gate conflates the SM's −0.996σ standing tension with the RS exclusion)? Is
   the 2σ fallback acceptable if you judge the loose-edge budget too permissive? **This is now my
   design call, not a deferred proposal — stress-test the source justification.**
6. **M6 defer.** Is deferring the bag-scale RG conversion acceptable given the kaon direction is
   conservative? Or must it ride with B3?
7. **Scope of M5/M2.** Are the structured-tag-hint (M5) and EW001 ΔT convention (M2) fixes in
   scope for this pass, or should they be separate items?
8. **Re-scan + floor-extraction method (now IN SCOPE, §8 — BLOCKER RESOLVED in v3).** The v2
   single-rescale cheap option was numerically wrong for the inclusive-floor gates (CR0xx ∝ 1/M_KK,
   EW001 ∝ 1/M_KK⁴) and biased for ΔF=2/T010 (running/ε drift + M_KK-dependent skip set); the cited
   tool cannot even read those ratios. **v3 takes ALL headline/inclusive/strict floors from the full
   §8.1 paired-1M grid** (cheap: ~30 min wall, ~115 core-hrs) and demotes §8.2 to an OPTIONAL,
   strict-ΔF2-ONLY, per-gate-exponent order-of-magnitude tripwire that sets NO headline number.
   Residual reviewer check: confirm option (a) (drop §8.2 as a floor source; full grid is the
   trustworthy headline) is the right call versus fully correcting §8.2, and that the corrected
   §8.2 recipe (ΔF2-only, `M_KK·ratio^(1/2)`, per-draw skip mask) is sound as a mere tripwire.

---

## 8. Re-scan + downstream restatement — IN SCOPE (v2 mandate)

**Scope change (v2):** the user mandate "actually run the scan" makes the re-scan and the corrected
minimum M_KK floors part of THIS program (no longer a deferred follow-up). The discipline is strict
sequencing: **all code fixes + M1 gate + full test suite GREEN first, THEN run the scan, THEN
restate** (§9). The budgets are exp/SM-anchored and were NOT calibrated to the buggy ⟨O⟩
(slice-4 I-1), so **no budget re-derivation is needed — only re-running the scans** with fixed code.

### 8.1 EXACT reproduction of the paired 1M runs (the SOLE floor source — v3)

Verified scan plumbing (paths absolute under repo root):
- Driver: `scripts/run_full_catalog_scan.py` (`--quark-only`, `--n-draws` per M_KK tile,
  `--m-kk-tev` CSV, `--ew-model {minimal_rs|custodial_rs_plr}`, `--quark-fit-r`, `--base-seed`,
  `--tile-seed-stride`, `--xi-kk` default 2.4487). Fully deterministic: `tile.seed = base_seed +
  tile_seed_stride·tile_id`, `draw_seed = tile.seed + draw_idx`, fresh `np.random.default_rng(seed)`
  per draw. Output: JSONL (`tile-NNNNN.jsonl` + `.summary.json`), atomic, auto-resume.
- Grid plan: `scripts/wq_quarkonly_1m_plan.py` — `R_GRID=(0.05,0.10,0.25,0.50,1.00)`,
  `M_KK_TEV=(1,2,3,5,7,10,15,20,30,50)`, `DRAWS_PER_MKK_PER_R=20_000`, `DRAW_SHARDS_PER_R=10`
  ⇒ 2000 draws/tile, **50 array tasks (5 r × 10 shards), 1,000,000 draws total** per run.
- SLURM wrappers (these ARE the prior-run launchers):
  - Minimal: `scripts/run_wq_quarkonly_1m_array.sbatch`
    (`--partition=serial_requeue --account=randall_lab --cpus-per-task=48 --time=04:00:00
    --mem=64G --array=0-49%50`); output root `scan_outputs/wq_quarkonly_1M_${SLURM_ARRAY_JOB_ID}`.
  - Custodial: `scripts/run_wq_quarkonly_1m_custodial_array.sbatch` (identical resources; the ONLY
    physics delta is `--ew-model custodial_rs_plr`); output root
    `scan_outputs/wq_quarkonly_1M_custodial_${SLURM_ARRAY_JOB_ID}`. Task 0 writes
    `custodial_scan_manifest.json` recording `baseline_minimal_root`.

**Exact commands (full reproduction):**
```bash
# 0. (pre-edit) update the custodial sbatch baseline_minimal_root literal to the NEW minimal job id
#    after step 1 returns it; or leave and patch the manifest post-hoc.
sbatch scripts/run_wq_quarkonly_1m_array.sbatch            # -> wq_quarkonly_1M_<NEWMIN>   (minimal)
sbatch scripts/run_wq_quarkonly_1m_custodial_array.sbatch  # -> wq_quarkonly_1M_custodial_<NEWCUST>
# analysis + comparison + explorer JSON:
python scripts/analyze_wq_quarkonly.py --input scan_outputs/wq_quarkonly_1M_<NEWMIN>  --output-dir scan_outputs/wq_quarkonly_1M_<NEWMIN>/analysis
python scripts/analyze_wq_quarkonly.py --input scan_outputs/wq_quarkonly_1M_custodial_<NEWCUST> --output-dir scan_outputs/wq_quarkonly_1M_custodial_<NEWCUST>/analysis
python scripts/build_wq_quarkonly_comparison.py --minimal scan_outputs/wq_quarkonly_1M_<NEWMIN> --custodial scan_outputs/wq_quarkonly_1M_custodial_<NEWCUST> --output scan_outputs/wq_quarkonly_1M_custodial_<NEWCUST>/comparison
# explorer builder hard-codes the OLD run ids at module scope (build_scan_explorer.py:41-42 MINIMAL_ROOT/CUSTODIAL_ROOT) -> must edit to the new ids before running:
python flavor_catalog/website/scripts/build_scan_explorer.py   # -> flavor_catalog/website/src/content/scan_explorer.json
```

**Compute cost (measured from the prior runs' 500 tile `.summary.json` timing blocks):**
- **~55 CPU-core-hours** (minimal) / **~60** (custodial) per 1M-draw run — **~115 core-hours for
  the pair**. ~0.20–0.22 s per evaluated point; ~878,700 of 1,000,000 draws evaluate (the rest
  skip on quark-fit-fail / non-perturbative, deterministically — the minimal/custodial skip sets
  match one-to-one by seed). A fixed ~40–46 s per-tile RS-EW cache build (×500 tiles ≈ 6 core-hrs).
- **Wall-clock:** with `--array=0-49%50` all 50 tasks run concurrently; each task does 20k draws in
  ~7–8 min, so **each 1M run finishes in well under its 4 h limit** (≈10–15 min wall + scheduling);
  the pair is **~30 min wall** if both arrays schedule promptly on serial_requeue.

### 8.2 Optional strict-ΔF2-ONLY order-of-magnitude cross-check — NOT a floor source

> **BLOCKER resolution (v3, option (a)).** PLAN v2 presented this section as the cheap *floor
> extractor* and claimed a single `M_KK·sqrt(max ratio)` rescale "valid for the FULL binding set."
> That claim is **STRUCK**: it is numerically WRONG for the gates that set the *inclusive* floor and
> biased for the strict floor (verified against production code — see the v2→v3 changelog). The
> headline/inclusive/strict floors are therefore taken **entirely from the full §8.1 paired-1M
> grid**, which is the trustworthy native source the explorer/comparison builders consume and is
> cheap anyway (~30 min wall, ~115 core-hrs). **§8.2 is dropped as a floor source.** What remains
> below is an OPTIONAL, strict-ΔF2-ONLY, order-of-magnitude sanity tripwire — it does NOT set ANY
> headline number, does NOT set the inclusive floor (CR0xx/EW001), and is NOT consumed by any
> artifact. Use it (if at all) only to catch a gross regression in the ΔF=2 sector before paying for
> the full run.

**Why the single-rescale shortcut fails (per-gate, verified in code).** The binding set is NOT just
ΔF=2; `constraint_veto_by_r_mkk.csv` shows it includes T010/Zbb, EW001/oblique, and 8
collider-resonance CR0xx bounds. The per-draw rescale `M_KK·sqrt(max_X ratio_X)` is only correct for
gates whose veto ratio scales as `1/M_KK²`. That is FALSE for the inclusive-floor drivers:

| Gate | `ratio` scaling in M_KK | Correct per-draw rescale | Source |
|---|---|---|---|
| ΔF=2 / ε_K (K001, B/D systems) | `1/M_KK²` (leading order only) | `M_KK · ratio^(1/2)` | `deltaf2.py` |
| Collider-resonance CR0xx (mass threshold) | **`1/M_KK` (linear)** | `M_KK · ratio` (≡ `m_limit`, a constant) | `collider_resonance.py:166–168`; `mass_tev ∝ M_KK` (:109–112,126) |
| EW001 oblique S/T (χ² quadratic form) | **`1/M_KK⁴`** | `M_KK · ratio^(1/4)` | `oblique_stu.py:252–267`; `chi2 = ds·Σ⁻¹·ds`, ds,dt ∝ 1/M_KK² |
| T010 / Zbb (R_b, A_b) | `~1/M_KK²` but NONLINEAR + ε-drift | (do not rescale — re-evaluate) | `T010.py`; δg carries ε=Λ/k f-factor drift, R_b is nonlinear in δg |

For CR0xx, `M_KK·sqrt(m_lim/M_KK) = sqrt(M_KK·m_lim)` is simply the wrong functional form — the
correct hard-mass-threshold floor is the constant `M_KK_min = m_limit`. For EW001,
`M_KK·sqrt(1/M_KK⁴) ∝ 1/M_KK` — also wrong. A blanket `sqrt` rescale therefore mis-locates the very
gates that set the inclusive floor.

**Additional biases even for the ΔF=2 gates** (why §8.2 is order-of-magnitude only, not trustworthy):
- The "1/M_KK²" of ΔF=2 is leading-order only: Wilsons are matched at `μ=M_KK` and QCD-run M_KK→2
  GeV (`deltaf2.py`/`qcd_running.py`), and `g_s²=4π·α_s(M_KK)` (`couplings.py`), so there is an
  `α_s(M_KK)`/running log drift; the overlap factors carry an `ε=Λ/k` drift. Log-level, but real.
- The quark-fit **perturbativity SKIP** (`run_full_catalog_scan.py:830–842`) depends on Yukawas
  fitted at the tile Λ_IR, so the **skip SET is M_KK-dependent**. A single-tile rescale evaluates
  the skip mask at ONE M_KK, biasing the survival-fraction crossing.
- The cited `scripts/rs_anarchy_mkk_min_hist*.py` rescales ONLY a `deltaf2_ratios` field that lives
  in the *separate* `run_rs_anarchy.py` pipeline. The full-catalog JSONL
  (`run_full_catalog_scan.py:891–901`) stores ONE opaque scalar `ratio` per process_id and NO
  per-system breakdown — so that tool **cannot read** the T010/EW001/CR0xx (or even per-system
  ΔF=2) ratios in any case.

**If you run the optional cross-check anyway, do it correctly:**
1. **Restrict to the strict ΔF=2 sector ONLY** (ε_K/K001 + B/D systems). Do NOT include CR0xx,
   EW001, or T010 — they set the inclusive floor and are not 1/M_KK² (per the table). The cross-check
   speaks ONLY to whether the corrected ΔF=2 sector lands in the expected order of magnitude.
2. **Use the per-gate rescale exponent** (`M_KK·ratio^(1/2)` for the leading-order ΔF=2 ratios). Do
   NOT apply a blanket exponent.
3. **Apply the per-draw perturbativity skip mask** evaluated at each draw's *rescaled* M_KK_min, not
   at the single tile Λ_IR (else the survival fraction is biased by the M_KK-dependent skip set).
4. Run a single tile (`--m-kk-tev 3`, `--n-draws ~2000`, 5 r values) per model — minutes of wall,
   a few core-hours — and read the resulting ΔF2-only `M_KK_min` histogram as an
   **order-of-magnitude** check on the corrected ε_K/B/D floor, NOT as a number to publish.

**Explicit non-status:** §8.2 output is NEVER the headline, NEVER the inclusive floor, and is NEVER
fed to the explorer/comparison/restatement artifacts. Those come from §8.1. §8.2 is an optional
tripwire that may be skipped entirely without affecting any deliverable.

**Recommendation:** go **straight to the full §8.1 paired-1M reproduction** for ALL floors — it is
~30 min wall / ~115 core-hrs, cheap enough that the §8.2 shortcut buys little, and it produces every
artifact in its native format (the explorer/comparison builders consume the full JSONL grid).
Optionally run the §8.2 strict-ΔF2-only cross-check (per the corrected recipe above) first as a fast
gross-regression tripwire; if skipped, nothing downstream changes.

### 8.3 Restatement targets — every artifact updated with the corrected floors (gates the deploy)

After the scan, restate ALL of:
- **`docs/STATE_OF_PROJECT.md`** — incl. the now-WRONG factor-of-2 claim near **line 29** ("do not
  fix ε_K by ×2" — the audit/review CONFIRM the ÷2-and-swap is required; this line is wrong and must
  be corrected), and the 25–30 / 2–3 / 7 TeV headline floors (line ~5) + the proxy-exclusion note
  (line ~142). Replace with grid-honest brackets (M4 below) and the corrected drivers.
- **`docs/quark_scan_methodology_note.tex`** — the prose floors and the 47.26/23.37 TeV crossings
  (lines ~587/654/703/995/1004), restated with **grid-honest brackets** (M4): the grid is
  `{1,2,3,5,7,10,15,20,30,50}` TeV (no 25 point), so the minimal strict floor is "in (20,30] TeV for
  r≥0.1, (15,20] TeV for r=0.05", NOT "25–30"; the "25–30" precision came from a literature
  10–12 TeV ×2.45 conversion, not the scan grid — say so. Custodial strict floor with the M3
  top-partner-loop caveat.
- **The minimal-vs-custodial comparison** — `build_wq_quarkonly_comparison.py` artifacts
  (`survival_by_r_mkk.csv`, `constraint_veto_by_r_mkk.csv`, `paired_*.parquet`, `manifest.json`,
  `README.md`) regenerated from the new paired roots.
- **The Scan Explorer data/JSON** — `flavor_catalog/website/src/content/scan_explorer.json` via
  `build_scan_explorer.py` (edit `MINIMAL_ROOT`/`CUSTODIAL_ROOT` literals to the new job ids first),
  and any `explore.astro` copy referencing the old numbers.
- **`.orchestration/PHASE2_PROGRAM_LEDGER.md`** — floor lines (~81/98/100), incl. the "dominated by
  the m_b² FULL-FLAVOR-SUM Casagrande admixture" note (now a corrected ~15% gauge correction).
- **The three `review_local/*.tex` notes** (the DEFENDANT, §0.1) — correct the O4/O5 coefficients,
  the (2/3)→(1/3) ⟨O1⟩, the "no 1/(2m_M)" claim, the Zbb `B(c,F)` bracket, and the "DO NOT fix the
  1/6 / a ×2 breaks literature" callout (half-wrong: clarify the ÷2 is on the MATRIX ELEMENT and is
  REQUIRED; the 1/6 *Wilson* is untouched).

**Deploy gate (load-bearing):** the Cloudflare site deploys from `main` (root
`flavor_catalog/website/`). Per the ledger the user does NOT want the Explorer tab deployed yet, so
the rebuilt `scan_explorer.json` and any website prose change **must NOT be pushed/merged to the
deploying `main` without explicit user go**. The restatement commit can land on `main` for the docs
but the website-deploy surface needs the user's explicit deploy approval — flag this to the
orchestrator/user at restatement time.

### 8.4 Expected post-fix picture (to CONFIRM by the scan, not assert)
Minimal: T010/Zbb collapses as the driver (B1 makes the fermion piece a correct ~15% of the gauge
piece; M1 gate removes the 0.004σ knife-edge → T010 floor O(few TeV)); corrected ε_K strengthens
~×1.31 and likely becomes the (now convention-stable, B2) dominant minimal strict floor; EW001
(proxy, M2) ~15–17 TeV but excluded from the strict lane. Custodial: strict floor driven by
corrected ΔF=2, with the 2–3 TeV claim suspended pending top-partner loops (M3). The
minimal-vs-custodial qualitative contrast survives; **every headline NUMBER changes** — read the
actual values off the scan.

---

## 9. Commit structure + sequencing (v2 — code first, scan, then restate)

Each commit message records plan approvers (codex+opus verdicts/SHAs) and impl approvers, per the
dual-signoff ledger. **All on a feature branch off `main` (not directly on `main`); the
website-deploy surface is gated on explicit user go (§8.3).**

**PHASE I — code fixes + tests GREEN (one commit per logical item):**
1. **`fix(quark): B2 — rephase SVD factors to PDG CKM convention (ε_K convention-stable)`**
   — `fit.py` rephasing (branch-free determinism recipe §1.2) + B2 absolute pins
   (rephasing-invariance, **determinism** bit-identical, PDG-convention, invariants).
2. **`fix(quark): B3 — ΔF=2 O4/O5 un-swap + 1/(2m_M) (all 6 ME sites, together)`**
   — `deltaf2.py` (sites 1–2) + `modern/phenomenology.py` (sites 3–6, incl. `_meson_matrix_
   elements_inline`) + B3 absolute pins (coeff ratio, singlet>crossed, SM-box anchor, abs ME) +
   re-pinned snapshots (ε_K 1.9286→≈3.3225; K001 line 197 unrun / 198 run; B003; C002 — all
   re-derived at the benchmark). ONE commit (swap + ÷2 inseparable, §0.2).
3. **`fix(rs-ew): B1 — retranslate CGHNP Zbb fermion-KK admixture (sign + flavour-sum index)`**
   — `rs_ew_couplings.py` bracket + flavour sums + B1 literature-anchored pins (replace `_manual_B`
   with the convention-dictionary check, refinement 4) + re-pinned T010/T011 snapshots.
4. **`fix(scan): M2 — EW001 ΔT uses physical-M_KK consistently with ΔS calibration`**
   — `oblique_stu.py` / `EW001.py` convention + a convention pin.
5. **`fix(scan): M5 — structured tag-class hint replaces brittle proxy/rigorous substring match`**
   — `run_full_catalog_scan.py` tag plumbing + explorer/comparison alignment (F7) + harness tests.
6. **`fix(rs-ew): M1 — T010 R_b veto uses the T011-style loose-edge NP-shift budget`**
   — `T010.py` veto scalar → loose-edge budget (factored from / mirroring `T011._build_budget`) +
   pin (SM-limit shift=0 ⇒ ratio 0 ⇒ SM passes by construction; known RS shift → known ratio) +
   re-pinned T010 snapshots. Lands AFTER B1 (the corrected couplings) since the gate operates on
   them; the two together turn the spurious 108 TeV into the O(few TeV) physics floor.

**GATE:** the FULL existing suite (after each item re-pins its own snapshots) **plus** all new
absolute pins must be GREEN before Phase II. No scan is launched against red tests.

**PHASE II — run the scan (§8):**
7. Run the **full paired 1M reproduction** (§8.1) — this is the **sole source of ALL headline,
   inclusive, and strict floors** (v3 BLOCKER resolution) and produces every artifact in native
   format. Optionally run the §8.2 strict-ΔF2-only order-of-magnitude tripwire FIRST per its
   corrected recipe (per-gate exponent + per-draw skip mask) as a fast gross-regression check; it
   sets NO headline number and may be skipped. This phase produces artifacts on disk; no
   production-code commit.

**PHASE III — scan-artifacts + restatement (ONE commit):**
8. **`docs(quark): restate corrected M_KK floors from re-scan + correct the DEFENDANT notes`**
   — the regenerated comparison/explorer artifacts + ALL §8.3 restatement targets
   (STATE_OF_PROJECT.md incl. line-29 factor-of-2 correction, methodology note with grid-honest
   brackets, ledger, the three review_local notes, scan_explorer.json). Website-deploy surface
   gated on explicit user go (§8.3) — do NOT push the deploying website change without it.

M6/M8/M3/M4: M3/M4 prose ride with commit 8; M6 (document-only caveat) rides with commit 2's docs
or its own tiny doc commit; M8 recorded as an open low-priority follow-up.

---

## 10. Summary of decisions

| Item | Decision | One-line justification |
|---|---|---|
| B2 ε_K SVD phase | FIX (approach A: PDG rephasing, branch-free determinism recipe) | prerequisite for any ε_K statement; fixes all phase-sensitive gates in one place |
| B3 O4/O5 + 1/(2m_M) | FIX TOGETHER (un-swap + ÷2, all 6 sites) | GGMS Eq.(8) confirms swap+double; sub-bugs cancel if split |
| B1 Zbb admixture | FIX (correct, don't drop) | genuine +15% correct term; correcting is honest and needed to state any bound |
| M1 T010 gate | **FIX (IN SCOPE): adopt T011-style loose-edge NP-shift budget on R_b** | R_b half of the same Z→bb̄ measurement must use the same gate as the A_b/A_FB half; CGHNP §6.4 constrains the RS *shift*, not the SM's standing −0.996σ tension; 2σ-cut fallback |
| M2 EW001 ΔT | FIX | cheap, anti-conservative convention bug; makes inclusive floor honest |
| M5 tag substring | FIX (structured tag-hint) | integrity/silent-non-veto risk; durable fix over prose match |
| M6 bag scales | DEFER (document only) | separate scheme question, conservative for kaon, needs proper NLO bag-matching |
| M3/M4 prose | DO at restatement (§8.3, in scope) | grid-honest brackets + custodial loop caveat, done when floors are restated |
| M8 qcd/ | DEFER | latent ≤2.2e-4, off the floor path |
| Re-scan + restatement | **IN SCOPE (v2): code→scan→restate**; **floors from full 1M grid ONLY (v3)** | full §8.1 paired 1M (~115 core-hrs, ~30 min wall) is the sole floor source; §8.2 demoted to optional ΔF2-only tripwire; site deploy gated on user go |
| §8.2 floor-extraction (v3 BLOCKER) | **DROP as floor source (option a)**; keep only as optional strict-ΔF2-only per-gate-exponent tripwire | blanket `M_KK·sqrt(ratio)` wrong for CR0xx (∝1/M_KK) & EW001 (∝1/M_KK⁴), biased for ΔF=2; full grid is cheap and trustworthy |

**PLAN v3 ready.**
