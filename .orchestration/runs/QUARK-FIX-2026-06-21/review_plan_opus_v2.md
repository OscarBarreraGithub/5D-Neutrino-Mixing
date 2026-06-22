# Independent review of PLAN v2 deltas — Opus reviewer lane

Scope: re-review of v2 deltas ONLY (M1 gate, re-scan recipe, refinement fold-in, restatement/
sequence). B1/B2/B3 algebra was verified in v1 (APPROVED) and is not re-derived here. READ-ONLY on
production code. HEAD `5328a22`.

**VERDICT: REVISE** — 1 blocking issue (the cheap floor-extraction recipe, §8.2, is numerically
wrong for the gates that set the *inclusive* floor and biased for the strict floor). M1 gate
decision is sound and source-defensible (APPROVE on that axis). Refinements folded correctly
(one minor literal imprecision, non-blocking). Restatement/sequence complete.

---

## 1. M1 gate decision — APPROVE (sound, source-defensible)

### (a) Is this how T011 gates today? YES — verified against code.
- `T011.py:_asymmetry_shift` (560–582): `np_shift = predicted − sm_prediction`;
  `ratio = |np_shift| / budget.hard_veto_budget`. `_build_budget` (481–515):
  `central = |exp − sm_fit|`, `combined = sqrt(σ_exp²+σ_sm²)`, `hard_veto_budget = central+combined`.
  `evaluate`: `selected = max(shifts, key ratio)`, `passes = selected.ratio ≤ 1` (830–831).
- T010 today: `_observable_pull` (T010.py:396–406) `pull=(predicted−experimental)/σ_combined`;
  `evaluate` `selected = max(pulls, key abs_pull)`, `ratio = abs_pull`, `passes = ratio ≤ 1`
  (648–655). T010's `_build_budget` (325–354) computes ONLY `combined_sigma` — it does NOT carry
  a `central_residual`/`hard_veto_budget`. The plan's "veto on max|pull| vs experiment" description
  is exactly correct, and its instruction to "factor out of T011 or mirror `_build_budget`" is
  correct: the implementer MUST add the `central = |exp − sm_fit|` term to T010's budget builder
  (it is not present today). Plan is aware of this. NOT a gap.

### (b) Is it physically right to gate the NP *shift* vs the EW-fit tolerance? YES.
This is the standard way RS Zbb reaches are quoted (CGHNP 0807.4937 §6.4 constrains the KK-induced
*shift* in the Zbb pseudo-observables against the EW fit's tolerance, not whether the SM itself
agrees with data). A `max|pred − exp|/σ` gate conflates the standing SM↔data tension (the
−0.996σ A_FB^b / R_b legacy anomaly, which RS does not author) with the NP exclusion. The
loose-edge budget `|exp−sm_fit| + sqrt(σ_exp²+σ_sm²)` cleanly absorbs the standing residual into the
tolerance so the SM-limit point passes by construction (shift=0 ⇒ ratio 0), while a genuine large RS
shift is still excluded. The decisive internal-consistency argument is airtight: **T010 (R_b) and
T011 (A_b/A_FB^b) are two halves of the SAME LEP/SLC Z→bb̄ measurement**; having them use two
different veto statistics is itself the bug. T011 is already dual-approved + rigorous-tagged with
this budget. Making T010 match is the correct unification, not a loosening for convenience.

### (c) Claimed floor behavior consistent? YES (modulo "confirm at scan", which the plan already
requires). The mechanism is correct: with the bare 1σ-vs-experiment gate and ~0.004σ headroom on the
low-R_b side, any standard-direction RS shift fails until 1/M_KK² decays it below 0.004σ → the 108 TeV
artifact. The ×16 swing to ~6.8 TeV under a 2σ gate is exactly the diagnostic that the number measures
anchor-snapshot tension, not RS. Under the loose-edge budget the SM-limit shift is 0 ⇒ passes by
construction ⇒ the 108 TeV collapses and T010 stops driving; the minimal strict floor is then set by
corrected ε_K (B3 ×1.31 up) + EW001. The plan correctly tags every NUMBER as "read off the scan, not
asserted." One caveat the reviewer flags as NON-BLOCKING: the plan's −0.996σ uses T010's
`sm_pull=(sm_prediction−experimental)/σ_combined` (line 401), where `sm_prediction` is the SM-LIMIT
formula R_b, not the YAML SM-fit anchor; these differ slightly. The 108-vs-6.8 TeV qualitative story
does not depend on the exact decimal, so this is immaterial to the decision — but the impl pin should
recompute the headroom number from code, not transcribe "0.004σ".

### (d) 2σ fallback sound? YES. A symmetric 2σ-vs-experiment cut still removes the knife-edge
artifact (floor ~6.8 TeV at the representative point) and is source-defensible. Correctly positioned
as the off-ramp, with the bare 1σ gate correctly ruled out post-B1.

**Not too permissive, not indefensible.** The budget is *wider* than 1σ_exp (adds the central
residual + σ_SM in quadrature) precisely so it does not manufacture exclusion from the SM's own
offset, yet still excludes a large RS shift. This recovers CGHNP-order R_b reaches rather than a
0.004σ artifact. APPROVE.

---

## 2. Scan recipe (cheap floor-extraction §8.2) — BLOCKING. The 1/M_KK² rescale assumption FAILS.

The plan asserts (§8.2:707–710) the per-draw rescale `M_KK_min = M_KK_tile · sqrt(max_X ratio_X)` is
"valid for the FULL binding set here" because "δg ∝ 1/M_KK²; ΔT ∝ 1/M_KK²". Verified against the
actual constraints that bind in the ensemble (`constraint_veto_by_r_mkk.csv` shows the binding set is
NOT just ΔF=2 — it includes T010/Zbb, EW001/oblique, and 8 collider-resonance CR0xx bounds). The
assumption breaks in three independent ways:

**(B-1) Collider-resonance CR0xx (set the INCLUSIVE floor) — ratio ∝ 1/M_KK, NOT 1/M_KK².**
`collider_resonance.py:166–168`: for a `MASS_LOWER_BOUND`, `ratio = limit.value / prediction.mass_tev`
and `mass_tev ∝ M_KK` (`:109–112`, `:126`). So the veto ratio is **linear** in 1/M_KK. Applying
`M_KK·sqrt(ratio)` gives `M_KK·sqrt(m_lim/M_KK) = sqrt(M_KK·m_lim)` — wrong. The correct per-draw floor
for a hard mass threshold is `M_KK_min = m_limit` (a constant), i.e. `M_KK · ratio`. The √ rescale
mis-locates these floors. CR001/002/003/004/008/010/012/013 all veto in the inclusive lane.

**(B-2) EW001 oblique T/S (inclusive floor) — ratio ∝ 1/M_KK⁴ (χ² quadratic form), NOT 1/M_KK².**
`oblique_stu.py:compare_oblique_st_to_fit` (252–267): `ratio = chi2/budget` with
`chi2 = ds·Σ⁻¹·ds + …` a quadratic form in (ds,dt), and ds,dt each ∝ 1/M_KK². So `ratio ∝ 1/M_KK⁴`
and `sqrt(ratio) ∝ 1/M_KK²` ⇒ `M_KK·sqrt(ratio) ∝ 1/M_KK` — wrong. The correct rescale for a
1/M_KK⁴ χ² ratio is `M_KK · ratio^(1/4)`. (Note: EW001's RS volume log is hardcoded
`DEFAULT_RS_VOLUME_LOG=35.0`, `oblique_stu.py:47`, so there is NO log(M_KK) drift in ΔT/ΔS — that
part of the plan's claim is fine. The break is the χ² nonlinearity, not a log.)

**(B-3) Tooling cannot consume the binding set; strict gates carry running/ε drift; skip set is
M_KK-dependent.**
- The cited tool `rs_anarchy_mkk_min_hist.py` rescales ONLY a `deltaf2_ratios` field (the 5 ΔF=2
  systems) that exists in the *separate* `run_rs_anarchy.py` pipeline. The full-catalog JSONL written
  by `run_full_catalog_scan.py:891–901` stores ONE opaque scalar `ratio` per process_id and NO
  per-system breakdown — so the cited tool cannot read T010/EW001/CR0xx ratios in any case.
- Even the strict ΔF=2 gates are only leading-order 1/M_KK²: Wilsons are matched at `μ=M_KK` and
  QCD-run M_KK→2 GeV (`deltaf2.py`/`qcd_running.py`), and `g_s²=4π·α_s(M_KK)` (`couplings.py`),
  giving an α_s(M_KK)/running log drift; T010 δg carries an ε=Λ/k f-factor drift and its veto ratio
  (R_b, A_b) is a nonlinear function of δg, not linear. These are log-level, not order-of-magnitude,
  but they bias a "trustworthy floor."
- The quark-fit perturbativity SKIP (`run_full_catalog_scan.py:830–842`) depends on Yukawas fitted at
  the tile Λ_IR, so the skip SET is M_KK-dependent. A single-tile rescale evaluates the skip mask at
  one M_KK, not at each draw's rescaled M_KK_min, biasing the survival-fraction crossing.

**Required change:** Do NOT present §8.2 as a trustworthy floor extractor for the headline. Either
(i) DROP the cheap option as a floor source and require the full paired 1M grid (§8.1) for ALL floors
(it is only ~30 min wall / ~115 core-hrs — cheap enough that the §8.2 shortcut buys little); OR
(ii) RESTRICT §8.2 explicitly to an order-of-magnitude cross-check of the ΔF=2-dominated STRICT floor
ONLY, with the rescale exponent chosen PER GATE (`M_KK·ratio^{1/2}` for ΔF=2 leading-order,
`M_KK·ratio` for CR0xx mass thresholds, `M_KK·ratio^{1/4}` for EW001 χ²), evaluate the skip mask at
each draw's rescaled M_KK, and STATE that it does not set the inclusive floor (CR0xx/EW001) and is not
the headline source. The plan's blanket "valid for the FULL binding set" must be struck. The headline
floors must come from §8.1, which the plan already runs — so this is a correction to the §8.2 framing
+ a guard against using its output as a floor, not a new compute burden.

---

## 3. Refinements folded correctly — APPROVE (one minor literal imprecision, non-blocking)

- **R2 K001 path/line:** CONFIRMED. `tests/constraints/primary/kaon/test_K001.py:197` =
  `unrun.ratio_to_budget == 0.26555011996248595`; `:198` = `run.ratio_to_budget == 2.844741580402456`.
  Plan's v2 line→value map (197=unrun, 198=run) is correct; v1's reversal is fixed.
- **R3 site-6 name:** CONFIRMED. `_meson_matrix_elements_inline` at phenomenology.py:566;
  `_m12_np_from_bridge_generic` (:584) calls it at :599/:636 (transitive, no 7th site). deltaf2.py
  site 2 `_meson_matrix_elements` (:908–910) has `(2/3)`, `(r_chi/6+0.25)`, `(r_chi/2+1/12)` — the
  exact wrong literals. All six sites enumerated correctly.
- **R1 B2 determinism:** the branch-free recipe (z/|z| not np.angle, `>=0.0` signed-zero tie-break,
  single `arg(V_ud)=0` anchor, documented tol fallback) + the NEW determinism bit-identical pin
  (§6.2) is sound for the paired-scan `(r,mkk,draw_seed)` join. APPROVE.
- **R4 in-test anchoring:** §6.1/§6.3 correctly mandate computing expected pins from the convention
  dictionary / GGMS rationals / CGHNP map, NOT off the production line (the M7 self-pin failure mode).
  APPROVE.
- **R5 D0/Bd/Bs re-derive:** §2.3 correctly mandates re-deriving C002/B003 literals at the benchmark.
  NON-BLOCKING imprecision: §6.4 lists "test_B003.py:151,298 (≈×1.05)" — but line 151
  (`budget == 2.635167648629676e-13`) and 298 (`budget == …`) are exp/SM-derived ANCHOR BUDGETS, which
  per slice-4 I-1 are NOT ⟨O⟩-dependent and should NOT move. Only `result.predicted` (4.76e-14) and
  `result.ratio` (0.18065) move. Same for C002: `anchor.budget`/`anchor.value`/`m12_budget_gev` are
  anchors (fixed); only the predicted/ratio move. Tighten §6.4 to re-pin the *predicted/ratio* lines,
  not the budget lines, else the impl may "re-derive" a literal that should be invariant. Cosmetic.

---

## 4. Restatement + sequence — APPROVE (complete)

- **STATE_OF_PROJECT "do not fix ε_K by ×2":** CONFIRMED present (the "Do not 'fix' epsilon_K by
  multiplying by 2" DEFENDANT line). §8.3 correctly targets it for correction (the ÷2-and-swap IS
  required per v1-approved B3). Good.
- **Methodology grid-honest brackets:** CONFIRMED the note carries `47.26 TeV / 127.13 TeV` crossings
  (`docs/quark_scan_methodology_note.tex:587`) and other prose floors; §8.3 correctly restates these
  as grid-honest `(20,30]`/`(15,20]` brackets given the grid `{1,2,3,5,7,10,15,20,30,50}` has no 25.
- **Scan Explorer JSON:** CONFIRMED `build_scan_explorer.py:41–42` hardcodes `MINIMAL_ROOT`/
  `CUSTODIAL_ROOT` job ids that must be edited; `_build_bare_floor` (:410), `FLOOR_THRESHOLD=0.5`
  (:47) as described. Deploy gate (do not push website surface without user go) correctly preserved.
- **Sequence:** PHASE I (all code fixes + M1 + tests GREEN) → PHASE II (scan) → PHASE III (restate),
  one restatement commit, M1 lands after B1. Order is correct and the GATE ("no scan against red
  tests") is explicit. APPROVE. (The §8.2 framing fix from §2 above does not change the sequence.)

---

## Blocking summary
1. **§8.2 cheap floor-extraction** — `M_KK·sqrt(max ratio)` is WRONG for CR0xx (ratio∝1/M_KK) and
   EW001 (ratio∝1/M_KK⁴), and biased for ΔF=2/T010 (running/ε drift + nonlinear ratio) and by the
   M_KK-dependent skip set; the cited tool cannot even read those ratios. Strike "valid for the FULL
   binding set", make §8.2 a per-gate-exponent, strict-ΔF2-only cross-check (or drop it), and source
   ALL headline floors from the full §8.1 grid.

## Non-blocking
- §5/M1: recompute the −0.996σ / 0.004σ headroom from code at impl, not transcribe; note
  sm_prediction (SM-limit formula) ≠ YAML SM-fit anchor.
- §6.4: re-pin C002/B003 *predicted/ratio* lines, not the exp/SM-derived anchor-budget lines
  (151/298, anchor.budget/value) which are ⟨O⟩-invariant.
- §5/M1 impl: T010's `_build_budget` lacks `central_residual`/`hard_veto_budget` today; the implementer
  must add them (plan already says "mirror T011" — just making the gap explicit).
