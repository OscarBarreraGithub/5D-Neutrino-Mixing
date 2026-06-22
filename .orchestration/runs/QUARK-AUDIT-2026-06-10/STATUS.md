# Quark-sector implementation audit — 2026-06-10

Goal: adversarial review of the quark-sector stack (constraints, model, scan
harness, QCD) to catch implementation mistakes/subtleties before the RS
"does it survive 2026 data" conclusions are trusted. User-requested; review
agents are READ-ONLY (no code changes without user sign-off).

## How to resume
Read this file, then the slice reports below. Launch any PENDING slice as a
general-purpose agent with the prompt stored in `prompts/` (verbatim), save
its final report to `slice-<n>-<name>.md`, update this table, and when all
slices are done write `SYNTHESIS.md` and report to the user.

## Slices

| # | Slice | Status | Report |
|---|-------|--------|--------|
| 1 | ΔF=2 Wilson matching, RG, ε_K normalization (deltaf2.py, qcd_running.py) | **DONE** | slice-1-deltaf2.md |
| 2 | MFV model + KK-gluon mass-basis couplings (model.py, couplings.py, scales.py) | **DONE** | slice-2-model-couplings.md |
| 3 | RS-EW / Zbb / custodial P_LR / oblique (rs_ew_couplings.py, T010/T011/T014, EW001) | **DONE** | slice-3-rsew.md |
| 4 | Catalog adapters + experimental anchors (K001, B002-4, C001-2, B011-12, YAMLs) | **DONE** | slice-4-adapters.md |
| 5 | Scan harness + M_KK/Λ_IR conventions (run_full_catalog_scan.py, explorer floors) | **DONE** | slice-5-harness.md |
| 6 | QCD running package (qcd/) | **DONE** | slice-6-qcd.md |

**ALL SLICES COMPLETE. Synthesis: SYNTHESIS.md (read that first on resume).**
Slice 6 summary: 0 BLOCKER, 2 MAJOR (wrong d3 3-loop mass-decoupling
constants, ≤2.2e-4 effect; pdg_quark_masses_at_scale crashes below 4.18 GeV
with PDG top input). alpha_s/γ_m running verified vs CRunDec to 1e-9,
1 GeV–50 TeV. Agent id af727ea41669273c3.

## Batch B summaries (full reports in slice files)

- **Slice 4 (0 BLOCKER, 2 MAJOR):** adapter/anchor layer sound — anchors
  current, B_K(2 GeV)=0.5503 correctly paired (no RGI mixup), ps⁻¹↔GeV exact,
  Δm→M12 halving once everywhere, B002/B004 phase logic correct. Quantified
  upstream interaction: kaon LR M12 understated 0.58× (K001 anti-conservative,
  ε_K floors ~1.3× low), VLL overstated 2.0×, B_s/B_d ≈1.06×, D0 ≈0.92×.
  Budgets exp/SM-anchored, NOT calibrated to buggy ⟨O⟩ ⇒ fix core + re-run,
  no budget re-derivation. MAJOR: bag-scale pairing (3 GeV bags at 2 GeV
  Wilson point, ~19% on kaon O4; B-meson bags at m_b ~7%); tests pin
  adapter≡core with no literature-anchored absolute pins (bugs structurally
  undetectable, 77/77 pass). Agent id a2c74aeeae398ff88.
- **Slice 5 (1 BLOCKER-interaction, 3 MAJOR):** harness mechanics sound; no
  additional Λ_IR/M_KK leaks beyond EW001; drops = perturbativity skips
  clustering at low r (not low M_KK), denominators exclude skips; pairing
  verified at draw level (2-ppm fit nondeterminism); floors threshold-
  insensitive; reports match artifacts. BUT: minimal floor set 100% by buggy
  T010 (694,123 vetoes, only rigorous veto above ~3 TeV); "25-30 TeV" is a
  literature 10-12 TeV ×2.45 conversion, not a scan measurement (grid only
  supports "(20,30]"); r=0.05 floor ≤20 TeV unstated. MAJOR: T001/T002
  evaluated+passing but mis-tagged partial→hard_not_evaluated on 100% of rows
  ("proxies" misses "proxy" substring) ⇒ coverage_complete=False everywhere,
  a failing partial-HARD would silently never veto. MAJOR: EW001/CR001/12/13/
  B012 veto counts bit-identical minimal vs custodial (EW001=352,491 both) ⇒
  "custodial inclusive 7 TeV" is a model-independent proxy wall inheriting
  the EW001 ×6 bug, not custodial physics. Agent id a3fe349291685f7f6.

## Batch A summaries (full reports in slice files)

- **Slice 3 (1 BLOCKER, 3 MAJOR):** Casagrande m_b² fermion-KK Zbb admixture
  mistranslated from CGHNP 0807.4937 (`rs_ew_couplings.py:1835-1847,957-971`):
  wrong per-generation 1/F² factor, wrong c-sign, missing F²=2f² factor ⇒
  δg_L^b wrong SIGN and ~190× too big; δg_R^b ~2000× off. This term dominates
  the minimal 25-30 TeV floor ⇒ **minimal floor is an artifact**; corrected
  minimal Zbb is gauge-dominated, ~25× weaker. MAJOR: T010 1σ gate leaves
  0.004σ headroom (SM R_b pull −0.996) — corrected floor swings 108 TeV (1σ)
  vs 6.8 TeV (2σ); EW001 ΔT uses Λ_IR-coefficient with physical M_KK ⇒ T
  ×6 underestimated (proxy lane); "custodial strict 2-3 TeV" violates the
  program's own top-partner-loop caveat (CORRECTED_PRESCRIPTION). Gauge-KK
  Zbb piece itself independently re-derived and verified exact.
  Agent id a171609616c2ebce9.
- **Slice 1 (1 BLOCKER, 1 MAJOR):** O4/O5 chiral-enhancement coefficients are
  BMU/SUSY basis-swapped (`deltaf2.py:719-724,908-913` + 3 vendored copies in
  modern/phenomenology.py): code (R/6+1/4)B4,(R/2+1/12)B5 vs correct
  (R/4+1/24)B4,(R/12+1/8)B5. MAJOR: all ⟨O⟩ are 2× M12-ready values (1/(2m_M)
  applied zero times) — contradicts STATE_OF_PROJECT.md:29 "no extra factor
  of 2". Net ε_K ×1.72 underestimate ⇒ ε_K floors ~31% too LOW
  (47.26→≈62 TeV p50; CFW-matched 23.4→≈30.7). B_d/B_s nearly cancel
  (1.02-1.08), D0 ≈0.9. Fix F1+F2 TOGETHER. RG/matching layer verified sound.
  Agent id ae4caca93ff50d30e.

## Headline findings so far

- **F1 (BLOCKER, high confidence, slice 2):** ε_K depends on arbitrary SVD
  column phases — `quarkConstraints/fit.py:245-258` never rephases the SVD
  factors to a physical CKM phase convention, so Im(M12^NP) (the ε_K input)
  is convention-dependent. Numerically demonstrated: a physically-equivalent
  rephasing changed ε_K ratio_to_bound 17.34 → 105.1 (×6) with masses, |CKM|,
  Jarlskog, and all |M12| observables bit-identical. ε_K binds ~99.5% of the
  ensemble ⇒ headline floors are percentiles of a convention-dependent number.
  Fix: rephase to PDG CKM convention after SVD, or use a rephasing-invariant
  ε_K construction relative to arg[(V_ts* V_td)²].
- F2 (documented convention): default g_s* = sqrt(4πα_s) ≈ 1.05 vs headline
  g*=3 vs CFW matching ~×2 higher; legacy scan.py/validation.py lanes hardcode
  the weak perturbative convention.
- F3 (latent): universal −g_s/√(2πkr_c)·𝟙 KK-gluon piece omitted; harmless for
  ΔF=2 off-diagonals, but diagonal couplings wrong in sign/magnitude (~10³)
  for any future width/interference consumer.
- F4 (minor): 1/ROOT² ratio-rescale in export_collaborator_5tev_points.py:313
  ignores α_s/RG scale shift (few-%).

## Context for the synthesizer

Headline claims under test: minimal RS physical-M_KK floor ~25-30 TeV
(T010/Zbb-driven); custodial strict floor ~2-3 TeV, inclusive ~7 TeV; scans =
two 1M quark-only runs (20128400 minimal, 20675555 custodial), grid
r×M_KK, ξ_KK = 2.4487, strict = rigorous-HARD-only. Key docs:
docs/STATE_OF_PROJECT.md, docs/quark_scan_methodology_note.tex,
.orchestration/PHASE2_PROGRAM_LEDGER.md.

Note: slice-2 agent (id aa4a9ba5107a5fc73) cited a methodology-note appendix
with "47.26 TeV p50 at g*=3" and ε_K binding 99.5% — reconcile that with the
25-30 TeV STATE_OF_PROJECT figure during synthesis (different lanes/conventions).
