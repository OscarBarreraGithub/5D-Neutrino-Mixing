# WQ — QUARK-ONLY (Bucket-1) MINI-SCAN MODE (codex author, gpt-5.x xhigh)

Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. You are the AUTHOR. First a SHORT plan, then implement code+tests. A separate codex reviewer AND an Opus reviewer will independently review; BOTH must APPROVE before commit. Build on the dual-approved harness `scripts/run_full_catalog_scan.py` (commit df947f3) and the W6a perf hook (RSEWSpectrum + RSEWOverlapSplineCache injection).

## GOAL
A fast, HONEST quark-sector-only mini scan: vary quark Yukawa anarchy + geometry (Lambda_IR / M_KK), DROP the lepton sector entirely (no lepton draws, no compute_all_yukawas, no lepton-perturbativity gate, no lepton extras), and evaluate ONLY the constraints that need ZERO swept-lepton data ("Bucket 1"). Purpose: see what the quark sector looks like — M_KK reach vs anarchic Yukawa structure — and later rank how constraining each process is. NOT claimed fully rigorous (lepton sector ignored); must be tagged as such.

## BUILD
Add a `--quark-only` mode (CLI flag + ScanConfig field; default OFF → existing full-catalog behavior BYTE-IDENTICAL). In quark-only mode the per-draw path:
1. Draw ONLY the quark anarchic Yukawa seed + geometry; NO `_draw_lepton_inputs`, NO `compute_all_yukawas`, NO `_require_perturbative_leptons`.
2. Run the existing quark fit (`fit_quark_sector(default_quark_targets(), …, fit_orientation=True)`) → real `QuarkFitResult`; keep the c=0.5 singularity guard on the QUARK c-values only.
3. Build the point via `build_from_rs_ew_inputs` with lepton pieces OFF and `lepton_yukawa_result=None`; set `include_charged_current` / `include_fermion_kk_mixing` / `include_higgs_yukawas` to whatever the ALLOWLIST below actually requires (see §ALLOWLIST — fermion-KK Zbb is QUARK-side b and IS needed for T010/T011; charged-current G_F touches muon decay → treat carefully). Inject the per-tile `spectrum=`/`rs_ew_cache=` exactly as the full path does.
4. Evaluate ONLY the allowlist constraints (filter `registry.evaluate_all`, or evaluate all then keep allowlist — your choice, but non-allowlist constraints must NOT appear as vetoes or as spurious stubs; they are simply OUT OF SCOPE for this mode and tagged `deferred_lepton_followup`).

## ALLOWLIST — YOU MUST VERIFY, NOT ASSUME
Candidate Bucket-1 (need NO swept-lepton data): ΔF=2 K001 K002 B001 B002 B003 B004 C001 C002; radiative B011 B012 B013 B014; hadronic CP B032 B033 B034 C003 K003 K013; top FCNC T001–T008; Z→qq̄ T010 T011 T012 T014; oblique/CKM EW001 EW002 EW003; quark/nuclear EDM E004 E006 E007 E008 E009.
For EACH candidate: OPEN its `evaluate` body + the extras it reads. CONFIRM it reads no swept-lepton coupling (z_delta_g_*_e, z_delta_g_L_nu, rs_higgs_yukawas, lepton_mass_basis_couplings, PMNS, lepton Yukawas). If a candidate genuinely needs a lepton coupling (likely suspects: EW002/EW003 first-row CKM & G_F via muon decay; any charged-current G_F normalization), you MUST either (a) keep it ONLY if the lepton side is a fixed SM input independent of the dropped sweep (document precisely why it is SM-input, not swept), or (b) DROP it from the allowlist and tag it deferred. Produce a table: constraint → extras read → lepton-dependent? Y/N → IN/OUT → one-line justification. The reviewers will recompute this independently; a single mis-classified constraint that silently uses an SM/zero default is the main failure mode — be conservative.

## OUTPUT WIRING (for "how constraining is each process")
The per-tile summary + the merged `run_summary.json` MUST include a per-constraint tally: for each allowlist id, count of points where it was active, evaluated, and FAILED (vetoed), plus its severity/tag — so we can rank constraining power. Keep the existing honest strict/inclusive HARD bookkeeping but restricted to the allowlist. Record `mode="quark_only"`, the allowlist, and `lepton_sector="dropped (not rigorous)"` in provenance.

## GATES / TESTS
- New unit test(s): quark-only mode runs, lepton stages are skipped, only allowlist ids appear, none degrade to silent SM defaults, determinism (same seed → same params/ratios), full-mode path unchanged.
- A SM/sanity check: an anarchic point with a high M_KK (e.g. 10–20 TeV) relaxes ΔF=2 vetoes vs a low M_KK (~1–2 TeV) — confirm the expected M_KK monotonicity (ε_K/Δm bite harder at low M_KK).
- Run a SMALL smoke (e.g. 1e3 points) NOW: report yield (expect ~100%), post-cache s/point, evaluated/active counts for the allowlist, and the per-constraint veto ranking. Do NOT run the 100k (the orchestrator launches that after the gate).
- `python -m pytest tests/ -q` stays green.

## CONSTRAINTS
Physics ONLY via existing adapters — do NOT change constraint physics or the builder's quark math. Touch: `scripts/run_full_catalog_scan.py` (+ its config/CLI) and tests only. Numeric outputs finite/sane.

## OUTPUT (≤18 lines)
Short plan; the allowlist table (id→extras→leptonDep→IN/OUT→why), flagging any dropped candidate; the mode/flag + which build-include flags you set and why; smoke results (yield, s/point, 1e3 extrapolated to 1e5, per-constraint veto ranking, M_KK monotonicity sanity); pytest counts. End with: WQ-AUTHOR-DONE.
