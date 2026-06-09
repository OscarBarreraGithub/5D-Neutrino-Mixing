# W7 plan REVISION — address dual-review fixes (codex + Opus, both NEEDS-FIXES)

Revise `.orchestration/runs/W7-MUEGAMMA/plan_codex.md` IN PLACE to resolve ALL items below. Both
reviewers verified the dipole reuse, the numeric oracle (BR_NP=1.5508e-10, BR/limit=1033.89), and
the quark-only byte-identity approach — keep those. The μ→eγ model stays LOCKED = LMFV
(Perez-Randall, spurion (Y_N Y_N†)_12). Re-verify each fix against the REAL repo. End with `PLAN-READY`.

## CRITICAL honesty framing (both reviewers, do this carefully)
The catalog μ→eγ veto is driven ENTIRELY by `BR_NP <= br_limit` (the predicted branching ratio vs
MEG II `1.5e-13`). The injected `C=0.02` (C_PAPER) only feeds the `dipole_rhs`/`dipole_ratio_to_bound`
DIAGNOSTIC and has NO effect on pass/fail. Note also: C=0.02 is paired with the paper-era limit
(1.2e-11); the MEG-II-consistent coefficient is ≈0.00194 (= the `C≈1.936e-3` already documented in
`flavor_catalog/processes/charged_lepton/L001.tex:40`). The plan MUST:
- State explicitly that C=0.02 is an inert diagnostic and does NOT drive the veto.
- Correct any plan text claiming the oracle "checks both the C=0.02 dipole RHS and the BR veto" as
  if both gate — only the BR gates.
- Make clear the three coexisting ratios (1033.89 BR-form / 3.11 paper-C dipole-form / 32.15
  consistent-C dipole-form) measure DIFFERENT things.
- Add a test pinning `passes` to the BR ratio INDEPENDENT of `c_lfv`.
- Do NOT silently overwrite the catalog L001 physics/text: keep the existing generic-bound
  description and the new LMFV-prediction path explicitly separate; if there is any genuine
  ambiguity about which the catalog should report, FLAG it (do not decide it).

## From codex review
1. Fix the files-touched list re L001.tex (`:40`): the catalog text documents `C≈1.936e-3` as the
   live default; W7 must NOT silently make catalog L001 use c_paper=0.02 — explicitly separate
   legacy `scanParams` behavior from catalog W7 behavior (ties to the honesty framing above).
2. Don't mix mass conventions in the oracle: plan says carrier uses physical
   `base_spectrum.kk_ew_mass_gev`, but the oracle uses `Lambda_IR=3000` AND `M_KK=3000`. Split into
   (a) a direct adapter/carrier oracle with explicit `M_KK=3000`, and (b) a separate builder test
   with scan-consistent `lambda_ir = mkk/xi` (so `build_from_rs_ew_inputs` first-KK mass is right).
3. The compatibility fallback marked `used_proxy=False` could be classified rigorous by the harness.
   Either remove the L001 fallback, or route it through the same validated LMFV carrier contract and
   record `extra_used`; otherwise tag it partial/transition, NOT rigorous.
4. Strengthen the quark-only regression to EXACT serialized-row comparison for a deterministic draw
   (not just config-hash + selected fields).

## From Opus review
5. (Covered by the honesty framing above — make the C=0.02-inert / BR-driven distinction explicit +
   the `passes`-pinned-to-BR test.)
6. Justify the NEW `LMFVLeptonParameters` carrier vs simply reusing the existing
   `RSLeptonMassBasisCouplings` (rs_ew_couplings.py:98-126) which ALREADY carries Y_N_bar/Y_N_matrix/
   pmns/lfv_dipole_spurion/kk_ew_mass_gev/c_L/c_E/c_N/M_N/params and lacks only `epsilon` (which the
   dipole does not use). Either document why a separate carrier is preferable (e.g. W8 isolation of
   rs_ew_couplings.py) OR reduce the duplication. Decide and state.
7. Add `"lepton_lmfv_parameters"` (or whatever the chosen extra key is) to
   `QUARK_ONLY_FORBIDDEN_EXTRAS` (run_full_catalog_scan.py:75-82) for documentation/consistency.
8. Require the carrier path to emit NO truthy "proxy"/"recast"-named diagnostic (keep
   `used_proxy=False`) so the `rigorous` tag is robust under `_proxy_flags`
   (run_full_catalog_scan.py:1045-1057).

End the revised plan with `PLAN-READY`.
