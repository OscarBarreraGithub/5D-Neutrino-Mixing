# W7 — IMPLEMENT the DUAL-APPROVED μ→eγ LMFV plan

The plan `.orchestration/runs/W7-MUEGAMMA/plan_codex.md` is DUAL-APPROVED (codex + Opus). Implement
it EXACTLY as written — production code + tests. Do NOT deviate from the approved plan; if you
discover the plan is infeasible, STOP and write the blocker to
`.orchestration/runs/W7-MUEGAMMA/impl_blocker.md` instead of improvising.

## Non-negotiables (from the approved plan)
- μ→eγ model LOCKED = LMFV; reuse `flavorConstraints/muToEGamma.py` unchanged.
- Veto driven ONLY by `BR_NP <= br_limit` (MEG II 1.5e-13); C=0.02 is an INERT diagnostic.
- Do NOT edit `flavor_catalog/processes/charged_lepton/L001.tex` or `L001.yaml`.
- Quark-only scans must stay BYTE-IDENTICAL (config hash unchanged; new extra key cannot alter it).
  Add the exact serialized-row regression the plan specifies.
- Add the `passes`-pinned-to-BR test (independent of c_lfv); the absent-extra graceful-degradation
  test; finite/perturbativity guards; the carrier-vs-builder mass-convention split oracle.
- Add `lepton_lmfv_parameters` (or the plan's chosen key) to `QUARK_ONLY_FORBIDDEN_EXTRAS`.
- Carrier path emits `used_proxy=False` and NO truthy proxy/recast diagnostic.

## After implementing
Run the FULL test suite and the new tests:
```
source ~/.bashrc && conda activate ising_bootstrap && export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
python -m pytest -q 2>&1 | tail -30
```
Write a terse implementation summary (files changed, tests added, pass count, the BR oracle number,
confirmation quark-only hash unchanged) to `.orchestration/runs/W7-MUEGAMMA/impl_summary.md` ending
with `IMPL-READY` if green, or `IMPL-BLOCKED` with the failure if not. Do NOT commit — the
orchestrator commits after dual review.
