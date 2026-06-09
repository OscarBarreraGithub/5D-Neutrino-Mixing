# W7 — μ→eγ LMFV implementation PLAN (codex author)

You are authoring an IMPLEMENTATION PLAN only. Do NOT write production code in this run.
Do NOT make physics/model decisions — the model is LOCKED. Output a written plan to
`.orchestration/runs/W7-MUEGAMMA/plan_codex.md` and end it with the literal line `PLAN-READY`.

## Locked decision (do not revisit)
μ→eγ model = **LMFV** (Perez-Randall arXiv:0805.4652): spurion `(Y_N Y_N†)_12`, NDA `C=0.02`.
The dipole/BR core is ALREADY coded in `flavorConstraints/muToEGamma.py`. MEG II 2025 limit
`br_limit = 1.5e-13`. Do not switch to anarchic/Chen-Yu/gauged — LMFV is final.

## Goal
Make μ→eγ a LIVE, evaluable constraint in the catalog scan by building the missing
**lepton-sector parameter carrier** and wiring the existing dipole into the constraint
that represents μ→eγ in `flavor_catalog_constraints/` (locate it: likely `primary/.../L001*`).

## What to investigate and specify (ground every claim in real files)
1. **Read** `flavorConstraints/muToEGamma.py`, the lepton Yukawa/seesaw code under `yukawa/`
   (`neutrino.py`, `compute_yukawas.py`) and `neutrinos/`, the existing quark `point_builder`
   / `ParameterPoint` and `get_extra(...)` pattern (see how `rs_ew_couplings` is attached and
   read by EW001), and the μ→eγ catalog constraint file(s).
2. Specify the **lepton-sector ParameterPoint extra**: what fields it must carry (Y_N, PMNS,
   M_KK / M_N, c_L, c_E, c_N, v, k, ε), how they are produced from a scan draw, and how this
   composes with the existing quark-only flow without breaking it.
3. Specify exactly how the catalog μ→eγ constraint reads that extra and calls the LMFV dipole
   with C=0.02 and br_limit=1.5e-13, including the rigorous|proxy|partial tag and graceful
   degradation when the lepton extra is absent (must NOT break quark-only scans).
4. **Tests:** numeric oracle for one benchmark BR; a regression that quark-only mode is byte-
   identical; an absent-extra graceful-degradation test; a perturbativity/finite guard test.
5. **Isolation:** list every file you will touch; flag any shared with W8 (custodial PR2:
   `quarkConstraints/rs_ew_couplings.py`, `oblique_stu.py`) or W9 (scan harness) so the
   orchestrator can serialize. Keep `minimal_rs` / quark-only config-hash unchanged.
6. State what stays NEEDS-HUMAN / deferred (loop-level, off-diagonal LFV) and why that is honest.

## Constraints
- Follow the dual-signoff gate: this plan will be independently reviewed by a second codex AND
  an Opus agent; write it to be checkable (cite file:line, give formulas, list tests).
- Determinism, graceful degradation, no fabricated physics. End with `PLAN-READY`.
