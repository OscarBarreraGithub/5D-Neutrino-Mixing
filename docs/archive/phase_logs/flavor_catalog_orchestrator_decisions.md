# Flavor Catalog — Orchestrator Decisions

**Date**: 2026-05-16
**Orchestrator**: Claude Code (acting under the PI's standing directive to
"make decisions on my behalf, but write them down").
**Plan under execution**: `docs/phase_logs/flavor_catalog_plan_v1.md` (commit
`744e70c`), approved by Opus at `docs/phase_logs/flavor_catalog_plan_v1_opus_signoff.md`
(commit `95f9f63`) with verdict **APPROVE-FOR-EXECUTION**.

This file records the orchestrator's resolution of Plan v1 Section I open
questions (which the planner left for the PI) and the disposition of the 5
NITs raised by Opus in the final approval. These decisions stand until the PI
overrides them in a future message.

---

## Section I — Open questions, with orchestrator decisions

### I.1 — `K1 → πγ` disambiguation
**Decision**: interpret the PI's seed `K1 → πγ` as **`K_L → π^0 γγ`** (process
id `K013`). Carry `K_S → π^0 γγ` (`K014`) as a secondary entry. Both belong to
the same PDG rare-radiative-kaon family and have published branching-fraction
measurements; this matches the "(off-diagonal ε')" annotation the PI placed
next to the seed. **Flag** for PI confirmation on the next PI touchpoint;
proceed in the meantime.

### I.2 — Scope
**Decision**: full default scope — quark + charged-lepton (LFV) + neutrino +
electroweak-pole + Higgs/top FCNC + EDM-adjacent. Pruning is cheaper after
PKA than before; the PI explicitly said "Don't pre-prune; I'll prune later."

### I.3 — Branch policy
**Decision**: execute on a new branch `flavor-catalog/2026q2`, cut from
`paper/quark-scan-2026q2` at the rc1.1 tag (`6300f82a` / peeled `9a13d164`).
**Do not** push catalog work onto the paper branch. The paper branch is
frozen at rc1.1.

### I.4 — PDG / publisher-PDF licensing
**Decision**: minimal text snapshots + URL + access date + sha256 only.
**Do not** track publisher PDFs even when access is technically permitted.
ArXiv preprints may be tracked as PDF or text snapshot at the agent's
discretion. License decisions are easy to upgrade later and expensive to
retract.

### I.5 — rc1.1 / rc2 relationship
**Decision**: companion artifact, **not** part of the rc1.1 / rc2 paper.
The quark-scan paper is sealed at rc1.1. The catalog informs the
next-paper roadmap (extension into the broader flavor sector).

### I.6 — Integration target if a process graduates from catalog to live constraint
**Decision**: target `quarkConstraints/modern/`. Treat
`quarkConstraints/deltaf2.py` as legacy-stable; new constraints should not
extend it. The graduation step is out of scope for the catalog effort
itself.

---

## Opus NITs — Disposition

### NIT 1 — Section C: oblique parameters and CKM unitarity not pre-seeded
**Disposition**: ACCEPT. DA-1 (the first Discovery Agent round) is
instructed to consider — at minimum — the following as proposed additions
in its first sweep:
- `EW001` S/T/U oblique electroweak parameters
- `EW002` first-row CKM unitarity (|V_ud| + |V_us| + |V_ub|, Cabibbo
  anomaly)
- `EW003` |V_cb| inclusive vs exclusive tension and |V_ub| equivalent
DA-1 may add or merge with EW process entries as the PKA evidence suggests.

### NIT 2 — Section C charm: D-mixing observables packed into one row
**Disposition**: ACCEPT. The PKA prompt for `C001` (charm D-mixing) will
require the YAML sidecar to emit a **separate** `pdg_or_equivalent` block
for each observable (`x_D`, `y_D`, `Δm_D`, `δy`). No plan edit; this is a
PKA-prompt clarification at dispatch time.

### NIT 3 — Section D Opus parallelism
**Disposition**: ACCEPT. The orchestrator will default to **one Opus
sign-off at a time** for the first two Opus rounds (to lock global
consistency), and only allow concurrent Opus sign-offs (max 2) once the
catalog index has passed a freeze-check.

### NIT 4 — Section E DA round 1 gating
**Disposition**: ACCEPT. The DA-1 spawn prompt will clarify that "50% CA
pass" is computed across the whole 128-entry initial process list, not
within a single family. Cosmetic clarification for the spawn prompt; no
plan edit.

### NIT 5 — Section H agent-hour total
**Disposition**: ACKNOWLEDGED. Budget expectation is roughly 466–768
agent-hours for Codex/Opus combined + 30–50 orchestrator-hours. Recorded
here so the PI sees it before the first PKA fires. The orchestrator will
flag if actual usage diverges materially from this band.

---

## Execution sequence (the orchestrator's actual plan)

The plan v1 prescribes the following execution sequence, which the
orchestrator commits to:

1. Create `flavor-catalog/2026q2` branch from `paper/quark-scan-2026q2`
   at rc1.1; push to origin.
2. Scaffold the `flavor_catalog/` subdirectory structure per Section A of
   the plan (directories, `README.md`, `catalog_index.yaml`,
   `catalog_index.tex`, `catalog_master.tex`, `latex/macros.tex`,
   `latex/process_template.tex`, empty `processes/<family>/`,
   `worklogs/`, `signoff/`, `references/`).
3. Dispatch first wave of PKAs: ~8–12 in parallel, one process each,
   batched by family. Each writes the per-process `.tex` + `.yaml`
   sidecar + `references/<process_id>/source_manifest.yaml`.
4. After ~50% PKA completion, dispatch WA batches (each WA polishes 3–5
   processes from the same family).
5. CA batches (separate from WA) verify each WA's batch.
6. DA-1: scan the catalog for missing processes, propose additions
   including the NIT-1 candidates.
7. Iterate WA/CA on additions; DA-2 if needed.
8. Opus per-process sign-off pass.
9. Master compile of `catalog_master.tex`.
10. Final orchestrator merge / branch hand-off to the PI.

Iteration caps from plan v1 Section G are binding: PKA micro 2, WA/CA 3,
DA 4, Opus corrective 1. Cap-hit transitions force `BLOCKED-PI` or
`DEFERRED-SCOPE` status with mandatory signoff log entries.

---

## What the orchestrator will NOT do

- Will not write per-process catalog content itself. All physics reading,
  PDG snapshot extraction, arXiv downloads, and TeX drafting goes to
  Codex/Opus.
- Will not skip the writer/checker separation rule. No agent will be
  both WA and CA on the same batch in the same round.
- Will not push catalog work onto the paper branch.
- Will not silently expand scope. If a process clearly falls outside the
  PI's intent, it goes to `DEFERRED-SCOPE` with a one-line note for the
  PI's eventual review.

---

This document is binding for the catalog workstream until the PI
explicitly overrides one of the above decisions.
