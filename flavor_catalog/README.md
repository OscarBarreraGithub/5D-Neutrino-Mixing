# Flavor Constraint Catalog

A discovery-mode catalog of low-energy flavor and direct-collider RS
constraints, maintained as a companion artifact to the quark-scan paper
at rc1.1. Not a live constraint backend; should stay separate from the
existing scan / solver / neutrino / Yukawa / paper-note code paths
until the PI approves a later integration step.

This work lives on branch `flavor-catalog/2026q2`. Do not merge to
`paper/quark-scan-2026q2` (frozen at rc1.1). Merge into `main` only on
PI sign-off.

## Layout

```text
flavor_catalog/
  README.md                       # this file
  AGENTIC_WORKFLOW.md             # reproducible playbook (multi-agent pipeline)
  PRIORITY_TIERS.md               # PRIMARY vs SECONDARY tier policy
  SESSION_NOTES.md                # final state and per-wave history
  CATALOG_METHODOLOGY.tex/.pdf    # 1-page collaborator-facing pitch
  HANDOFF_PROMPT.md               # cold-boot prompt for a new orchestrator
  catalog_master.tex              # master compile (179 pp at v0.4)
  catalog_master.pdf              # current compiled output
  catalog_index.tex / .yaml       # auto-generated TOC + status ledger
  master_compile_v0N_report.md    # per-tag compile reports

  latex/
    macros.tex
    process_template.tex

  processes/                      # PRIMARY family entries
    kaon/                         # 13 entries
    charm/                        # 8 entries
    beauty/                       # 22 entries
    top_higgs_ew/                 # 19 entries
    charged_lepton/               # 11 entries (LFV)
    edm_neutrino/                 # 7 entries
    collider_rs/                  # 14 entries (Wave-9: KK gauge resonances,
                                  #             custodial top/bottom partners,
                                  #             VLQ, KK-graviton, EW-tail)
    secondary/                    # SECONDARY (Wave-8+ deferred re-promotions)
      kaon/                       # 3 entries (K019/K020/K021 LFV trio)
      beauty/                     # 4 entries (B007/B008/B013/B014)
      top_higgs_ew/               # 1 entry  (T014)

  references/                     # per-process snapshots, manifests, sha256
    <id>/                         # FLAT; IDs are globally unique

  worklogs/
    pka/                          # per-process PKA worklogs
    writer/                       # per-batch WA worklogs
    checker/                      # per-batch CA worklogs
    discovery/                    # DA-1..DA-4 + addendum
    orchestration/                # wave_NNN_runbook.md (live state per wave)

  signoff/
    by_process/                   # per-arbitration docs (only when needed)
      L001.md, B021_B023.md, B001_B003.md
    round_NNN_index.md            # per-round Opus verdicts (round_001..round_005)

  audits/
    factcheck_<family>.md         # per-family fact-check tables (WebFetch)
    factcheck_status.md           # master status matrix
```

Process drafts belong under `processes/<family>/` (PRIMARY) or
`processes/secondary/<family>/` (SECONDARY). Raw source snapshots and
manifests belong under `references/<id>/`; use minimal text snapshots
with URL, access date, and sha256 for PDG, HFLAV, FLAG, CKMfitter,
UTfit, and experiment pages. Do not track publisher PDFs. ArXiv PDFs or
text extracts may be deposited only when useful and compatible with the
catalog size policy.

## Build

`catalog_master.pdf` at the current tag `flavor-catalog-v0.4` is built
from **102 OPUS-APPROVED entries** across seven families.

- **9 waves** of additions:
  - Waves 1–3: 26 PI-seed + early-discovery entries
  - Waves 4–6: 50 DA-1/2/3 additions
  - Wave-7: 5 entries from external GPT Deep Research reconciliation
  - Wave-8: 8 SECONDARY DA-deferred re-promotions
  - Wave-9: 14 PRIMARY new-family `collider_rs` entries
- **5 Opus sign-off rounds**: round 1 (50 APPROVE), round 2 (21 APPROVE),
  round 3 (5 + subtleties + addendum APPROVE), round 4 (8 APPROVE),
  round 5 (14 APPROVE). Cumulative: 102/102 APPROVE, 0 RETURN-TO-WA,
  0 ESCALATE-TO-PI.
- **5 individual arbitrations** completed: L001, B021, B023, B001, B003
  (all APPROVE-OVERRIDE under the CHK-1 carve-out precedent).
- **Fact-check status**: **100 VERIFIED / 2 PARTIAL** (E009 INSPIRE
  JS-only, APS cross-check confirms; K020 NA62 author convention
  cleared post-v0.3) / **0 MISMATCH / 0 FAILED** across all 102
  entries.
- Current master build: **179 pages** (913 KB PDF).
- Tags: `flavor-catalog-v0` (75) → `v0.1` (fact-checked) → `v0.2` (80,
  +Wave-7) → `v0.3` (88, +Wave-8 SECONDARY) → `v0.4` (102, +Wave-9
  collider_rs PRIMARY).
- Last update: 2026-05-17.

## Status

Per-process metadata lives in YAML sidecars; the TeX files contain
prose only. The current status is the `to` value of the final
`status_history` entry; there is no separate scalar `status` field.

The normal approval path for writer/checker work is:

```text
DRAFT -> WRITER-INITIATED -> WRITER-DONE -> CHECKER-DONE -> OPUS-APPROVED -> FACT-CHECKED
```

Any non-terminal state may transition to `BLOCKED-PI` or
`DEFERRED-SCOPE` on an iteration-cap hit. Every `status_history` entry
should include `from`, `to`, `at`, `agent_id`, and `reason`.

## Tier policy

`flavor_catalog/PRIORITY_TIERS.md` is the binding policy doc. In short:

- **PRIMARY** (the default, no `priority_tier` field): the Waves 1–7
  flavor + Wave-9 collider_rs entries (94 total). These are the
  canonical RS-flavor + RS-collider constraints; downstream code
  implementation should start here.
- **SECONDARY** (`priority_tier: SECONDARY` + `promoted_in_wave`):
  Wave-8+ DA-deferred re-promotions (8 entries so far). Lower
  implementation priority. Located at
  `processes/secondary/<family>/<id>.{tex,yaml}`.

Two extension patterns for future waves are codified in
`PRIORITY_TIERS.md §7` ("How to Extend"):

- **Pattern A**: SECONDARY re-promotion of deferred entries (Wave-8
  pattern).
- **Pattern B**: PRIMARY new-family directory for a new scope class
  opened by PI directive (Wave-9 pattern).

## Governance + reproducibility inputs

- [Reproducible playbook](AGENTIC_WORKFLOW.md) — multi-agent pipeline
  mechanics (roles, 8-item CA checklist, CHK-1 carve-out policy,
  arbitration, fact-check round, concrete examples through Wave-9)
- [Tier policy](PRIORITY_TIERS.md)
- [Session notes](SESSION_NOTES.md) — final state, per-wave history,
  cluster quirks, candidate next steps
- [Cold-boot orchestrator handoff](HANDOFF_PROMPT.md)
- [Collaborator-facing pitch (1 page)](CATALOG_METHODOLOGY.tex)
- [Plan v1 (original design)](../docs/phase_logs/flavor_catalog_plan_v1.md)
- [Plan-v1 orchestrator decisions](../docs/phase_logs/flavor_catalog_orchestrator_decisions.md)
- [Plan-v1 Opus approval](../docs/phase_logs/flavor_catalog_plan_v1_opus_signoff.md)
- Per-wave runbooks: `worklogs/orchestration/wave_NNN_runbook.md`
  (Wave-8, Wave-9; future waves follow the same template)
