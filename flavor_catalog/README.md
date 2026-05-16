# Flavor Constraint Catalog

This subdirectory is the discovery-mode flavor-constraint catalog, maintained as a companion artifact to the quark-scan paper at rc1.1. It is not a live constraint backend and should stay separate from the existing scan, solver, neutrino, Yukawa, and paper-note code paths until the PI approves a later integration step.

This work lives on branch `flavor-catalog/2026q2`. Do not merge it to `paper/quark-scan-2026q2` until PI sign-off is recorded.

## Layout

```text
flavor_catalog/
  README.md
  catalog_index.tex
  catalog_index.yaml
  catalog_master.tex
  latex/
    macros.tex
    process_template.tex
  processes/
    kaon/
    charm/
    beauty/
    top_higgs_ew/
    charged_lepton/
    edm_neutrino/
  worklogs/
    writer/
    checker/
    discovery/
    pka/
  signoff/
    by_process/
    round_index/
  references/
```

Process drafts belong under `processes/<family>/`. Raw source snapshots and manifests belong under `references/`; use minimal text snapshots with URL, access date, and sha256 for PDG, HFLAV, FLAG, CKMfitter, UTfit, and experiment pages. Do not track publisher PDFs. ArXiv PDFs or text extracts may be deposited only when useful and compatible with the catalog size policy.

## Status

Per-process metadata lives in YAML sidecars; the TeX files contain prose only. The current status is the `to` value of the final `status_history` entry, and there is no separate scalar `status` field.

The normal approval path for writer/checker work is:

```text
DRAFT -> WRITER-INITIATED -> WRITER-DONE -> CHECKER-DONE -> OPUS-APPROVED
```

Any non-terminal state may transition to `BLOCKED-PI` or `DEFERRED-SCOPE` on an iteration-cap hit. Every `status_history` entry should include `from`, `to`, `at`, `agent_id`, and `reason`.

## Governance Inputs

- [Flavor Catalog Planner v1](../docs/phase_logs/flavor_catalog_plan_v1.md)
- [Flavor Catalog Orchestrator Decisions](../docs/phase_logs/flavor_catalog_orchestrator_decisions.md)
- [Opus Final Approval](../docs/phase_logs/flavor_catalog_plan_v1_opus_signoff.md)
