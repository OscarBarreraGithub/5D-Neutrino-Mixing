# WA Worklog: wa_w6_edm_v2

**Date**: 2026-05-16
**Family**: edm_neutrino
**Batch ID**: wa_w6_edm_v2
**Cycle**: 2
**Process IDs**: E009

## Required Reading
- Read `docs/phase_logs/flavor_catalog_plan_v1.md` Sections B and D.
- Read CA worklog `flavor_catalog/worklogs/checker/ca_w6_edm.md`.
- Read current E009 `.tex`, `.yaml`, `flavor_catalog/references/E009/source_manifest.yaml`, and `flavor_catalog/worklogs/pka/E009.md`.

## CA Findings Addressed
- E009 CHK-1: promoted the TeX line 53-60 benchmark numerical claims into value-bearing `pdg_or_equivalent.values` entries with `year`, `value`, `uncertainty`, `units`, `source_url`, `access_date`, and `sha256` metadata.
- Moved the relevant Pospelov-Ritz and Haisch-Hala benchmark content out of `paper_era_reference` / `auxiliary_theory_inputs` so the CA-listed claims no longer live only in auxiliary blocks.

## Source Decisions
- Used existing tracked snapshots under `flavor_catalog/references/E009/`.
- Used `source_manifest.yaml` URLs and access dates for PDG Live 2026, Pospelov-Ritz 2005, and Haisch-Hala 2019.
- Recomputed/checked snapshot hashes with `sha256sum` for the PDG, Pospelov-Ritz, and Haisch-Hala snapshots.

## TeX Changes
- No TeX prose changes were made; the CA finding was handled in sidecar metadata only.

## Bibliography Changes
- None.

## Status Transition
- Appended E009 `WRITER-DONE` with `agent: WA-v2`, `cycle: 2`, and batch ID `wa_w6_edm_v2`.

## Unresolved Items
- No new `\textbf{CHECK}` items introduced.
