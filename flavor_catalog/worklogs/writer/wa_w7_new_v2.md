# WA Worklog: wa_w7_new_v2

Date: 2026-05-16
Family label: mixed
Processes: `T003`, `T004`, `B012`
Agent: WA-v2
Cycle: 2

## Scope

Second-pass rework for the CHK-1 failures listed in
`flavor_catalog/worklogs/checker/ca_w7_new_processes.md`.  Edits were limited
to the three process TeX/YAML files and this new worklog.

## T003

- Removed the exact SM central value `4.6e-14` from `T003.tex`.
- Kept only the allowed order-of-magnitude SM statement at the `10^-14` level.
- Appended a cycle-2 `WRITER-DONE` status-history entry and updated
  `last_updated_at`.

## T004

- Removed the flagged exact `139 fb^-1`, `138 fb^-1`, `10^-5`, and `10^-6`
  TeX claims instead of promoting projection/context values.
- Appended a cycle-2 `WRITER-DONE` status-history entry and updated
  `last_updated_at`.

## B012

- Added explicit `year` keys to each `pdg_or_equivalent` observable block.
- Renamed `sha256_of_text_snapshot` metadata fields to `sha256`.
- Removed the flagged exact `3.1 sigma`, `5.2 sigma`, and `2019--2022` TeX
  claims instead of promoting post-2008 context into `pdg_or_equivalent`.
- Appended a cycle-2 `WRITER-DONE` status-history entry and updated
  `last_updated_at`.

## Source Checks

- Recomputed the B012 HFLAV, PDG, Belle, LHCb, and Belle II snapshot hashes
  with `sha256sum`; the values match the sidecar/source-manifest hashes.
- No reference snapshots, manifests, bibliography files, index files, or
  non-batch process files were modified.
