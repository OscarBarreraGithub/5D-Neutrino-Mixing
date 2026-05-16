# WA Worklog: wa_wave1_top_higgs_ew_v2

**Date**: 2026-05-16
**Family**: top_higgs_ew
**Cycle**: 2
**Processes**: T001, T010
**Checker log**: `flavor_catalog/worklogs/checker/ca_wave1_top_higgs_ew.md`

## Scope

Addressed only the CA CHK-1 metadata-placement findings for T001 and T010.
No TeX narrative was rewritten.

## Changes

- T001: promoted all CA-flagged contextual numerical entries from
  `auxiliary_values` into `pdg_or_equivalent.values`; added explicit `year`
  fields while preserving value, uncertainty, units, source URL, access date,
  normalized value where present, and local snapshot sha256.
- T010: promoted the LEP/SLC `2.8 sigma` pull and the FCC-ee order-`0.01%`
  projected relative uncertainty from `additional_numerical_context` into
  `pdg_or_equivalent`; preserved the existing metadata fields.
- T001 and T010: appended `WRITER-DONE` status-history transitions for
  `WA-v2`, cycle 2, and updated `last_updated_at`.

## Source checks

Local snapshot hashes were recomputed with `sha256sum` for both process
reference directories. The hashes matched the sidecar/source-manifest values
used in the promoted entries.

## Open issues

No new `\textbf{CHECK}` items or unresolved bibliography changes were introduced.
