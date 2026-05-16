# WA Worklog: wa_w7_new_processes

Date: 2026-05-16
Family label: mixed
Processes: `T003`, `T004`, `T008`, `T012`, `B012`
Agent: WA

## Scope

Polished the Wave-7 follow-up PKA drafts for the mixed top/Higgs/electroweak
and beauty batch.  Edits were limited to the five process TeX/YAML sidecars and
this writer worklog.  No PKA worklogs, reference snapshots, index files,
macros, templates, or bibliography files were modified.

## T003

- Added explicit standard notation and normalized the Section-B headings.
- Tightened the CMS-vs-ATLAS prose and added source keys for the CMS dataset,
  ATLAS comparison, PDG generic row, and SM benchmark numerical claims.
- Appended a `WRITER-DONE` status-history entry and updated `last_updated_at`.

## T004

- Added explicit standard notation and normalized the Section-B headings.
- Kept the PDG combined row separate from the flavor-resolved ATLAS/CMS limits,
  with all numerical limits tied to sidecar value ids.
- Appended a `WRITER-DONE` status-history entry and updated `last_updated_at`.

## T008

- Normalized the Section-B headings and tightened the PDG/CMS/ATLAS limit prose.
- Added source keys to the ATLAS and CMS companion numerical limits and the
  95% CL model-dependence statement.
- Appended a `WRITER-DONE` status-history entry and updated `last_updated_at`.

## T012

- Added explicit standard notation and normalized the Section-B headings.
- Added source keys for the PDG charm Z-pole values, the LEP/SLC cross-check,
  and the post-2008 RS/SM/FCC-ee context.
- Appended a `WRITER-DONE` status-history entry and updated `last_updated_at`.

## B012

- Normalized the Section-B headings and kept the HFLAV/PDG numerical averages
  clearly separated.
- Added source keys for the branching fractions, isospin/direct-CP averages,
  time-dependent CP averages, and post-2008 BaBar/Belle/LHCb/Belle-II claims.
- Appended a `WRITER-DONE` status-history entry and updated `last_updated_at`.

## Bibliography And Source Notes

- No bibliography patch was made because this dispatch forbids edits outside
  `flavor_catalog/processes/<family>/` and `flavor_catalog/worklogs/writer/`.
- All accepted PKA source snapshots remain represented in each sidecar's
  `source_shas` block or equivalent source metadata.
- No `\textbf{CHECK}` markers were introduced.

## Open Issues For CA

- `T003`, `T004`, `T008`, and `T012` inherited sidecars whose initial PKA state
  is `WRITER-INITIATED` rather than `PKA-DONE`; WA appended `WRITER-DONE` from
  the existing state without rewriting PKA-owned history.
- `T003`: CA should confirm that the CMS charm-specific limit is the desired
  headline rather than an ATLAS-only or PDG-generic convention.
- `T004`: CA should confirm whether the ATLAS left-handed benchmark can remain
  the headline while the right-handed benchmark is kept in prose.
- `T008`: CA should confirm that the PDG headline CMS diphoton row remains the
  prose headline despite the numerically identical CMS 2025 observed
  combination.
- `T012`: CA/orchestrator should decide whether the combined T012 entry also
  closes the plan's separate T013 charm-asymmetry row.
