# WA Worklog: wa_w23_top_higgs_ew

Family: `top_higgs_ew`
Processes: `T002`, `T007`, `T018`
Writer timestamp: 2026-05-16T12:15:15-04:00

## Source and Scope Notes

- Read plan v1 Section B/D, the orchestrator decisions, each PKA TeX/YAML
  sidecar, each PKA worklog, and the process-local source manifests and
  snapshots.
- Kept writes scoped to `flavor_catalog/processes/top_higgs_ew/` and this
  writer worklog. No reference snapshots, PKA worklogs, catalog indexes,
  macros, or other families were modified.
- Did not create a `catalog.bib` patch because this batch prompt's hard rules
  prohibit writes outside the process files and writer worklog. Process-local
  source keys remain listed in each entry.

## Process Updates

### T002 -- \(t \to u Z\)

- Added the explicit standard-notation line and normalized the
  `PDG-or-equivalent value` heading to the plan-v1 template.
- Tightened the PKA prose and added sidecar-resolving IDs for the PDG, ATLAS,
  CMS, SM-estimate, and RS-suppression numerical claims.
- Appended `WRITER-DONE` to `T002.yaml` and updated `last_updated_at`.

Open for CA: decide whether both ATLAS left- and right-handed tensor-SMEFT
benchmarks should remain equally prominent, and retain the implementation note
that a live scan needs an RS-to-SMEFT or direct `Ztu` convention.

### T007 -- \(t \to H c\)

- Added the explicit standard-notation line and normalized the
  `PDG-or-equivalent value` heading.
- Smoothed the PDG/ATLAS/CMS discussion while keeping the PDG ATLAS-combination
  headline and the CMS comparison value tied to sidecar IDs.
- Appended `WRITER-DONE` to `T007.yaml` and updated `last_updated_at`.

Open for CA: confirm that the PDG/ATLAS `3.4e-4` limit should remain the
headline over the newer but weaker CMS `3.7e-4` combination.

### T018 -- \(h \to \mu\tau\)

- Added the explicit standard-notation line and normalized the
  `PDG-or-equivalent value` heading.
- Tightened the Run-2 and historical CMS-hint paragraphs, with sidecar IDs on
  the PDG, CMS, ATLAS, Harnik-Kopp-Zupan, and CFW numerical/source claims.
- Appended `WRITER-DONE` to `T018.yaml` and updated `last_updated_at`.

Open for CA: decide whether the final display should headline the strongest
CMS `0.15%` limit alone or keep the PDG ATLAS/CMS pair as the canonical value.

## Batch-Level CA Notes

- No `\textbf{CHECK}` markers were introduced in the TeX.
- The PKA sidecars arrived with `WRITER-INITIATED` as the latest status rather
  than `PKA-DONE`; the writer transition was appended from that latest state
  without rewriting PKA history.
- The PKA-owned `pdg_or_equivalent`, `code_coverage`, and
  `implementation_difficulty` blocks were not modified.
