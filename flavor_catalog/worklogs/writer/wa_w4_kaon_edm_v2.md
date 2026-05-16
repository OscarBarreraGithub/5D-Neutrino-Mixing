# Writer Worklog: wa_w4_kaon_edm_v2

Family: mixed (`kaon`, `edm_neutrino`)
Processes: K017, E006
Writer: WA-v2
Cycle: 2
Date: 2026-05-16

## Scope and Source Decisions

Read plan v1 Sections B and D, the checker worklog
`flavor_catalog/worklogs/checker/ca_w4_kaon_edm.md`, the current TeX/YAML
sidecars, process-local source manifests and snapshots, and PKA worklogs for
K017 and E006.

This pass addressed only the CA-listed CHK-1 findings.  No TeX narrative,
reference snapshots, catalog indexes, shared LaTeX files, or other-family files
were modified.  Snapshot hashes used for the metadata fixes were verified with
`sha256sum`.

No bibliography changes were made.  No new `\textbf{CHECK}` markers were
introduced.

## Per-Process Changes

### K017 -- kaon \(R_K\) LFU

- Added explicit `year` fields to the existing `pdg_or_equivalent` blocks for
  the PDG average, NA62 input, KLOE input, and Cirigliano--Rosell SM
  prediction.
- Appended the cycle-2 `WRITER-DONE` status transition and updated
  `last_updated_at`.

### E006 -- mercury-199 EDM

- Added a value-bearing `pdg_or_equivalent` block for the Graner 2016 claim
  that the \(^{199}\mathrm{Hg}\) EDM bound improved the previous limit by a
  factor of 4, including year, value, uncertainty, units, source URL, access
  date, snapshot path, and SHA-256 metadata.
- Appended the cycle-2 `WRITER-DONE` status transition and updated
  `last_updated_at`.
