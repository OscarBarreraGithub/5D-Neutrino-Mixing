# Writer Worklog: wa_w23_kaon_charm_edm_v2

Family: mixed (`kaon`, `charm`, `edm_neutrino`)
Processes: K002, C001, E001
Writer: WA-v2
Cycle: 2
Date: 2026-05-16

## Scope

Read plan v1 Sections B and D, the checker worklog
`flavor_catalog/worklogs/checker/ca_w23_kaon_charm_edm.md`, the current
process `.tex`/`.yaml` files, process-local source manifests, reference
snapshots, and PKA worklogs for K002, C001, and E001.

This pass addressed only the CA CHK-1 findings.  No TeX narrative, reference
snapshots, PKA worklogs, shared LaTeX files, catalog indexes, or other process
families were modified.

## Per-Process Changes

### K002 -- `Delta m_K`

- Promoted the PDG no-CPT companion fit, the Bai 2014
  `3.19(41)(96) x 10^-12 MeV` lattice claim, and the Wang 2023 near-9%
  statistical-uncertainty claim from `supporting_values` to keyed
  `pdg_or_equivalent` blocks.
- Added explicit `year`, structured value/uncertainty fields where applicable,
  source URL, access date, snapshot path, and local sha256 metadata to the
  promoted blocks.
- Appended the cycle-2 `WRITER-DONE` status transition.

### C001 -- D-mixing

- Added two CFW contextual `pdg_or_equivalent` numerical blocks for the
  approximately 21 TeV generic RS KK-gluon scale and approximately 33 TeV
  composite pseudo-Goldstone scale.
- Used the existing C001 CFW manifest URL and local snapshot sha256.
- Appended the cycle-2 `WRITER-DONE` status transition.

### E001 -- electron EDM

- Added explicit `year` fields to the PDG 2026 limit, Roussy 2023 primary
  measurement, and ACME/Andreev 2018 benchmark blocks.
- Split the Roussy and ACME numerical measurements/limits out into structured
  `value`, `uncertainty`, component uncertainty, units, and confidence-level
  fields while retaining the existing summary strings.
- Added direct `sha256` fields alongside the existing local-snapshot hash
  fields for the touched numerical blocks.
- Appended the cycle-2 `WRITER-DONE` status transition.

## Remaining Items for CA

No new open issues were introduced in this v2 pass.  The next checker pass
should re-run CHK-1 against the updated YAML placement and structured metadata.
