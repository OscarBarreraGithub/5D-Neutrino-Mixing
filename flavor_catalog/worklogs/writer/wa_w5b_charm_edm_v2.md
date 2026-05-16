# WA Worklog: wa_w5b_charm_edm_v2

Family: charm_edm / `edm_neutrino`
Processes: `E007`
Writer timestamp: 2026-05-16T15:27:52-04:00
Cycle: 2

## Scope

- Read plan v1 Sections B and D, `flavor_catalog/worklogs/checker/ca_w5b_charm_edm.md`, `E007.tex`, `E007.yaml`, the E007 source manifest and reference snapshots, and `flavor_catalog/worklogs/pka/E007.md`.
- Addressed only the CA CHK-1 finding for E007.
- Did not modify `latex/*`, `catalog_index.*`, reference snapshots, PKA worklogs, checker worklogs, or any other family's process files.

## E007 -- Radium-225 and xenon-129 EDMs

- Resolved the CHK-1 finding by softening the future/context numerical claims in the Post-2008 developments paragraph rather than promoting proposal-context blocks into `pdg_or_equivalent`.
- Left the measured direct Ra and Xe numerical limits unchanged because CA already found those traced through `pdg_or_equivalent`.
- Appended `WRITER-DONE` with `cycle: 2` to `E007.yaml` and updated `last_updated_at`.

Checker-facing notes:

- No new source snapshots or reference entries were added.
- No `\textbf{CHECK}` markers were introduced.
- The remaining numerical TeX claims are the direct measured limits already covered by `pdg_or_equivalent`, plus the existing Ra improvement factor and Xe direct-limit scale already supported by the accepted direct-measurement blocks.
