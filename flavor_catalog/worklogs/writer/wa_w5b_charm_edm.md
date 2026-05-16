# WA Worklog: wa_w5b_charm_edm

Family: mixed charm/EDM batch
Processes: `C006`, `C008`, `E007`
Writer timestamp: 2026-05-16T14:46:21-04:00

## Scope

- Read plan v1 Section B/D, the orchestrator decisions, the three PKA TeX/YAML files, PKA worklogs, process-local source manifests, and local reference snapshots.
- Edited only the assigned process files in `flavor_catalog/processes/charm/` and `flavor_catalog/processes/edm_neutrino/`, plus this writer worklog.
- Did not modify `catalog_index.*`, `latex/*`, reference snapshots, PKA worklogs, or unassigned process files.

## C006 -- \(D^0\to e^\pm\mu^\mp\)

- Tightened the process, PDG/equivalent, post-2008, and implementation prose.
- Added sidecar source-key tags for the PDG/LHCb current limit, Belle and BaBar predecessor limits, and CFW 2008 baseline.
- Appended `WRITER-DONE` to `C006.yaml` and updated `last_updated_at`.

Checker-facing issues:

- The sidecar has no explicit `PKA-DONE` transition before WA; I did not rewrite PKA history and appended `WRITER-DONE` from the existing terminal state.

## C008 -- \(D^+\to\pi^+e^\pm\mu^\mp\)

- Tightened the PDG/equivalent and post-2008 prose while preserving the two PDG charge-mode limits.
- Added sidecar source-key tags for the two PDG entries, the LHCb 2021 current limits, the BaBar predecessor limits, and the CFW 2008 baseline.
- Appended `WRITER-DONE` to `C008.yaml` and updated `last_updated_at`.

Checker-facing issues:

- The sidecar has no explicit `PKA-DONE` transition before WA; I did not rewrite PKA history and appended `WRITER-DONE` from the existing terminal state.
- The PKA open issue remains: decide whether future recasts should combine the two charge modes or preserve the PDG S031.110/S031.111 split.

## E007 -- Radium-225 and xenon-129 EDMs

- Tightened the PDG-equivalent, post-2008, and code-coverage prose without adding new physics claims beyond the PKA deposit.
- Added sidecar snapshot/source-key tags for the Ra and Xe direct limits, the Ra and Xe future-context numbers, and the CFW 2008 baseline.
- Appended `WRITER-DONE` to `E007.yaml` and updated `last_updated_at`.

Checker-facing issues:

- The PKA open issue remains: verify whether a current PDG Review text table includes standalone `225Ra` or `129Xe` atom-EDM rows.
- If promoted to live code, the PI must specify the CP-odd low-energy operator basis and acceptable nuclear/atomic Schiff-moment inputs for Ra and Xe.
- Coordinate with E004, E006, E008, and E009 to avoid double-counting neutron, Hg, chromo-EDM, and Weinberg-operator translations as independent measurements.

## Source and Bibliography Notes

- Checked the numerical values cited in the TeX against the sidecar/source-manifest records and local snapshots.
- Reused existing PKA sources only; no new source snapshots or reference files were added.
- No `\textbf{CHECK}` markers were introduced.
- No `catalog.bib` patch was made because this WA batch was restricted to process files and the batch writer worklog.
