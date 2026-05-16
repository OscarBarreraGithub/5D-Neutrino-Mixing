# WA Worklog: wa_w6_edm

Family: `edm_neutrino`
Processes: `E009`
Writer timestamp: 2026-05-16T15:53:32-04:00

## Scope

- Read plan v1 Section B/D, the orchestrator decisions, the E009 PKA TeX/YAML, the PKA worklog, the process-local source manifest, and the local reference snapshots.
- Edited only `flavor_catalog/processes/edm_neutrino/E009.tex`, `flavor_catalog/processes/edm_neutrino/E009.yaml`, and this writer worklog.
- Did not modify `catalog_index.*`, `latex/*`, reference snapshots, PKA worklogs, or other families' files.

## E009 -- Weinberg three-gluon operator

- Tightened the PDG/equivalent and post-2008 prose while preserving the PKA's distinction between the measured neutron EDM anchor and convention-dependent Weinberg-coefficient translations.
- Added explicit source-key anchors for the PDG 2026 neutron limit, Abel 2020 measurement, Pospelov--Ritz 22 MeV conversion, Haisch--Hala 74 MeV convention, and lattice/composite-context claims already present in the sidecar.
- Appended `WRITER-DONE` to `E009.yaml` and updated `last_updated_at`.

Checker-facing issues:

- The two displayed coefficient bounds remain central, one-source-at-a-time benchmark translations, not correlated global EDM-fit intervals.
- CA should verify whether the final catalog should keep both Pospelov--Ritz and Haisch--Hala conventions or privilege one convention in global EDM prose.
- Coordinate with E004/E008 so neutron-EDM, qCEDM, and Weinberg-operator translations are not double-counted as independent experimental measurements.

## Source and Bibliography Notes

- Checked the numerical values cited in the TeX against the sidecar/source-manifest records and local snapshots.
- Reused existing PKA sources only; no new source snapshots or reference files were added.
- No `\textbf{CHECK}` markers were introduced.
- No `catalog.bib` patch was made because this WA batch was restricted to process files and the batch writer worklog.
