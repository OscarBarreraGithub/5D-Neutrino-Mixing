# Writer Worklog: wa_w6_charged_lepton_top

Batch: `wa_w6_charged_lepton_top`
Family label: `charged_lepton_top` mixed batch
Processes: `L023`, `T020`
Writer timestamp: 2026-05-16T15:53:39-04:00

## L023 -- neutrino trident

- Tightened the process text and normalized the PDG-or-equivalent heading while preserving the PKA's historical trident-ratio inputs.
- Added explicit process-local source keys for the CHARM-II/CCFR/NuTeV ratios, CCFR event counts, Altmannshofer 2014 reinterpretation, Belle-II auxiliary reach, and DUNE projection numbers.
- Appended a `WRITER-DONE` status-history entry and updated `last_updated_at`.

Open issues for CA:

- CHARM-II uncertainty convention remains source-dependent: `Altmannshofer2014TridentProbe` quotes `1.58 +/- 0.57`, while the DUNE compilation used for the catalog value quotes `1.58 +/- 0.64`.
- DUNE projection wording remains split between about `25%` accuracy in the abstract and a `40%` baseline plus `25%` improved contour in the body; both are documented in the sidecar.

## T020 -- h to e mu

- Tightened the collider-limit prose and kept the PDG-rounded Run-2 CMS and ATLAS limits as the headline values.
- Kept the ATLAS abstract's `6.1e-5` observed / `5.8e-5` expected value as a primary-source cross-check, while making the PDG `6.2e-5` / `5.9e-5` rounding the catalog value.
- Added source-key/value-ID anchors for the `125 GeV`, `13 TeV`, luminosity, branching-ratio, mass-scan, and significance numbers; appended `WRITER-DONE` and updated `last_updated_at`.

Open issues for CA:

- The T020 sidecar entered this WA pass with `WRITER-INITIATED` as the only prior status-history state, not an explicit `PKA-DONE`; I appended `WRITER-DONE` from the latest existing state rather than mutating PKA history.
- A future live implementation still needs the PKA-noted mapping from RS lepton-sector parameters to off-diagonal physical Higgs Yukawa couplings and a decision on non-SM Higgs production or total-width assumptions.

## Bibliography / source decisions

- No `catalog.bib` patch was made because this batch's hard rules restrict edits to the two process entries and the writer worklog.
- Process-local source keys from `flavor_catalog/references/L023/source_manifest.yaml` and `flavor_catalog/references/T020/source_manifest.yaml` were retained in the TeX pending bibliography consolidation.
