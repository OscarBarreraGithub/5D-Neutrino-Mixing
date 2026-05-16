# WA Worklog: wa_w5a_charged_lepton_edm

Family: `charged_lepton_edm` mixed batch
Processes: `L005`, `L006`, `E002`
Writer timestamp: 2026-05-16T14:23:54-04:00

## Scope

- Read plan v1 Section B/D, the orchestrator decisions, the three PKA TeX/YAML files, PKA worklogs, process-local source manifests, and local reference snapshots.
- Edited only the assigned process files in `flavor_catalog/processes/charged_lepton/` and `flavor_catalog/processes/edm_neutrino/`, plus this writer worklog.
- Did not modify `catalog_index.*`, `latex/*`, reference snapshots, PKA worklogs, or unassigned process files.

## L005 -- \(\mu-e\) conversion in titanium

- Tightened the process and post-2008 wording to keep the Ti SINDRUM II bound distinct from Al/Au conversion programs.
- Added explicit sidecar-block/source-key tags for the PDG Ti limit, Mu2e projection context, COMET Phase-I numbers, and CFW 2008 baseline.
- Appended `WRITER-DONE` to `L005.yaml` and updated `last_updated_at`.

Checker-facing issues:

- The PKA open issue remains: confirm whether the PDG live `NODE=S004` stream should remain the canonical source rather than a static annual listing.
- A live implementation still needs an operator basis and target-specific nuclear-overlap convention before comparing Ti, Au, and Al limits.

## L006 -- Muonium-antimuonium conversion

- Added sidecar-block tags for the PDG \(G_C/G_F<0.0030\) limit and the MACS/PSI \(P_{M\bar M}<8.3\times10^{-11}\) probability limit.
- Tightened the post-2008 section to distinguish unchanged current bounds from MACE/Snowmass prospects and EFT context.
- Appended `WRITER-DONE` to `L006.yaml` and updated `last_updated_at`.

Checker-facing issues:

- The sidecar had no explicit `PKA-DONE` transition before WA; I did not rewrite PKA history and appended `WRITER-DONE` from the existing terminal state.
- The PKA open issue remains: confirm the PDG/final-publication \(8.3\times10^{-11}\) probability value versus the arXiv abstract's \(8.2\times10^{-11}\) wording.
- A live implementation still needs Lorentz-operator, magnetic-field, and spin-state conventions.

## E002 -- Muon electric dipole moment

- Added sidecar-block/source-key tags for the PDG direct limit, PDG combined measurement, Bennett 2009 direct result, PSI frozen-spin projection, and 2024 PSI status numbers.
- Tightened the post-2008 wording while preserving the PKA distinction between current direct bounds and prospective PSI reach.
- Appended `WRITER-DONE` to `E002.yaml` and updated `last_updated_at`.

Checker-facing issues:

- The PKA open issue remains: decide how much PSI projection detail belongs in final prose rather than auxiliary metadata.
- If promoted to live code, the PI must specify lepton-sector matching assumptions and whether muon EDM is a hard constraint or roadmap diagnostic.

## Source and Bibliography Notes

- Checked local snapshots for the numerical values cited in the TeX against the sidecar/source manifests.
- No `\textbf{CHECK}` markers were added.
- No `catalog.bib` patch was made because this WA batch was restricted to process files and the batch writer worklog.
