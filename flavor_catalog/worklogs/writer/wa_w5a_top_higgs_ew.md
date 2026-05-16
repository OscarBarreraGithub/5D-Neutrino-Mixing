# WA Worklog: wa_w5a_top_higgs_ew

Family: `top_higgs_ew`
Processes: `T005`, `T015`, `T019`
Writer timestamp: 2026-05-16T14:21:56-04:00

## Scope

- Read plan v1 Section B/D, the orchestrator decisions, the three PKA TeX/YAML files, PKA worklogs, process-local source manifests, and local reference snapshots.
- Edited only `flavor_catalog/processes/top_higgs_ew/` and this writer worklog.
- Did not modify `catalog_index.*`, `latex/*`, reference snapshots, PKA worklogs, or other-family files.

## T005 -- \(t \to c g\)

- Added explicit sidecar/source-key tags for the PDG/ATLAS headline limit, CMS comparison limit, SM benchmark, and CFW 2008 RS scale context.
- Tightened the PDG/equivalent and post-2008 wording while preserving the PKA-selected ATLAS Run-2 headline value.
- Appended `WRITER-DONE` to `T005.yaml` and updated `last_updated_at`.

Checker-facing issues:

- The PKA open issue remains: decide whether the older ATLAS 8 TeV PDG row should be mentioned in final prose.
- A live implementation still needs a chiral top chromomagnetic convention and RS-to-SMEFT/collider reinterpretation policy.

## T015 -- \(Z \to e\mu\)

- Replaced loose "latest" wording with sidecar-scoped headline wording for the CMS 2025 experiment-equivalent value.
- Added source-key tags for the CMS 2025, PDG 2025, ATLAS Run-2, ATLAS Run-1, and Tera-\(Z\) numerical claims.
- Appended `WRITER-DONE` to `T015.yaml` and updated `last_updated_at`.

Checker-facing issues:

- The PKA open issue remains: verify whether a future PDG listing has incorporated the CMS 2025 value before Opus approval.
- A live implementation still needs a left/right off-diagonal \(Z e\mu\) coupling convention and a combination policy with low-energy \(\mu\)-flavor observables.

## T019 -- \(h \to e\tau\)

- Added source-key tags for the PDG, ATLAS, CMS, early ATLAS Run-2, and Harnik-Kopp-Zupan numerical claims.
- Kept the sidecar-selected presentation of the PDG ATLAS/CMS pair while making the ATLAS \(2.0\times10^{-3}\) headline trace explicit.
- Appended `WRITER-DONE` to `T019.yaml` and updated `last_updated_at`.

Checker-facing issues:

- The PKA open issue remains: confirm whether final prose should foreground the strongest ATLAS observed limit alone or the PDG ATLAS/CMS pair.
- A live implementation still needs the RS mapping from lepton-sector parameters to off-diagonal physical Higgs Yukawa couplings and a Higgs-width/rate policy.

## Source and Bibliography Notes

- Checked local snapshots for the numerical values cited in the TeX against the sidecar/source manifests.
- No `\textbf{CHECK}` markers were added.
- No `catalog.bib` patch was made because this WA batch was restricted to process files and the batch writer worklog.
