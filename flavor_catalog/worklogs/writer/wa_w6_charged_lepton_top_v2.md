# WA-v2 Worklog: wa_w6_charged_lepton_top_v2

Date: 2026-05-16
Family: charged_lepton_top
Processes: `L023`, `T020`
Cycle: 2
Agent: WA-v2

## Scope

Addressed only the CA-listed CHK-1 metadata-completeness findings from
`flavor_catalog/worklogs/checker/ca_w6_charged_lepton_top.md`.  TeX narrative
was left unchanged.

## L023

- Added `pdg_or_equivalent.values[]` metadata entries for the CA-listed
  auxiliary numerical claims retained in TeX: Altmannshofer 2014 `95%` CL
  contours, Belle-II `10.58 GeV` and `50 ab^-1`, and the DUNE `25%`, `about 3`
  years per beam mode, `40%`, and `25%` projection wording.
- Used source URLs and years from `flavor_catalog/references/L023/source_manifest.yaml`.
- Used local snapshot SHA-256 values from the existing sidecar and verified
  with `sha256sum`.

## T020

- Added `PDG2025:T020:higgs_mass_hypothesis` under
  `pdg_or_equivalent.values[]` for the `125 GeV` Higgs mass hypothesis quoted in
  TeX.
- Used the PDG 2025 Higgs review source URL and local snapshot SHA-256 from the
  existing sidecar and `source_manifest.yaml`.

## Status

- Appended `WRITER-DONE` status-history entries with `agent: WA-v2` and
  `cycle: 2` to both sidecars.
- No bibliography changes.
- No unresolved `\textbf{CHECK}` markers introduced.
