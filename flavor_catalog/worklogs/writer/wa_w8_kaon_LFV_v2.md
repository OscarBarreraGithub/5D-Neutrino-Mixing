# WA Worklog: wa_w8_kaon_LFV_v2

Date: 2026-05-17
Family label: kaon LFV, Wave-8 SECONDARY
Processes: `K020`
Agent: WA-v2
Cycle: 2

## Scope

Second-pass rework for the K020 CHK-1 failure listed in
`flavor_catalog/worklogs/checker/ca_w8_kaon_LFV.md`. Edits were limited to
`flavor_catalog/processes/secondary/kaon/K020.yaml` and this new worklog.
K019 and K021 were not modified.

## K020 Rework

- Promoted the Sher/BNL E865
  `BR(K+ -> pi+ mu+ e-) < 2.1 x 10^-11` 90% CL measured upper limit into
  `pdg_or_equivalent.values`.
- Populated the promoted value block with value, uncertainty, limit type, CL,
  year, units, source URL, access date, snapshot path, and local snapshot
  sha256 metadata from the existing K020 source manifest.
- Appended the required cycle-2 `WRITER-DONE` status-history entry for
  `WA-v2` and updated `last_updated_at`.
- Left `K020.tex` unchanged because the cycle-1 prose already accurately
  distinguishes the E865-only limit from the PDG combined E865/E777 limit.

## Source Checks

- Confirmed the Sher/E865 source snapshot
  `flavor_catalog/references/K020/sher_2005_arxiv_hep_ex_0502020.txt`
  contains the `2.1 x 10^-11` 90% CL upper limit in the arXiv abstract text.
- Spot-checked the local snapshot hash with `sha256sum`:
  `0ae93baf10b71d30e476330661aee67a6a3c9786fa6436d5caa015446f5ad551`.
- The hash matches `flavor_catalog/references/K020/sha256sums.txt`,
  `flavor_catalog/references/K020/source_manifest.yaml`, and the existing
  `K020.yaml` `source_shas` entry.

## Open Issues

- None expected.
