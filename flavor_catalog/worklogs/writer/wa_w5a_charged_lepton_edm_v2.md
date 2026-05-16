# WA Worklog: wa_w5a_charged_lepton_edm_v2

**Date**: 2026-05-16
**Family**: charged_lepton_edm
**Cycle**: 2
**Agent**: WA-v2
**Input CA log**: `flavor_catalog/worklogs/checker/ca_w5a_charged_lepton_edm.md`
**Processes**: L005, L006, E002

## Scope

Addressed only the CA-listed CHK-1 metadata-placement findings.  The TeX
narratives were left unchanged; the quoted projection/context numbers now have
complete `pdg_or_equivalent.values` traceability blocks in the sidecars.

## Fixes

- L005: added `pdg_or_equivalent.values` entries for the Mu2e roughly
  four-orders-of-magnitude program context and the COMET Phase-I
  `3.1e-15` single-event sensitivity plus `7.0e-15` expected 90% C.L. upper
  limit.  Metadata uses the existing Mu2e and COMET snapshots and sha256s.
- L006: added `pdg_or_equivalent.values` entries for the Snowmass/MACE
  prospective improvement by more than two orders of magnitude and the 2024
  MACE reach beyond `1.0e-13`.  Metadata uses the existing Bai 2022 and Bai
  2024 snapshots and sha256s.
- E002: added `pdg_or_equivalent.values` entries for the PSI projection/status
  numbers quoted in the TeX: `125 MeV/c`, `3 T`, `1 GV/m`,
  `6.0e-23 e cm`, first-phase target year `2026`, improvement factor `100`,
  and the `early 2030s` ultimate-goal period.  Metadata uses the existing
  Adelmann 2021 and Renga 2024 snapshots and sha256s.
- All three sidecars received a cycle-2 `WRITER-DONE` status-history entry
  with `agent: WA-v2`; no `CHECKER-DONE` entry was added.

## Verification

- Recomputed sha256 values for the six existing source snapshots used in the
  new blocks.
- Parsed the three edited YAML sidecars after patching.
- Confirmed no changes were made to `latex/*`, `catalog_index.*`, or other
  families for this batch.
