# WA Worklog: wa_w5b_top_higgs_ew_v2

Family: `top_higgs_ew`
Processes: `T006`
Writer: `WA-v2`
Cycle: 2
Timestamp: `2026-05-16T15:27:46-04:00`
Checker log: `flavor_catalog/worklogs/checker/ca_w5b_top_higgs_ew.md`

## Required Reading

- Read plan v1 Section B/D, CA worklog `ca_w5b_top_higgs_ew.md`, the current
  `T006.tex`, `T006.yaml`, `references/T006/source_manifest.yaml`, and
  `worklogs/pka/T006.md`.

## CA Finding Addressed

- T006 CHK-1: the TeX-cited SM
  `B(t -> u g) ~= 3.6e-14` context value was present only under
  `theory_context`, not under `pdg_or_equivalent.values`.

## Changes

- Promoted `AguilarSaavedra2004:T006:t_ug_SM` into
  `pdg_or_equivalent.values`, preserving the existing value, uncertainty,
  units, source URL, access date, and sha256.
- Added the local snapshot path
  `flavor_catalog/references/T006/aguilar_saavedra_hepph0409342_tug_sm.txt`;
  `sha256sum` confirmed
  `ed810970e049bd20488740565ba46036ab0c29892a9991db7819a6c3805f2098`.
- Left `T006.tex` unchanged because the CA allowed promotion as the fix and the
  existing prose already cites `AguilarSaavedra2004:T006:t_ug_SM`.
- Appended a cycle-2 `WRITER-DONE` status-history transition to `T006.yaml` and
  updated `last_updated_at`.

## Scope Notes

- No `CHECKER-DONE` transition was added.
- No bibliography, index, LaTeX build files, or other family files were
  modified.
