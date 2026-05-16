# WA Worklog: wa_w5a_kaon_charm

**Date**: 2026-05-16
**Family**: mixed kaon/charm
**Writer agent**: WA
**Process IDs**: K008, C005

## Scope
Polished the assigned PKA drafts only: `K008` and `C005`.  Reference snapshots
and PKA worklogs were read but not modified.

## Changes
- K008: tightened the process, RS relevance, post-2008, and validity prose;
  added `\texorpdfstring` to the section heading; tied numerical claims to
  sidecar-listed local snapshots and process-local source keys.
- K008: appended a `WRITER-DONE` status-history entry and updated
  `last_updated_at`.
- C005: tightened prose, normalized source-key references around the PDG,
  Belle, BaBar, CFW, and rare-charm theory claims, and clarified the PDG
  record's absence of an LHCb `D^0 -> e^+e^-` entry.
- C005: appended a `WRITER-DONE` status-history entry and updated
  `last_updated_at`.

## Source Decisions
- Reused existing PKA sources only; no new source snapshots or reference files
  were added.
- K008 main measured/supporting values are cited in TeX by the local snapshot
  filenames listed under `source_shas`, since the PKA-owned measured-value YAML
  blocks do not carry explicit `source_key` fields.
- C005 numerical claims are cited by sidecar `source_key` entries and their
  local snapshot filenames.

## Bibliography
No bibliography changes.  No `catalog.bib` merge was attempted in this batch.

## Open Checks For CA
- K008: re-query PDG `S013.20` before freeze, as requested in the PKA sidecar,
  in case the 2026 pdgLive entry changes after the snapshot.
- K008: confirm that K010 carries `K_S -> pi^0 e^+ e^-` as its own primary
  process; K008 uses it only as an indirect-CP supporting input.
- C005: reconcile the dispatch note's "LHCb upper limit" wording with the
  accessed PDG/API record, which lists Belle 2010 as the current used
  `D^0 -> e^+e^-` limit and no LHCb entry.

No new `\textbf{CHECK}` markers were introduced.
