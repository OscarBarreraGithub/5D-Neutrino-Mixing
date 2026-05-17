# Writer Worklog: WA-w8-B-rare-leptonic

Date: 2026-05-17
Batch ID: `WA-w8-B-rare-leptonic`
Agent: `WA-WA-w8-B-rare-leptonic`
Cycle: 1
Processes: `B007` (`B_s,d -> e+e-`), `B008` (`B_s,d -> tau+tau-`)
Referenced CA log path for cycle 1: `flavor_catalog/worklogs/checker/ca_w8_B_rare_leptonic.md`

## Scope

`B007`: polished the prose, normalized the validity-section heading, removed a loose order-of-magnitude prose comparison, and expanded code-coverage evidence to exact `quarkConstraints/deltaf2.py` file-line references.  The SECONDARY-tier note was already present in the Relevance section.

`B007.yaml`: kept the existing SECONDARY fields, moved measured experimental limits into `pdg_or_equivalent.values`, left dataset metadata in `recent_experimental_inputs`, left CFW theory reference scales in `paper_era_reference`, and replaced historical numeric-limit duplication with value IDs.  Removed stale PKA open issues.

`B008`: removed precise theory-prediction/enhancement numerals from prose while retaining source traceability in auxiliary theory inputs, tightened the code-coverage paragraph, and kept the SECONDARY-tier Relevance note.

`B008.yaml`: kept the existing SECONDARY fields, moved the canonical and historical measured limits into `pdg_or_equivalent.values`, moved luminosity and sqrt(s) metadata into `supporting_measurements`, and removed noncanonical related-mode and theory-bound numerals from the active sidecar claims.

## Source Checks

Re-ran checksum validation against the tracked snapshot checksum files:

```bash
(cd flavor_catalog/references/B007 && sha256sum -c sha256sums.txt)
sha256sum -c flavor_catalog/references/B008/sha256sums.txt
```

Result: all 7 B007 text snapshots and all 9 B008 text snapshots returned `OK`.

Manifest audit:

```bash
git ls-files flavor_catalog/references/B007/*.txt flavor_catalog/references/B008/*.txt \
  flavor_catalog/references/B007/source_manifest.yaml \
  flavor_catalog/references/B008/source_manifest.yaml
```

All source-manifest entries have non-empty `snapshot_path` and
`sha256_of_local_snapshot` fields, and all referenced text snapshots are tracked.

## Bibliography

No process-local `refs.bib` or other `.bib` file exists under
`flavor_catalog/references/B007` or `flavor_catalog/references/B008`.  This repo
checkout also has no `flavor_catalog/references/catalog.bib`; shared
bibliography edits were outside this batch's allowed write scope.

Unresolved process-local keys retained pending catalog-wide bibliography
consolidation:

`B007`: `PDG2026_BsBdEe`, `LHCb2020_BsBdEe`,
`BobethEtAl2013_BqllSM`, `CDF2009_BsBdEe`,
`FleischerJaarsmaTetlalmatziXolocotzi2017_BqllNP`,
`BaBar2008_BdEe`, `CsakiFalkowskiWeiler2008_RSFlavor`.

`B008`: `PDG2026_BsTauTau`, `PDG2026_BdTauTau`,
`AaijEtAl2017_BqTauTau`, `AubertEtAl2006_BdTauTau`,
`LeesEtAl2017_BKTauTauRelated`, `BobethEtAl2014_BqTauTauSM`,
`CapdevilaEtAl2018_BSTauTauNP`, `BordoneNavarro2023_BSTauTauNP`,
`CsakiFalkowskiWeiler2008_RSFlavor`.

## Open Issues

No `\textbf{CHECK}` markers were introduced.

`B007`: no open issues.

`B008`: future implementation must choose between independent one-mode limits
and a two-dimensional likelihood for simultaneous `B_s`/`B_d` tauonic decays.

## Status Transitions

Appended `WRITER-DONE` with `agent_id: "WA"`, `cycle: 1`, and
`at: "2026-05-17T02:37:14-04:00"` to both `B007.yaml` and `B008.yaml`; updated
both `last_updated_at` fields to the same timestamp.
