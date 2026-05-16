# WA Worklog: wa_wave1_beauty_v2

Date: 2026-05-16
Family: beauty
Cycle: 2
Processes: B009, B015
Checker log: `flavor_catalog/worklogs/checker/ca_wave1_beauty.md`

## Scope

Second-pass writer rework addressing only the CA CHK-1 metadata findings for
B009 and B015.  No TeX narrative rewrite was performed.

## Process fixes

### B009

- Added `source_url` and `access_date` to each
  `pdg_or_equivalent.experimental_inputs` branching-fraction block for the
  five post-2008 inputs listed in `B009.tex`: Belle 2013, Belle 2015,
  Belle II 2025, BaBar 2013, and BaBar 2010.
- Source URLs were copied from
  `flavor_catalog/references/B009/source_manifest.yaml`.
- Local snapshot SHA-256 values were left unchanged and verified with
  `sha256sum`.
- Appended a cycle-2 `WRITER-DONE` status-history transition from
  `WRITER-REWORK`.

### B015

- Added top-level `pdg_or_equivalent.source_url` for the HFLAV December 2024
  source.
- Added per-observable `year`, `source_url`, `access_date`, and `sha256`
  metadata for the HFLAV average, PDG-listed value, BaBar input, Belle input,
  and low/high-\(q^2\) weighted averages flagged by CA.
- Added the Belle II projection luminosity observable with value `50`,
  uncertainty `not quoted`, units `ab^-1`, and the existing
  Huber--Hurth--Lunghi snapshot metadata.
- Added the same metadata fields to the existing Huber--Hurth--Lunghi SM-bin
  observable blocks so every numerical observable block carries the required
  metadata fields.
- Source URLs were copied from
  `flavor_catalog/references/B015/source_manifest.yaml`.
- Local snapshot SHA-256 values were left unchanged and verified with
  `sha256sum`.
- Appended a cycle-2 `WRITER-DONE` status-history transition from
  `WRITER-REWORK`.

## Bibliography and open issues

No bibliography changes.  No new `\textbf{CHECK}` markers or YAML open issues
were introduced.
