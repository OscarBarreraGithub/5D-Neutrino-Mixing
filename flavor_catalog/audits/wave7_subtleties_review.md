# Wave 7 Subtleties Peer Review

Reviewer: subtlety-checker-wave7
Writer commit: `837e709` (`flavor-catalog(wave7): thread 6 cross-cutting subtleties from external review into existing process entries`)
Baseline parent: `dbf2b960d3854804385afadb53bef2796a48e7ff`

## Verdict

RETURN-TO-WRITER

The writer touched exactly the declared 22 process `.tex` files and their 22
matching `.yaml` sidecars, plus the writer worklog.  All TeX additions are in
the "Constraint validity and model dependence" subsection or the local
near-equivalent spelling, and the sidecars all received a `SUBTLETY-ADDED`
history entry from `subtlety-writer-wave7`.

The only return item is CHK-C for T007: the inserted kinematics note adds
`125 GeV` to the TeX.  That is a true numerical physics claim, even though it
comes from the external-review wording.  Wave-7 instructions said these
subtlety edits should not introduce new numerical values.  Rephrase as "the
observed Higgs" or otherwise remove the number.

## Per-file table

Legend: CHK-A phrasing present; CHK-B placement/order; CHK-C no new numerical
claim; CHK-D sidecar history; CHK-E no value/uncertainty/sha/citation change;
CHK-F no unintended TeX scope.

| file path | subtlety_id | CHK-A | CHK-B | CHK-C | CHK-D | CHK-E | CHK-F |
|---|---:|---|---|---|---|---|---|
| `flavor_catalog/processes/beauty/B011.tex` | D | OK | OK | OK | OK | OK | OK |
| `flavor_catalog/processes/beauty/B015.tex` | A | OK | OK | OK | OK | OK | OK |
| `flavor_catalog/processes/beauty/B017.tex` | A+C | OK | OK | OK | OK | OK | OK |
| `flavor_catalog/processes/beauty/B018.tex` | A+C | OK | OK | OK | OK | OK | OK |
| `flavor_catalog/processes/beauty/B019.tex` | A+C | OK | OK | OK | OK | OK | OK |
| `flavor_catalog/processes/charm/C001.tex` | B | OK | OK | OK | OK | OK | OK |
| `flavor_catalog/processes/charm/C002.tex` | B | OK | OK | OK | OK | OK | OK |
| `flavor_catalog/processes/charm/C003.tex` | B | OK | OK | OK | OK | OK | OK |
| `flavor_catalog/processes/charm/C004.tex` | B | OK | OK | OK | OK | OK | OK |
| `flavor_catalog/processes/top_higgs_ew/T001.tex` | B | OK | OK | OK | OK | OK | OK |
| `flavor_catalog/processes/top_higgs_ew/T002.tex` | B | OK | OK | OK | OK | OK | OK |
| `flavor_catalog/processes/top_higgs_ew/T005.tex` | B | OK | OK | OK | OK | OK | OK |
| `flavor_catalog/processes/top_higgs_ew/T006.tex` | B | OK | OK | OK | OK | OK | OK |
| `flavor_catalog/processes/top_higgs_ew/T007.tex` | B+F | OK | OK | FAIL | OK | OK | OK |
| `flavor_catalog/processes/top_higgs_ew/T010.tex` | A | OK | OK | OK | OK | OK | OK |
| `flavor_catalog/processes/top_higgs_ew/EW001.tex` | A | OK | OK | OK | OK | OK | OK |
| `flavor_catalog/processes/top_higgs_ew/EW002.tex` | A | OK | OK | OK | OK | OK | OK |
| `flavor_catalog/processes/edm_neutrino/E001.tex` | E | OK | OK | OK | OK | OK | OK |
| `flavor_catalog/processes/edm_neutrino/E004.tex` | E | OK | OK | OK | OK | OK | OK |
| `flavor_catalog/processes/edm_neutrino/E006.tex` | E | OK | OK | OK | OK | OK | OK |
| `flavor_catalog/processes/edm_neutrino/E008.tex` | E | OK | OK | OK | OK | OK | OK |
| `flavor_catalog/processes/edm_neutrino/E009.tex` | E | OK | OK | OK | OK | OK | OK |

## Issues

1. `flavor_catalog/processes/top_higgs_ew/T007.tex`: CHK-C fails because the
   added note says "for a 125 GeV Higgs".  This was absent from the parent TeX
   and is not accompanied by new sidecar value metadata.  Since the assignment
   forbade new numerical claims, the clean fix is to say "for the observed
   Higgs" while preserving the on-shell `h -> tc` caveat.

No unintended process TeX files were modified.  A diff grep found no edits to
`pdg_or_equivalent`, uncertainties, sha256 fields, source snapshots, or citation
keys.  Spot checks covered one file from each touched family: B011, C001, T010,
and E001.

## Recommendation

Return to the writer for the one T007 wording fix.  After that, the batch should
be approvable without further physics or sidecar changes.
