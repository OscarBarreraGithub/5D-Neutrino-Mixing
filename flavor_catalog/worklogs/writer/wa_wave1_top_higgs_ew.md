# WA Worklog: wa_wave1_top_higgs_ew

Family: `top_higgs_ew`
Processes: `T001`, `T010`
Writer timestamp: 2026-05-16T11:32:50-04:00

## T001 -- \(t \to c Z\)

- Tightened the process and model-dependence language while preserving the PKA claims.
- Added sidecar-resolving source/value tags beside the PDG, ATLAS, CMS, SM, and HL-LHC numerical claims.
- Appended `WRITER-DONE` to `T001.yaml` and updated `last_updated_at`.

Checker-facing issues:

- The PKA open issue remains: decide whether the headline value should quote both ATLAS left- and right-handed tensor-SMEFT benchmarks or choose one.
- A live implementation still needs an RS-to-SMEFT or direct `tcZ` convention before applying this limit to scan points.

## T010 -- \(Z \to b\bar b\) pole observables

- Normalized the section headings and tightened the combined `R_b^0`, `A_{\rm FB}^{0,b}`, and `A_b` discussion.
- Added sidecar/snapshot tags for the PDG values, LEP/SLC pull, Freitas two-loop context, and FCC-ee projection.
- Appended `WRITER-DONE` to `T010.yaml` and updated `last_updated_at`.

Checker-facing issues:

- The assignment combines plan rows `T010` and `T011`; checker/orchestrator should confirm whether `T011` should be marked covered by this combined entry or kept as a separate process.

## Source and Bibliography Notes

- Checked local snapshot hashes against the sidecar `source_shas`; all process-local text snapshots matched the recorded hashes.
- No `\textbf{CHECK}` markers were added.
- No `catalog.bib` patch was made because this WA batch was restricted to `flavor_catalog/processes/top_higgs_ew/` and `flavor_catalog/worklogs/writer/`.
