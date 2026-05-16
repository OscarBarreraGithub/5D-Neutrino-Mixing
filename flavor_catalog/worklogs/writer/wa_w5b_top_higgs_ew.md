# WA Worklog: wa_w5b_top_higgs_ew

Family: `top_higgs_ew`
Processes: `T006`, `T016`, `T017`
Writer: `WA`
Timestamp: `2026-05-16T14:44:48-04:00`

## Summary

### T006 `t -> u g`

- Added the explicit standard-notation field in the Process section and
  normalized branching-fraction notation.
- Tightened the PDG/equivalent section and attached sidecar value IDs to the
  ATLAS, CMS, and SM-context numerical claims.
- Kept the PKA open issues: CA should confirm the headline choice between the
  current PDG/ATLAS SMEFT convention and the numerically stronger older CMS
  comparison limit.

### T016 `Z -> e tau`

- Preserved the PDG/ATLAS headline limit and CMS 2025 cross-check while adding
  explicit sidecar traces for the repeated PDG limit and Tera-Z scale.
- Trimmed repeated post-2008 language without adding new physics claims.
- Open issue for CA: future implementation still needs a left/right
  off-diagonal `Z e tau` convention and a combination policy with tau LFV.

### T017 `Z -> mu tau`

- Removed filler from the PDG/equivalent section and added an explicit sidecar
  trace for the post-2008 ATLAS `10^-6`-level statement.
- Preserved the PDG/ATLAS headline observed limit and CMS 2025 independent
  cross-check.
- Open issue for CA: verify future PDG treatment of CMS 2025 and the eventual
  `Z_mu_tau` coupling convention if promoted to code.

## Status And Sources

- Appended `WRITER-DONE` status-history entries to `T006.yaml`, `T016.yaml`,
  and `T017.yaml`; updated only `last_updated_at` outside the status history.
- No checker placeholders were introduced.
- No bibliography files were modified; all references used here remain
  process-local keys or sidecar value IDs.
