# Writer Worklog: wa_w4_kaon_edm

Family: mixed (`kaon`, `edm_neutrino`)
Processes: K017, E004, E006, E008
Writer: WA
Date: 2026-05-16

## Scope and Source Decisions

Read plan v1 Section B/D, the orchestrator decisions, each PKA TeX/YAML
sidecar, each PKA worklog, and each process-local source manifest and snapshot
set.  The requested `kaon_edm` family directory does not exist in this checkout,
so I treated the explicit process IDs as the batch scope: K017 under
`processes/kaon/` and E004/E006/E008 under `processes/edm_neutrino/`.

Writes were limited to the four process TeX/YAML files and this writer
worklog.  I did not modify reference snapshots, PKA worklogs, catalog indexes,
shared LaTeX files, or `catalog.bib`.  Because the batch hard rules prohibit
shared bibliography edits, accepted references remain as process-local source
keys.  Local snapshot hashes were checked against the YAML `source_shas`.

No explicit checker markers were introduced.  All quantitative physics values
left in the TeX now have adjacent process-local source keys or are code
file-line evidence already recorded by the PKA sidecar.

## Per-Process Changes

### K017 -- kaon \(R_K\) LFU

- Tightened the process, RS-relevance, and post-2008 prose while preserving the
  charged-current/lepton-extension classification.
- Added adjacent source keys for the PDG API average, NA62 and KLOE inputs, the
  Cirigliano--Rosell SM prediction, FlaviaNet context, and the CFW baseline.
- Appended `WRITER-DONE` to `status_history` and updated `last_updated_at`.

Open for CA:
- Decide whether to keep the process-local `PDG2025:K017RK` key if PDGLive
  publishes a 2026 API edition before final signoff.
- If promoted later, PI input is still needed on the BSM parameterization
  scalar, right-handed charged-current, sterile-neutrino, or EFT ratio wrapper.

### E004 -- neutron EDM

- Tightened the process and PDG-value prose and added source keys next to the
  direct limit, Abel measurement, and superseded Pendlebury comparison.
- Rephrased the lattice-development sentence to avoid over-specific operator
  counting while keeping the PKA-documented theta/Weinberg/qCEDM scope.
- Appended `WRITER-DONE` to `status_history` and updated `last_updated_at`.

Open for CA:
- Decide how much theta-term, qCEDM, and Weinberg-operator interpretation should
  remain in E004 versus E008/E009.
- A live implementation still needs PI choices for the CP-odd low-energy basis,
  QCD running, thresholds, and neutron matrix-element inputs.

### E006 -- mercury-199 EDM

- Made the Graner direct limit the displayed catalog value and retained the
  central-value/erratum caveat in prose rather than displaying a signed central
  value.
- Added source keys for the PDG Hg-derived neutron rows and Sahoo
  interpretation-level limits.
- Appended `WRITER-DONE` to `status_history` and updated `last_updated_at`.

Open for CA:
- Verify the signed central value against the PRL erratum before any final
  display beyond the absolute \(^{199}\mathrm{Hg}\) limit.
- Coordinate with E004/E008/E009 so neutron/proton/qCEDM/Weinberg translations
  are not double-counted as independent measurements.

### E008 -- quark chromo-EDM bounds

- Recast the entry as benchmark EFT translations rather than direct PDG
  measurements.
- Added adjacent source keys for the neutron and mercury anchors, the
  Pospelov--Ritz neutron qCEDM translation, the Olive--Pospelov--Ritz--Santoso
  mercury estimate, CFW, Koenig--Neubert--Straub, and Bhattacharya 2024.
- Appended `WRITER-DONE` to `status_history` and updated `last_updated_at`.

Open for CA:
- Confirm whether the qCEDM numbers should remain headline benchmark bounds or
  be demoted to translation examples in the final catalog.
- A live implementation still requires PI choices on the CP-odd basis,
  Peccei--Quinn/\(\bar\theta\) treatment, and hadronic/nuclear defaults.
