# Writer Worklog: wa_w23_kaon_charm_edm

Family: mixed (`kaon`, `charm`, `edm_neutrino`)
Processes: K001, K002, C001, E001
Writer: WA
Date: 2026-05-16

## Scope and Source Decisions

Read the required plan-v1 Section B/Section D writer instructions, the
orchestrator decisions, each PKA draft, each process sidecar, each PKA worklog,
and each process-local source manifest.  I did not modify reference snapshots,
PKA worklogs, catalog indexes, shared LaTeX files, or `catalog.bib`.

All quantitative physics values left in the TeX now carry an adjacent
process-local source key or are part of the code-coverage evidence already
recorded in the sidecar.  No `\textbf{CHECK}` markers remain in this batch.
Because the prompt forbade shared-file edits, bibliography consolidation is
left as process-local source-key usage rather than a `catalog.bib` patch.

## Per-Process Changes

### K001 -- `epsilon_K`

- Tightened the process, RS-relevance, and post-2008 prose; normalized notation
  to `\varepsilon_K` in the TeX.
- Added explicit source keys next to the PDG value, BGS SM prediction, FLAG bag
  inputs, and CFW RS-scale numbers.
- Appended `WRITER-DONE` to `status_history` and updated `last_updated_at`.

Open issues for CA:
- The sidecar entered WA with latest status `WRITER-INITIATED`, not
  `PKA-DONE`; I appended the requested `WRITER-DONE` transition without
  rewriting PKA-owned history.
- CA should review the existing PKA issue about PDG review extract versus
  pdgLive JSON as the canonical source path.

### K002 -- `Delta m_K`

- Tightened wording around the CPT-assuming PDG value and conservative
  long-distance-limited interpretation.
- Added source keys for PDG, KTeV, lattice long-distance developments, FLAG,
  and the 2024 BSM kaon-mixing bag-input source.
- Appended `WRITER-DONE` to `status_history` and updated `last_updated_at`.

Open issues for CA:
- The sidecar entered WA with latest status `WRITER-INITIATED`, not
  `PKA-DONE`; I appended the requested `WRITER-DONE` transition without
  rewriting PKA-owned history.
- CA should decide whether the no-CPT PDG fit should remain in the final
  displayed value or stay as a sidecar companion value.
- CA should decide whether the modern bridge evidence changes only prose or
  eventually warrants a coverage-label update by a PKA/orchestrator owner.

### C001 -- D-mixing

- Added explicit source keys for the HFLAV `x_D`, `y_D`, interval, `delta_y`,
  derived `DeltaGamma/Gamma`, PDG `Delta m_D`, LHCb 2021, and CFW RS numbers.
- Removed exact significance numbers from the prose because they were not
  structured values in the sidecar; kept the qualitative post-2008 conclusion.
- Appended `WRITER-DONE` to `status_history` and updated `last_updated_at`.

Open issues for CA:
- CA should choose whether `Delta y` or `DeltaGamma_D/Gamma_D = 2y_D` is the
  preferred headline convention for the orchestrator-requested `delta_y` block.

### E001 -- electron EDM

- Added explicit source keys for the PDG limit, Roussy measurement, ACME 2018
  benchmark, CFW baseline, and Panico/Pomarol/Riembau EFT interpretation.
- Tightened the EDM relevance and post-2008 text while preserving the PKA's
  caveat that this is not implemented in the quark-only `Delta F=2` pipeline.
- Appended `WRITER-DONE` to `status_history` and updated `last_updated_at`.

Open issues for CA:
- CA should decide how prominently the final entry should caveat paramagnetic
  molecule sensitivity to semileptonic CP-odd electron-nucleon operators.
- If this graduates to live code, PI/orchestrator input is still needed on
  lepton-sector matching assumptions and whether EDMs are hard constraints or
  contextual diagnostics.
