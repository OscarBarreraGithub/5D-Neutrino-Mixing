# Writer Worklog: wa_w23_charged_lepton

Family: `charged_lepton`
Processes: `L001`, `L002`, `L007`
Writer agent: `WA`
Timestamp: `2026-05-16T12:15:26-04:00`

## L001 -- \(\mu\to e\gamma\)

- Normalized the process, PDG-or-equivalent, RS relevance, validity, coverage, implementation-difficulty, and reference sections.
- Checked the MEG II 2025 limit, PDG 2025 lagging listing value, MEG/MEG II prior limits, Perez--Randall paper-era coefficient, and live repo default against `L001.yaml` and matching local snapshot hashes.
- Appended `WRITER-DONE` to `L001.yaml`, updated `last_updated_at`, and moved the sidecar version to `0.2.0`.

Open issues for CA:
- Confirm the catalog display policy for the newer MEG II 2025 primary result versus the lagging PDG 2025 muon-listing value.
- Recheck the Perez--Randall NDA mass-scaling convention before any later paper-quality derivation; WA preserved the PKA's model-dependence caveat.

## L002 -- \(\mu\to 3e\)

- Tightened the value/prospect split so the SINDRUM/PDG limit is the current bound and Mu3e sensitivities are explicitly prospects.
- Checked the PDG 2026 S004R4 limit, SINDRUM metadata value, Mu3e phase-I and upgrade sensitivity scales, and EFT/RG source against `L002.yaml` and matching local snapshot hashes.
- Appended `WRITER-DONE` to `L002.yaml`, updated `last_updated_at`, and moved the sidecar version to `0.2.0`.

Open issues for CA:
- `L002.yaml` did not contain a `PKA-DONE` transition before WA work began; WA appended `WRITER-DONE` from the sidecar's current final state and did not invent a PKA transition.
- Future implementation should decide whether a tree-level contact approximation is enough or whether the first repo version should include correlated EFT/RG mixing with \(\mu\to e\gamma\) and \(\mu\)-to-\(e\) conversion.

## L007 -- \(\tau\to\mu\gamma\)

- Normalized the entry to the template fields and made the PDG/Belle source trail explicit for the \(4.2\times10^{-8}\) bound.
- Checked the PDG value, Belle 988 fb\(^{-1}\) source, BaBar predecessor, and Belle II 50 ab\(^{-1}\) projection against `L007.yaml` and matching local snapshot hashes.
- Appended `WRITER-DONE` to `L007.yaml`, updated `last_updated_at`, and moved the sidecar version to `0.2.0`.

Open issues for CA:
- Confirm whether any Belle II \(\tau\to\mu\gamma\) measurement supersedes the Belle 2021 result; the PKA material records only Belle II projections.
- Decide whether a future repo implementation should generalize `muToEGamma.py` or create a separate tau radiative LFV module.

## Batch Notes

- No unresolved `\textbf{CHECK}` markers were left in the TeX files.
- No bibliography files were changed; process-local reference keys remain listed in each entry pending catalog-wide bibliography consolidation.
- No PKA reference snapshots or PKA worklogs were modified.
