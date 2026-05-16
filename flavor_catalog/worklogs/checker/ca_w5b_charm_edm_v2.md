# CA Worklog: ca_w5b_charm_edm_v2
**Date**: 2026-05-16
**Family**: charm_edm / `edm_neutrino`
**Process IDs**: E007

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| E007 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| E007 | \(|d(^{225}\mathrm{Ra})|\) current direct limit | \(<1.4\times10^{-23}\ e\,\mathrm{cm}\), 95% CL; improvement factor 36 | Bishof et al. 2016, arXiv:1606.04931, access 2026-05-16 | a9042e38b7b3 |
| E007 | \(|d(^{225}\mathrm{Ra})|\) first measurement | \(<5.0\times10^{-22}\ e\,\mathrm{cm}\), 95% CL | Parker et al. 2015, arXiv:1504.07477, access 2026-05-16 | ff5adb5cdce9 |
| E007 | \(d_A(^{129}\mathrm{Xe})\) best direct result and limit | \((1.4 \pm 6.6_{\rm stat} \pm 2.0_{\rm syst})\times10^{-28}\ e\,\mathrm{cm}\); \(<1.4\times10^{-27}\ e\,\mathrm{cm}\), 95% CL | Sachdeva et al. 2019, arXiv:1909.12800, access 2026-05-16 | abce59966ada |
| E007 | \(d_{\mathrm{Xe}}\) independent Heidelberg result and limit | \((-4.7\pm6.4)\times10^{-28}\ e\,\mathrm{cm}\); \(<1.5\times10^{-27}\ e\,\mathrm{cm}\), 95% CL | Allmendinger et al. 2019, arXiv:1904.12295, access 2026-05-16 | 94e4abc0ed9f |

## Evidence notes
- CHK-1: Scanned TeX numeric tokens. The load-bearing measured values at lines 27-34, 39-48, the 95% CL labels, the Ra improvement factor 36, and the \(10^{-27}\ e\,\mathrm{cm}\)-level Xe statement all trace to `pdg_or_equivalent` blocks with year, value/uncertainty, source URL, access date, snapshot path, and sha256. Remaining numbers are process/source identifiers, isotope labels, reference years, or code line citations rather than standalone physics values.
- CHK-2: All TeX source keys and process-local snapshot names resolve to `flavor_catalog/references/E007/source_manifest.yaml`; each manifest entry has a non-empty tracked `snapshot_path` under `flavor_catalog/references/E007/`.
- CHK-3: `flavor_catalog/references/E007/` contains only `.txt` snapshots and `source_manifest.yaml`; no `.pdf` files are tracked.
- CHK-4: `status_history` contains `WRITER-INITIATED` from PKA and the WA-v2 `WRITER-DONE` in order with ISO 8601 timestamps; before this CA update the most recent entry was `WRITER-DONE`.
- CHK-5: The literal requested `rg -l -E` form is not accepted by this ripgrep because `-E` is parsed as an encoding flag. Equivalent focused `rg -l -e` searches for radium/xenon/Schiff/diamagnetic/chromo-EDM/Weinberg terms over the required implementation directories returned no source files after excluding notebooks; a generic `EDM|dipole` search found only `flavorConstraints/muToEGamma.py` and its README docs, matching the TeX `NO` code-coverage claim.
- CHK-6: `HIGH` is consistent with the rubric: E007 needs new CP-odd quark/gluon and semileptonic matching, RG to hadronic scales, and nuclear/atomic Schiff-response inputs, not the existing \(\Delta F=2\) SLL/SLR/VLL/VRR/LR basis.
- CHK-7: No rc1.1 load-bearing number is revised; the 2008 CFW source is cited only as RS-flavor baseline context.
- CHK-8: No unresolved `\textbf{CHECK}`, `\cite`, `\ref`, or `\label` markers found; math-mode/bracing inspection passed.

## Issues (if any)
- None.
