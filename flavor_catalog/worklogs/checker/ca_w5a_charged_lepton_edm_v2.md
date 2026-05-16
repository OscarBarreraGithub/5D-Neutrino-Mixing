# CA Worklog: ca_w5a_charged_lepton_edm_v2
**Date**: 2026-05-16
**Family**: charged_lepton_edm
**Process IDs**: L005 L006 E002

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| L005 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| L006 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| E002 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| L005 | \(R_{\mu e}^{\rm Ti}\) | \(<4.3\times10^{-12}\), 90% C.L. | PDG Live 2026 NODE=S004 / DOHMEN 93 | 41061ed0e0e6 |
| L005 | Mu2e program context | roughly four orders beyond previous conversion limits | Mu2e TDR, arXiv:1501.05241 | da321f22d972 |
| L005 | COMET Phase-I Al sensitivity | \(3.1\times10^{-15}\) SES; \(7.0\times10^{-15}\) expected 90% C.L. upper limit | COMET Phase-I TDR, arXiv:1812.09018 | f5521799ba53 |
| L006 | \(G_C/G_F\) | \(<0.0030\), 90% C.L. | PDG Live 2026 S004MC | ba46a431af4e |
| L006 | \(P_{M\bar M}\) | \(<8.3\times10^{-11}\), 90% C.L., 0.1 T field | PDG Live 2026 S004MC / Willmann 1999 | ba46a431af4e |
| L006 | MACE prospective improvement | more than two orders of magnitude | Bai 2022 Snowmass, arXiv:2203.11406 | 04c55fc736d4 |
| L006 | MACE prospective reach | beyond \(10^{-13}\) conversion probability | Bai 2024 MACE CDR, arXiv:2410.18817 | 6acea5b25a8a |
| E002 | \(|d_\mu|\) direct limit | \(<1.8\times10^{-19}\ e\,\mathrm{cm}\), 95% CL | PDG Live 2026 S004EDM | 2918e0eb7819 |
| E002 | PDG combined measurement | \((0.0\pm0.9)\times10^{-19}\ e\,\mathrm{cm}\) | PDG Live 2026 S004EDM | 2918e0eb7819 |
| E002 | Bennett 2009 BNL measurement | \(-0.1(0.9)\times10^{-19}\ e\,\mathrm{cm}\); \(|d|<1.9\times10^{-19}\ e\,\mathrm{cm}\), 95% CL | Bennett et al., arXiv:0811.1207 | 16ce7a84e7ba |
| E002 | PSI frozen-spin projection | \(p=125\,\mathrm{MeV}/c\), \(|B|=3\,\mathrm{T}\), \(1\,\mathrm{GV/m}\), \(\sigma(d_\mu)\le6\times10^{-23}\ e\,\mathrm{cm}\) | Adelmann et al., arXiv:2102.08838 | bc86e027b5ee |
| E002 | PSI status projection | first phase by 2026; factor 100 improvement goal by early 2030s | Renga 2024, arXiv:2409.20050 | c9c8a8cb2c32 |

## Issues (if any)
- None. WA-v2 added `pdg_or_equivalent.values` traceability for all projection/context numbers flagged in `ca_w5a_charged_lepton_edm.md`.
- CHK-1: Numeric physics/prospect claims in the TeX now resolve to `pdg_or_equivalent` blocks with source URL, access date, snapshot path, and sha256. Upper limits carry CL metadata rather than uncertainty values; prospective values use `uncertainty: null` where applicable.
- CHK-2: All process-local source keys in the TeX key-reference sections resolve to `source_manifest.yaml` entries. Each manifest entry has a non-empty `snapshot_path`, and `git ls-files flavor_catalog/references/{L005,L006,E002}` shows the snapshots are tracked.
- CHK-3: `ls flavor_catalog/references/{L005,L006,E002}/` showed only `.txt` snapshots, `source_manifest.yaml`, and L005 `sha256sums.txt`; no publisher PDFs are tracked.
- CHK-4: Each sidecar has `WRITER-INITIATED` followed by `WRITER-DONE` with ISO 8601 timestamps, and the most recent pre-check state was the WA-v2 `WRITER-DONE` entry. This re-check appends `CHECKER-DONE`.
- CHK-5: The requested `rg -l -E "<process keyword>" ...` form fails in this environment because ripgrep treats `-E` as an encoding flag. The semantic reruns with `-e` found no L005 or L006 implementation hits. For E002, notebook-only EDM mentions appear unless `-g '!*.ipynb'` is used; with notebooks excluded there is no muon-EDM implementation. The cited adjacent `muToEGamma` and scan file lines exist and are off-diagonal \(\mu\to e\gamma\), not these observables.
- CHK-6: Difficulty labels are consistent. L005 is HIGH because coherent conversion needs new lepton-quark operators plus target-specific nuclear/capture inputs. L006 is MEDIUM because it needs a new four-lepton/bound-state observable but no lattice or hadronic long-distance inputs. E002 is HIGH because a flavor-diagonal CP-odd lepton dipole needs new matching and likely running.
- CHK-7: Searches outside `flavor_catalog/` found only plan/prospect mentions for these processes, plus unrelated `mu_had` variable names; no rc1.1 load-bearing process number is contradicted.
- CHK-8: No `\textbf{CHECK}`, unresolved `\cite`/`\ref`, TODO/FIXME marker, missing-label marker, or obvious math-mode/bracing issue was found in the three TeX files.
