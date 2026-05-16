# CA Worklog: ca_w5a_charged_lepton_edm
**Date**: 2026-05-16
**Family**: charged_lepton_edm
**Process IDs**: L005 L006 E002

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| L005 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| L006 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| E002 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| L005 | \(R_{\mu e}^{\rm Ti}\) | \(<4.3\times10^{-12}\), 90% C.L. | PDG Live 2026 NODE=S004 / DOHMEN 93 | 41061ed0e0e6 |
| L005 | COMET Phase-I Al sensitivity / UL | \(3.1\times10^{-15}\) SES; \(7.0\times10^{-15}\) 90% UL | COMET Phase-I TDR, arXiv:1812.09018 | f5521799ba53 |
| L005 | Mu2e program context | roughly four orders beyond previous conversion limits | Mu2e TDR, arXiv:1501.05241 | da321f22d972 |
| L006 | \(G_C/G_F\) | \(<0.0030\), 90% C.L. | PDG Live 2026 S004MC | ba46a431af4e |
| L006 | \(P_{M\bar M}\) | \(<8.3\times10^{-11}\), 90% C.L., 0.1 T field | PDG Live 2026 S004MC / Willmann 1999 | ba46a431af4e / 3e22390c63c8 |
| L006 | MACE prospective reach | beyond \(10^{-13}\); more than two orders of magnitude improvement | Bai 2024 MACE CDR / Bai 2022 Snowmass | 6acea5b25a8a / 04c55fc736d4 |
| E002 | \(|d_\mu|\) direct limit | \(<1.8\times10^{-19}\ e\,\mathrm{cm}\), 95% CL | PDG Live 2026 S004EDM | 2918e0eb7819 |
| E002 | PDG combined measurement | \((0.0\pm0.9)\times10^{-19}\ e\,\mathrm{cm}\) | PDG Live 2026 S004EDM | 2918e0eb7819 |
| E002 | Bennett 2009 BNL measurement | \(-0.1(0.9)\times10^{-19}\ e\,\mathrm{cm}\); \(|d|<1.9\times10^{-19}\ e\,\mathrm{cm}\), 95% CL | Bennett et al., arXiv:0811.1207 | 16ce7a84e7ba |
| E002 | PSI frozen-spin projection | \(p=125\,\mathrm{MeV}/c\), \(|B|=3\,\mathrm{T}\), \(1\,\mathrm{GV/m}\), \(\sigma(d_\mu)\le6\times10^{-23}\ e\,\mathrm{cm}\) | Adelmann et al., arXiv:2102.08838 | bc86e027b5ee |
| E002 | PSI status projection | first phase by 2026; factor 100 improvement goal by early 2030s | Renga 2024, arXiv:2409.20050 | c9c8a8cb2c32 |

## Issues (if any)
- L005: CHK-1 fails under the strict sidecar-placement rule. The headline PDG Ti limit is in `pdg_or_equivalent.primary_current_limit` with year, value, CL, source URL, access date, snapshot path, and sha256. However, `L005.tex` also quotes Mu2e four-orders context and COMET Phase-I \(3.1\times10^{-15}\) / \(7.0\times10^{-15}\) / 90% prospect numbers; those are only under `post_2008_developments`, not `pdg_or_equivalent`.
- L006: CHK-1 fails under the same rule. The PDG \(G_C/G_F<0.0030\) and \(P_{M\bar M}<8.3\times10^{-11}\) values are in `pdg_or_equivalent`, but `L006.tex` quotes MACE/Snowmass prospective numbers ("more than two orders of magnitude" and \(10^{-13}\)) that live only under `prospects`.
- E002: CHK-1 fails under the same rule. The current PDG/Bennett muon-EDM values are in `pdg_or_equivalent`, but `E002.tex` quotes PSI projection/status numbers \(125\,\mathrm{MeV}/c\), \(3\,\mathrm{T}\), \(1\,\mathrm{GV/m}\), \(6\times10^{-23}\ e\,\mathrm{cm}\), 2026, factor 100, and early 2030s that live only under `auxiliary_values`.
- CHK-2: All process-local source keys cited in the TeX key-reference sections resolve to `source_manifest.yaml` entries. Each manifest entry has a non-empty `snapshot_path`; `git ls-files flavor_catalog/references/{L005,L006,E002}` shows the snapshot files are tracked.
- CHK-3: `ls flavor_catalog/references/{L005,L006,E002}/` showed only `.txt` snapshots, `source_manifest.yaml`, and L005 `sha256sums.txt`; no publisher PDFs are present.
- CHK-4: Each sidecar had `WRITER-INITIATED` before `WRITER-DONE` with ISO 8601 timestamps, and the most recent pre-check state was `WRITER-DONE`. This checker pass appends `WRITER-REWORK`.
- CHK-5: The requested `rg -l -E "<process keyword>" ...` form fails in this environment because ripgrep treats `-E` as an encoding flag; I reran the semantic searches with `-e`. Focused L005 and L006 greps returned no files. Focused E002 greps for `muEDM|muon.?EDM|muon electric dipole`, exact `d_mu`/`dmu`, and `electric dipole|EDM` returned no non-notebook implementation. The adjacent `muToEGamma` line references in the TeX exist at `flavorConstraints/muToEGamma.py:3`, `:21`, `:75`, and `:81`; `scanParams/scan.py:33` and `:524` also exist.
- CHK-6: Difficulty labels are consistent. L005 is HIGH because coherent conversion needs new lepton-quark operators plus target-specific nuclear/capture inputs. L006 is MEDIUM because it needs a new four-lepton/bound-state observable but no lattice or hadronic long-distance inputs. E002 is HIGH because a flavor-diagonal CP-odd lepton dipole needs new matching and likely running.
- CHK-7: No rc1.1 load-bearing number is contradicted. Searching outside `flavor_catalog/` found only catalog plan rows for these process IDs/prospects, not paper claims.
- CHK-8: No `\textbf{CHECK}`, unresolved `\cite`/`\ref`, TODO, missing-label, or obvious math-mode issue was found in the three TeX files.
