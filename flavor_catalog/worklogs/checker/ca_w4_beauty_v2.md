# CA Worklog: ca_w4_beauty_v2
**Date**: 2026-05-16
**Family**: beauty
**Process IDs**: B006

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| B006 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| B006 | PDG canonical BR(B0 -> mu+ mu-) | <1.5e-10 at 90% CL | PDG live/API S042.7, accessed 2026-05-16, `pdg_2026_bdmumu.txt` | 45ef5f9721c6 |
| B006 | CMS input BR(B0 -> mu+ mu-) | <1.5e-10 at 90% CL; <1.9e-10 at 95% CL | CMS PLB 842 (2023) 137955 / arXiv:2212.10311, accessed 2026-05-16, `cms_2212_10311_bdmumu.txt` | 3c83d78acb15 |
| B006 | CMS Run 2 dataset | 140 fb^-1 at sqrt(s)=13 TeV | CMS PLB 842 (2023) 137955 / arXiv:2212.10311, accessed 2026-05-16, `cms_2212_10311_bdmumu.txt` | 3c83d78acb15 |
| B006 | HFLAV BR(B0 -> mu+ mu-) / BR(B_s0 -> mu+ mu-) | <8.1e-2 at 90% CL | HFLAV Apr. 2023 rare-decay table, accessed 2026-05-16, `hflav_2023_bd_over_bs_mumu.txt` | 6f4876eae175 |
| B006 | ATLAS combined BR(B0 -> mu+ mu-) | <2.1e-10 at 95% CL | ATLAS JHEP 04 (2019) 098 / arXiv:1812.03017, accessed 2026-05-16, `atlas_1812_03017_bdmumu.txt` | 4df38ad17985 |
| B006 | LHCb BR(B0 -> mu+ mu-) | <2.6e-10 at 95% CL | LHCb PRL 128 (2022) 041801 / arXiv:2108.09283, accessed 2026-05-16, `lhcb_2108_09283_bdmumu.txt` | 7b5442cafe6e |
| B006 | SM BR(B_d -> mu+ mu-) | (1.06 +/- 0.09)e-10 | Bobeth et al. PRL 112 (2014) 101801 / arXiv:1311.0903, accessed 2026-05-16, `bobeth_1311_0903_bsdll_sm.txt` | 77747a869720 |

## Issues (if any)
- None.

## Verification notes
- CHK-1: grep of `B006.tex` found the numerical physics/source claims: PDG S042.7 accessed 2026-05-16, <1.5e-10 at 90% CL, CMS Run 2 <1.9e-10 at 95% CL, HFLAV Apr. 2023 <8.1e-2 at 90% CL, Bobeth SM (1.06 +/- 0.09)e-10, ATLAS <2.1e-10 at 95% CL, LHCb <2.6e-10 at 95% CL, and CMS dataset 140 fb^-1 at sqrt(s)=13 TeV. The cycle-2 YAML now includes the CMS dataset block under `pdg_or_equivalent.cms_run2_dataset` with year, value, uncertainty null, source URL, access date, and sha256 for both dataset numerics. The remaining value claims are represented under `pdg_or_equivalent` and agree with the source snapshots.
- CHK-2: every TeX source key (`PDG2026:BdMuMu`, `HFLAV2023:BdOverBsMuMu`, `CMS2023:BdMuMu`, `LHCb2022:BdMuMu`, `ATLAS2019:BdMuMu`, `BobethEtAl2013:BdMuMuSM`, `CsakiFalkowskiWeiler2008:CompositeFlavor`) resolves to `references/B006/source_manifest.yaml`. Each manifest entry has a non-empty `snapshot_path`; `git ls-files flavor_catalog/references/B006` confirms the snapshots are tracked.
- CHK-3: `ls flavor_catalog/references/B006/` shows only `.txt` snapshots plus `source_manifest.yaml`; no publisher PDFs are tracked.
- CHK-4: before this CA mutation, `B006.yaml` contained `WRITER-INITIATED` followed by `WRITER-DONE`, `WRITER-REWORK`, and the cycle-2 `WRITER-DONE`, with ISO-8601 timestamps. The most recent pre-CA entry was `WRITER-DONE`.
- CHK-5: this ripgrep build treats `-E` as an encoding flag, so the literal prompt form `rg -l -E "<process keyword>" ...` errors. I reran the equivalent `rg -l -e` searches over the required directories. Broad process-keyword searches hit notebooks/base64 and generic `B_d/B_s` Delta-F=2 mixing code only; focused non-notebook searches for `mu+ mu-`, `B0 -> mu`, `dimuon`, and related B006 strings returned no rare-leptonic implementation. The TeX `NO` coverage claim is consistent, and the cited `quarkConstraints/deltaf2.py:225`, `quarkConstraints/deltaf2.py:903`, and `quarkConstraints/modern/phenomenology.py:646` lines exist as B_d mixing, not B0 -> mu+ mu-.
- CHK-6: `MEDIUM` is consistent with the rubric. B006 needs a new Delta-B=1 leptonic rare-decay observable and C10/scalar/pseudoscalar normalization, but uses standard decay-constant inputs and does not require a full angular likelihood or new Delta-F=2 operator basis.
- CHK-7: no rc1.1 contradiction found. This is a companion-mode catalog entry; the local snapshots and TeX do not revise any load-bearing rc1.1 quark-scan number.
- CHK-8: no `\textbf{CHECK}`, `TODO`, `FIXME`, `\ref{}`, `\cite{}`, or `??` tokens were found in `B006.tex`; brace balance delta is zero.
