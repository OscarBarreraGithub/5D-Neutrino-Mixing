# CA Worklog: ca_w6_kaon
**Date**: 2026-05-16
**Family**: kaon
**Process IDs**: K012 K018

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| K012 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| K018 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| K012 | BR(K_S^0 -> mu+ mu-) headline limit | <2.1e-10 at 90% CL | PDG 2025 K0S listing | 5bed3ef38dfd |
| K012 | LHCb 2016-2018 standalone BR(K_S^0 -> mu+ mu-) | <2.2e-10 at 90% CL | LHCb 2020 arXiv:2001.10354 | d5a95b6ce20 |
| K012 | LHCb combined BR(K_S^0 -> mu+ mu-) | <2.1e-10 at 90% CL | LHCb 2020 arXiv:2001.10354 | d5a95b6ce20 |
| K012 | LHCb Run-1 BR(K_S^0 -> mu+ mu-) | <0.8(1.0)e-9 at 90%(95%) CL | LHCb 2017 arXiv:1706.00758 | a35a20eef31d |
| K012 | SM-scale BR(K_S^0 -> mu+ mu-) | approximately 5e-12 | Chobanova et al. 2018 arXiv:1711.11030 | 7b25e1bf225f |
| K012 | time-dependent ell=0 hadronic uncertainty | <1% | Dery et al. 2021 arXiv:2104.06427 | ecf9094e5782 |
| K018 | \|V_us\| f_+(0) K_l3 average | 0.21656(35) | PDG 2025 Vud/Vus review | e900e36b74e3 |
| K018 | \|V_us\| f_+(0) K+- e3 | 0.2169(6) | PDG 2025 Table 67.1 | e900e36b74e3 |
| K018 | \|V_us\| f_+(0) K+- mu3 | 0.2168(10) | PDG 2025 Table 67.1 | e900e36b74e3 |
| K018 | \|V_us\| f_+(0) K_L e3 | 0.2162(5) | PDG 2025 Table 67.1 | e900e36b74e3 |
| K018 | \|V_us\| f_+(0) K_L mu3 | 0.2165(6) | PDG 2025 Table 67.1 | e900e36b74e3 |
| K018 | \|V_us\| f_+(0) K_S e3 | 0.2169(8) | PDG 2025 Table 67.1 | e900e36b74e3 |
| K018 | \|V_us\| f_+(0) K_S mu3 | 0.2125(47) | PDG 2025 Table 67.1 | e900e36b74e3 |
| K018 | f_+(0), Nf=2+1+1 | 0.9698(17) | FLAG 2024 arXiv:2411.04268 | 7a1062d3a82b |
| K018 | \|V_us\| derived from K_l3 and FLAG f_+(0) | 0.22330(53) | PDG 2025 Eq. 67.13 | e900e36b74e3 |
| K018 | FlaviaNet K_l3 average | 0.2163(5) | FlaviaNet 2010 arXiv:1005.2323 | c8855c1e4539 |
| K018 | r_mu e = K_mu3/K_e3 | 1.002(5) | FlaviaNet 2010 arXiv:1005.2323 | c8855c1e4539 |
| K018 | \|V_us\| derived from FLAG f_+(0) | 0.22328(58) | FLAG 2024 arXiv:2411.04268 | 7a1062d3a82b |

## Check evidence
- K012 CHK-1: the headline PDG limit is in `pdg_or_equivalent`, but the TeX also displays LHCb standalone/Run-1 limits, the SM-scale estimate, and the sub-1% hadronic-uncertainty number. Those trace to `supporting_numeric_values`, not `pdg_or_equivalent`, so this fails the batch prompt's strict `pdg_or_equivalent` requirement.
- K018 CHK-1: the PDG K_l3 average, mode values, FlaviaNet average, and `r_mu e` are in `pdg_or_equivalent`. The TeX also displays `f_+(0)=0.9698(17)`, `|V_us|=0.22330(53)`, and `|V_us|=0.22328(58)`, which trace to `auxiliary_inputs`, not `pdg_or_equivalent`; the `PDG2025:K018:Vus_from_Kl3` auxiliary block lacks `source_url` and `access_date`.
- CHK-2: all key-reference entries resolve in the process-local source manifests: K012 keys `PDG2025:K012`, `LHCb2020:KSmumu`, `LHCb2017:KSmumu`, `Chobanova2018:KSmumuSUSY`, and `DeryGhoshGrossmanSchacht2021:KMuMuClean`; K018 keys `PDG2025:K018Kl3`, `FLAG2024:K018KaonSemileptonic`, `FlaviaNet:K018Kaon2010`, and `CsakiFalkowskiWeiler:RSFlavor2008`. Each manifest entry has a non-empty tracked `snapshot_path` under `flavor_catalog/references/<process_id>/`.
- CHK-3: `ls flavor_catalog/references/K012/` and `ls flavor_catalog/references/K018/` show only `.txt`, `source_manifest.yaml`, and K018 `sha256sums.txt`; no `.pdf` files are tracked.
- CHK-4: both YAML sidecars contain `WRITER-INITIATED` before `WRITER-DONE` with ISO 8601 timestamps; before this CA update the most recent entry was `WRITER-DONE`.
- CHK-5: K012 `NO` coverage matches the exact no-hit grep for `K0S|K0_S|K_S^0|mumu|dimuon|mu+ mu|s->d.*mu`; the broad kaon grep only finds Delta F=2 mixing, including `quarkConstraints/deltaf2.py:615`. K018 `PARTIAL` coverage is supported by existing `Vus` targets at `quarkConstraints/scan.py:43`, `quarkConstraints/scan.py:247`, and `tests/test_quark_fit.py:587`, while the K_l3-specific grep has no source-level K_l3 likelihood hits.
- CHK-6: K012 and K018 `HIGH` difficulty is consistent with the rubric because both require a new mode calculation or likelihood beyond the existing Delta F=2 SLL/SLR/VLL/VRR/LR basis.
- CHK-7: neither entry silently revises load-bearing rc1.1 paper numbers; these are companion-catalog entries and K018 uses the 2008 RS reference only as context.
- CHK-8: no `\textbf{CHECK}`, unresolved `\ref{...}`, or unresolved `\cite{...}` markers found; visual LaTeX inspection found no blocking math-mode or brace issues.

## Issues (if any)
- K012: Move the displayed supporting numerical values into `pdg_or_equivalent` entries or remove them from TeX; currently they live only in `supporting_numeric_values`.
- K018: Move the displayed FLAG/derived auxiliary numerical values into `pdg_or_equivalent` entries or remove them from TeX. Also add `source_url` and `access_date` to `PDG2025:K018:Vus_from_Kl3` if that displayed value remains.
