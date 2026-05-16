# CA Worklog: ca_w6_kaon_v2
**Date**: 2026-05-16
**Family**: kaon
**Process IDs**: K012 K018

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| K012 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| K018 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| K012 | BR(K_S^0 -> mu+ mu-) headline limit | <2.1e-10, 90% CL | PDG 2025 K0S listing | 5bed3ef38dfd |
| K012 | LHCb 2016-2018 standalone BR(K_S^0 -> mu+ mu-) | <2.2e-10, 90% CL | LHCb 2020 arXiv:2001.10354 | d5a95b6ce20e |
| K012 | LHCb combined BR(K_S^0 -> mu+ mu-) | <2.1e-10, 90% CL | LHCb 2020 arXiv:2001.10354 | d5a95b6ce20e |
| K012 | LHCb Run-1 BR(K_S^0 -> mu+ mu-) | <0.8e-9 at 90%; <1.0e-9 at 95% | LHCb 2017 arXiv:1706.00758 | a35a20eef31d |
| K012 | SM estimate for BR(K_S^0 -> mu+ mu-) | approximately 5e-12 | Chobanova et al. 2018 arXiv:1711.11030 | 7b25e1bf225f |
| K012 | Time-dependent ell=0 hadronic uncertainty | <1% | Dery, Ghosh, Grossman, Schacht 2021 arXiv:2104.06427 | ecf9094e5782 |
| K018 | \|V_us\| f_+(0) K_l3 average | 0.21656 +/- 0.00035 | PDG 2025 Vud/Vus review | e900e36b74e3 |
| K018 | \|V_us\| f_+(0) from K+- e3 | 0.2169 +/- 0.0006 | PDG 2025 Table 67.1 | e900e36b74e3 |
| K018 | \|V_us\| f_+(0) from K+- mu3 | 0.2168 +/- 0.0010 | PDG 2025 Table 67.1 | e900e36b74e3 |
| K018 | \|V_us\| f_+(0) from K_L e3 | 0.2162 +/- 0.0005 | PDG 2025 Table 67.1 | e900e36b74e3 |
| K018 | \|V_us\| f_+(0) from K_L mu3 | 0.2165 +/- 0.0006 | PDG 2025 Table 67.1 | e900e36b74e3 |
| K018 | \|V_us\| f_+(0) from K_S e3 | 0.2169 +/- 0.0008 | PDG 2025 Table 67.1 | e900e36b74e3 |
| K018 | \|V_us\| f_+(0) from K_S mu3 | 0.2125 +/- 0.0047 | PDG 2025 Table 67.1 | e900e36b74e3 |
| K018 | FlaviaNet \|V_us\| f_+(0) average | 0.2163 +/- 0.0005 | FlaviaNet 2010 arXiv:1005.2323 | c8855c1e4539 |
| K018 | K_mu3 / K_e3 universality ratio r_mu e | 1.002 +/- 0.005 | FlaviaNet 2010 arXiv:1005.2323 | c8855c1e4539 |
| K018 | \|V_us\| from K_l3 and FLAG f_+(0) | 0.22330 +/- 0.00053 | PDG 2025 Eq. 67.13 | e900e36b74e3 |
| K018 | f_+(0), direct lattice average, Nf=2+1+1 | 0.9698 +/- 0.0017 | FLAG 2024 arXiv:2411.04268 | 7a1062d3a82b |
| K018 | \|V_us\| derived from f_+(0), Nf=2+1+1 | 0.22328 +/- 0.00058 | FLAG 2024 arXiv:2411.04268 | 7a1062d3a82b |

## Check evidence notes
- K012 CHK-1: numerical grep found the displayed limits <2.1e-10, <2.2e-10, <0.8(1.0)e-9, 5e-12, and <1%; each now has a `pdg_or_equivalent` entry with year, value, uncertainty or null uncertainty, source URL, access date, snapshot path, and sha256. Local snapshots contain the same values and `sha256sum` matches the sidecar/source manifest.
- K018 CHK-1: numerical grep found the displayed K_l3 products, FlaviaNet values, FLAG f_+(0), and derived |V_us| values; each has a `pdg_or_equivalent` entry with complete metadata and matching local snapshot hash.
- CHK-2: all key-reference items in the TeX resolve to process-local manifest entries. Manifest snapshot paths are non-empty and tracked under `flavor_catalog/references/K012/` or `flavor_catalog/references/K018/`.
- CHK-3: `ls` and `find ... -name '*.pdf'` showed no tracked PDFs under either process reference directory.
- CHK-4: both sidecars contain `WRITER-INITIATED` followed by cycle-2 `WRITER-DONE` with ISO 8601 timestamps; before this CA update, the most recent status was `WRITER-DONE`.
- CHK-5: K012 `NO` is consistent. The K012 process grep has no source-code implementation hits when notebooks are excluded; the only unfiltered hit was base64 notebook payload in `warpConfig/wavefuncsTest.ipynb`. The cited nearest line `quarkConstraints/deltaf2.py:615` exists and is Delta F=2 mixing, not this decay. K018 `PARTIAL` is consistent: the cited generic `Vus` targets exist in `quarkConstraints/scan.py:43`, `quarkConstraints/scan.py:247`, and `tests/test_quark_fit.py:587`; K_l3/f_+(0)/semileptonic-kaon specific grep returned no implementation hits, and `quarkConstraints/modern/phenomenology.py:23` lists only the neutral-meson systems.
- CHK-6: `HIGH` is consistent for both. K012 needs a new rare-kaon Delta S=1 dimuon mode treatment with long-distance/interference handling. K018 needs a charged-current semileptonic likelihood/EFT wrapper plus radiative/isospin and lattice-provenance bookkeeping; neither is covered by the existing Delta F=2 SLL/SLR/VLL/VRR/LR basis.
- CHK-7: repo-wide non-catalog grep for K012/K018 process labels and displayed headline values found no rc1.1 paper text to contradict.
- CHK-8: no `\textbf{CHECK}`, unresolved `\cite{}`, unresolved `\ref{}`, `TODO`, `FIXME`, or `??` markers were found; inspected math/texttt usage is syntactically acceptable.

## Issues (if any)
- None.
