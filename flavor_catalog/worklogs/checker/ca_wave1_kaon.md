# CA Worklog: ca_wave1_kaon
**Date**: 2026-05-16
**Family**: kaon
**Process IDs**: K003 K004 K005 K006 K013

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| K003 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| K004 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| K005 | PASS | PASS | PASS | PASS | PASS | FAIL | PASS | PASS | WRITER-REWORK |
| K006 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| K013 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| K003 | Re(epsilon'/epsilon) | 0.00166 +/- 0.00023 | PDG pdgLive S013 row | 230285c9b8cb |
| K003 | PDG DataBlock OUR FIT | 1.67 +/- 0.23 in 10^-3 | PDG DataBlock S013EPS | 283335664646 |
| K003 | PDG DataBlock OUR AVERAGE | 1.68 +/- 0.20 in 10^-3 | PDG DataBlock S013EPS | 283335664646 |
| K003 | KTeV Re(epsilon'/epsilon) | (19.2 +/- 2.1) x 10^-4 | KTeV arXiv:1011.0127 | 8257d37c79b0 |
| K003 | NA48 Re(epsilon'/epsilon) | (14.7 +/- 2.2) x 10^-4 | NA48 arXiv:hep-ex/0208009 | 72baed88144b |
| K003 | RBC/UKQCD Re(epsilon'/epsilon) | 21.7(2.6)(6.2)(5.0) x 10^-4 | RBC/UKQCD arXiv:2004.09440 | ed6a62e309e1 |
| K003 | SM octet/nonet benchmarks | (17.4 +/- 6.1), (13.9 +/- 5.2) x 10^-4 | Aebischer-Bobeth-Buras arXiv:2005.05978 | c133e14c8a13 |
| K004 | BR(K+ -> pi+ nu nubar) | (9.6 +1.9/-1.8) x 10^-11 | NA62 2026 arXiv:2604.12649 | db5484ce8b53 |
| K004 | BR(K+ -> pi+ nu nubar) published anchor | (13.0 +3.3/-3.0) x 10^-11 | NA62 2025 arXiv:2412.12015 | 4bdb16bba0ed |
| K004 | Signal/background counts | 51 candidates, 18 +3/-2 background | NA62 2025 arXiv:2412.12015 | 4bdb16bba0ed |
| K004 | SM BR(K+ -> pi+ nu nubar) | (8.60 +/- 0.42) x 10^-11 | Buras-Venturini arXiv:2203.10099 | a50d05504c22 |
| K004 | SM BR(K+ -> pi+ nu nubar) | 7.73(61) x 10^-11 | Brod-Gorbahn-Stamou arXiv:2105.02868 | 128d73845482 |
| K004 | CFW RS context | about 21 TeV and about 33 TeV | Csaki-Falkowski-Weiler arXiv:0804.1954 | 88e65fabef3b |
| K005 | BR(KL -> pi0 nu nubar) | <2.2 x 10^-9 at 90% CL | PDG pdgLive S013R40 | 977270390246 |
| K005 | KOTO background and SES | 0.252 +/- 0.055 +0.052/-0.067; (9.33 +/- 0.06 +/- 0.84) x 10^-10 | KOTO arXiv:2411.11237 | 97178b5f26d1 |
| K005 | Observed signal-region events | 0 | KOTO arXiv:2411.11237 | 97178b5f26d1 |
| K005 | SM BR(KL -> pi0 nu nubar) | 2.59(29) x 10^-11 | Brod-Gorbahn-Stamou arXiv:2105.02868 | 864ac68f32c5 |
| K005 | SM BR(KL -> pi0 nu nubar) | (2.94 +/- 0.15) x 10^-11 | Buras-Venturini arXiv:2203.10099 | 088ac9757e7a |
| K005 | Related BR(K+ -> pi+ nu nubar) | (9.6 +1.9/-1.8) x 10^-11 | NA62 2026 arXiv:2604.12649 | a9c940a6313f |
| K006 | BR(KL -> mu+ mu-) | (6.84 +/- 0.11) x 10^-9 | PDG 2024 KL listing | 7b8d82405b9c |
| K006 | E871 BR(KL -> mu+ mu-) | (7.18 +/- 0.17) x 10^-9 | Ambrose et al. PubMed snapshot | d38cf3e1320b |
| K006 | E871 normalization ratio | (3.474 +/- 0.057) x 10^-6 | Ambrose et al. PubMed snapshot | d38cf3e1320b |
| K006 | Short-distance bound | <2.5 x 10^-9 | Isidori-Unterdorfer arXiv:hep-ph/0311084 | 20e19ea8bbf0 |
| K006 | Clean-probe hadronic uncertainty | below 1% | Dery-Ghosh-Grossman-Schacht arXiv:2104.06427 | b8fe01acff7a |
| K006 | Lattice target accuracy | 10% | Chao-Christ arXiv:2406.07447 | 0f6e418c75fa |
| K013 | BR(KL -> pi0 gamma gamma) | (1.273 +/- 0.033) x 10^-6 | PDG 2025 KL listing | 48bdda1e34bb |
| K013 | PDG-rescaled KTeV input | 1.28 +/- 0.06 +/- 0.01 in 10^-6 | PDG 2025 KL listing / KTeV | 48bdda1e34bb |
| K013 | PDG-rescaled NA48 input | 1.27 +/- 0.04 +/- 0.01 in 10^-6 | PDG 2025 KL listing / NA48 | 48bdda1e34bb |
| K013 | KTeV original BR and aV | (1.29 +/- 0.03 +/- 0.05) x 10^-6; aV=-0.31 +/- 0.05 +/- 0.07 | KTeV arXiv:0805.0031 | 5f872c5a8480 |
| K013 | NA48 original BR and aV | (1.36 +/- 0.03 +/- 0.03 +/- 0.03) x 10^-6; aV=-0.46 +/- 0.03 +/- 0.04 | NA48 repository snapshot | b8d0f2ce120f |
| K013 | Related charged-mode z cut and BRs | z>0.2; 0.877 +/- 0.089 and 0.910 +/- 0.075 in 10^-6 | NA48/2 arXiv:1310.5499 | 5727e0ee9052 |

## Check Evidence
- CHK-1: Grepped each TeX for numeric claims and matched all load-bearing values to YAML numeric blocks plus local snapshots. Snapshot hashes were recomputed with `sha256sum` and matched sidecar/source-manifest shas. Some supporting K006/K013 numeric blocks carry access provenance through the process source manifest and snapshot text rather than repeating `access_date` on every nested value.
- CHK-2: All TeX local source keys resolve to entries in the process `source_manifest.yaml`; each manifest entry has a non-empty `snapshot_path` under `flavor_catalog/references/<process_id>/`, and `git ls-files` shows those snapshots are tracked.
- CHK-3: `find flavor_catalog/references/{K003,K004,K005,K006,K013} -maxdepth 1 -name '*.pdf'` returned no publisher PDFs.
- CHK-4: Each sidecar had `WRITER-INITIATED` before `WRITER-DONE`, with ISO-like timestamps, and `WRITER-DONE` was the terminal pre-CA state.
- CHK-5: Installed `rg` treats `-E` as an encoding flag, so the literal prompt form errors. Equivalent `rg -l -e` process-keyword sweeps found no direct rare-kaon implementations; broad false positives were notebooks, generic warp `epsilon`, or old `K_LR` Delta-F=2 variable names. Cited lines `quarkConstraints/deltaf2.py:1`, `:209`, `:615`, `:729`, and `quarkConstraints/modern/phenomenology.py:23` exist and are Delta-F=2/mixing-only surfaces.
- CHK-6: K003 HIGH is consistent because it needs Delta-S=1 matching, RG, and hadronic matrix elements. K004 MEDIUM is consistent because it needs a new semileptonic operator and standard clean hadronic input. K006 and K013 HIGH are consistent because they need long-distance two-photon/radiative/chiral treatment. K005 HIGH is not consistent with K004 under the rubric; it should be MEDIUM unless the writer explicitly requires a full RS tower/neutrino-sector matching implementation.
- CHK-7: No load-bearing rc1.1 paper number is silently revised. The CFW contextual 21 TeV / 33 TeV statements are preserved as context only and agree with the local CFW snapshots.
- CHK-8: No `\textbf{CHECK}`, TODO/FIXME, `\cite`, `\ref`, or unresolved-reference markers were found in the batch TeX. Spot inspection found no blocking math-mode or brace errors.

## Issues (if any)
- K005: CHK-6 fails. The TeX/YAML classify implementation difficulty as HIGH, but the rubric and K004 comparator support MEDIUM: a new clean semileptonic `s -> d nu nubar` operator/formula is needed, with standard hadronic inputs and no identified new lattice/RG/long-distance calculation.
