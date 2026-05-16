# CA Worklog: ca_wave1_kaon_v2
**Date**: 2026-05-16
**Family**: kaon
**Process IDs**: K005

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| K005 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| K005 | BR(KL -> pi0 nu nubar) | <2.2 x 10^-9 at 90% CL | PDG pdgLive S013R40 | 977270390246 |
| K005 | KOTO background and SES | 0.252 +/- 0.055 +0.052/-0.067; (9.33 +/- 0.06 +/- 0.84) x 10^-10 | KOTO arXiv:2411.11237 / PRL 134 081802 | 97178b5f26d1 |
| K005 | Observed signal-region events | 0 | KOTO arXiv:2411.11237 / PRL 134 081802 | 97178b5f26d1 |
| K005 | SM BR(KL -> pi0 nu nubar) | 2.59(29) x 10^-11 | Brod-Gorbahn-Stamou arXiv:2105.02868 | 864ac68f32c |
| K005 | SM BR(KL -> pi0 nu nubar) | (2.94 +/- 0.15) x 10^-11 | Buras-Venturini arXiv:2203.10099 | 088ac9757e7a |
| K005 | Related BR(K+ -> pi+ nu nubar) | (9.6 +1.9/-1.8) x 10^-11 | NA62 2026 arXiv:2604.12649 | a9c940a6313 |

## Check evidence
- CHK-1: Grepped `K005.tex` for numeric claims. Load-bearing values match K005 sidecar numeric blocks (`pdg_or_equivalent`, `experimental_inputs`, and `theory_inputs`) with year, value or supporting number, uncertainty/CL where applicable, source URL, access date, snapshot path, and sha256. Recomputed `sha256sum flavor_catalog/references/K005/*.txt`; all six hashes matched the sidecar and source manifest.
- CHK-2: The TeX keys `PDG2026:K005`, `KOTO2025:KLPi0NuNu`, `BrodGorbahnStamou2021:KPiNuNu`, `BurasVenturini2022:RareK`, `NA622026:KPlusPiNuNu`, and `BuchallaBuras1996:KPiNuNu` all resolve in `flavor_catalog/references/K005/source_manifest.yaml`. Each manifest entry has a non-empty tracked `.txt` snapshot under `flavor_catalog/references/K005/`.
- CHK-3: `ls flavor_catalog/references/K005/` shows only `.txt` snapshots plus `source_manifest.yaml`; `find flavor_catalog/references/K005 -maxdepth 1 -type f -name '*.pdf' -print` returned no PDFs.
- CHK-4: `status_history` contains `WRITER-INITIATED` before `WRITER-DONE`; the terminal pre-CA entry was the WA-v2 `WRITER-DONE` at `2026-05-16T11:54:07-04:00`, with ISO 8601 timestamps.
- CHK-5: The literal requested `rg -l -E` form errors on this installed `rg` because `-E` is the encoding flag. Equivalent `rg -l -e` and `rg -n --glob '!*.ipynb' -e` sweeps found no KOTO, NA62, `nu nubar`, `nunu`, `pi0 nu`, `pi^0 nu`, rare-kaon, or `s->d` rare-decay implementation. The broad `K_L` hit `quarkConstraints/deltaf2.py:615` exists and is only the `K_L - K_S` mass-difference input; other hits are Delta-F=2 `K_LR` identifiers or one `drift` false positive in a test comment.
- CHK-6: WA-v2 changed K005 to `MEDIUM`. This matches the rubric: K005 needs a new Delta-S=1 semileptonic `s -> d nu nubar` observable/operator path and CP-odd `K_L` projection, while using standard clean hadronic inputs and no identified new lattice, RG, or long-distance calculation.
- CHK-7: No rc1.1 load-bearing number is revised. A non-catalog search outside `flavor_catalog/` found only seed/planning mentions, generic code strings, and unrelated numeric false positives after excluding notebooks.
- CHK-8: No `\textbf{CHECK}`, TODO/FIXME, `\cite`, `\ref`, unresolved marker, or blocking math/bracing issue was found in `K005.tex`.

## Issues (if any)
- None.
