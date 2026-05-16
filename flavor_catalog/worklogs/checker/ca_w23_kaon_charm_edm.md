# CA Worklog: ca_w23_kaon_charm_edm
**Date**: 2026-05-16
**Family**: mixed
**Process IDs**: K001 K002 C001 E001

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| K001 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| K002 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| C001 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| E001 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| K001 | \(|\varepsilon|\) | \((2.228 \pm 0.011)\times10^{-3}\) | PDG 2026 kaon CP-violation review | 9ada2a952ab9 |
| K001 | \(\phi_\varepsilon\) | \((43.5 \pm 0.5)^\circ\) | PDG 2026 kaon CP-violation review | 9ada2a952ab9 |
| K001 | \(|\varepsilon_K|_{\rm SM}\) | \((2.161 \pm0.153 \pm0.076 \pm0.065)\times10^{-3}=2.16(18)\times10^{-3}\) | Brod-Gorbahn-Stamou 2020 | 56faa9c2fed8 |
| K001 | \(\hat B_K\), \(B_K^{\overline{\rm MS}}(2\,{\rm GeV})\) | \(0.7533(91)\), \(0.5503(66)\) | FLAG Review 2024 | 513caa2edf3d |
| K001 | CFW RS context | about 21 TeV and about 33 TeV | Csaki-Falkowski-Weiler 2008 | e6330846f148 |
| K002 | \(\Delta m_K\), CPT-assuming fit | \((0.5293\pm0.0009)\times10^{10}\hbar\,{\rm s}^{-1}\), scale factor 1.3 | PDG pdgLive S013D 2026 | b19ce40faf7d |
| K002 | \(\Delta m_K\), no-CPT companion fit | \((0.5289\pm0.0010)\times10^{10}\hbar\,{\rm s}^{-1}\) | PDG pdgLive S013D 2026 | b19ce40faf7d |
| K002 | exploratory lattice \(\Delta M_K\) | \(3.19(41)(96)\times10^{-12}\) MeV | Bai et al. 2014 | 94522a76f8bf |
| K002 | physical-mass lattice statistical uncertainty | near 9% | Wang 2023 | 26d8bab47c57 |
| C001 | \(x_D\) | \((0.405\pm0.043)\%\), 95% C.L. [0.320, 0.489]% | HFLAV CKM25 all-CPV fit | 85254267f434 |
| C001 | \(y_D\) | \((0.636\pm0.024)\%\), 95% C.L. [0.590, 0.682]% | HFLAV CKM25 all-CPV fit | 85254267f434 |
| C001 | \(\Delta m_D\) | \((0.997\pm0.116)\times10^{10}\hbar\,{\rm s}^{-1}=(6.56\pm0.76)\times10^{-15}\) GeV | PDG Live S032D / HFLAV | 398d3420f232 |
| C001 | \(\Delta y\) | \((0.031\pm0.035\pm0.013)\%\) | HFLAV CKM25 input table | 024569bebcc2 |
| C001 | \(\Delta\Gamma_D/\Gamma_D\) | \(2y_D=(1.272\pm0.048)\%\) | HFLAV CKM25 all-CPV fit | 85254267f434 |
| C001 | CFW RS context | about 21 TeV and about 33 TeV | Csaki-Falkowski-Weiler 2008 | 224e28f70091 |
| E001 | \(|d_e|\) current limit | \(<4.1\times10^{-30}\ e\,\mathrm{cm}\) at 90% CL | PDG Live S003EDM 2026 | aab24570becf |
| E001 | Roussy et al. measurement | \((-1.3\pm2.0_{\rm stat}\pm0.6_{\rm syst})\times10^{-30}\ e\,\mathrm{cm}\) | Roussy et al. 2023 | 85dfc60b010c |
| E001 | ACME/Andreev 2018 limit and note | \(<1.1\times10^{-29}\ e\,\mathrm{cm}\) at 90% CL; \((4.3\pm3.1\pm2.6)\times10^{-30}\ e\,\mathrm{cm}\) | ACME / Andreev et al. 2018 | c81e0a1936ec |

## Check Evidence
- CHK-1: Grepped each TeX for numeric claims and recomputed snapshot hashes with `sha256sum`. K001 passes: all load-bearing numerical values are represented in `pdg_or_equivalent` blocks with year/provenance and matching local snapshot hashes. K002, C001, and E001 fail the strict placement/structured-metadata rule listed below.
- CHK-2: All TeX source keys resolve to entries in `flavor_catalog/references/<process_id>/source_manifest.yaml`. Each manifest entry has a non-empty `snapshot_path`; `git ls-files flavor_catalog/references/{K001,K002,C001,E001}/` shows all referenced snapshots are tracked.
- CHK-3: `ls flavor_catalog/references/{K001,K002,C001,E001}/` showed only `.txt` snapshots, manifests, and C001 `sha256sums.txt`; `find ... -iname '*.pdf'` returned no PDF files.
- CHK-4: The pre-CA sidecars all contained `WRITER-INITIATED` before `WRITER-DONE` with ISO 8601 timestamps, and `WRITER-DONE` was the most recent pre-check state.
- CHK-5: The literal prompt command form `rg -l -E "<pattern>" ...` fails with this installed ripgrep because `-E` is parsed as an encoding flag. Equivalent `rg -l -e` sweeps confirmed the coverage claims: K001 hits `epsilon_K/epsilon_k` in `quarkConstraints/deltaf2.py`, modern policy files, and tests; K002 hits `DELTA_M_K/evaluate_delta_mk`; C001 hits `D0/d_mix/evaluate_d0_mixing`; the focused E001 electron-EDM grep returned no hits in the required implementation/test directories. Cited file:line locations exist.
- CHK-6: The implementation-difficulty labels are consistent. K001, K002, and C001 use the existing \(\Delta F=2\) neutral-meson operator basis and are LOW for catalog/live conservative bounds. E001 needs a new flavor-diagonal CP-odd lepton dipole observable plus matching/running and is HIGH.
- CHK-7: No load-bearing rc1.1 paper number is contradicted. The CFW 21/33 TeV values agree with local CFW snapshots and existing audit text; current HFLAV/PDG values are companion-catalog inputs, not silent paper revisions.
- CHK-8: No `\textbf{CHECK}`, TODO/FIXME, `\cite`, `\ref`, unresolved-marker hits, or brace imbalances were found in the four TeX files. `git diff --check` over the TeX paths was clean.

## Issues (if any)
- K002: CHK-1 fails under the strict instruction that every numerical claim in the TeX trace to a `pdg_or_equivalent` entry with year, value, uncertainty/CL where applicable, source URL, access date, snapshot path, and sha256. The headline PDG CPT-assuming value satisfies this, but the no-CPT companion fit, Bai 2014 \(3.19(41)(96)\times10^{-12}\) MeV claim, and Wang 2023 near-9% statistical-uncertainty claim are in `supporting_values`, not `pdg_or_equivalent`.
- C001: CHK-1 fails because `C001.tex` quotes CFW contextual KK-gluon scales of about 21 TeV and 33 TeV, and the snapshot verifies them, but `C001.yaml` has no corresponding `pdg_or_equivalent` numerical block for those two values.
- E001: CHK-1 fails because the `pdg_or_equivalent` blocks for the PDG limit, Roussy 2023 measurement, and ACME 2018 benchmark lack explicit `year` fields; the Roussy/ACME numerical values and uncertainties are also stored only in summary strings rather than structured value/uncertainty/CL fields.
