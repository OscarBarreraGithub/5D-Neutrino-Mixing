# CA Worklog: ca_w4_kaon_edm_v2
**Date**: 2026-05-16
**Family**: kaon_edm
**Process IDs**: K017 E006

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| K017 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| E006 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| K017 | \(R_K\) PDG average | \((2.488\pm0.009)\times10^{-5}\) | PDG API S010R28 / PDG2025:K017RK | c00ed1ff90dd |
| K017 | \(R_K\) NA62 input | \((2.488\pm0.007_{\rm stat}\pm0.007_{\rm syst})\times10^{-5}\) | NA62 2013 | b6389e671653 |
| K017 | \(R_K\) KLOE input | \((2.493\pm0.025_{\rm stat}\pm0.019_{\rm syst})\times10^{-5}\) | KLOE 2009 | c473017a6d6d |
| K017 | \(R_K^{\rm SM}\) | \((2.477\pm0.001)\times10^{-5}\) | Cirigliano-Rosell 2007 | 9178029253ec |
| E006 | \(|d_{\mathrm{Hg}}|\) direct limit | \(<7.4\times10^{-30}\ e\,\mathrm{cm}\) at 95% CL | Graner et al. 2016 | 44522bf9473b |
| E006 | Graner Hg-bound improvement | factor of \(4\) | Graner et al. 2016 | 44522bf9473b |
| E006 | PDG Hg-derived neutron rows | Graner \(<0.16\), Sahoo \(<0.22\) in \(10^{-25}\ e\,\mathrm{cm}\), both 95% CL | PDG Live 2026 S017EDM cross-reference | 15f7baddbe77 |
| E006 | Sahoo Hg-derived translation limits | \(d_n<2.2\times10^{-26}\), \(d_p<2.1\times10^{-25}\), \(|\bar\theta|<1.1\times10^{-10}\), \(|\tilde d_u-\tilde d_d|<5.5\times10^{-27}\) | Sahoo 2017 | 88e10fc1c214 |

## Check Evidence
- CHK-1: `rg -n '[0-9]'` over both TeX files found the physical numerical claims listed above. K017 now has value-bearing `pdg_or_equivalent` blocks with `year`, value, uncertainty, source URL, access date, snapshot path, and sha256 for the PDG average, NA62 input, KLOE input, and SM prediction. E006 has strict metadata for the direct Hg limit, the factor-4 Graner improvement, the PDG Hg-derived neutron rows, and the Sahoo interpretation limits.
- CHK-2: Every process-local key cited in the TeX resolves in `flavor_catalog/references/<process_id>/source_manifest.yaml`. Each manifest entry has a non-empty `snapshot_path`; `git ls-files --error-unmatch` confirmed the snapshot files are tracked under the process reference directory.
- CHK-3: `ls flavor_catalog/references/K017/` and `ls flavor_catalog/references/E006/` show only `.txt` snapshots plus `source_manifest.yaml`; `find ... -iname '*.pdf'` returned no PDFs.
- CHK-4: Before the CA append, both sidecars had `WRITER-INITIATED` followed by `WRITER-DONE` with ISO 8601 timestamps, and the latest pre-check entry was the cycle-2 `WRITER-DONE`.
- CHK-5: The installed `rg` treats `-E` as an encoding flag, so the exact requested `rg -l -E "<process keyword>" ...` form fails with `unknown encoding`; I ran the semantic equivalent with `-e`. K017 has no focused `K+ e nu`, `K+ mu nu`, `Kell2`, `K_l2`, `Kmu2`, or `Ke2` matches; raw `R_K` matches are unrelated Python names such as `VAR_KEYWORD` and `M12_K_NP`. The cited near-miss lines `quarkConstraints/deltaf2.py:209` and `quarkConstraints/modern/phenomenology.py:23` exist. E006 focused Hg/Schiff/diamagnetic/EDM greps found no Hg implementation; only the unrelated `mu->e gamma` dipole helper appears at `flavorConstraints/muToEGamma.py:3`, `:21`, and `:81`, matching the TeX.
- CHK-6: The difficulty labels are consistent with the rubric. K017 is `MEDIUM` because it needs a new charged-current leptonic observable but no new lattice inputs. E006 is `HIGH` because a live constraint would need new CP-odd matching, RG/threshold treatment, and atomic/nuclear translation.
- CHK-7: No load-bearing rc1.1 number is contradicted. The TeX entries cite CFW 2008 only as baseline context and do not quote or revise the rc1.1 21/33 TeV comparison numbers.
- CHK-8: No `\textbf{CHECK}`, TODO/FIXME, `\cite`, `\ref`, unresolved-marker hits, or simple brace imbalance was found in either TeX file.

## Issues (if any)
- None.
