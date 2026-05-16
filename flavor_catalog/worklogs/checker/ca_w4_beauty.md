# CA Worklog: ca_w4_beauty
**Date**: 2026-05-16
**Family**: beauty
**Process IDs**: B004 B006 B019 B026

## Per-process verdict
| process_id | CHK-1 | CHK-2 | CHK-3 | CHK-4 | CHK-5 | CHK-6 | CHK-7 | CHK-8 | overall |
|---|---|---|---|---|---|---|---|---|---|
| B004 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| B006 | FAIL | PASS | PASS | PASS | PASS | PASS | PASS | PASS | WRITER-REWORK |
| B019 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |
| B026 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | CHECKER-DONE |

## Verified source/value table
| process_id | observable | value | source | sha256 (first 12) |
|---|---|---|---|---|
| B004 | HFLAV all-combined \(\phi_s^{c\bar c s}\) | \(-0.041\pm0.016\) rad; \(\Delta\Gamma_s=0.079\pm0.006~{\rm ps}^{-1}\) | HFLAV PDG 2025, accessed 2026-05-16, `hflav_pdg2025_phis_inputs.txt` | ae8171d0b32e |
| B004 | HFLAV \(B_s^0\to J/\psi K^+K^-\) combination | \(-0.050\pm0.017\) rad; \(\Delta\Gamma_s=0.077\pm0.006~{\rm ps}^{-1}\) | HFLAV PDG 2025, accessed 2026-05-16, `hflav_pdg2025_phis_inputs.txt` | ae8171d0b32e |
| B004 | LHCb \(B_s^0\to J/\psi K^+K^-\) input | \(-0.039\pm0.022_{\rm stat}\pm0.006_{\rm syst}\) rad, 6 fb\(^{-1}\), ~349000 signal decays | LHCb PRL 132 (2024) 051802 / arXiv:2308.01468, accessed 2026-05-16, `lhcb_2024_jpsikk_phis_arxiv.txt` | d94a1b07a02b |
| B004 | SM \(-2\beta_s\) comparison | \(-0.0368^{+0.0009}_{-0.0006}\) rad | LHCb PRL 132 (2024) 051802 / arXiv:2308.01468, accessed 2026-05-16, `lhcb_2024_jpsikk_phis_arxiv.txt` | d94a1b07a02b |
| B006 | PDG canonical \(\mathcal{B}(B^0\to\mu^+\mu^-)\) | \(<1.5\times10^{-10}\) at 90% CL | PDG live/API S042.7, accessed 2026-05-16, `pdg_2026_bdmumu.txt` | 45ef5f9721c6 |
| B006 | CMS input \(\mathcal{B}(B^0\to\mu^+\mu^-)\) | \(<1.5\times10^{-10}\) at 90% CL; \(<1.9\times10^{-10}\) at 95% CL | CMS PLB 842 (2023) 137955 / arXiv:2212.10311, accessed 2026-05-16, `cms_2212_10311_bdmumu.txt` | 3c83d78acb15 |
| B006 | HFLAV \(\mathcal{B}(B^0\to\mu^+\mu^-)/\mathcal{B}(B_s^0\to\mu^+\mu^-)\) | \(<8.1\times10^{-2}\) at 90% CL | HFLAV Apr. 2023 rare-decay table, accessed 2026-05-16, `hflav_2023_bd_over_bs_mumu.txt` | 6f4876eae175 |
| B006 | SM \(\mathcal{B}(B_d\to\mu^+\mu^-)\) | \((1.06\pm0.09)\times10^{-10}\) | Bobeth et al. PRL 112 (2014) 101801 / arXiv:1311.0903, accessed 2026-05-16, `bobeth_1311_0903_bsdll_sm.txt` | 77747a869720 |
| B019 | HFLAV low-\(q^2\) \(R_{K^*}\) | \(0.927\pm0.097\), \(0.1<q^2<1.1~{\rm GeV}^2/c^4\) | HFLAV Dec. 2025, accessed 2026-05-16, `hflav_dec2025_rkstar_lowq2_lhcb.txt` | 7fbe070228eb |
| B019 | HFLAV central-\(q^2\) \(R_{K^*}\) | \(1.028\pm0.074\), \(1.1<q^2<6.0~{\rm GeV}^2/c^4\) | HFLAV Dec. 2025, accessed 2026-05-16, `hflav_dec2025_rkstar_centralq2.txt` | 585fdc5f096d |
| B019 | LHCb 2017 low/central \(R_{K^*}\) | \(0.66^{+0.11}_{-0.07}\pm0.03\); \(0.69^{+0.11}_{-0.07}\pm0.05\) | LHCb JHEP 08 (2017) 055 / arXiv:1705.05802, accessed 2026-05-16, `lhcb_2017_rkstar_arxiv1705_05802.txt` | 0a77136a5169 |
| B019 | LHCb 2023 combined LFU SM compatibility | \(0.2\sigma\) | LHCb PRD 108 (2023) 032002 / arXiv:2212.09153, accessed 2026-05-16, `lhcb_2023_rk_rkstar_arxiv2212_09153.txt` | c8bba3b0498f |
| B019 | BIP central-bin SM/QED reference | \(1.00\pm0.01_{\rm QED}\) | Bordone--Isidori--Pattori EPJC 76 (2016) 440 / arXiv:1605.07633, accessed 2026-05-16, `bordone_isidori_pattori_2016_rkstar_sm_arxiv1605_07633.txt` | 9cbbd94af6b2 |
| B026 | HFLAV CKM 2025 \(R_{D^*}\) | \(0.281\pm0.011\) | HFLAV CKM 2025, updated 2025-09-28, accessed 2026-05-16, `hflav_ckm2025_rdrds.txt` | 4da37c4b6e64 |
| B026 | HFLAV joint \(R_D\)-\(R_{D^*}\) context | \(R_D=0.358\pm0.024\), \(\rho=-0.374\), \(\chi^2/{\rm dof}=16.683/14\) | HFLAV CKM 2025, accessed 2026-05-16, `hflav_ckm2025_rdrds.txt` | 4da37c4b6e64 |
| B026 | HFLAV displayed SM reference and pulls | \(R_D=0.296\pm0.004\), \(R_{D^*}=0.254\pm0.005\), combined \(3.8\sigma\), \(p=1.48\times10^{-4}\), \(R_{D^*}\)-only \(2.3\sigma\) | HFLAV CKM 2025, accessed 2026-05-16, `hflav_ckm2025_rdrds.txt` | 4da37c4b6e64 |

## Issues (if any)
- B006: CHK-1 failed. `B006.tex` states that the CMS Run 2 input used \(140\,{\rm fb}^{-1}\) at \(\sqrt{s}=13~{\rm TeV}\). The local CMS snapshot confirms this dataset statement, but `B006.yaml` `pdg_or_equivalent` records only the branching-fraction limits and not those dataset numerical claims with source URL/access-date/sha256 provenance. Rework can either add a YAML provenance block for the dataset numerics or remove the numerics from the TeX prose.

## Verification notes
- Source manifests: every process-local source key cited in the TeX appears in the corresponding `source_manifest.yaml`. Each manifest entry has a non-empty `snapshot_path`; `git ls-files flavor_catalog/references/B004 B006 B019 B026` confirms the pointed snapshots are tracked text files.
- Checksums: recomputed `sha256sum *.txt` for all four reference directories. The snapshot checksums match the YAML/source-manifest values used above.
- PDF policy: `ls flavor_catalog/references/B004`, `B006`, `B019`, and `B026` shows only `.txt`, `.yaml`, and checksum files; no `.pdf` files are tracked.
- Status lineage: before CA mutation, all four sidecars contained `WRITER-INITIATED` followed by `WRITER-DONE` with ISO-8601 timestamps, and `WRITER-DONE` was the latest status.
- Code coverage: this `rg` treats `-E` as an encoding flag, so the literal prompt form `rg -l -E "<process keyword>" ...` fails. I used the equivalent `rg -l -e` searches over the required directories. B004 direct \(\phi_s\) searches returned no direct observable while all cited partial-coverage file:lines exist. B006's hits are generic pseudoscalar/Delta-F=2 code, not \(B^0\to\mu^+\mu^-\). B019 and B026 targeted searches returned no implementation hits after excluding notebooks, matching the TeX claims.
- Implementation difficulty: B004 LOW matches existing \(\Delta F=2\) operator machinery; B006 MEDIUM matches a new \(\Delta B=1\) leptonic operator normalization with standard decay-constant inputs; B019 HIGH matches a new semileptonic mode/bin/QED/likelihood calculation; B026 HIGH matches new charged-current semileptonic matching, form-factor integration, and correlated \(R_D\)-\(R_{D^*}\) likelihood.
- rc1.1 consistency: these catalog processes are companion-mode additions. A repo search outside `flavor_catalog/` found only generic \(B_d/B_s\) bag/mixing audit references, not load-bearing rc1.1 numerical statements for B004/B006/B019/B026 that the catalog contradicts.
- Cosmetic LaTeX: no `\textbf{CHECK}`, `TODO`, `\ref{}`, or `\cite{}` tokens were present in the four TeX files. A brace-balance pass returned zero for all four files.
