# Charged-lepton fact-check report

Date: 2026-05-16
Agent: factcheck-codex-charged_lepton
Scope: L001 L002 L003 L004 L005 L006 L007 L008 L009 L010 L023

Overall verdict: VERIFIED for 11/11 processes.
Mismatches: 0.
Unresolvable sources: 0.

Method: I fetched the catalog `source_url` targets where possible, including PDG live pages/PDF listings, arXiv abstract pages, arXiv e-prints for table/body-only values, Crossref/CDS/PubMed metadata, and the local snapshots under `flavor_catalog/references/<id>/`. I skipped theory normalization constants and repo-derived NDA prefactors where the signoff policy classifies them as theory/policy inputs rather than measured observables.

## L001

Verdict: VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| MEG II 2025 `BR(mu -> e gamma) < 1.5e-13` at 90% CL; sensitivity `2.2e-13` | https://arxiv.org/abs/2504.15711 | Y | ArXiv title/authors match MEG II; abstract contains both values. |
| PDG 2025 listing `BR < 3.1e-13` at 90% CL | https://pdg.lbl.gov/encoder_listings/s004.pdf | Y | PDG muon listing contains `NODE=S004R3` and `<3.1 x 10^-13`; catalog notes PDG lag relative to MEG II 2025. |
| MEG II 2024 first-dataset `7.5e-13` and MEG+MEG II combination `3.1e-13` | https://arxiv.org/abs/2310.12614 | Y | Local snapshot contains both limits. |
| MEG 2016 full-dataset `4.2e-13`; `7.5e14` stopped muons | https://arxiv.org/abs/1605.05081 | Y | Local snapshot contains the limit and exposure. |
| Perez-Randall paper-era input `1.2e-11` | https://arxiv.org/abs/0805.4652 | Y | Local snapshot contains the paper-era bound context. |
| Perez-Randall/repo prefactor `4.0e-8`, `c=0.02`, derived `C=0.001936...` | https://arxiv.org/abs/0805.4652 | N/A | Theory/repo normalization and derived values; skipped under policy carve-out. |

No mismatch or unresolvable source.

## L002

Verdict: VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| PDG `BR(mu -> 3e) < 1.0e-12` at 90% CL | https://pdgprod.lbl.gov/pdgprod/pdgLive/DataBlock.action?node=S004R4 | Y | PDG live page row is `<1.0` in `10^-12` units, BELLGARDT 1988/SINDRUM. |
| SINDRUM original limit `Gamma(mu -> 3e)/Gamma(mu -> e 2nu) < 1.0e-12` | https://inspirehep.net/api/literature/251865 | Y | Inspire/local snapshot records the SINDRUM title and `1.0e-12` limit. |
| Mu3e phase-I single-event sensitivity `2.0e-15`; rate up to `1.0e8` muon decays/s | https://arxiv.org/abs/2009.11690 | Y | ArXiv TDR abstract contains both claims. |
| Mu3e status: phase-I `O(1e-15)`, upgrade `O(1e-16)`, SM-with-massive-neutrino rate `O(1e-54)` | https://arxiv.org/abs/2501.14667 | Y | Local snapshot contains the stated sensitivity scales and SM estimate. |
| RS/EFT context references | https://arxiv.org/abs/1702.03020, https://arxiv.org/abs/hep-ph/0606021, https://arxiv.org/abs/0804.1954 | N/A | Bibliographic/key-reference checks only; no measured observable value used. |

No mismatch or unresolvable source.

## L003

Verdict: VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| Published benchmark limit `R_mue(Au) < 7.0e-13` at 90% CL, not an Al measurement | https://juser.fz-juelich.de/record/55087 | Y | JuSER/local snapshot matches SINDRUM II gold paper; catalog explicitly labels it as a cross-target benchmark. |
| No published direct aluminum-target limit found | Mu2e/COMET projection sources below | Y | Supported by catalog synthesis and projection-only source set. |
| Mu2e Al projection: SES `2.8e-17`; expected 90% UL `6.7e-17` | https://arxiv.org/abs/1606.05559 | Y | ArXiv abstract contains both numbers. |
| Mu2e Run I projection: expected UL `6.2e-16`; 5 sigma discovery `1.2e-15`; background `0.11 +- 0.03` | https://arxiv.org/abs/2210.11380 | Y | Local snapshot contains all values. |
| COMET Phase-I Al projection: SES `3.1e-15`; expected 90% UL `7.0e-15`; background `0.032` | https://arxiv.org/abs/1812.09018 | Y | ArXiv abstract/local snapshot contains all values. |
| COMET Phase-II strategy SES `2.6e-17` after `2.0e7` s | https://arxiv.org/abs/1812.07824 | Y | Local snapshot contains the Phase-II sensitivity and runtime. |
| Mu2e-II goal: at least one additional order beyond Mu2e | https://arxiv.org/abs/2203.07569 | Y | ArXiv abstract contains the qualitative sensitivity goal. |
| CFW RS baseline | https://arxiv.org/abs/0804.1954 | N/A | Theory context only. |

No mismatch or unresolvable source.

## L004

Verdict: VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| PDG/SINDRUM II gold limit `<7.0e-13` at 90% CL | https://pdg.lbl.gov/encoder_listings/s004.pdf | Y | Fetched PDG muon listing contains BERTL 06/SINDRUM II `<7 x 10^-13`. |
| SINDRUM II gold paper title/metadata | https://api.crossref.org/works/10.1140/epjc/s2006-02582-x | Y | Crossref/local snapshot matches `A search for mu-e conversion in muonic gold`. |
| Mu2e TDR four-orders improvement context | https://arxiv.org/abs/1501.05241 | Y | ArXiv/local snapshot confirms the program context. |
| COMET Phase-I `3.1e-15`, `7.0e-15`, factor `100`, background `0.032` | https://arxiv.org/abs/1812.09018 | Y | ArXiv abstract/local snapshot contains all values. |
| CFW RS baseline | https://arxiv.org/abs/0804.1954 | N/A | Theory context only. |

No mismatch or unresolvable source.

## L005

Verdict: VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| PDG/SINDRUM II titanium limit `<4.3e-12` at 90% CL | https://pdg.lbl.gov/encoder_listings/s004.pdf | Y | Fetched PDG listing contains DOHMEN 93/SINDRUM II `<4.3 x 10^-12` and ground-state note. |
| SINDRUM II titanium paper title/metadata | https://api.crossref.org/works/10.1016/0370-2693(93)91383-X | Y | Crossref/local snapshot matches the titanium conversion paper. |
| Mu2e TDR four-orders improvement context | https://arxiv.org/abs/1501.05241 | Y | ArXiv/local snapshot confirms the program context. |
| COMET Phase-I SES `3.1e-15` and expected 90% UL `7.0e-15` | https://arxiv.org/abs/1812.09018 | Y | ArXiv abstract/local snapshot contains both values. |
| CFW RS baseline | https://arxiv.org/abs/0804.1954 | N/A | Theory context only. |

No mismatch or unresolvable source.

## L006

Verdict: VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| PDG `P(M -> anti-M) < 8.3e-11` at 90% CL in 0.1 T | https://pdgprod.lbl.gov/pdgprod/pdgLive/DataBlock.action?node=S004MC | Y | PDG live footnote contains `P_M Mbar < 8.3 x 10^-11` at 90% CL. |
| PDG effective coupling `G_C/G_F < 0.0030` | https://pdgprod.lbl.gov/pdgprod/pdgLive/DataBlock.action?node=S004MC | Y | PDG live row contains `<0.0030` for `R_g = G_C/G_F`. |
| Willmann/MACS original experiment `P <= 8.3e-11` | https://arxiv.org/abs/hep-ex/9807011 | Y | Article/PDG canonical value is `8.3e-11`; arXiv abstract says `8.2e-11`, already noted in YAML. This is not a mismatch. |
| MACE/Snowmass improvement `>2` orders of magnitude | https://arxiv.org/abs/2203.11406 | Y | Local snapshot contains the improvement claim. |
| MACE CDR probability reach beyond `1.0e-13` | https://arxiv.org/abs/2410.18817 | Y | ArXiv abstract/local snapshot contains the `10^-13` reach. |
| RS/operator context references | https://arxiv.org/abs/hep-ph/0606021, https://arxiv.org/abs/0804.1954, https://arxiv.org/abs/2401.09580 | N/A | Theory/operator context only. |

No mismatch or unresolvable source.

## L007

Verdict: VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| PDG `BR(tau -> mu gamma) < 4.2e-8` at 90% CL | https://pdglive.lbl.gov/view/S035.31 | Y | PDG live page contains the `4.2 x 10^-8` limit. |
| Belle full-data result `<=4.2e-8` with `988 fb^-1` | https://arxiv.org/abs/2103.12994 | Y | ArXiv abstract contains the mu-gamma limit and luminosity. |
| BaBar prior limit `<4.4e-8`; `(963 +- 7)e6` tau decays | https://arxiv.org/abs/0908.2381 | Y | Local snapshot contains both values. |
| Belle II Physics Book projection `<5.0e-9` at `50 ab^-1` | https://arxiv.org/abs/1808.10567 | Y | Local snapshot contains the projection. |
| CFW/Perez-Randall context references | https://arxiv.org/abs/0804.1954, https://arxiv.org/abs/0805.4652 | N/A | Theory context only. |

No mismatch or unresolvable source.

## L008

Verdict: VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| PDG `BR(tau -> e gamma) < 3.3e-8` at 90% CL | https://pdglive.lbl.gov/view/S035.32 | Y | PDG live page contains the `3.3 x 10^-8` limit. |
| BaBar primary result `<3.3e-8`; `(963 +- 7)e6` tau decays | https://arxiv.org/abs/0908.2381 | Y | Local snapshot contains the e-gamma limit and sample size. |
| Belle cross-check `<=5.6e-8` with `988 fb^-1` | https://arxiv.org/abs/2103.12994 | Y | ArXiv abstract contains the e-gamma limit and luminosity. |
| Belle II Physics Book projection `<1.2e-8` at `50 ab^-1` | https://arxiv.org/abs/1808.10567 | Y | Local snapshot contains the projection. |
| CFW/Perez-Randall context references | https://arxiv.org/abs/0804.1954, https://arxiv.org/abs/0805.4652 | N/A | Theory context only. |

No mismatch or unresolvable source.

## L009

Verdict: VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| PDG tau listing `BR(tau -> 3mu) < 1.9e-8` at 90% CL | https://pdg.lbl.gov/2025/listings/rpp2025-list-tau.pdf | Y | Local PDG snapshot contains ADACHI 24R/Belle II and `<1.9e-8`. |
| Belle II primary result `<1.9e-8` with `424 fb^-1` | https://arxiv.org/abs/2405.07386 | Y | ArXiv abstract contains the limit and luminosity. |
| LHCb Run 2 preprint `1.9e-8` at 90% CL and `2.3e-8` at 95% CL; `5.4 fb^-1` | https://arxiv.org/abs/2601.20785 | Y | ArXiv abstract contains both limits and sample. |
| Belle 2010 historical limit `<2.1e-8`; `782 fb^-1` | https://arxiv.org/abs/1001.3221 | Y | Local snapshot contains the limit and luminosity. |
| BaBar 2010 historical limit `<3.3e-8`; `468 fb^-1` | https://arxiv.org/abs/1002.4550 | Y | ArXiv e-print table contains the `mmm` mode limit and sample. |
| LHCb 2015 historical limit `<4.6e-8` | https://arxiv.org/abs/1409.8548 | Y | Local snapshot contains the limit. |
| CFW RS baseline | https://arxiv.org/abs/0804.1954 | N/A | Theory context only. |

No mismatch or unresolvable source.

## L010

Verdict: VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| PDG live `BR(tau -> 3e) < 2.7e-8` at 90% CL | https://pdglive.lbl.gov/BranchingRatio.action?desig=38&parCode=S035&home= | Y | PDG live page contains the `2.7 x 10^-8` limit and Belle reference. |
| Belle primary result `<2.7e-8`; `782 fb^-1`, `719M` tau pairs, efficiency `6.0%`, background `0.21 +- 0.15`, 0 observed | https://arxiv.org/abs/1001.3221 | Y | Local snapshot contains all table values. |
| BaBar support result `<2.9e-8`; `468 fb^-1`, efficiency `8.6 +- 0.2`, background `0.12 +- 0.02`, expected limit `3.4e-8`, 0 observed | https://arxiv.org/abs/1002.4550 | Y | ArXiv e-print table contains all values. |
| Belle II Physics Book prospect `O(10^-9)-O(10^-10)` at `50 ab^-1` | https://arxiv.org/abs/1808.10567 | Y | Local snapshot contains the prospect statement. |
| LHCb tau->3mu context only | https://arxiv.org/abs/2601.20785 | N/A | Companion-mode context, not a L010 numerical input. |
| Agashe/CFW RS context references | https://arxiv.org/abs/hep-ph/0606021, https://arxiv.org/abs/0804.1954 | N/A | Theory context only. |

No mismatch or unresolvable source.

## L023

Verdict: VERIFIED

| Claim | Source URL | Verified? | Note |
|---|---|---:|---|
| CHARM-II ratio `sigma_exp/sigma_SM = 1.58 +- 0.64` | https://arxiv.org/abs/1902.06765 | Y | DUNE trident e-print contains the 1.58 +- 0.64 compiled ratio; CDS record confirms CHARM-II paper metadata. |
| CHARM-II primary paper metadata | https://cds.cern.ch/record/208653 | Y | CDS/local snapshot confirms the first-observation paper. |
| CCFR event counts `37.0 +- 12.4` observed, `45.3 +- 2.3` SM; ratio `0.82 +- 0.28` | https://pubmed.ncbi.nlm.nih.gov/10043703/ and https://arxiv.org/abs/1902.06765 | Y | PubMed/local snapshot contains the CCFR paper metadata; DUNE e-print and local snapshot contain counts/ratio. |
| NuTeV data `17`, background `9.8 +- 2.9`, SM `10.8 +- 0.3`, ratio `0.72 +1.73 -0.72` | https://arxiv.org/abs/hep-ex/9811012 and https://arxiv.org/abs/1902.06765 | Y | NuTeV snapshot and DUNE e-print contain the event counts and ratio. |
| Altmannshofer 2014 CCFR exclusion contours at `95%` CL | https://arxiv.org/abs/1406.2332 | Y | ArXiv e-print contains the 95% CL contour statement. |
| Altmannshofer 2014 CHARM-II convention `1.58 +- 0.57` | https://arxiv.org/abs/1406.2332 | Y | E-print contains the alternate uncertainty convention; YAML note covers the difference from DUNE 2019. |
| Kaneta-Shimomura Belle II auxiliary inputs: `sqrt(s)=10.58 GeV`, `50 ab^-1` | https://arxiv.org/abs/1701.00156 | Y | ArXiv/local snapshot contains both auxiliary inputs. |
| DUNE trident projections: abstract `25%`; about `3` years per beam mode; body `40%` baseline and `25%` improved contours | https://arxiv.org/abs/1902.06765 | Y | ArXiv e-print contains all projection values. |
| CFW RS baseline | https://arxiv.org/abs/0804.1954 | N/A | Theory context only. |

No mismatch or unresolvable source.

## Summary

All charged-lepton process entries passed independent content checks. The only noteworthy source nuance is L006: the Willmann arXiv abstract reports `8.2e-11`, while PDG live and the catalog use `8.3e-11`; the YAML already documents this and the PDG canonical value is verified.
