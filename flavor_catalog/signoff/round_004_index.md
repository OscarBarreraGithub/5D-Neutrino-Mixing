# Flavor Catalog Opus Sign-Off — Round 004

**Date**: 2026-05-17
**Round ID**: `opus_round_004`
**Reviewer**: Opus round 4 sign-off agent
**Branch**: `flavor-catalog/2026q2`
**Input commit (head before this round)**: `7e121c3a1764553cdb5a76da6546da3705dfb82f`
**Scope**: 8 Wave-8 SECONDARY entries — K019, K020, K021 (kaon LFV trio), B007, B008 (rare leptonic B -> e e / tau tau), B013, B014 (exclusive radiative b -> s gamma / b -> d gamma), T014 (Z FCNC). All tagged `priority_tier: SECONDARY` per `flavor_catalog/PRIORITY_TIERS.md` and located at `flavor_catalog/processes/secondary/<family>/<id>.{tex,yaml}`. All reached `CHECKER-DONE` and ran the family fact-check addenda.

**Plan and precedent references**: `flavor_catalog/AGENTIC_WORKFLOW.md` (Sections 2, 4, 5 — roles, 8-item checklist, L001 CHK-1 carve-out); `flavor_catalog/PRIORITY_TIERS.md` (Wave-8 SECONDARY policy); `flavor_catalog/worklogs/orchestration/wave_008_runbook.md` (live wave-8 state); `flavor_catalog/signoff/round_003_index.md` (Wave-7 sign-off template); `flavor_catalog/signoff/by_process/L001.md`, `B021_B023.md`, `B001_B003.md` (CHK-1 carve-out precedents).

**Method**: For each Wave-8 SECONDARY process I (i) confirmed the YAML sidecar carries `priority_tier: SECONDARY`, a non-empty `priority_rationale`, and `promoted_in_wave: 8`, and lives under `processes/secondary/<family>/`; (ii) confirmed `status_history` reached `CHECKER-DONE` after the terminal WA edit (cycle 1 for K019, K021, B007, B008, B013, B014, T014; cycle 2 for K020 after the L001-aligned CHK-1 rework); (iii) read the terminal CA worklog and verified CHK-1..8 plus CHK-W8 all PASS, with the L001 / B001_B003 carve-out applied uniformly for theory normalization scales (CFW 21/33 TeV, Perez-Randall 3 TeV), dataset descriptors (LHCb 3 / 9 fb^-1, 7/8/13 TeV, event yields), FLAG/HFLAV lattice averages, FCC-ee projections, SM EFT/RG predictions, and coupling translations; (iv) confirmed the Wave-8 fact-check addendum row VERIFIED for 7 of 8 entries and PARTIAL for K020. The K020 PARTIAL is metadata-only (NA62 2021 first-author convention: manifest "E. Cortina Gil et al.", live arXiv "NA62 Collaboration", PDG "R. Aliberti et al."); all numerical values, year, DOI, and URL are correct. This is the same class as the v0.2 E009 PARTIAL — author-line convention drift on a collaboration paper — and is ACCEPTed per that precedent. The bias is `APPROVE` unless a substantive physics error is found.

## Summary statistics

| Verdict | Count |
|---|---|
| `APPROVE` | **8** |
| `RETURN-TO-WA` | 0 |
| `ESCALATE-TO-PI` | 0 |
| **Total** | **8** |

| Fact-check | Count |
|---|---|
| VERIFIED | 7 |
| PARTIAL | 1 (K020: metadata-only, accepted per v0.2 E009 precedent) |
| MISMATCH | 0 |
| FAILED | 0 |

## Per-process verdicts

| process_id | observable | terminal CA worklog | fact-check | verdict | justification |
|---|---|---|---|---|---|
| K019 | `BR(K_L -> e^+/- mu^-/+)` PDG/E871 | `ca_w8_kaon_LFV.md` (cycle-1 PASS) | VERIFIED (0) | APPROVE | PDG 2025 `< 4.7e-12` at 90% CL; BNL E871 (Ambrose 1998) `< 4.7e-12` at 90% CL. Both in `pdg_or_equivalent.values` with full strict metadata. MEDIUM difficulty correct. CHK-W8 PASS (SECONDARY sidecar, secondary/kaon/ path, `.tex` notes PRIORITY_TIERS.md). |
| K020 | `BR(K^+ -> pi^+ e^± mu^∓)` PDG/Sher/NA62 | `ca_w8_kaon_LFV_v2.md` (cycle-2 PASS) | PARTIAL (1: NA62 author-line metadata) | APPROVE | PDG/API `< 1.3e-11` (S010.29) and `< 6.6e-11` (S010.25); cycle-2 promoted Sher/E865-only `< 2.1e-11` to `pdg_or_equivalent.values` (measured observable, not theory normalization — L001 precedent applied in the intended direction). PARTIAL fact-check is NA62 2021 first-author convention drift (manifest "E. Cortina Gil et al." vs. live arXiv "NA62 Collaboration" vs. PDG "R. Aliberti et al."); values/year/DOI/URL all correct; same metadata-only class as v0.2 E009 PARTIAL — ACCEPT. HIGH difficulty correct. |
| K021 | `BR(K_L -> pi^0 e^± mu^∓)` PDG/KTeV + NA62 charged-companion limits | `ca_w8_kaon_LFV.md` (cycle-1 PASS) | VERIFIED (0) | APPROVE | PDG 2025 `< 7.6e-11` at 90% CL summed; KTeV exact `< 7.56e-11`; NA62 2021 companion `< 6.6e-11`, `< 4.2e-11`, `< 3.2e-10`. All in `pdg_or_equivalent.values` with full metadata. KTeV zero-event context kept under `supporting_measurements` per L001. HIGH difficulty correct. |
| B007 | `BR(B_{s,d} -> e^+ e^-)` PDG/LHCb/CDF/BABAR/Belle | `ca_w8_B_rare_leptonic.md` (cycle-1 PASS) | VERIFIED (0) | APPROVE | PDG 2026 `BR(B_s -> ee) < 9.4e-9` and `BR(B_d -> ee) < 2.5e-9` at 90% CL; LHCb 2020 95% CL companions `< 11.2e-9`, `< 3.0e-9`; historical CDF / BABAR / Belle limits also in `pdg_or_equivalent.values`. Bobeth/Fleischer SM predictions and CFW reference scales kept in supporting/auxiliary blocks per L001/B001_B003. MEDIUM difficulty correct. |
| B008 | `BR(B_{s,d} -> tau^+ tau^-)` PDG/LHCb/BABAR | `ca_w8_B_rare_leptonic.md` (cycle-1 PASS) | VERIFIED (0) | APPROVE | PDG 2026 `BR(B_s -> tautau) < 6.8e-3` and `BR(B_d -> tautau) < 2.1e-3` at 95% CL; historical BABAR `< 4.1e-3`; LHCb 3 fb^-1 dataset descriptors and Bobeth SM (B_s 7.73e-7, B_d 2.22e-8) kept in supporting/theory blocks per L001. MEDIUM difficulty correct. |
| B013 | `B_s -> phi gamma` BR / time-dep S, C / A_Delta | `ca_w8_B_radiative.md` (cycle-1 PASS) | VERIFIED (0) | APPROVE | PDG 2025 `(3.4 +/- 0.4)e-5`; HFLAV Dec-2024 `(34.0 +/- 3.2)e-6`; LHCb 2019 `S = 0.43 +/- 0.30(stat) +/- 0.11(syst)`, `C = 0.11 +/- 0.29 +/- 0.11`, `A_Delta = -0.67 +0.37/-0.41 +/- 0.17`. All in `pdg_or_equivalent` with full metadata. HIGH difficulty correct. |
| B014 | `B -> rho gamma, B -> omega gamma` PDG / HFLAV / Belle II / LHCb 2025 | `ca_w8_B_radiative.md` (cycle-1 PASS) | VERIFIED (0) | APPROVE | PDG 2025 `(9.8 +/- 2.5)e-7`, `(8.6 +/- 1.5)e-7`, `(4.4 +1.8/-1.6)e-7`; HFLAV Dec-2024 `B -> rho gamma = (1.40 +/- 0.22)e-6`; isospin `-0.46 +/- 0.17` (HFLAV naive caveat preserved); semi-inclusive `B -> X_d gamma = (9.2 +/- 3.0)e-6`. HIGH difficulty correct. |
| T014 | `Z -> bs, bd, sd` ECFA 2025 / Kamenik 2024 / Abu-Ajamieh 2026 | `ca_w8_T014.md` (cycle-1 PASS) | VERIFIED (0) | APPROVE | Three direct non-standard hadronic-Z width limits `< 2.9e-3` at 95% CL (bs, bd, sd) in `pdg_or_equivalent.values` with full ECFA/Kamenik/Abu-Ajamieh provenance; PDG `(69.911 +/- 0.056)%` inclusive hadronic context kept as context only; SM predictions, FCC-ee projections, and coupling translations kept in supporting/theory blocks per L001/B001_B003. HIGH difficulty correct. |

All eight are APPROVE. None requires RETURN-TO-WA; none requires ESCALATE-TO-PI.

## Cross-cutting observations

- L001 / B001_B003 / B021 precedent uniformly applied across the wave: theoretical reference scales (CFW 21/33 TeV, Perez-Randall 3 TeV, KK-gluon scales), dataset metadata (LHCb 3 / 9 fb^-1, 7/8/13 TeV, event yields, KTeV/E871 zero-event statements), FLAG/HFLAV lattice averages, FCC-ee projections, SM EFT predictions, and coupling translations live under `paper_era_reference`, `supporting_measurements`, `auxiliary_theory_inputs`, or `theory_context` — not in `pdg_or_equivalent`. No Wave-8 entry triggered Opus arbitration.
- The single cycle-2 rework (K020 CHK-1: Sher/E865 `< 2.1e-11` promotion to `pdg_or_equivalent.values`) is policy-aligned in the **intended** direction of the L001 carve-out: a measured experimental upper limit was promoted, while theory normalization and dataset descriptors stayed in supporting blocks. This is the correct application of the precedent and not a CHK-1 false positive.
- K020 fact-check PARTIAL is metadata-only: the NA62 2021 first-author convention disagrees between the manifest (`E. Cortina Gil et al.`), the live arXiv landing page (`NA62 Collaboration`), and the PDG/journal convention (`R. Aliberti et al.`). The numerical limit, year, DOI, URL, and snapshot sha256 all verify. This is the same class as the v0.2 E009 PARTIAL (collaboration author-line convention drift on a published-paper citation) and is ACCEPTed per that precedent without TeX or YAML edits.
- SECONDARY tier policy uniformly satisfied: all 8 sidecars carry `priority_tier: SECONDARY`, `priority_rationale: <one-line: why deferred + why now>`, and `promoted_in_wave: 8`; all 8 file pairs live under `flavor_catalog/processes/secondary/<family>/`; every `.tex` carries a SECONDARY-tier note pointing at `flavor_catalog/PRIORITY_TIERS.md`. CHK-W8 PASS on all eight.
- Iteration caps respected: K020 used WA/CA cycle 2 inside the 3-cycle cap; the other seven entries passed cycle 1.
- Total catalog after Wave-8: 80 PRIMARY (Waves 1-7) + 8 SECONDARY (Wave-8) = **88 entries**.

## What is intentionally NOT changed

No `.tex` edits in this round. Sidecar edits are only the two `OPUS-APPROVED` + `FACT-CHECKED` `status_history` transitions on each of the 8 SECONDARY YAMLs, plus matching `opus_approved_at` / `fact_checked_at` / `fact_check_verdict` / `last_updated_at` field updates. No `pdg_or_equivalent`, `paper_era_reference`, `supporting_measurements`, `auxiliary_*`, `source_shas`, `code_coverage`, or `theory_context` blocks modified. No PRIMARY (Waves 1-7) sidecar touched. No `catalog_master.tex` edits.

===OPUS_ROUND_004_END===
