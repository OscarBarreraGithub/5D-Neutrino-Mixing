# Audit Fix Ledger — full_repo_audit_2026-07

Tracks disposition of **every** finding in [`AUDIT_COMPENDIUM.md`](AUDIT_COMPENDIUM.md).
Process per code finding (per user directive): **Codex researches → Codex fixes → Claude audits**;
serial (never parallel) to stay within session limits. A finding advances to `VERIFIED`
only after a Claude audit subagent signs off on the Codex fix. Findings that require **no code
change** (already-correct, revised-away, documentation-only, or deliberately-deferred) are marked
with an explicit rationale — nothing is silently dropped.

Disposition legend:
- `PENDING` — not yet started
- `IN-PROGRESS` — Codex researching/fixing or Claude auditing
- `FIXED` — Codex fix applied, awaiting Claude audit
- `VERIFIED` — Claude audit signed off
- `DOC-ONLY` — resolved by documentation/banner change, no physics-code change
- `WONT-FIX` — deliberately not changed; rationale recorded
- `ALREADY-CORRECT` — audit re-checked; code was right (revision in §8.3 or on re-verification)

Column `Grp` groups findings fixed together in one Codex/Claude cycle.

---

## Priority 1 — Lane C (`paper_0710_1869/`) bundle (§8.6 item 1)

| ID | Grp | Severity | Title | Location | Disposition |
|----|-----|----------|-------|----------|-------------|
| C-1 | L | CRIT | Lane C RG inverted + untransposed ADM | `paper_0710_1869/eft_deltaf2/rg.py:288-313,1361-1365,1516-1519` | VERIFIED |
| C-2 | L | CRIT | Lane C default path zeroes C4 (C4≡0) | `paper_0710_1869/eft_deltaf2/matching_kkgluon.py:191-198,504-527` | VERIFIED (capability fixed; RESIDUAL: default alignment model pending paper) |
| C-5 | L | CRIT | Lane C missing √(2πkr_c) volume factor | `paper_0710_1869/.../couplings.py:241-243` | VERIFIED |
| C-6 | L | CRIT | Lane C seed→profile inverts RS localization | `paper_0710_1869/.../fit.py` (physical branch) | VERIFIED (orientation fixed; RESIDUAL: exact Table-I coeffs pending paper) |
| M-8 | L | MAJ | Lane C Fierz sign in paper→BMU map (+2 vs −2) | `paper_0710_1869/eft_deltaf2/rg.py:110-125` | VERIFIED |
| M-27 | L | MAJ | Lane C matrix elements ×4 (P_L/P_R vs (1∓γ5)) | `paper_0710_1869/.../hadronic.py` | VERIFIED |
| M-28 | L | MAJ | Lane C circular validation + strict_paper Table I fail | `paper_0710_1869/` tests/verifier | VERIFIED (RG test de-circularized; benchmark/model self-formula tests remain regression-guards only, acceptable for quarantined lane) |
| M-13(C) | L | MAJ | Lane C ξ_g=1 KK-gluon default | `paper_0710_1869/.../kkgluon.py` | VERIFIED |

## Priority 2 — LFV normalizations (§8.6 item 2)

| ID | Grp | Severity | Title | Location | Disposition |
|----|-----|----------|-------|----------|-------------|
| M-29 | LFV | MAJ | μ→3e interference coeff 2√2·e vs KO 8e (+sign) | `lfv_three_body.py:778` | VERIFIED |
| M-30 | LFV | MAJ | τ→3ℓ / τ→ℓγ missing BR(τ→ℓνν̄)≈0.177 (×5.65) | `lfv_three_body.py`, L007/L008 prefactor | VERIFIED |

## Priority 3 — Mechanical ×2 / kernel / proxy (§8.6 item 3)

| ID | Grp | Severity | Title | Location | Disposition |
|----|-----|----------|-------|----------|-------------|
| M-9 | TOP | MAJ | t→qγ/qg dipole widths ×2 too large | `top_fcnc.py:225,241` | VERIFIED |
| M-31 | TOP | MAJ | M-9 propagates to T003–T006 + tests hardcode buggy BR | `processes/top_higgs_ew/`, tests | VERIFIED |
| M-32 | BK | MAJ | B→K*μμ kernel: missing q² + transverse ×2 + primed-Wilson sign | `rare_b_kstar_dilepton.py`; `test_B019` | VERIFIED |
| M-33 | RD | MAJ | R(D)/R(D*) stress proxy Wilson dimensionful (GeV²) | `semileptonic_lfu.py:192` | VERIFIED |

## Priority 4 — compare_2007_vs_modern (§8.6 item 4)

| ID | Grp | Severity | Title | Location | Disposition |
|----|-----|----------|-------|----------|-------------|
| C-7 | CMP | CRIT | Rescales scan ratios by bounds scan never used (~1e6 units) | `compare_2007_vs_modern.py` | VERIFIED (invalid all-system comparison retracted; superseded/corrected banners; no dangling refs) |

## Priority 5 — everything else (§7)

### εK / ΔF=2 budget coherence (M-1…M-7)
| ID | Grp | Severity | Title | Location | Disposition |
|----|-----|----------|-------|----------|-------------|
| M-1 | EPSK | MAJ | εK budget split-brain ×4.5 core vs catalog | `deltaf2.py:794-795`, `K001.py:206-231` | VERIFIED (unified policy, mandate-strict band) |
| M-2 | EPSK | MAJ | K001 budget sign-blind | `K001.py:211-218` | VERIFIED (sign-aware, K004 pattern) |
| M-3 | EPSK | MAJ | Opposite budget philosophy for Δm | `deltaf2.py:956-968` | PENDING |
| M-4 | EPSK | MAJ | Mixed CLs across HARD vetoes catalog-wide | catalog-wide | VERIFIED (policy-id + CL labeled; global-CL disclaimed) |
| M-5 | EPSK | MAJ | B022 (B→Kνν) HARD-vetoes SM limit | `primary/beauty/B022.py` | PENDING |
| M-6 | RG | MAJ | B-meson Wilsons over-run below m_b (~25% C4 inflation) | `deltaf2.py:464`, `modern/evaluation.py:177` | PENDING |
| M-7 | RG | MAJ | LO-only ΔF=2 running biases εK floor low ~10-15% | `deltaf2.py` | PENDING |

### Normalization / convention (M-10..M-14)
| ID | Grp | Severity | Title | Location | Disposition |
|----|-----|----------|-------|----------|-------------|
| M-10 | HYUK | MAJ | rs_higgs_yukawas.py unfixed twin of B1 Zbb bug | `rs_higgs_yukawas.py:292-304` | PENDING |
| M-11 | LANEA | MAJ | Lane A Bauer bridge factor-2 (Y_max≈1.5 not 3) | `anarchic_bauer_s1.py:24-33,191-201` | PENDING |
| M-12 | LANEA | MAJ | Lane A drops Wolfenstein A + hash() seeds | `anarchic_bauer_s1.py:194,408` | PENDING |
| M-13 | GSSTAR | MAJ | Lane B/C KK-gluon g_s* normalization ~13× | `quark_model_core` scales/matching | PENDING |
| M-14 | LFVCONV | MAJ | μ→eγ M_KK convention split ×2.45⁴ | `scanParams/scan.py:226-227`, `rs_ew_builder.py:131` | PENDING |

### Fit / scan logic (M-15..M-18)
| ID | Grp | Severity | Title | Location | Disposition |
|----|-----|----------|-------|----------|-------------|
| M-15 | FIT | MAJ | Reported fit seed not gauge-equivalent | `fit.py:627-655,843` | PENDING |
| M-16 | FIT | MAJ | Scan seed chaining double-applies overall_scale | `scan.py:360-373`, `validation.py:242-271` | PENDING |
| M-17 | HARDP | MAJ | Evaluated hard-partial never vetoes/miscounted | `run_full_catalog_scan.py:913-921` | PENDING |
| M-18 | HARDP | MAJ | Catalog silent-pass on invalid input | `TEMPLATE.py` design | PENDING |

### QCD running (M-19..M-21)
| ID | Grp | Severity | Title | Location | Disposition |
|----|-----|----------|-------|----------|-------------|
| M-19 | QCD | MAJ | 3-loop mass-decoupling d3 wrong (value+slope sign) | `qcd/decoupling.py:117-128` | PENDING |
| M-20 | QCD | MAJ | Spurious 6→6 threshold matching (m_t 163.5 vs 162.5) | `qcd/mass_running.py:187-191` | PENDING |
| M-21 | QCD | MAJ | PDG m_t(MS̄) σ understated ~2.5× + load-bearing test | `pdg_quark_masses.py:110-116`, test | PENDING |

### Collider / EW existence floor (C-3, M-26)
| ID | Grp | Severity | Title | Location | Disposition |
|----|-----|----------|-------|----------|-------------|
| C-3 | EW | CRIT | EW001 oblique anchors + 15.96 vs 18-20 TeV floor | `EW001.yaml:82-118` | PENDING |
| M-26 | EW | MAJ | Collider SSM benchmarks as HARD vetoes to ×L-suppressed RS | CR005/CR006, CR001 | PENDING |

### Bessel solver (C-4)
| ID | Grp | Severity | Title | Location | Disposition |
|----|-----|----------|-------|----------|-------------|
| C-4 | BESSEL | CRIT | KK Bessel solver skips/misorders roots n≥5 | `solvers/bessel.py:153-179,264-285` | PENDING |

### Lane A modern verification (M-25, M-34, M-35)
| ID | Grp | Severity | Title | Location | Disposition |
|----|-----|----------|-------|----------|-------------|
| M-25 | MODERN | MAJ | modern/ circular verification + wrong bound fields | `modern/` verifiers, `matching.py:517-543` | PENDING |
| M-34 | MODERN | MAJ | Self-referential validation scripts + PR tension | `benchmark_quark_0710_1869.py`, `audit_wilson_rg.py` | PENDING |
| M-35 | WEB | MAJ | Website envelope floor statistically optimistic (max ignores joint) | website Scan Explorer | PENDING |

### Documentation / reproducibility (M-22, M-23, M-24, M-36)
| ID | Grp | Severity | Title | Location | Disposition |
|----|-----|----------|-------|----------|-------------|
| M-22 | DOC | MAJ | review_local/ certifies pre-B3 physics as CONFIRMED | `review_local/*.tex`, `docs/audits/*` | PENDING |
| M-23 | DOC | MAJ | Headline floor parquet absent from tree | `FLOOR_SUMMARY.md` refs | PENDING |
| M-24 | DOC | MAJ | Un-bannered stale Zbb 25-30 TeV artifact | `scan_outputs/wq_quarkonly_1M_.../` | PENDING |
| M-36 | DOC | MAJ | Un-bannered pre-B1 notebook (Zbb 25-30 TeV) | `notebooks/wq_quarkonly_explore.ipynb` | PENDING |

## §8.3 Revisions (first-pass claims corrected — track for documentation)

| ID | Disposition | Note |
|----|-------------|------|
| C-3 anchors | ALREADY-CORRECT (anchors) | S,T anchors plausibly legit recent fit; only 15.96-vs-18-20 TeV floor mismatch stands (kept under C-3/EW) |
| CR001 anchor | ALREADY-CORRECT | CMS-B2G-25-009 real (5.5 TeV observed); provenance concern withdrawn |
| K004 anchor | ALREADY-CORRECT | NA62 Moriond-2026 9.6e-11 real; residual is policy only (see M-4) |
| Seesaw factor-2 | ALREADY-CORRECT | self-consistent (g₀²=2f²); documentation gap only |
| NuFIT 6.1 | DOC | plausible; pin table values before publication |

## Minor findings (§4, §8.4) — triaged in a dedicated pass after majors

Enumerated and dispositioned in [`MINOR_FINDINGS_LEDGER.md`](MINOR_FINDINGS_LEDGER.md) (created during the minors pass).

---

## Cycle log

### Cycle 1 — Lane C bundle (C-1,C-2,C-5,C-6,M-8,M-27,M-28,M-13C) — VERIFIED
- Codex research (scoping): all 8 confirmed with file:line + corrected constants; noted numeric
  cancellation between C-5 (×8.48), M-27 (÷4), M-13 (×1/6) → must fix coordinately.
- Codex fix L1 (commit `26b181a`): RG direction+transpose, Fierz −2, ME ÷4, √(2L) volume + g_s* unfreeze,
  ξ_g=2.4487, C-6 negative-slope orientation, C-2 LR/C4 capability. Sanity: C1=0.729, C4=3.543.
- Codex fix L2 (commit `676c131`): de-circularize RG test (independent BBL Eq-3.94 oracle),
  reconcile convention/contract-ID test pins, regenerate goldens/results, quarantine banner + residual markers.
- Claude audit: **APPROVE** — all 7 checks PASS, numerics executed independently (C1=0.729, C4=3.543,
  Fierz round-trip, ME ratios =4.0, √(2L)=8.48, ξ_g=2.4487); RG oracle genuinely independent; isolated to Lane C; 66/66 green.
- Residuals (documented, pending arXiv:0710.1869 not in-repo): C-6 exact Table-I affine coefficients;
  C-2 default RH-down alignment model; M-28 strict-paper Eq.(3) re-extraction. Lane C stays NON-PRODUCTION-quarantined.

### Cycle 2 — LFV normalizations (M-29, M-30) — VERIFIED (commit `23c95e8`)
- Codex research+fix: M-29 mu->3e interference 2sqrt(2)e -> -8e (KO RMP 73 Eq 2.14);
  M-30 tau->3l/tau->l-gamma width->BR via PDG tau leptonic BR (0.1782 e / 0.1739 mu). Muon unchanged.
- Claude audit: **APPROVE** — verified -8e magnitude (8-vs-16 resolved via contact (2,1) structure),
  subtle sqrt(leptonic_br) envelope scaling correct, muon bit-identical, test pins literature-anchored. 137 tests green.

### Cycle 3 — Mechanical (M-9, M-31, M-32, M-33) — VERIFIED (commit `75bc84d`)
- Codex research+fix: M-9 top dipole 0.5->0.25; M-31 T003-T006 BR=old/2 + fixtures; M-32 B->K* kernel
  q^2+factor-2 + (C±C') split, independent Gauss-Legendre oracle; M-33 R(D) proxy dimensionless.
- Claude audit: **APPROVE** — M-9 confirmed by auditor's own Dirac trace (ratio exactly 2), M-31 bit-exact
  halves, M-32 transversity structure correct (SM proxy 6.87e-7 = ~73% PDG, acceptable), M-33 dimensionless,
  B025 unaffected. 93 tests green.

### Cycle 4 — C-7 compare_2007_vs_modern retraction — VERIFIED (commit pending)
- Codex research+fix: chose option (b) RETRACT. Modern B/D ratios are hadronic |M12|/GeV budgets; the only
  legacy "2007" B/D numbers are dimensionless operator-weights (~1e6 apart) so a same-convention all-system
  comparison is NOT reconstructable. Removed RESCALE/ratios_2007/accepted_2007, added SUPERSEDED/CORRECTED
  banners, guarded plot_publication_figures use_2007=True to raise, removed publication BOUND_RATIOS.
  9.4x D0 "tightening" retracted as a units artifact.
- Claude audit (direct): retraction complete, no dangling refs to removed symbols anywhere, banners honest,
  guard raises, both scripts compile. VERIFIED.

### Cycle 5 — Lane C straggler tests — VERIFIED (commit `911b879`)
- Cycle-1's pytest -k filter missed 9 test_paper_* files -> 16 tests failed on the already-audited Lane C fixes.
- Codex reconciled: (A) re-pins to corrected ME formula IDs / sqrt(2L) couplings / affine coefficients;
  (B) test_paper_observables converted from asserting the OLD C-2 bug (reject LR) to asserting corrected
  behavior (LR flows + independent hand-check M12 = <Q4>*C4/(2 m_K)). Claude spot-audit: category-B genuine
  corrected-physics, not re-pin. 226 paper tests green.

### Cycle 6 — epsilon_K budget M-1/M-2/M-4 — VERIFIED (commit `dbae398`)
- Codex research proposal (orchestrator steered: mandate-strict sigma = sqrt(sigma_BGS^2+sigma_exp^2), SM-choice
  diagnostic-only) then implement. One shared EpsilonKBudgetPolicy consumed by core deltaf2 + catalog K001.
- Claude audit: **APPROVE** (executed) - mandate fidelity (1.80336e-4, no SM-choice inflation), split-brain
  killed (same policy object, YAML cross-validated), sign-aware per K004, re-pins mechanical (20.67/5.60=3.69),
  floor band arithmetic correct. Production eps_K floor now a BAND (raise-edge ~3.3-3.6, central ~6.3-7,
  lower-edge ~4.8-5.4 TeV). 107 tests green. NOTE: 2 Lane-A repro scripts still use own central budget (fold into P5b).
