# DA-1 Discovery Worklog: round_001_full_scope
**Date**: 2026-05-16
**Discovery agent**: DA-1
**Scope**: full catalog

## Current catalog inventory
- Total drafted: 27
- Per family: kaon 7, charm 1, beauty 10, top_higgs_ew 5, charged_lepton 3, edm_neutrino 1
- Total CHECKER-DONE so far: 17
- Plan v1 Section C coverage by family: kaon 7/22, charm 1/12, beauty 10/37, top_higgs_ew 5/22, charged_lepton 3/25, edm_neutrino 1/10.
- Drafted IDs: K001, K002, K003, K004, K005, K006, K013; C001; B002, B005, B009, B011, B015, B017, B018, B025, B032, B033; T001, T002, T007, T010, T018; L001, L002, L007; E001.

I cross-checked the current process tree against plan v1 Section C and the orchestrator NIT-1 decision. The most visible imbalances are charm, charged-lepton LFV, and EDM-adjacent CP. Beauty has the largest absolute number of drafts, but several canonical paired observables are still missing. Kaon coverage has the core CP/mixing and rare-neutrino entries, but leaves lepton-universality and muon rare modes open. PI-seeded or seed-adjacent entries already drafted include K013, B005, B009, B011, and B033; I did not find another unambiguous PI seed that should jump ahead of the imbalance fixes, except for K014 as the alternate `K1 -> pi gamma` interpretation, which should wait for PI confirmation.

## Pre-approved additions (NIT-1)
| proposed_id | family | one-line rationale | plan v1 row? |
|---|---|---|---|
| EW001 | top_higgs_ew / EW-global | S/T/U oblique parameters are broad electroweak precision constraints that can bind new vectorlike or KK sectors even when flavor observables are indirect. | No Section C row; explicitly accepted by NIT-1. |
| EW002 | top_higgs_ew / EW-global | First-row CKM unitarity and the Cabibbo anomaly are high-impact weak-current consistency tests tied to V_ud and V_us inputs. | No Section C row; explicitly accepted by NIT-1. |
| EW003 | top_higgs_ew / EW-global | Inclusive-vs-exclusive \|V_cb\| and \|V_ub\| tensions are central semileptonic flavor inputs and should not be hidden inside individual B rows only. | No Section C row; explicitly accepted by NIT-1. |

## Proposed Wave-4 additions (15-25 entries)
| proposed_id | family | one-line rationale | plan v1 row |
|---|---|---|---|
| C002 | charm | CP violation in D mixing completes the current C001 mixing entry with \|q/p\| and phi_D sensitivity to new up-sector phases. | C002 |
| C003 | charm | Delta A_CP in D0 -> K+K-, pi+pi- is the flagship direct charm-CP observable and fixes the charm imbalance. | C003 |
| C004 | charm | D0 -> mu+ mu- is a clean rare c -> u dilepton upper limit with direct FCNC relevance. | C004 |
| C007 | charm | D+ -> pi+ mu+ mu- adds semileptonic charm rare-decay coverage beyond the purely leptonic C004 channel. | C007 |
| L003 | charged_lepton | Mu-e conversion in Al is the target for Mu2e/COMET and essential for coherent conversion coverage. | L003 |
| L004 | charged_lepton | Mu-e conversion in Au carries the classic strong SINDRUM II bound and anchors present-day conversion limits. | L004 |
| L005 | charged_lepton | Mu-e conversion in Ti provides complementary nuclear-target dependence and historical limits. | L005 |
| L008 | charged_lepton | Tau -> e gamma completes the radiative tau dipole pair with drafted L007 tau -> mu gamma. | L008 |
| L009 | charged_lepton | Tau -> 3mu is the leading four-lepton tau-LFV mode with active LHCb/Belle/Belle II relevance. | L009 |
| L010 | charged_lepton | Tau -> 3e is a cheap companion to L009 and separates e-sector contact LFV from the muonic mode. | L010 |
| L023 | charged_lepton | Neutrino trident production is not LFV, but it is an independently strong muon-gauge/contact-operator constraint. | L023 |
| E004 | edm_neutrino | Neutron EDM is the standard hadronic CP benchmark for quark EDM, chromo-EDM, and Weinberg operators. | E004 |
| E006 | edm_neutrino | Mercury EDM is the main diamagnetic-atom complement to neutron EDM and probes different nuclear/CP combinations. | E006 |
| E008 | edm_neutrino | Quark chromo-EDM bounds are directly useful for translating loop-level RS or KK CP phases into EDM constraints. | E008 |
| B004 | beauty | phi_s in Bs -> J/psi phi is the canonical Bs CP-phase partner to B_s mixing and high priority for Delta B = 2 phases. | B004 |
| B006 | beauty | Bd -> mu+ mu- complements drafted B005 and gives the b -> d leptonic rare-decay analog. | B006 |
| B019 | beauty | R_K* should be paired with drafted R_K and B -> K*ll because it is a central b -> s lepton-universality observable. | B019 |
| B026 | beauty | R_D* is the indispensable companion to drafted R_D in charged-current b -> c tau nu tests. | B026 |
| K012 | kaon | K_S -> mu+ mu- is a rare short-lived neutral-kaon dimuon limit with direct FCNC/long-distance caveats. | K012 |
| K017 | kaon | K -> l nu / R_K adds a clean kaon lepton-universality ratio that is easier to interpret than many rare modes. | K017 |
| K018 | kaon | K_l3 provides V_us and semileptonic kaon input needed to contextualize first-row CKM unitarity. | K018 |

## Near-duplicates / merge candidates

- T001 and T002 are separate t -> cZ and t -> uZ rows, but many LHC searches report tqZ limits in a shared analysis. Keep separate if the catalog wants flavor-index-specific couplings; otherwise a future merged `t -> qZ` entry could reduce duplicated source work.
- B017 and B018 are tightly coupled through B -> K(*)ll phenomenology, but they should stay separate: B017 is angular/branching-fraction heavy, while B018 is a lepton-universality ratio.
- B002 and B033 are both time-dependent CP asymmetry entries, but not duplicates. B002 is the golden sin2beta input; B033 is a penguin-phase control/context observable and should remain contextual unless the PI wants it promoted.
- K013 and possible K014 are the real duplicate risk around the PI's ambiguous `K1 -> pi gamma` seed. Current catalog has only K013, so no merge is needed now.
- Proposed EW002 overlaps with proposed K018, and proposed EW003 overlaps conceptually with B029/B030/B031. Treat the EW entries as global overview constraints and the Section C entries as source-level constituent rows if both are pursued.

## PI escalations (if any)

- Confirm whether K013 alone satisfies the PI's `K1 -> pi gamma` seed, or whether K014 should also be drafted as an alternate radiative-kaon entry.
- Decide whether EW002 should be a standalone global CKM-unitarity entry in addition to K018, or whether K018 should be drafted first and EW002 deferred to synthesis.
- Decide whether EW003 should remain a single inclusive/exclusive tension overview or trigger separate Wave-4 PKAs for B029, B030, and B031.
- For EDMs, confirm whether E004/E006/E008 are intended as catalog-only constraints or likely future hard cuts; this affects how much PKA effort should go into EFT translation rather than value capture.

## Summary recommendation

Dispatch the three NIT-1 EW additions first, then the 21 Section C rows above as Wave-4 PKAs. This wave repairs the largest family imbalances while preserving high-impact flavor physics: charm CPV/rare decays, mu-e conversion, tau LFV, hadronic/atomic EDMs, and key paired B/kaon observables. Defer lower-priority long tails such as broad hadronic tau LFV, rare radiative charm, and additional B invisible/LFV modes until after this wave unless the PI explicitly redirects.

===DA1_DISCOVERY_END===

---

## Closure Addendum (added by cleanup-C11, 2026-05-26)

**Closes**: R11-I2 INFO (unresolved PI-escalation thread) and R12-I4 INFO (Wave-4 partial resolution not back-annotated).

### Resolution log for the four PI escalations (lines 55-60)

| # | Escalation | Resolution | Resolver |
|---|---|---|---|
| (a) | K013 alone vs. K014 alternate `K1 -> pi gamma` | RESOLVED — K013 retained as the active radiative-kaon row; K014 (`K_S -> pi^0 gamma gamma`) recorded as DEFERRED-SCOPE pending an explicit PI choice over K013. | `worklogs/discovery/round_004_addendum_deferred_scope.md:43` (DA-4 Wave-7 addendum). |
| (b) | EW002 standalone CKM-unitarity vs. defer | RESOLVED — EW002 drafted as a standalone NIT-1 global electroweak row in Wave-4 (`flavor_catalog/processes/top_higgs_ew/EW002.{yaml,tex}`). | Wave-4 PKA + WA-v2/CA-v2 chain: `worklogs/writer/wa_w4_ew_v2.md`, `worklogs/checker/ca_w4_ew_v2.md` (commits incl. `52acd5e`, `d6e78b6`). |
| (c) | EW003 single overview vs. split into B029/B030/B031 PKAs | RESOLVED — EW003 retained as a single inclusive/exclusive `|V_cb|`/`|V_ub|` overview row in Wave-4 (`flavor_catalog/processes/top_higgs_ew/EW003.{yaml,tex}`); B029/B030/B031 explicitly recorded as DEFERRED-SCOPE "folded into EW003" by DA-4. | Wave-4 PKA + WA-v2/CA-v2 chain (above); deferred-scope assignment at `worklogs/discovery/round_004_addendum_deferred_scope.md:38`. |
| (d) | E004/E006/E008 catalog-only vs. future hard cuts | PARTIALLY RESOLVED — all three drafted in Wave-4 as catalog rows (`flavor_catalog/processes/edm_neutrino/E004.{yaml,tex}` etc.); hard-cut intent NOT recorded in the yamls and remains open for the PI/synthesis stage. | Wave-4 PKA + `worklogs/writer/wa_w4_kaon_edm{,_v2}.md`, `worklogs/checker/ca_w4_kaon_edm{,_v2}.md`; hard-cut decision deferred to post-paper synthesis. |

### Status

DOC CLOSED for discovery-side bookkeeping. The hard-cut intent for atomic/hadronic EDMs (item (d)) is tracked at the synthesis layer, not on the DA-1 inventory; no DA-side action remains.

Cross-link to dependent worklog: see DA-2 (`round_002_followup.md`) closure addendum for the residual K013/K014 carry-forward and DA-4 (`round_004_convergence.md`) for the final deferred-scope ruling.

