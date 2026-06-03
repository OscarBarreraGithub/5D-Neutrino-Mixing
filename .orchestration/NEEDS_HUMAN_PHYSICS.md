# Needs Human Physics — structured list of flagged gaps

Running list of items where a constraint's implementation could NOT be fully
grounded in established physics with the inputs currently available, and uses a
documented proxy / placeholder pending human review. Each constraint is
self-contained (one isolated class), so these gaps do NOT break other
constraints — they are contained to the listed constraint's `evaluate()` and
recorded in its docstring + result diagnostics.

Maintained per-constraint as the rebuild proceeds; presented as the final
deliverable. See `.orchestration/runs/<ID>/` for the full agent trail of each.

---

## ★ EXECUTIVE SUMMARY — what needs human physics input (read this first)

The rebuild is COMPLETE: **103 constraints (95 PRIMARY + 8 SECONDARY), full suite 1054 passed, 0 registry import failures.** Every constraint computes a rigorous, validated SM (or SM≈0 for forbidden processes) side and a defensible experimental budget. The human-input items below are confined to the **new-physics (RS) matching** and a handful of **incalculable-SM observables** — they are documented proxies, never fabricated numbers, and each is flagged in the constraint's docstring + result `diagnostics["needs_human_physics"]`. Of 103 constraints, **98 carry a proxy caveat; 5 are fully rigorous end-to-end** (K001 ε_K, K002 Δm_K, B001 Δm_d, B003 Δm_s, C001 D⁰ mixing — the ΔF=2 sector, because the KK-gluon coupling IS on `ParameterPoint`).

The gaps cluster into **6 categories**:

**G1 — Cross-cutting RS electroweak-coupling matching (the dominant gap; ~majority of the 98).**
`ParameterPoint` exposes quark mass-basis couplings + the KK-gluon mass, but NOT the full RS EW sector (KK W/Z masses & profiles, Z-fermion coupling shifts, Z′ couplings, lepton/neutrino bulk profiles). So every EW-mediated NP contribution (rare semileptonic K/D/B decays, b→sνν, b→sℓℓ, Z-pole, dipoles) uses a documented coupling proxy for the **NP part only** (SM part is rigorous). *Decision needed:* build the shared RS-EW coupling machinery + extend `point_builder`/`ParameterPoint` (unlocks rigorous NP for ~all of these at once), or accept the bounded proxies. Affects: most of beauty, kaon, charm, top_higgs_ew rare/Z families.

**G2 — Lepton-flavor-violating / off-diagonal lepton couplings not on `ParameterPoint`.**
All LFV constraints supply the off-diagonal charged-lepton coupling as an explicit proxy: L001/L002/L005/L006/L008/L009/L010 (μ→eγ, μ→3e, μ–e conv, muonium, τ→eγ, τ→3μ, …), K019/K020/K021 (K→(π)eμ), Z→eμ/eτ/μτ LFV (T015–T017), Higgs LFV (T018–T020). *Decision needed:* RS lepton-sector flavor structure on `ParameterPoint`.

**G3 — Collider σ×BR / recast gap (all of collider_rs, CR001–CR014).**
Mass-vs-published-limit recasts use a single-number mass proxy; the full σ×BR, branching surface, width dependence, production mix, interference, and detector-acceptance/limit-curve recast are flagged. *Decision needed:* proper signal recasts (or accept mass-proxy bounds).

**G4 — CKM-phase machinery absent in the core.**
B002 (S_ψKs / sin2β) and B004 (φ_s) take the SM reference phase from the catalog yaml — no core CKM-phase computation exists. *Decision needed:* add CKM-phase machinery, then derive these in-core.

**G5 — Incalculable-/dirty-SM honest stubs (INFO/SOFT, NON-VETOING) — dual NEEDS-HUMAN (SM side AND RS side).**
These intentionally do not veto; they record the catalogued bound + flag both the missing SM calculation and the RS matching: **K003** (ε′/ε), **K013** (radiative kaon), **B032/B033/B034** (nonleptonic / φ_s^sss), **C003** (charm direct CPV), **E004/E006/E007/E008/E009** (EDMs: atomic/hadronic CP-odd matrix elements + Weinberg 3-gluon), **CR011** (needs full recast). Plus charged-current/G_F-matching pulls **EW002** (CKM 1st-row unitarity), **EW003** (|Vcb|/|Vub|), **K018** (|Vus| K_l3). *Decision needed:* genuine SM hadronic/ChPT/lattice inputs (ε′/ε, EDM matrix elements) or RS charged-current matching.

**G6 — EDM CP-odd loop machinery + exclusive form factors / exclusive-SM predictions.**
E001 (e-EDM) needs a genuine RS one-loop CP-odd dipole; E004–E009 need hadronic CP-odd matrix elements (overlaps G5). Exclusive radiative B013 (Bs→φγ) and B014 (B→ργ/ωγ) need real exclusive form factors, photon-helicity observables (A_Δ/S_φγ), and a genuine exclusive-SM theory BR prediction (none in catalog → their HARD budgets are honest *measurement-consistency* bands, not theory-vs-data room).

**Bottom line / highest-leverage action:** building the shared **RS electroweak-coupling sector** on `ParameterPoint`/`point_builder` (G1) would convert the single largest block of proxies into rigorous NP in one move, and G2 is the lepton-sector extension of the same work. G4 (CKM phase) is a small, self-contained core addition. G5/G6 are genuinely hard SM-theory inputs (lattice/ChPT/loop) best sourced from the literature per observable. None of these block the catalog from running or affect the 5 fully-rigorous ΔF=2 anchors or any constraint's SM side.

---

## ★★ DOES EACH GAP *ACTUALLY* NEED HUMAN INPUT? (honest triage — added 2026-06-02)

The "98 constraints carry a proxy" figure is misleading: most of it is **work we do (build G1)**, not human intervention. Triage:

**(I) NOT human input — we build it (in the active plan):**
- **G1 RS-EW couplings + G2 lepton couplings** — the bulk. Builds the RS Z-coupling shift / KK gauge sector from the existing `f_IR/f_UV`+`bessel` machinery onto `ParameterPoint`. Closes ~24–26 constraints to fully rigorous, ~13–16 to partial. (See `rs_ew_sector_design.md`.)
- **G4 CKM phase (B002 sin2β, B004 φ_s)** — trivial arithmetic from the CKM matrix already in-repo. Small core addition. *(Reclassified: NOT a human decision.)*
- **Exclusive form factors (B013 Bs→φγ, B014 B→ργ/ωγ)** — standard literature values (LCSR/lattice, e.g. Bharucha–Straub–Zwicky). Wire in with a cited source; human only sign-off on the source, no calculation. *(Reclassified: NOT a human decision, modulo source approval.)*

**(II) GENUINELY needs human input — with the reason (TRICKY-PHYSICS vs YOUR-DECISION):**
- **Collider σ×BR recasts — all CR* (14). → YOUR DECISION (scope/effort), not unknown physics.** Rigorous recast needs detector-level signal sim + experiment efficiency/acceptance maps; standard scan practice uses the published *mass* limit as a conservative bound (what we do). Decide: accept conservative mass bounds (recommended) vs invest in full recasts (large infra, marginal gain for a scan).
- **ε′/ε — K003. → TRICKY PHYSICS (unsettled SM).** Lattice (RBC-UKQCD) vs large-Nc dual-QCD disagree on the SM central value AND error; making it vetoing forces endorsing one disputed value. Decide: keep advisory/non-vetoing, or adopt a specific reference + uncertainty.
- **Charm direct CPV / nonleptonic — C003, B032–B034. → TRICKY PHYSICS.** Long-distance hadronic matrix elements, no controlled expansion ("theoretically dirty"); SM not reliably calculable. Same shape as ε′/ε.
- **EDMs — E001, E004, E006–E009 (~6). → TRICKY PHYSICS *and* YOUR MODELING CHOICE — the one set G1 does NOT close.** (a) RS contribution = multiloop CP-odd dipole needing the **complex 5D-Yukawa phase structure**, which is a model input *you specify*; (b) quark/gluon CP-odd operators → atomic/neutron EDMs needs nuclear matrix elements (Schiff moments) with large literature uncertainties.
- **(Minor, mostly closed by G1)** Charged-current/G_F RS matching for the data-pull constraints **EW002, EW003, K018**: the SM-vs-data tension is already rigorous (that IS the constraint); only the RS *shift* is a proxy, closed by G1's charged-current piece. Residual: whether the leading RS W-KK correction suffices.

**Net after G1:** genuine human-input residual ≈ **EDMs (~6, hard physics + your CP-phase assumptions), ε′/ε + charm/nonleptonic CPV (~5, unsettled SM — adopt-or-keep-advisory), collider recasts (14, scope decision — recommend accept).** All are honest non-vetoing/advisory today, so they never silently exclude points.

---

The per-constraint detail (with validated SM numbers, budgets, and exact proxy descriptions) follows below, wave by wave.

---

## Cross-cutting infrastructure gap (affects many EW/rare/LFV constraints)

**RS electroweak coupling inputs are not on `ParameterPoint`.** The scaffold's
`ParameterPoint` currently exposes quark mass-basis couplings + KK gluon mass.
Rigorous RS *new-physics* matching for electroweak-mediated processes (rare
semileptonic/neutrino decays, LFV, Z-pole, EDMs) needs the full RS EW sector:
KK W/Z masses & profiles, Z-fermion coupling shifts, Z′ couplings, and lepton/
neutrino bulk profiles. Until `point_builder` is extended to produce these,
such constraints implement a rigorous **SM** part + a **documented proxy** for
the RS-NP part.
**Human decision needed:** build the shared RS-EW coupling machinery +
extend `ParameterPoint`/`point_builder` (enables rigorous NP for these
families), or accept bounded proxies. (Independent of per-constraint correctness.)

---

## Per-constraint flags

### K004 — BR(K⁺→π⁺νν̄)   [status: committed @ 4ff15a3]
- **Rigorous:** SM short-distance BR (Buras/BGS parametrization; validated
  BR_SM ≈ 8.50×10⁻¹¹ vs literature ≈ 8.6×10⁻¹¹).
- **Proxy / needs human input:** RS NP contribution uses a Z-like effective
  coupling proxy `X_NP=(Δ_sd^L+Δ_sd^R)·Δ_ν/(g²M_KK²)`. The full RS electroweak
  KK/Z/Z′ tower + neutrino-coupling matching is not available on
  `ParameterPoint` (see cross-cutting gap above).
- **Trail:** `.orchestration/runs/K004/`.

### K005 — BR(K_L→π⁰νν̄)   [status: committed @ ed736c3]
- **Rigorous:** SM short-distance BR (purely CP-violating top term; validated 2.95×10⁻¹¹).
- **Proxy / needs human input:** RS NP = same Z-like proxy as K004 (Im part); full RS-EW matching absent (cross-cutting gap).
- **Trail:** `.orchestration/runs/K005/`.

### B002 — S_ψKs / sin2β   [status: committed @ 80d4224]
- **Rigorous:** NP mixing-phase shift from the existing running ΔF=2 complex M₁₂.
- **Needs human input:** the SM reference phase 2β is taken from B002.yaml (β=22.63°) — no CKM-phase
  computation exists in the core. Substitute a core-derived 2β when CKM-phase machinery is added.
- **Trail:** `.orchestration/runs/B002/`.

### C002 — CP violation in neutral charm mixing   [status: committed @ b85f482]
- **Rigorous:** NP CP amplitude |Im(M₁₂^NP)| from running ΔF=2 D⁰ mixing; budget from HFLAV φ_M interval.
- **Needs human input:** no grounded SM long-distance charm phase / Γ₁₂; CP-phase verdict omits SM-side
  long-distance contribution (charm is theoretically dirty).
- **Trail:** `.orchestration/runs/C002/`.

### B004 — φ_s (B_s→J/ψφ)   [committed @ 6034c96]
- **Needs human input:** SM φ_s = −2β_s taken from B004.yaml — no CKM-phase computation in the core.

### EW002 — CKM first-row unitarity   [committed @ d640d2e]
- **Needs human input:** RS charged-current / G_F matching not on ParameterPoint; only the SM-vs-data tension is rigorous (NP shift = 0).

### EW003 — |V_cb|/|V_ub| incl-vs-excl   [committed @ 4947894]
- **Needs human input:** charged-current EW/WET RS matching absent; rigorous part is the data/theory pull only.

### L001 — μ→eγ   [committed @ 99063a2]
- **Needs human input:** lepton-sector RS dipole couplings not on the quark-only ParameterPoint → dipole coefficient is a documented proxy; missing-input case returns explicitly `evaluated=False` (not a pass).

### B005 / B011 / B022 — B-meson rare decays   [committed be82e38 / f746d56 / e4c64d2]
- **Rigorous:** SM rates validated (B_s→μμ 3.65e-9; B→X_sγ 3.40e-4 with LL RG running; B⁺→K⁺νν̄ 5.58e-6 with correct long-distance split).
- **Needs human input:** RS NP is a documented proxy in each shared module
  (`rare_b_dilepton` C9/C10 penguin, `bsgamma` C7/C8 dipole, `rare_b_nunu` Z-like) —
  full b→s EW-penguin/dipole RS matching needs EW KK couplings not on ParameterPoint (cross-cutting gap).
- These shared modules are reused by the downstream fan-out (B006/B012/B015-B019/B021/B023); the NP-proxy caveat propagates to all of them.

### Wave-6 family pioneers (K006, C004, T010, T001, B016)   [committed df0eb34/cf65c4b/f8ad10d/72a22e4/3cafb41]
- **Rigorous:** SM sides validated (K_L→μμ SD 0.82e-9 vs Gorbahn-Haisch; Z→bb̄ R_b 0.2156/A_b 0.935; B⁺→K⁺μμ 5.75e-7; t→cZ SM negligible; D⁰→μμ SD≈0).
- **Needs human input (NP proxies, all from the cross-cutting EW-coupling gap):**
  `rare_kaon_dilepton` s→dℓℓ Z/penguin; `rare_charm_dilepton` c→uℓℓ; `zpole` Zbb-coupling-shift (classic RS Z→bb);
  `top_fcnc` top-Z FCNC coupling; exclusive `rare_b_dilepton` C9/C10 (B016 budget carries a 30% proxy-theory inflation).
  These propagate to all downstream consumers of each module.

### Wave-7 (K008, C005, T012, T002, B015)   [committed — see git log]
- **Rigorous:** SM sides validated (K_L→π0ee direct-CP y7V/y7A vs KTeV; Z→cc̄ R_c 0.1721/A_c 0.668; B→X_sℓℓ[1,6] 1.62e-6; t→uZ negligible; D⁰→ee SD≈0).
- **Needs human input:** RS NP proxies inherited from the shared modules (s→dℓℓ y7V/y7A, c→uℓℓ, Zcc-shift, top-Z FCNC, C9/C10) — all the cross-cutting EW-coupling gap.

### Wave-8 (K009, C007, T015, T003, B017)   [committed — see git log]
- **Rigorous:** K_L→π0μμ y7V/y7A SM; B→K(*)ℓℓ R_K proxy SM=1; t→cγ pure-NP; Z→eμ LFV pure-NP; D⁺→π⁺μμ SD.
- **Needs human input:** RS NP proxies (s→dℓℓ y7V/y7A, c→uℓℓ semileptonic, off-diagonal Z-eμ LFV coupling,
  top-photon dipole, C9/C10) — cross-cutting EW-coupling gap. C007 is a FULL-q² proxy (resonance/LD windows
  not applied) + 8.96% form-factor normalization uncertainty.

### Wave-9 (K010, C006, T004, T016, B018)   [committed — see git log]
- RS NP proxies inherited from shared modules (s→dℓℓ K_S a_S w/ sign-envelope, c→uℓℓ LFV eμ coupling,
  top-photon dipole t→uγ, off-diagonal Z-eτ LFV, lepton-non-universal C9/C10 for R_K) — cross-cutting EW-coupling gap.

### Wave-11 pioneers (B009, K018, CR001, E001, L002)   [committed — see git log]
- **Rigorous:** B⁺→τν tree SM 8.63e-5; |V_us| K_l3 extraction + unitarity pull; KK-gluon mass vs tt̄ limit; e-EDM SM≈0; μ→3e dipole(L001)+contact+interference.
- **Needs human input:** charged-current RS proxy (B009); CR σ×BR recast proxy (CR001); RS one-loop CP-odd dipole for EDMs (E001 — real loop machinery); μ→3e Z/box + dipole-contact relative PHASE (L002 uses conservative constructive envelope); K018 RS CKM shift ~0.

> NOTE: waves 12–19 per-constraint flags are recorded in their commit messages + `.orchestration/runs/<ID>/`; a consolidated end-of-project pass will fold them into this file when 95/95 is reached.

### Wave-20 (L006, CR010, E009, K019, B007)   [committed a6f0d74 / 75afcfc / 45e9477 / 6e329c5 / ca580ca]
- **L006 — muonium→antimuonium P(M→M̄):** Rigorous pure-NP proxy P=P_limit·|G_C/G_F|²/(G_C/G_F)²_limit vs MACS/PSI bound P<8.3e-11. **Needs human input:** full RS ΔL_μ=−ΔL_e=2 four-lepton matching not on ParameterPoint (effective coupling supplied as documented proxy).
- **CR010 — VLQ (T,B) pair production:** Rigorous mass-vs-limit recast (m_T=m_B=M_KK proxy) vs ATLAS/PDG 2025 ≥1.37 TeV. **Needs human input:** σ×BR / branching-surface / T,B spectrum / widths / acceptance / full limit-surface recast (cross-cutting collider recast gap).
- **E009 — Weinberg three-gluon CP-odd operator:** INFO/non-vetoing; records neutron-EDM-derived |C_6|<1.2e-11 GeV⁻² (w<4.1e-11 GeV⁻², n-EDM 1.8e-26 e·cm). **Needs human input (dual):** hadronic CP-odd gluonic matrix element AND RS CP-odd gluonic matching on ParameterPoint.
- **K019 — K_L→eμ LFV (secondary):** Rigorous pure-NP bound (SM≈0) BR<4.7e-12 (BNL E871) via rare_kaon_dilepton κ_μ adapted to the unequal-lepton (e,μ) two-body rate. **Needs human input:** off-diagonal charged-lepton neutral-current RS matching (eμ coupling supplied as documented proxy; cross-cutting EW-coupling gap).
- **B007 — B0/Bs→ee (secondary):** Rigorous SM (Bs→ee 8.54e-14, B0→ee 2.45e-15, m_e² helicity suppression) vs limits BR(Bs→ee)<9.4e-9, BR(B0→ee)<2.5e-9. **Needs human input:** RS C9/C10 b→s(e e) NP matching (documented proxy; cross-cutting EW-coupling gap).

### Wave-21 (CR012, B008, K020, T014)   [committed 69648e2 / 6ece7d8 / 966af3b / 0183016]
- **CR012 — diboson high-mass resonance (spin-1):** Rigorous mass-vs-limit recast vs HVT-B W'→WZ ≥4.4 TeV. **Needs human input:** σ×BR / branching-surface / width / production-mix / interference / acceptance / limit-curve recast (cross-cutting collider recast gap; documented single-number mass proxy).
- **B008 — Bs/B0→ττ (secondary):** Rigorous SM (Bs→ττ 7.74e-7, B0→ττ 2.19e-8; heavy-τ, β_τ=0.749 phase space) vs BR(Bs→ττ)<6.8e-3, BR(B0→ττ)<2.1e-3. **Needs human input:** RS C9/C10 b→s(ττ) NP matching (documented proxy; cross-cutting EW-coupling gap).
- **K020 — K⁺→π⁺eμ LFV semileptonic (secondary):** Rigorous pure-NP three-body rate (SM≈0) with K→π vector form factor f₊(q²) + genuine q²+angular phase-space integral vs BR<1.3e-11 (E865/PDG). **Needs human input:** off-diagonal e-μ lepton coupling not on ParameterPoint (documented proxy; cross-cutting EW-coupling gap).
- **T014 — FCNC Z→bs/bd/sd (down-sector, secondary):** Rigorous pure-NP off-diagonal Z width; B=Γ_FCNC/(Γ_Z,total^SM+Γ_FCNC) vs B<2.9e-3 (ECFA 2025). **Needs human input:** RS FCNC-Z off-diagonal down-sector coupling matching (documented (m_Z/M_KK)² overlap proxy; cross-cutting EW-coupling gap).

### Wave-22 (CR013, K021, B013)   [committed cc5e457 / 481369c / 9a6db0b]
- **CR013 — diphoton high-mass resonance (spin-0/2):** Rigorous spin-2 KK-graviton mass (m_G=3.8317·Λ_IR) vs CMS 2024 RSG γγ ≥4.8 TeV (k̃=0.1). **Needs human input:** σ×BR / branching-surface / k̃ / width / acceptance / spin-0 interpretation recast (cross-cutting collider recast gap).
- **K021 — K_L→π⁰eμ LFV semileptonic (secondary):** Rigorous pure-NP neutral-mode three-body rate (SM≈0; K_L/π⁰ kinematics, neutral K_l3 f₊(0), charge-sum 2) vs BR<7.6e-11 (KTeV/PDG). **Needs human input:** off-diagonal e-μ lepton coupling proxy AND K_L CP/charge-orientation matching (cross-cutting EW-coupling gap).
- **B013 — Bs→φγ exclusive radiative (secondary):** Rigorous C7/C7' dipole-power NP with LL RG running (zero-NP recovers measured 3.4e-5); budget is an HONEST measurement-consistency band (HFLAV vs PDG). **Needs human input (multiple):** RS C7/C7' dipole proxy; exclusive form factors; A_Δ / S_φγ photon-helicity observables; AND a genuine exclusive-SM theory/form-factor BR prediction (none exists in the catalog → the budget is measurement-consistency, not theory-vs-data room).

### Wave-23 (CR014, B014)   [committed 6910282 / 67b2f46]   *** FINAL WAVE ***
- **CR014 — four-top top-philic vector (PRIMARY, 95th):** Rigorous mass-vs-limit recast vs CMS four-top Z' excluded-up-to 850 GeV (Γ/m=50%). **Needs human input:** σ×BR / width-dependence (Γ/m=50%) / top-philic couplings / four-top acceptance / SM-four-top background recast (cross-cutting collider recast gap; documented low-mass width-dependent mass proxy).
- **B014 — B→ργ/ωγ exclusive b→dγ (secondary):** Rigorous C7/C7' dipole-power NP with LL running, CKM-suppressed |V_td/V_ts|²=0.0457; zero-NP recovers measured BRs; honest measurement-consistency band. **Needs human input:** RS C7/C7' dipole proxy; exclusive B→ρ/ω down-sector form factors; weak-annihilation/spectator/helicity; missing exclusive-SM theory/form-factor prediction.
