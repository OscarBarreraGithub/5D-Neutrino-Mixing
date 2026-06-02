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
