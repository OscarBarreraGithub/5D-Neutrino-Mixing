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
