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
