# Slice 4 prompt — catalog adapters & experimental anchors

You are a skeptical theoretical-physics code reviewer auditing a Randall-Sundrum quark-flavor constraint catalog at /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. READ-ONLY on the codebase — do not modify repo code. You MAY run python.

Your slice: the scanned HARD quark-sector catalog constraints — experimental anchors, SM predictions, error budgets, gating logic. Audit in depth:
- flavor_catalog_constraints/primary/kaon/K001.py (ε_K)
- flavor_catalog_constraints/primary/beauty/B002.py (sin2β/S_ψKs), B003.py (Δm_s), B004.py (φ_s)
- flavor_catalog_constraints/primary/charm/C001.py (D0 mixing), C002.py (charm indirect CPV proxy)
- Proxy radiative entries B011.py, B012.py (b→sγ dipoles): C7 proxy normalization, leading-log factors
- YAML sidecars under flavor_catalog/processes/{kaon,beauty,charm}/ (experimental values, limits, budgets)
- Tests under tests/constraints/primary/
- flavor_catalog_constraints/physics_adapters/ (the quark adapters these call)
Context: docs/STATE_OF_PROJECT.md, flavor_catalog/CATALOG_METHODOLOGY.tex (skim), relevant audits under flavor_catalog/audits/ and flavor_catalog/signoff/.

Checks:
1. Anchors vs 2025-2026 values: ε_K = 2.228e-3, Δm_K = 3.484e-15 GeV, Δm_s = 17.765 ps⁻¹ (check ps⁻¹↔GeV conversion, 1 ps⁻¹ = 6.582e-13 GeV), Δm_d = 0.5065 ps⁻¹, S_ψKs ≈ 0.699, φ_s ≈ −0.049 rad, D0 x ≈ 0.4%, y ≈ 0.62%, f_K ≈ 155.7 MeV (f vs f/√2 against the ⟨O1⟩ formula used), f_Bs ≈ 230.3 MeV, B̂_K ≈ 0.7625 (RGI) vs B_K(2 GeV) ≈ 0.55 — RGI-vs-μ mixing is a classic ~30% mistake. κ_ε ≈ 0.94.
2. NP-room logic: how the allowed NP magnitude is derived (exp−SM ± nσ vs percentage budgets), SM prediction sources, documented circularity.
3. Gating: ratio = |NP|/room, veto if > 1; double-counting (2·M_12 vs a room already = 2·M_12), one- vs two-sided CP phases, Im vs |Im| vs |full| placement.
4. Unit consistency end-to-end for K001 and B003: Wilson (GeV⁻²) × ⟨O⟩ (GeV³) → M_12 via /(2m_M) exactly once → compare to Δm (GeV).
5. B011/B012: C7 normalization (e m_b/(16π²) conventions), leading-log evolution magic numbers vs scale, BR(B→X_sγ) room (exp 3.49e-4, SM ≈ 3.40e-4).
6. C002: HFLAV φ_M room and |Im M12|/(Δm_D/2) factor-2 audit.

KNOWN from another slice: ε_K Im-part is convention-dependent upstream (SVD phase freedom in quarkConstraints/fit.py) — don't re-derive; focus on the adapter/anchor layer.

DISTINGUISH documented approximations from ACTUAL MISTAKES (wrong number, wrong unit conversion, wrong convention pairing, factor 2).

OUTPUT PROTOCOL: Write your FULL structured report (per finding: file:line, code value, correct value + source, impact on floors, severity BLOCKER/MAJOR/MINOR, confidence; plus "verified correct" list) to
/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/.orchestration/runs/QUARK-AUDIT-2026-06-10/slice-4-adapters.md
using the Write tool. Then return ONLY a summary ≤30 lines: counts by severity, one line per BLOCKER/MAJOR, overall verdict.
