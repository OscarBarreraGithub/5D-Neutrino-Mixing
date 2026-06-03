# Independent design task (codex, gpt-5.x xhigh). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. Produce your OWN independent design — do NOT look for or read any existing `rs_ew_sector_design*.md` (a separate design exists; we want an INDEPENDENT second opinion to cross-review). Do NOT modify production code. WRITE your design to `.orchestration/runs/G1-DESIGN/codex_design.md` and print a tight summary as your final message.

GOAL: design (not implement) the "RS electroweak-coupling sector" so the flavor-catalog constraints that currently use NEW-PHYSICS proxies (the G1/G2 gap) can be re-wired to RIGOROUS RS new-physics, and so a 100M-point parameter scan rigorously constrains the model.

GROUND IN (read these):
- `CLAUDE.md` (conventions: c=M_5/k; f_IR/f_UV overlap formulas; warp params; seesaw; PMNS).
- `warpConfig/baseParams.py` (get_warp_params) and `warpConfig/wavefuncs.py` (f_IR, f_UV) — the existing zero-mode overlap machinery to build on.
- `solvers/bessel.py` (KK tower mass solver) for KK W/Z masses + profiles.
- `derivations/` (LaTeX — the authoritative physics; find RS gauge-KK couplings, Z-coupling shifts, dipoles).
- `flavor_catalog_constraints/point_builder.py` + the `ParameterPoint` definition (note KNOWN_EXTRA_KEYS already lists RS-EW/lepton placeholders no builder fills) + `build_from_quark_couplings`.
- The proxy sites: grep `flavor_catalog_constraints/physics_adapters/` and `quarkConstraints/` for `NEEDS-HUMAN-PHYSICS`/`proxy` — esp. zpole.py (Z-coupling shifts), bsgamma.py + lepton dipole, rare_b_dilepton/rare_kaon_dilepton/rare_charm (C9/C10 + s→dℓℓ Z-penguin), rare_b_nunu, the LFV modules (mu_e_conversion, lfv_three_body, zpole_lfv, higgs_lfv), edm.py. Read 2–3 affected constraints' docstrings for exactly what proxy input they expect.

DELIVER IN `codex_design.md`:
1. **Physics inventory** — the RS-EW quantities the constraints need, each with a computing formula in terms of the EXISTING machinery (c-values, f_IR/f_UV, epsilon, KK masses) and a `derivations/` reference: KK gauge masses (W'/Z'); Z-fermion coupling shifts δg_L^f/δg_R^f (quarks AND leptons; the classic RS Zbb shift); flavor off-diagonal Z couplings (FCNC-Z, Z→ℓ_iℓ_j LFV, s→dℓℓ/b→sℓℓ Z-penguin); charged-current/W shifts (G_F, |Vij| for EW002/EW003/K018); lepton & neutrino bulk profiles (c_L,c_e,c_N→f's) for LFV dipoles / μ-e conversion / EDMs; dipole coefficients (C7, lepton dipole) where tractable — HONESTLY mark genuinely loop-hard pieces.
2. **ParameterPoint schema** — exact new fields/extras (names, types, units, real-vs-complex, complex→diagnostics), reusing KNOWN_EXTRA_KEYS; keep the frozen-container + fail-loud contract.
3. **point_builder logic** — how a new build path computes each field from sweep inputs (c_L, c_E[3], c_N, M_KK/Λ_IR, Yukawa structure) via warpConfig/wavefuncs/bessel; note NEW sweep inputs the current quark-only builder lacks.
4. **Matching map** — table: for EACH affected adapter/constraint family, CURRENT proxy → NEW rigorous input → rigorous formula. Group: Z-pole (T010–T013), Z-LFV (T015–T017), FCNC-Z (T014), Higgs-LFV (T018–T020), b→sℓℓ/s→dℓℓ/c→uℓℓ, b→sνν, dipoles (μ→eγ/τ→ℓγ, b→sγ), μ-e conversion + LFV 3-body, EDMs, charged-current (EW002/EW003/K018). Mark FULLY-rigorous vs still-partial.
5. **Phasing & risk** — implementation order (which sub-block unlocks the most constraints first); what stays NEEDS-HUMAN even after this (genuine loop/lattice items); top correctness risks (coupling normalization vs the existing validated zpole/dipole/rare_* cores; double-counting; basis/PMNS rotation for leptons).
6. **Validation strategy** — how to check against known RS results (Zbb shift sign/magnitude, KK mass spacing) and against the constraints' already-validated SM sides.

Final message (<=25 lines): the inventory headline, proposed ParameterPoint additions, count of constraints becoming fully-rigorous vs partial vs still-human, recommended implementation order, top-3 correctness risks. Confirm the doc path.
