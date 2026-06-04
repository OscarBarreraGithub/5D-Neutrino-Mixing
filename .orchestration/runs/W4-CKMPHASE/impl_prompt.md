W4 — CKM PHASE in-core (codex, gpt-5.x xhigh). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. Implement W4: compute the SM CKM phases (2β for B002, φ_s=-2β_s for B004) IN-CORE from the CKM matrix, instead of reading them from the yaml. First a SHORT plan, then implement code+tests. Codex + Opus dual-review (both must APPROVE).

CONTEXT: B002 (S_ψKs/sin2β) and B004 (φ_s, Bs->J/psi phi) currently take the SM reference phase from their yaml (B002: β=22.63deg; B004: φ_s^SM=-2β_s) and flag NEEDS-HUMAN "no CKM-phase computation in the core". The repo HAS the CKM matrix (the quark fit produces V_CKM = U_L_u^dag U_L_d; there are also standard Wolfenstein/CKM inputs in the constants). Compute the unitarity-triangle angles in-core.

BUILD:
- A small, well-tested core function (e.g. in quarkConstraints/ckm_extraction.py or a new ckm_phases.py) computing:
  - `beta = arg(-(V_cd V_cb*)/(V_td V_tb*))` ; `sin(2 beta)` ; report 2β.
  - `beta_s = arg(-(V_ts V_tb*)/(V_cs V_cb*))` ; `phi_s = -2 beta_s`.
  from a CKM matrix input (use the repo's standard CKM / Wolfenstein values; if the constraint has access to the fit V_CKM, prefer that, else the PDG Wolfenstein anchor — be explicit about which).
- Validate against PDG: 2β ≈ 45.5deg (sin2β ≈ 0.70), φ_s^SM ≈ -0.0368 rad (-2β_s, β_s≈0.018-0.022 rad). State the computed numbers.
- Wire B002/B004 to use the in-core SM phase instead of the yaml literal; REMOVE the "no CKM-phase in core" NEEDS-HUMAN flag (the SM phase is now rigorous). The NP mixing-phase shift (from the running DeltaF2 complex M12) is UNCHANGED. Keep anchors via load_anchor; the in-core phase replaces only the SM reference.

TESTS: the computed 2β/sin2β and φ_s match PDG within stated tolerance; B002/B004 SM-limit (no NP) reproduces the committed SM S_ψKs/φ_s; an independent recompute of the angle from the CKM matrix; the NEEDS-HUMAN-for-CKM-phase flag is removed (but any OTHER partial flags kept). `python -m pytest tests/ -q` stays green.

CONSTRAINTS: physics via the core/adapter; numeric fields real floats; touch ONLY the new CKM-phase core + B002/B004 + tests. Do NOT alter the ΔF=2 M12 machinery (separately confirmed correct).

OUTPUT (<=14 lines): short plan; new function + the computed 2β/sin2β/φ_s vs PDG; B002/B004 SM-limit; NEEDS-HUMAN-CKM-phase removed; test-count; pytest counts. End with: W4-DONE.
