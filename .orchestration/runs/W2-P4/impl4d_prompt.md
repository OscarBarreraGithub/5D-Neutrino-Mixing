# W2 PHASE 4 — SUB-STEP 4d: b→sνν / s→dνν rewire (codex, gpt-5.x xhigh). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. Implement sub-step 4d of the DUAL-APPROVED plan `.orchestration/runs/W2-P4/plan.md` (step 14), using the 4a builder (committed 711b56d: `quarkConstraints/rs_semileptonic_wilsons.b_to_s_nunu` / `.s_to_d_nunu`) and the established rewire+degradation pattern. First a SHORT plan, then implement. Codex + Opus dual-review (both must APPROVE). SCOPE = neutrino-pair channels ONLY: B022, B023, K004, K005. This is the LAST Phase-4 sub-step.

REWIRE (quark-FCNC × active-neutrino):
- **B022, B023** (B→Kνν̄ / B→K*νν̄, b→sνν̄): replace the one-Z-like proxy `Δ_q Δ_ν/(g_SM² M_KK²)` with `rs_semileptonic_wilsons.b_to_s_nunu` mapped as `X_NP = C/g_SM²` (LH active-ν), into the existing `rare_b_nunu` core's X_NP input — NO `_wilson_prefactor`, NO second 1/M_KK².
- **K004, K005** (K⁺→π⁺νν̄ / K_L→π⁰νν̄, s→dνν̄): replace the proxy with `rs_semileptonic_wilsons.s_to_d_nunu` mapped `X_NP=C/g_SM²` into the existing `rare_kaon_snd` core. PRESERVE the SM short-distance (Buras/BGS) + charm/long-distance (P_c, δP_{c,u}) treatment and the NA62/KOTO/HPQCD anchors EXACTLY.
- **v1 PHYSICS (note — these DO bite)**: unlike lepton-LFV, the νν channels get NONZERO NP in v1 because the QUARK FCNC δg (b→s, s→d) is nonzero (from 3a/3c, c_Q non-degeneracy) × the active-ν coupling. So these become RIGOROUS and VETOING. The old proxy NEEDS-HUMAN (RS EW coupling) is RESOLVED. (Active-ν is LH, uses c_L; Majorana=Dirac for light νν.)
- SM rates / anchors / long-distance UNCHANGED.

GRACEFUL DEGRADATION: rigorous `rs_semileptonic_wilsons` when present; ABSENT ⇒ non-vetoing evaluated=False + missing_extra (no crash/fake pass).

TESTS: rigorous-path (point via `build_from_rs_ew_inputs` ⇒ nonzero NP X_NP; the BR shifts from SM; cross-check the X_NP from the quark FCNC × ν coupling independently of the adapter); SM-limit (universal-c ⇒ quark FCNC=0 ⇒ X_NP=0 ⇒ recover committed SM BR); absent-path; Majorana=Dirac. Replace proxy-only tests — enumerate, no silent coverage loss. `python -m pytest tests/ -q` stays green.

CONSTRAINTS: physics ONLY via adapters (rare_b_nunu, rare_kaon_snd); ConstraintResult numeric fields real finite floats (complex in diagnostics); touch ONLY B022/B023/K004/K005 + their adapters/tests.

OUTPUT (<=14 lines): short plan; rewired files + IDs; one rigorous nonzero X_NP→BR example (e.g. K004 or B022) + SM-limit recovery + absent-path + Majorana=Dirac; old proxy NEEDS-HUMAN resolved; test-count change; pytest counts. End with: P4D-DONE.
