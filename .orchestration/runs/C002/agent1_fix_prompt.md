You are agent1 fixing C002 (CP violation in charm mixing). Reviews: .orchestration/runs/C002/agent2_out.md (PHYSICS-NEEDS-FIXES) and agent3_out.md (CODE-NEEDS-FIXES). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing.
FIX (physics):
1. Budget is too loose. Use the HFLAV φ_M = Arg(M12) 95% interval (~[-1.48, 1.35] deg) to set the CP budget → |Im M12^NP| <= ~8.5e-17 GeV, instead of the looser |q/p|/φ_D rooms. Source the interval from C002.yaml / its snapshot. Keep severity HARD.
2. The diagnostic q_over_p_phase_proxy_degrees = -arg(M12^NP) is misleading (a real-negative NP amplitude reports -180 deg while the CP-odd Im(M12^NP) is zero). Either ground it properly or clearly relabel it as a raw-argument diagnostic; the pass/fail VERDICT must keep using |Im(M12^NP)| (that part is correct).
FIX (code):
3. Make the test CP oracle INDEPENDENT: compute expected complex M12 via quarkConstraints.deltaf2 core directly (_evolve_wilsons + compute_m12_np; evaluate_d0_mixing_with_running for abs), not via the adapter C002 calls. Report numbers.
4. (SHOULD-FIX) Route q/p, phi_D anchors through load_anchor instead of load_pdg_block/find_block.
Keep the NEEDS-HUMAN-PHYSICS caveat (no grounded SM long-distance charm phase / Gamma12) — accepted.
IGNORE the "isolation/worktree has other files" finding — parallel-wave artifact. C002 must only touch C002.py, its test, and the append-only deltaf2 adapter D0 helper.
Run `python -m pytest tests/constraints/primary/charm/test_C002.py -q` and `tests/constraints/ -q` green. OUTPUT <=12 lines: changes, new budget value, cross-check numbers, pytest counts.
