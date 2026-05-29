1. NIT: Correct amplitude: `B002.py:413-419` uses complex `arg(1 + M12_NP/M12_SM)` for the CP phase, not `abs(M12_NP)`; `B002.py:445-449` records Re/Im and `phi_d_np`. Test point: `M12_NP = 2.985e-15 i GeV`, `phi_d^NP = 0.950 deg`, `S = 0.721883`.

2. NIT: QCD running is actually used: verdict path calls `bd_mixing_m12_np_from_wilsons_with_running(..., mu_had=2.0)` at `B002.py:404-407`; non-phase magnitude path at `B002.py:408-411` is diagnostic only. Running effect for the test point: `|M12|` goes from `4.094e-15` to `2.985e-15 GeV`, factor `0.729`.

3. NIT: Budget is defensible and uncertainty-aware: `B002.py:299-309` builds `sqrt(0.011^2 + 0.01106^2 + 0.008354^2) = 0.0176929`; central SM-exp residual is `3.08e-4`. This is moderately loose but justified by the YAML `beta` and penguin systematics policy.

4. NIT: Anchors match YAML/snapshots: YAML has `sin2beta=0.710 +/- 0.011`, `S(J/psi K_S)=0.712 +/- 0.011`, `beta=22.63^{+0.45}_{-0.44} deg`, `|Delta phi_d| <= 0.68 deg` at `B002.yaml:71-149`, matching snapshots `hflav...txt:10-21` and `frings...txt:14-15`.

5. NIT: Units/severity/notes OK: `M12_SM = DELTA_M_BD_SM/2 = 1.8e-13 GeV` from `B002.py:291-294` and `deltaf2.py:672-673`; `HARD` severity at `B002.py:375-376` is appropriate for a percent-level CP observable. `NEEDS-HUMAN-PHYSICS` for missing core CKM `2 beta` is present at `B002.py:76-79,497-498`.

PHYSICS-OK