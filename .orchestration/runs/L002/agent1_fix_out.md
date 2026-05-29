Implemented L002 fixes.

- Interference: `lfv_three_body_from_components` now includes `I_dc ∝ Re[A_R(2G_LL+G_LR)^* + A_L(2G_RR+G_RL)^*]`. Since L001 only gives parent BR, HARD uses constructive sign-envelope upper bound; lower envelope is diagnostic. Flagged `NEEDS-HUMAN-PHYSICS`.
- Flavor pin: `mu_to_3e_from_lepton_input` rejects explicit non-`mu/e` labels; added tau->mu regression.
- Cross-check: predicted/core `2.011338841566640003e-14`, diff `0`, ratio `2.011338841566640001e-02`.
- Components: dipole `6.941358611456529606e-16`, contact `1.884227536782848111e-14`, dipole-contact envelope `+5.769771866922660645e-16`.
- Safe/excluded probes: `6.758796580084927110e-15` pass, `2.000000000000000331e-10` fail.
- Tests: `test_L002.py` 12 passed; `tests/constraints/ -q` 514 passed.