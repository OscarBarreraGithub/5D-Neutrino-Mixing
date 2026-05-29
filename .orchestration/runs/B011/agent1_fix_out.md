Changes: fixed chirality so left overlap feeds `C7/C8`, right overlap feeds `C7'/C8'`; added LL `C7-C8` running to `mu_b=4.8 GeV`; kept `NEEDS-HUMAN-PHYSICS`; added shared-module recomputation test.

SM BR after running: `3.400000e-4`.

Running coefficients at `M_KK=3 TeV -> mu_b=4.8 GeV`: `C7=0.511998*C7(MKK)+0.119166*C8(MKK)`, `C8=0.556685*C8(MKK)`.

Pytest: `tests/constraints/primary/beauty/test_B011.py -q` -> `11 passed`; `tests/constraints/ -q` -> `158 passed, 1 failed` due unrelated `B022.load_pdg_block` monkeypatch attribute failure.