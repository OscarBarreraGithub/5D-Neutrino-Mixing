Changes:
- `B002.py`: beta and penguin numeric anchors now load via `load_anchor(..., value_key=...)`; 2β `NEEDS-HUMAN-PHYSICS` flag kept.
- `test_B002.py`: independent oracle now calls core `_evolve_wilsons + compute_m12_np` and `evaluate_bd_mixing_with_running`, not the B002 phase adapter.
- Adapter contract: B002 remains on append-only complex phase helper; magnitude wrapper unchanged.

Cross-check:
- safe `(1e-3+1e-3j)`: `M12=0+2.9852989857064115e-15j GeV`, `phi=0.950163068823887 deg`, `S=0.7218829884163276`, `ratio=0.6716246168494265`, pass, agreement `0`.
- excluded `(5e-3+5e-3j)`: `M12=0+7.463247464266028e-14j GeV`, `phi=22.52011154516739 deg`, `S=0.925739373309186`, `ratio=12.19355509419042`, fail, agreement `0`.

Pytest:
- `python -m pytest tests/constraints/primary/beauty/test_B002.py -q`: `8 passed`
- `python -m pytest tests/constraints/ -q`: `89 passed`