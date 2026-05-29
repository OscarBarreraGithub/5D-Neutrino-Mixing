Plan: studied scaffold/K001/K004, added STU core+adapter, implemented EW001, tested.

Physics: PDG 2025 U-fixed `S,T,rho(S,T)` from `EW001.yaml`; RS proxy `ΔS=cS v²/M_KK²`, `ΔT=πL/(2cW²) v²/M_KK²`, `U=0`.

Built: `quarkConstraints/oblique_stu.py`; adapter `physics_adapters/oblique_stu.py`; constraint `primary/top_higgs_ew/EW001.py`; test `test_EW001.py`.

Budget: `χ²_2,95 = 5.991464547`, from PDG S/T covariance plus standard 95% 2D fit contour.

Validation: SM reference `S=T=U=0` gives `χ²=0.962793485`; `M_KK=6 TeV` gives `χ²=3.631258796`, ratio `0.606`, passes; `3 TeV` fails.

NEEDS-HUMAN-PHYSICS: full EW KK/Z/Z′/custodial/brane-term matching is absent on `ParameterPoint`; proxy is explicitly flagged in docs/diagnostics.

Tests: `python -m pytest tests/constraints/ -q` -> `851 passed`; focused EW001 -> `10 passed`. `ruff` unavailable; `py_compile` passed.