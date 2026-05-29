Plan executed: read scaffold/K001/K004/collider examples, added CR011 adapter + constraint + tests, ran the suite.

- Physics/source: records ATLAS 2025 longitudinal same-sign WW VBS fiducial upper limit from `CR011.yaml`; CMS 2020 retained as comparator.
- Parametrization: INFO/non-vetoing hard stub. No SM/EFT or RS recast; only optional human-supplied `cr011_human_fiducial_sigma_fb` gets advisory `sigma_fid / limit`.
- Built: new small adapter `flavor_catalog_constraints/physics_adapters/vbs_longitudinal.py`; reused `load_anchor`, registry, `ConstraintResult`.
- Budget: `0.45 fb` observed 95% CL, expected `0.70 fb`, ATLAS arXiv:2503.11317; CMS comparator `1.17 fb`.
- Validation: independent test recomputes `0.225 / 0.45 = 0.5`; no SM prediction by design.
- NEEDS-HUMAN-PHYSICS: diagnostics separately flag SM/EFT likelihood/unitarity/form factors and RS strong-EWSB matching/widths/acceptance.
- Files changed: `CR011.py`, `vbs_longitudinal.py`, `test_CR011.py`.
- Tests: CR011 focused `8 passed`; full `python -m pytest tests/constraints/ -q` -> `786 passed`.

Unrelated untracked orchestration/other-agent files are still present and untouched.