Changed `test_L004.py`: expected Au conversion now comes from `quarkConstraints.mu_e_conversion.mu_e_conversion_from_components(...)` with `gold_nuclear_inputs()`, using independently computed L001 dipole BR.

`python -m pytest tests/constraints/primary/charged_lepton/test_L004.py -q`: 9 passed  
`python -m pytest tests/constraints/ -q`: 680 passed