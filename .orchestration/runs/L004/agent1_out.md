Plan: read scaffold/YAML/core, add L004 only, run targeted and full constraint tests.

Implemented L004 as a HARD observed Au bound using the existing `mu_e_conversion` adapter/core and gold KKO inputs: `D=0.189`, `V_p=0.0974`, `V_n=0.146`, `S_p=0.0614`, `S_n=0.0918`, capture `13.07e6 s^-1`.

Budget: `7.0e-13` from `L004.yaml` `primary_current_limit` / PDG live S004 SINDRUM II Au. SM prediction is `0.0`, pure NP.

Validation: independent Au-overlap recomputation test matches sample `CR_lower=2.797698061e-18`, `CR_upper=3.096905868e-18`, ratio `3.996711516e-6`.

NEEDS-HUMAN-PHYSICS remains flagged for scalar/vector RS matching and unknown dipole-contact phase; verdict uses the conservative lower envelope.

Files changed: [L004.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/charged_lepton/L004.py:259), [test_L004.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/charged_lepton/test_L004.py:118). No L003/Al edits.

Tests: `test_L004.py` `9 passed`; `python -m pytest tests/constraints/ -q` `664 passed`.