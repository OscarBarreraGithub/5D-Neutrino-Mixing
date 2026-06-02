No theory-only SM prediction block existed in `B013.yaml`.

Budget now represents a conservative measurement-consistency band using HFLAV 2024 vs PDG 2025 measurement rows, not SM-theory-vs-data room. Value: `5.122499389946279e-06`.

Added explicit NEEDS-HUMAN-PHYSICS diagnostic for the missing exclusive SM theory/form-factor prediction. C7/C7' dipole-power path unchanged.

NIT fixed: `A_Delta` diagnostics now preserve stat `+0.37/-0.41` and combined plus/minus errors.

Pytest: `tests/constraints/secondary/beauty/test_B013.py -q` -> 12 passed. `tests/constraints/ -q` -> 1033 passed.

FIX-DONE.