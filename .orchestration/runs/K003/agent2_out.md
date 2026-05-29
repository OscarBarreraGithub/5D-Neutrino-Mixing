1. NIT: Check 1 OK. No amplitude misuse: K003 computes no ΔS=1 NP amplitude (`predicted=None`) and no ΔF=2 `M12` proxy; a real CP calculation would need imaginary ΔS=1 penguin/chromomagnetic amplitudes. `K003.py:392`, `K003.py:402`, `K003.py:420`, `kaon_direct_cp.py:82`.

2. NIT: Check 2 OK/N/A for hard stub. No non-running ΔF=2 path is used for the verdict; no ΔS=1 RG exists, and that absence is explicitly flagged. Running effect not quantifiable without the missing ΔS=1 matching/RG. `K003.py:410`, `K003.py:435`, `kaon_direct_cp.py:41`.

3. NIT: Check 3 OK. Budget is non-vetoing bookkeeping, not an SM residual: `|Re(epsilon'/epsilon)| = 0.00166`, ratio `1.0`, `passes=True`, severity INFO. Experimental uncertainty `0.00023` is recorded; acceptable only because this is not a veto. `K003.py:121`, `K003.py:395`, `K003.py:401`, `K003.py:444`.

4. NIT: Check 4 OK. Anchors match YAML/snapshots: PDG `0.00166 +/- 0.00023`; KTeV `19.2 +/- 2.1 x10^-4`; NA48 `14.7 +/- 2.2 x10^-4`; RBC/UKQCD `21.7(2.6)(6.2)(5.0) x10^-4`. `K003.yaml:66`, `K003.yaml:97`, `K003.yaml:107`, `K003.yaml:117`.

5. NIT: Check 5 OK. Units are dimensionless throughout; `sm_prediction=0.00217` is clearly labeled RBC/UKQCD context, not veto-grade SM subtraction, with combined uncertainty `0.00083785`; both SM and NP NEEDS-HUMAN-PHYSICS flags present. `K003.py:403`, `K003.py:429`, `K003.py:426`.

6. NIT: Verification OK. `pytest -q tests/constraints/primary/kaon/test_K003.py` passes: `9 passed`.

PHYSICS-OK