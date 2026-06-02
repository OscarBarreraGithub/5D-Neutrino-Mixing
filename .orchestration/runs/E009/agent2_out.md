1. NIT: ΔF=2 amplitude check is N/A for E009; no `M_12`, Δm, or CP-mixing observable is computed. Correct for this stub: `predicted=None`, `ratio=None`, non-vetoing `passes=True`. `flavor_catalog_constraints/primary/edm_neutrino/E009.py:812`, `:824-832`.

2. NIT: QCD-running check is N/A; E009 does not use any ΔF=2 running/non-running path, and no RS Weinberg Wilson is available to evolve. Missing RG/operator convention is correctly flagged as human physics. `physics_adapters/weinberg_operator.py:8-12`, `E009.py:837-856`.

3. NIT: Bookkeeping budget is acceptable because it is INFO-only/non-vetoing, not a hard residual. Default is Haisch-Hala `|C_6| < 1.2e-11 GeV^-2`; uncertainty is tracked as 50% and assumptions say central one-source benchmark. `E009.py:806`, `:827-840`, `:896-902`; `E009.yaml:217-232`.

4. NIT: Anchor numbers match YAML/snapshots: PDG `|d_n|<1.8e-26 e cm`, table `<0.18 x 10^-25 e cm`; Abel `0.0 +/- 1.1_stat +/- 0.2_sys`. `E009.yaml:107-146`; snapshots `pdg2026...txt:18-23`, `abel2020...txt:15-18`.

5. NIT: Weinberg translations check numerically: `9.12e-13/0.022 = 4.15e-11 -> 4.1e-11 GeV^-2`; `9.12e-13/0.074 = 1.23e-11 -> 1.2e-11 GeV^-2`. `E009.yaml:148-236`, `:255-259`; audit `factcheck_edm_neutrino.md:113-117`.

6. NIT: Dual `NEEDS-HUMAN-PHYSICS` coverage is correct: one for hadronic CP-odd gluonic matrix elements, one for RS CP-odd gluonic matching absent from `ParameterPoint`. `physics_adapters/weinberg_operator.py:41-52`; `E009.py:851-859`.

7. NIT: Units and severity are physically consistent: EDM anchor in `e cm`, conversion to `GeV^-1`, Weinberg bounds in `GeV^-2`, `Severity.INFO`. Test file passes: `9 passed`. `E009.py:658-682`, `:801-832`.

PHYSICS-OK