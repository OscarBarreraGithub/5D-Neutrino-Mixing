1. BLOCKER — T014 compares `B(Z -> qq') < 2.9e-3` to a hadronic-width fraction, not a total-Z branching fraction: `flavor_catalog_constraints/physics_adapters/zpole.py:273` sums only `u,d,s,c,b`, then `:289-290` uses `Gamma_FCNC/(Gamma_had^SM+Gamma_FCNC)`. Correct physics is `Gamma_FCNC/(Gamma_Z^SM+Gamma_FCNC)` for the catalog `B` limit; see `flavor_catalog/references/T014/ecfa2025_fv_z_decays_extract.txt:29-37` and YAML rows `flavor_catalog/processes/secondary/top_higgs_ew/T014.yaml:76-120`. This makes the HARD veto too tight by about `1/0.699 = 1.43` in BR; code threshold corresponds to true `B ~= 2.0e-3`, not `2.9e-3`.

2. SHOULD-FIX — The effective-coupling diagnostic inherits that normalization problem: `flavor_catalog_constraints/secondary/top_higgs_ew/T014.py:308` / `flavor_catalog_constraints/physics_adapters/zpole.py:330-352` gives `sqrt(|dgL|^2+|dgR|^2)=0.03497` for `2.9e-3`, while a total-Z denominator in the same charge-summed zpole convention gives about `0.0421`; the source convention translation is `4.96e-2` at `flavor_catalog/references/T014/abu_ajamieh_2026_fv_z_quark_extract.txt:30-32`.

3. NIT — Correct amplitude structure: for this decay-width constraint the code uses `|delta_g_L|^2 + |delta_g_R|^2` at `flavor_catalog_constraints/physics_adapters/zpole.py:267-270`; no Δm/CP real-vs-imaginary issue applies.

4. NIT — QCD running checklist is not applicable here: T014 is an on-shell EW Z-width bound, not a ΔF=2 Wilson/matrix-element observable; the evaluated path calls zpole proxy helpers at `flavor_catalog_constraints/secondary/top_higgs_ew/T014.py:426-432`, not a non-running ΔF=2 path.

5. NIT — Anchors, SM treatment, and severity are otherwise physics-consistent: all three rows load `2.9e-3` at `95% CL`; SM `bs=4.2e-8`, `bd=1.8e-9` in `T014.yaml:153-160` are negligible vs the limit; HARD and NEEDS-HUMAN-PHYSICS labeling is appropriate at `T014.py:98-103` and `:486-500`.

PHYSICS-NEEDS-FIXES