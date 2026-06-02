# B013 FIX PASS (agent1) — you implemented B013 (Bs→φγ); the physics reviewer raised a SHOULD-FIX about the budget's SM side + 1 NIT. Code review was CODE-OK. Fix honestly, re-run tests. Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. Your files: flavor_catalog_constraints/secondary/beauty/B013.py, flavor_catalog_constraints/physics_adapters/bsgamma.py, tests/constraints/secondary/beauty/test_B013.py.

## SHOULD-FIX (agent2, PHYSICS) — the "SM" side of the budget is a PDG MEASUREMENT, not a theory prediction
Currently B013.py:72 selects `branching_fraction_pdg_2025` (PDG (3.4±0.4)e-5) as the "SM" value and builds the HARD budget band as |exp − SM| + sqrt(σ_exp² + σ_SM²) with σ_SM = 4.0e-6 — but that 4.0e-6 is a MEASUREMENT error, so the band is NOT defensibly "SM-theory-vs-measurement room."

**Do this:**
1. Inspect B013.yaml for a genuine THEORY / form-factor SM prediction block for BR(Bs→φγ) (e.g. an SCET/QCD-factorization value with its theory uncertainty, distinct from the PDG measured BR). 
   - IF a theory SM prediction block exists: use IT (value + theory uncertainty) as the SM side, and build the uncertainty-aware band from |exp − SM_theory| + sqrt(σ_exp² + σ_SM,theory²). Load it via load_anchor (loud fail).
   - IF NO theory-only SM block exists in the yaml: do NOT fabricate one. Instead (a) RELABEL the budget honestly in code + diagnostics as a "measurement-consistency band" (NP must not push BR outside the measured interval), NOT as theory room; (b) add an explicit NEEDS-HUMAN-PHYSICS note that a rigorous exclusive Bs→φγ theory/form-factor SM prediction is required to set a true theory-vs-data budget; (c) keep the band conservative (it already is). The exclusive form factors are already NEEDS-HUMAN — make the SM-prediction gap explicit too.
2. Keep the constraint HARD and the C7/C7' dipole-power NP path unchanged (agent2 confirmed |C7_SM+C7_NP|²+|C7p_SM+C7p_NP|² with LL running is correct).

## NIT (do if trivial)
A_Δ is loaded with only the positive statistical error as Anchor.uncertainty (B013.py:263-269) while the yaml has asymmetric +0.37/−0.41 + syst (B013.yaml:152-160). Helicity observables are diagnostics-only (not evaluated) so the verdict is unaffected, but preserve BOTH error sides in the diagnostic if easy.

## Rules
- Numeric ConstraintResult fields real floats; complex C7 in diagnostics. Anchors via load_anchor (loud fail). Keep NEEDS-HUMAN-PHYSICS flags (RS C7/C7' proxy; exclusive form factors; A_Δ/S_φγ helicity; and now the SM-theory-prediction gap if applicable).
- Tests must cross-check INDEPENDENTLY of the adapter (recompute the LL-run C7 from quarkConstraints/bsgamma core). Update pinned budget numbers that legitimately change; a safe point passes, an excluded point fails.
- Re-run: `python -m pytest tests/constraints/secondary/beauty/test_B013.py -q` AND `python -m pytest tests/constraints/ -q`. Report counts.

OUTPUT (≤14 lines): whether a theory SM block existed, what the budget now represents + its new value, the NIT outcome, pytest counts. End with: FIX-DONE.
