1. SHOULD-FIX: The “SM” anchor is a PDG experimental listing, not a theory/form-factor SM prediction: `B013.py:72` selects `branching_fraction_pdg_2025`, and `B013.py:446` admits no theory-only block. YAML `B013.yaml:77` gives PDG `(3.4 +/- 0.4)e-5`; HFLAV exp is `34.0 +/- 3.2e-6` at `B013.yaml:89`. Current budget `5.122e-6` is uncertainty-aware, but `sigma_SM=4.0e-6` is a measurement error, so the HARD budget is not defensibly SM-vs-measurement theory/form-factor room.

2. NIT: `A_Delta` is loaded with only the positive statistical error as `Anchor.uncertainty` via `B013.py:263-269`, while YAML has asymmetric `+0.37/-0.41` plus syst at `B013.yaml:152-160`. Verdict is unaffected because helicity observables are not evaluated and diagnostics preserve both sides at `B013.py:497-507`.

3. NIT: ΔF=2 M12 real/imag checks are inapplicable to B013; the BR path correctly uses dipole power `|C7_SM+C7_NP|^2 + |C7p_SM+C7p_NP|^2` at `B013.py:384-390` through the b→sγ adapter. LL running is actually used: `quarkConstraints/bsgamma.py:346-356,489-493`; for 3 TeV→4.8 GeV I get `u77=0.512`, `u78=0.119`, `C7_low=0.631*C7_MKK` when `C7_MKK=C8_MKK`.

PHYSICS-NEEDS-FIXES