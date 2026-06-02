1. NIT: DeltaF2 amplitude check is N/A. CR010 is a collider mass lower bound, not a Δm/CP observable; it imports only `collider_resonance` helpers at `flavor_catalog_constraints/primary/collider_rs/CR010.py:54` and compares masses in `quarkConstraints/collider_resonance.py:166`.

2. NIT: QCD running check is N/A. No Wilson coefficients or hadronic-scale matrix elements enter; verdict uses `evaluate_collider_resonance_limit` at `CR010.py:448`. Running effect is undefined for this direct LHC mass exclusion.

3. NIT: Budget is physics-defensible for the requested proxy. Active bound is `max(m_T,m_B)=1.37 TeV`, from `CR010.yaml:95` and `CR010.yaml:114`; code uses the max at `CR010.py:172`. Correct for the documented `m_T=m_B=M_KK` proxy.

4. NIT: Anchor numbers match snapshots/audit. PDG snapshot has `t'` doublet `>1370 GeV` and `b'` doublet `>1370 GeV`; ATLAS snapshot states `(T,B)` masses below `1.37 TeV` excluded. Diagnostic endpoints `1.70/1.57/1.56/1.54/1.48 TeV` also match.

5. NIT: Severity/units/notes are appropriate. `HARD` is set at `CR010.py:393`; input `M_KK` GeV is converted to TeV by the collider core, and `sigma*BR`/branching-surface limitations are explicitly marked `NEEDS-HUMAN-PHYSICS` at `CR010.py:500`.

PHYSICS-OK