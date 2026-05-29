1. NIT: ΔF=2 amplitude/CP-part check is N/A for K017; correct rate physics is used: `|1+δ_e|^2/|1+δ_mu|^2` at `quarkConstraints/leptonic_tree.py:484-490`, and K017 compares `|np_shift_ratio|/budget` at `flavor_catalog_constraints/primary/kaon/K017.py:434-437`.

2. NIT: QCD running is N/A; K017 calls the leptonic-tree adapter, not Delta-F=2 Wilson/M12 code (`K017.py:38-42,407-410`; `physics_adapters/leptonic_tree.py:156-169`). No non-running Delta-F=2 path is used; running effect = N/A/0 Wilson coefficients.

3. NIT: Budget is defensible and uncertainty-aware: YAML gives `R_exp=(2.488±0.009)e-5`, `R_SM=(2.477±0.001)e-5`; code uses `|Δ|+sqrt(σ_exp^2+σ_SM^2)=1.10e-7+9.055e-8=2.006e-7` (`K017.py:285-304`; `K017.yaml:96-154`).

4. NIT: Anchor numbers match YAML/snapshots: PDG display `(2.488±0.009)e-5`, NA62 `2.488±0.007±0.007`, KLOE `2.493±0.025±0.019`, SM `2.477±0.001` (`pdg2025_rk_api_20260516.txt:14-17,30-45`; `cirigliano_rosell_2007_arxiv0707_4464.txt:10-12`).

5. NIT: Severity/units/notes OK: observable and budget are dimensionless, masses are in GeV, tree ratio `2.56896e-5` with radiative multiplier `0.964202` reproduces `2.477e-5`; HARD proxy is explicitly flagged `NEEDS-HUMAN-PHYSICS` (`K017.py:360,497-503`; `quarkConstraints/leptonic_tree.py:136-140,241-260`).

PHYSICS-OK