1. NIT: Correct-amplitude check OK/N.A. K013 is a Delta S=1 branching-ratio stub, not a Delta F=2 Delta m or CP-mixing observable; it computes no M12 real/imag projection. `flavor_catalog_constraints/primary/kaon/K013.py:297`, `:305-313`.

2. NIT: QCD-running check OK/N.A. No Wilsons are evaluated and no non-running DeltaF2 verdict path is used; K013 only calls the radiative stub adapter. Running effect is undefined, not omitted. `K013.py:299-303`; `physics_adapters/radiative_kaon.py:72-112`.

3. NIT: Budget/severity OK. Non-vetoing bookkeeping room is the full measured BR, not an SM residual: BR=1.273e-6, sigma=0.033e-6, budget=1.273e-6, ratio=1.0, INFO, passes=True. `K013.py:290-319`, `:321-350`.

4. NIT: Anchor values OK. YAML matches snapshots: PDG average 1.273 +/- 0.033 in units 1e-6; KTeV 1.28 +/- 0.06 +/- 0.01; NA48 1.27 +/- 0.04 +/- 0.01. `flavor_catalog/processes/kaon/K013.yaml:69-97`; `references/K013/pdg2025_kl_pi0gammagamma.txt:12-14`.

5. NIT: SM/NP limitation wording OK. Code explicitly flags SM ChPT O(p^4)/O(p^6) loops/vector exchange and NP Delta S=1 plus chiral-counterterm matching as NEEDS-HUMAN-PHYSICS, with no fake ChPT/RS amplitude. `K013.py:5-15`, `:333-339`; `radiative_kaon.py:32-43`.

PHYSICS-OK