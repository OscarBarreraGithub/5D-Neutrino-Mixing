R2 verdict OK: no findings. Wrappers delegate to core `*_with_running` once (e.g. [deltaf2.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/physics_adapters/deltaf2.py:207), :235, :263, :332, :401). Recomputed 3 TeV -> 2 GeV RG: `C1=0.729130912171`, `C4=3.538163974859`, `C5 -> (C4=0.894757448992,C5=0.853891627884)`.

R3 verdict OK: no findings. Complex phase path uses evolved complex `M12`, not `abs(M12)`: B002 [B002.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/beauty/B002.py:383) -> `arg(1+M12NP/M12SM)`; C002 [C002.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/charm/C002.py:493) -> `abs(Im M12NP)`. Recomputed B002 `S=0.721882988416`, `phi=0.950163068824 deg`; C002 `|ImM12|/budget=0.0234065081741`, ratio `0.906245459724`.

1. R5 BLOCKER: [quarkConstraints/deltaf2.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/deltaf2.py:722) and :910 use Q1 `M12` normalization half of the paper-mode aligned path ([hadronic.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/paper_0710_1869/eft_deltaf2/hadronic.py:762), [observables.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/paper_0710_1869/eft_deltaf2/observables.py:561)); recomputed K Q1 `0.004425624659` vs `0.008851249318` GeV^3, same factor 2 for Bd/D.

2. R5 SHOULD-FIX: [quarkConstraints/deltaf2.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/deltaf2.py:643) documents K LR `B4/B5` are FLAG 3 GeV while Wilsons are evolved to `mu_had=2 GeV` at [deltaf2.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/deltaf2.py:611); this violates scale/scheme alignment for LR-dominated K001/K002.

R5 scales verdict OK: generic roots are correct (`2.448687135269161`, `3.8317059702075125`) in [scales.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/scales.py:21); paper-mode `m_g1_GeV` is explicit and used as propagator mass.

pytest `tests/constraints/ -q`: 1072 passed.

RBATCH-NEEDS-FIXES: R5