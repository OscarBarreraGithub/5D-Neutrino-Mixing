Option (a).

Added YAML anchors for `f_B=0.1900±0.0013 GeV` and `|V_ub|=0.00368`, with FLAG/PDG provenance: https://arxiv.org/abs/2411.04268 and https://pdg.lbl.gov/2024/reviews/rpp2024-rev-vcb-vub.pdf.

B009 now loads both via `load_anchor` into `quarkConstraints/leptonic_tree.py`; RS proxy remains `NEEDS-HUMAN-PHYSICS`. SM BR stays `8.6308e-5` vs UTfit `8.65e-5`.

Pytest: B009 focused `10 passed`; `tests/constraints/ -q` `514 passed`.