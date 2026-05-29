Changes: [B003.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/beauty/B003.py) now validates anchor `units == "ps^-1"` and `observable == "Delta m_s"` before budget construction; auxiliary-anchor error is now auxiliary-specific.  
[test_B003.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/beauty/test_B003.py) adds probes for bad units/observable and missing auxiliary key.

Top-level now: `sm_prediction=1.170000000000000063e-11 GeV`, `experimental=1.169379362628539876e-11 GeV`.  
Diagnostics keep half-M12: `sm_m12_gev=5.850000000000000317e-12`, `experimental_m12_gev=5.846896813142699378e-12`; `predicted` remains `|M12^NP|`.

Pytest: `test_B003.py -q` -> `11 passed`; `tests/constraints/ -q` -> `55 passed`.  
`git diff --stat` tracked output: `2 files changed, 257 insertions(+), 5 deletions(-)` from pre-existing tracked changes; B003 files are currently untracked, so standard `git diff --stat` does not include them.