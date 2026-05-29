1. QCD running: K001 now uses `evaluate_epsilon_k_with_running` via an additive adapter wrapper.
2. Budget: HARD veto uses the documented uncertainty-aware loose edge, not central residual.
3. Uncertainties: band uses BGS grouped SM uncertainties, PDG exp uncertainty, and SM-choice sensitivity; kappa is folded through BGS non-pert group.
4. Wording: K001 docstring/notes now say uncertainty-aware band and QCD-evolved Wilsons.
5. Test: `tests/constraints/primary/kaon/test_K001.py` is included via intent-to-add.

Budget: central `6.700000e-05`; combined sigma `2.367086817e-04`; band `1.1e-05` to `3.037086817e-04`; HARD budget `3.037086817e-04`. Construction cites `docs/audits/epsilon_k_sm_decision.md:34-37,71-100`.

Test couplings central-budget ratios: unrun `0.265550119962486`; run `2.84474158040246`. K001 hard-budget run ratio is `0.627567459743646`.

Pytest: `tests/constraints/ -q` -> `24 passed`; `tests/test_epsilon_k_physics.py -q` -> `20 passed`.

Isolation: `git diff --stat` shows only adapter, `K001.py`, and `tests/constraints/primary/kaon/test_K001.py` (`3 files changed, 634 insertions, 40 deletions`).