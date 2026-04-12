# paper_0710_1869 Current Status

This is the authoritative handoff note for the next Codex instance working on
the paper-facing `quarkConstraints.paper_0710_1869` path.

## Current Repo State

- branch state: local `main` tracks `origin/main`
- sync state: clean worktree and `main == origin/main`
- checkpoint commit: `eb9089d48858d79b525cb86a2c1f624b63e74e0a`

## Current Completed Milestone

- completed through `D0-Q1-CUSTOM-1`
- completed milestone: `LR-RCHI-FREEZE-1`

`LR-RCHI-FREEZE-1` froze only the standalone kaon LR
`R_chi(mu_had = 2.0 GeV)` semantics and provenance under the exact BV 2004
definition
`R_chi(mu) = [m_K / (m_s(mu) + m_d(mu))]^2`
with the PDG 2024 `N_L = 4` mass-source chain and explicit no-hidden-conversion
policy.

## Preserved Boundary

- default/exported kaon bundles remain Q1-only
- kaon LR-only and custom combined `Q1 + LR` surfaces remain custom-input-only
- default LR `B4/B5` hadronic inputs are not yet frozen
- artifacts and the standalone verifier remain unchanged
- no paper-facing claim was widened beyond the completed `LR-RCHI-FREEZE-1`
  slice

## Verification Already Run

The completed `LR-RCHI-FREEZE-1` checkpoint was verified with:

- `pytest -q tests/test_paper_*.py`
- `python scripts/benchmark_quark_0710_1869.py --require-package`
- `python scripts/export_quark_0710_1869_artifacts.py`

Those checks passed before this handoff checkpoint was synced to `main`.

## Next Honest Milestone

The next honest milestone is `LR-DEFAULT-HAD-1`.

Its narrow scope is:

- freeze sourced default kaon LR `B4/B5` inputs
- build the default LR hadronic bundle around the already-frozen
  `R_chi(mu_had = 2.0 GeV)` contract
- keep the default/exported kaon observable path, artifacts, and verifier
  unchanged unless a later explicitly reviewed slice widens them

Before coding starts again, re-run the required workflow:

- planner convergence with two `gpt-5.4` / `xhigh` planners
- disjoint ownership split across implementation, tests/acceptance, and docs
- separate non-author logic, numerical, and physics reviews
- second-pass non-author re-review for any blocking-finding fix
