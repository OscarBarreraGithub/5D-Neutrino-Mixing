# paper_0710_1869 Current Status

This is the authoritative handoff note for the next Codex instance working on
the paper-facing `quarkConstraints.paper_0710_1869` path.

Updated on 2026-04-12.

## Current Repo State

- branch state: the authoritative handoff is the exact published commit that
  carries this note for this paper-facing module, not any later-moving branch
  tip
- the verified `LR-DEFAULT-HAD-1` closure patchset was committed and pushed on
  2026-04-12 as `df67be1dd5d08a27563c4bb31d62e94bd2e82e1b`
  (`paper_0710_1869: close LR-DEFAULT-HAD-1 patchset`)
- this note describes the durable published handoff state, not the earlier
  pre-release local/dirty worktree
- the published closure patchset includes these owned package, acceptance, and
  paper-artifact/golden updates:
  `quarkConstraints/paper_0710_1869/__init__.py`,
  `quarkConstraints/paper_0710_1869/eft_deltaf2/__init__.py`,
  `quarkConstraints/paper_0710_1869/eft_deltaf2/hadronic.py`,
  `results/paper_0710_1869/default_kaon/hadronic.json`,
  `results/paper_0710_1869/default_kaon/provenance.json`,
  `scripts/benchmark_quark_0710_1869.py`,
  `tests/golden/paper_0710_1869/default_kaon_np_only/hadronic.json`,
  `tests/golden/paper_0710_1869/default_kaon_np_only/provenance.json`,
  `tests/test_paper_benchmarks.py`, and
  `tests/test_paper_hadronic_lr_inputs.py`
- expendable untracked SLURM stdout/stderr residue from the pre-release runs
  was removed and is not part of the handoff state

## Current Completed Milestone

- completed milestone: `LR-DEFAULT-HAD-1`
- completed milestone chain:
  `LR-MAP-1 -> LR-RG-1 -> LR-HAD-1 -> LR-OBS-1 -> LR-TOTAL-1 -> BS-Q1-CUSTOM-1 -> D0-Q1-CUSTOM-1 -> LR-RCHI-FREEZE-1 -> LR-DEFAULT-HAD-1`

`LR-DEFAULT-HAD-1` closed the default-kaon LR hadronic slice by freezing the
already-reviewed ETM Collaboration, JHEP 03 (2013) 089, arXiv:1207.1287,
Table 1, MS(Buras) `2 GeV` default source package with `B4(2 GeV) = 0.78(3)`
and `B5(2 GeV) = 0.57(4)` around the already-frozen kaon LR
`R_chi(mu_had = 2.0 GeV)` contract. The frozen ids are
`bundle_id = hadronic.kaon.lr.default.etm2013_ms_2gev.v1`,
`source_id = hadronic.kaon.lr.default.etm2013_ms_2gev.aggregate.v1`, and
`input_policy_id = default_source.etm2013.table1.ms_2gev.no_hidden_conversion.v1`.

## Preserved Boundary

- default/exported kaon observable path, benchmark numerics, paper artifacts,
  and standalone verifier remain Q1-only
- the default kaon LR hadronic bundle is now frozen in the paper package under
  the ETM 2013 Table 1 MS(Buras) `2 GeV` source gate, but it is not
  auto-consumed by the kaon LR-only or combined `Q1 + LR` custom surfaces
- kaon LR-only and custom combined `Q1 + LR` surfaces remain custom-input-only
- custom `B_d`, `B_s`, and `D0` Q1 surfaces remain separate custom-input-only
  entrypoints
- the already-frozen BMU-linked MS/NDR projector/operator-normalization
  contract stays unchanged, and no hidden `3 GeV -> 2 GeV` running or
  `RI-MOM -> MS` conversion is allowed in the default LR hadronic path
- artifact/verifier widening, conservative-bound `D0` interpretation, and
  `epsilon_K` remain out of scope

## Verification State

Second-pass review state on 2026-04-12:

- logic second pass: no blocking or medium logic, scope, or contract finding
  remained
- numerical second pass: no blocking or medium numerical, determinism, or
  verification finding remained
- physics second pass: no blocking or medium finding remained; the custom LR
  policy id and note text were judged honest for the closed
  `LR-DEFAULT-HAD-1` boundary
- post-ruff-fix logic re-review: clean and semantics-preserving; the adjacent
  string-literal split preserved the searched text and surrounding boolean
  logic at the benchmark/test sites

Fresh clean verification jobs on 2026-04-12 (verification root:
`/n/holylabs/randall_lab/Lab/obarrera/slurm-verification/lr-default-had-post-ruff-fix-20260412T075512`):

- `5228468` `ruff`: `COMPLETED`, exit `0:0`, `All checks passed!`
- `5228469` `pytest`: `COMPLETED`, exit `0:0`,
  `231 passed, 1 skipped in 699.20s`
- `5228470` `benchmark`: `COMPLETED`, exit `0:0`,
  `"tracked_default_exports_match_current_export": true`,
  `"writer_outputs_are_deterministic": true`
- `5228471` `export`: `COMPLETED`, exit `0:0`, `manifest_stable=yes`; the same
  four tracked artifact SHA256 values were present before and after export

## Next Honest Milestone

- next milestone: `none`
- there is no remaining scientific/code milestone relative to the current
  frozen paper-facing claim boundary
- the paper-facing module is publishable on its current claim boundary, with no
  remaining milestone or release-hygiene follow-up inside that boundary

Any future widening beyond this boundary, including artifact/verifier
widening, conservative-bound `D0` interpretation, `epsilon_K`, or automatic
consumption of the default LR hadronic bundle by custom surfaces, must be
opened as a new later milestone rather than inferred from
`LR-DEFAULT-HAD-1`.
