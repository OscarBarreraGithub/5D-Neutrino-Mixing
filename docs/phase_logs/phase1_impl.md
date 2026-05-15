# Phase 1 Implementation Report

Status: DONE

Branch: `paper/quark-scan-2026q2`
Base: `a775dc3cb7224dfd3b3e00db594dda4e36d7f2dc` (`origin/main`)
Date: 2026-05-15

## Recovery And Preflight

- Verified the retry began on `main` with `HEAD == origin/main`.
- Preserved the existing `snapshots/` directory and its prior provenance files.
- Appended the preflight local-state ignore rules for `snapshots/` and `.claude/` before re-running the tar snapshot.
- Rewrote:
  - `snapshots/pre-paperwork-status.txt`
  - `snapshots/pre-paperwork-tracked.diff`
  - `snapshots/pre-paperwork-staged.diff`
- Created `snapshots/pre-paperwork-untracked.tgz` successfully after `snapshots/` was ignored.
- Created and verified `snapshots/pre-paperwork-20260515T041424.bundle`.
- Verified `git branch --list paper/quark-scan-2026q2` returned no collision.
- Stashed only `.gitignore` as `phase1-preflight-gitignore`, switched to the new paper branch, and popped the stash there.

## Commit List

1. `chore(gitignore): track paper docs + phase logs, ignore .claude/ and snapshots/`
2. `feat(qcd): add PDG 2024 MS-bar quark mass running`
3. `feat(quark): integrate PDG targets into scan gates`
4. `feat(scan): add RS-anarchy driver and dispatch scripts`
5. `feat(analysis): add quark-scan follow-up utilities`
6. `docs(paper): add quark-scan methodology and notebooks`

The final docs commit was amended to include this report.

## Verification

- `pytest -q tests/test_mass_running.py tests/test_pdg_quark_masses.py tests/test_diagnostics.py tests/test_rs_anarchy_priors.py tests/test_modern_input_registry.py tests/test_quark_benchmarks.py tests/test_quark_fit.py tests/test_quark_scan.py tests/test_quark_target_regression.py`
  - Result: `75 passed, 2 xfailed in 166.99s`.
- `pdflatex -interaction=nonstopmode -halt-on-error -file-line-error quark_scan_assumptions_compact.tex`
  - Result: success; rebuilt `docs/quark_scan_assumptions_compact.pdf`.
  - Remaining warnings are overfull boxes from long monospaced path/label strings.
- `pdfinfo docs/quark_scan_assumptions_compact.pdf`
  - Result: 7 pages, 204920 bytes.

## Notes

- `.claude/`, `data/`, and `snapshots/` are ignored as local state.
- The two `tests/test_quark_fit.py` xfails are intentionally preserved and remain Phase 2 work per the roadmap.
- No scan-output trees were committed.
