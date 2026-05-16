# Phase 3 Hole #10 Implementation Report

## Figure Inventory Summary

- Methodology-note figure references: 9 PDF files.
- Referenced and present at cited paths: 9.
- Referenced but missing: 0.
- Referenced but untracked before this pass: 1 (`results/figures/quark/yukawa_size_envelope_vs_anarchic.pdf`); now tracked and whitelisted.
- Unreferenced figures moved aside: 49 total.
  - 39 top-level unreferenced PNG/PDF files moved to `results/figures/quark/exploratory/`.
  - 10 unreferenced `runA/` scratch PNG/PDF files moved to `results/figures/quark/exploratory/runA/`.
- Figures regenerated: 0.
- Historical snapshots left in place: `results/figures/quark_pre_audit_constants/` and `results/figures/quark_baseline_800k/`.

Detailed disposition inventory: `docs/audits/figure_prune_inventory.md`.

## PDF Rebuild

- Command sequence: removed LaTeX auxiliaries, ran `pdflatex -interaction=nonstopmode quark_scan_methodology_note.tex` twice from `docs/`, then removed auxiliaries again.
- Second pass unresolved references/citations: none detected.
- Final methodology-note page count: 19.

## Test Suite

Command:

```bash
/n/home09/obarrera/.conda/envs/ising_bootstrap/bin/python -m pytest -q
```

Result: 543 passed, 1 skipped, 0 failed, 0 xfailed in 889.51 seconds.

## Lint

Command:

```bash
/n/home09/obarrera/.conda/envs/ising_bootstrap/bin/ruff check . 2>&1 || echo "ruff not installed or violations found"
```

Result: `/n/home09/obarrera/.conda/envs/ising_bootstrap/bin/ruff` is not installed; no ruff violations were evaluated.

## Tracking And Document Hygiene

- `git status --short` was empty after the figure-move, inventory, and PDF commits.
- `git ls-files docs/quark_scan_methodology_note.{tex,pdf}` confirms both methodology-note files are tracked.
- `git ls-files | grep -E "phase_logs/|audits/|artifact_manifest" | wc -l` returned 40 after adding the figure inventory.
- README and `CLAUDE.md` still reference `paper/quark-scan-2026q2` and `docs/quark_scan_methodology_note.tex`/`.pdf`.
- README figure-tracking wording now matches the pruned methodology-note figure set.
- No methodology-note figure references were added or removed.
- No physics code, constants, scan outputs, or methodology-note source text were changed.

Hole #10 ready for peer review.
