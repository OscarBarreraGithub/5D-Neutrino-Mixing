# Figure Prune Inventory

Phase 3 hole #10 figure-pruning audit for `docs/quark_scan_methodology_note.tex`.

## Source Reference Scan

Command:

```bash
grep -nE "includegraphics" docs/quark_scan_methodology_note.tex
```

Output:

```text
355:  \includegraphics[width=0.92\linewidth]{yukawa_size_envelope_vs_anarchic.pdf}
567:  \includegraphics[width=0.95\linewidth]{rs_anarchy_mkk_min_hist_gsstar.pdf}
591:  \includegraphics[width=0.78\linewidth]{rs_anarchy_max_ratio_vs_mkk.pdf}
601:  \includegraphics[width=0.78\linewidth]{rs_anarchy_per_system_pass_rate.pdf}
628:  \includegraphics[width=0.95\linewidth]{rs_anarchy_gate_sensitivity.pdf}
736:  \includegraphics[width=0.78\linewidth]{rs_anarchy_mkk_min_hist_by_cvals.pdf}
793:  \includegraphics[width=0.95\linewidth]{rs_anarchy_mkk_min_hist_by_yprior_gsstar.pdf}
826:  \includegraphics[width=0.78\linewidth]{rs_anarchy_mkk_min_hist_by_pdg_tightness.pdf}
875:  \includegraphics[width=0.78\linewidth]{rs_anarchy_cfw_comparison.pdf}
```

## Summary

- Referenced methodology-note figures: 9 PDF files.
- Referenced and present at `results/figures/quark/<filename>.pdf`: 9.
- Referenced but missing: 0.
- Referenced but untracked before this pass: 1 (`yukawa_size_envelope_vs_anarchic.pdf`); it is now whitelisted in `.gitignore` and tracked.
- Unreferenced top-level PNG/PDF files moved to `results/figures/quark/exploratory/`: 39.
- Unreferenced `results/figures/quark/runA/` PNG/PDF files moved to `results/figures/quark/exploratory/runA/`: 10.
- Figures regenerated during this pass: 0.
- Historical snapshot directories left untouched: `results/figures/quark_pre_audit_constants/` and `results/figures/quark_baseline_800k/`.

After pruning, the top-level `results/figures/quark/` submission figure set contains only the nine PDF files referenced by the methodology note.

## Referenced Figures

| Line | Figure | Status | Regeneration source |
|---:|---|---|---|
| 355 | `results/figures/quark/yukawa_size_envelope_vs_anarchic.pdf` | Present; newly tracked | `scripts/yukawa_envelope_vs_anarchic.py` |
| 567 | `results/figures/quark/rs_anarchy_mkk_min_hist_gsstar.pdf` | Present; tracked | `scripts/rs_anarchy_mkk_min_hist_gsstar.py --draws scan_outputs/rs_anarchy_runA_20260515T085316/draws.jsonl` |
| 591 | `results/figures/quark/rs_anarchy_max_ratio_vs_mkk.pdf` | Present; tracked | `scripts/plot_rs_anarchy_summary.py --summary scan_outputs/rs_anarchy_runA_20260515T085316/tile_summary.json --draws scan_outputs/rs_anarchy_runA_20260515T085316/draws.jsonl` |
| 601 | `results/figures/quark/rs_anarchy_per_system_pass_rate.pdf` | Present; tracked | `scripts/plot_rs_anarchy_summary.py --summary scan_outputs/rs_anarchy_runA_20260515T085316/tile_summary.json --draws scan_outputs/rs_anarchy_runA_20260515T085316/draws.jsonl` |
| 628 | `results/figures/quark/rs_anarchy_gate_sensitivity.pdf` | Present; tracked | `scripts/rs_anarchy_gate_sensitivity.py --draws scan_outputs/rs_anarchy_runA_20260515T085316/draws.jsonl --summary scan_outputs/rs_anarchy_runA_20260515T085316/tile_summary.json` |
| 736 | `results/figures/quark/rs_anarchy_mkk_min_hist_by_cvals.pdf` | Present; tracked | `scripts/rs_anarchy_mkk_min_hist_by_cvals.py` with the four Run3 directories listed in `docs/artifact_manifest.md` |
| 793 | `results/figures/quark/rs_anarchy_mkk_min_hist_by_yprior_gsstar.pdf` | Present; tracked | `scripts/rs_anarchy_mkk_min_hist_by_yprior_gsstar.py` with RUNA baseline and the three RunB directories listed in `docs/artifact_manifest.md` |
| 826 | `results/figures/quark/rs_anarchy_mkk_min_hist_by_pdg_tightness.pdf` | Present; tracked | `scripts/rs_anarchy_mkk_min_hist_by_pdg_tightness.py --draws scan_outputs/rs_anarchy_runA_20260515T085316/draws.jsonl` |
| 875 | `results/figures/quark/rs_anarchy_cfw_comparison.pdf` | Present; tracked | `scripts/rs_anarchy_cfw_comparison.py --run scan_outputs/rs_anarchy_runA_20260515T085316` |

## Unreferenced Figure Dispositions

| Original path | Disposition | Rationale |
|---|---|---|
| `results/figures/quark/c_values_2sigma_vs_1sigma.pdf` | Moved to `results/figures/quark/exploratory/c_values_2sigma_vs_1sigma.pdf` | Legacy envelope diagnostic; not cited by the methodology note. |
| `results/figures/quark/c_values_2sigma_vs_1sigma.png` | Moved to `results/figures/quark/exploratory/c_values_2sigma_vs_1sigma.png` | Legacy envelope diagnostic; not cited by the methodology note. |
| `results/figures/quark/fig1_exclusion_boundaries.png` | Moved to `results/figures/quark/exploratory/fig1_exclusion_boundaries.png` | Older envelope/publication snapshot; superseded for the methodology-note figure set. |
| `results/figures/quark/fig1_exclusion_boundaries_1sigma.pdf` | Moved to `results/figures/quark/exploratory/fig1_exclusion_boundaries_1sigma.pdf` | Older envelope/publication variant; not cited by the methodology note. |
| `results/figures/quark/fig1_exclusion_boundaries_1sigma.png` | Moved to `results/figures/quark/exploratory/fig1_exclusion_boundaries_1sigma.png` | Older envelope/publication variant; not cited by the methodology note. |
| `results/figures/quark/fig1_exclusion_boundaries_2sigma.pdf` | Moved to `results/figures/quark/exploratory/fig1_exclusion_boundaries_2sigma.pdf` | Older envelope/publication variant; not cited by the methodology note. |
| `results/figures/quark/fig1_exclusion_boundaries_2sigma.png` | Moved to `results/figures/quark/exploratory/fig1_exclusion_boundaries_2sigma.png` | Older envelope/publication variant; not cited by the methodology note. |
| `results/figures/quark/fig1_exclusion_boundaries_pdg2024.pdf` | Moved to `results/figures/quark/exploratory/fig1_exclusion_boundaries_pdg2024.pdf` | Older envelope/publication variant; not cited by the methodology note. |
| `results/figures/quark/fig1_exclusion_boundaries_pdg2024.png` | Moved to `results/figures/quark/exploratory/fig1_exclusion_boundaries_pdg2024.png` | Older envelope/publication variant; not cited by the methodology note. |
| `results/figures/quark/fig2_mkk_bound_1sigma.pdf` | Moved to `results/figures/quark/exploratory/fig2_mkk_bound_1sigma.pdf` | Older bound-plot variant; not cited by the methodology note. |
| `results/figures/quark/fig2_mkk_bound_1sigma.png` | Moved to `results/figures/quark/exploratory/fig2_mkk_bound_1sigma.png` | Older bound-plot variant; not cited by the methodology note. |
| `results/figures/quark/fig2_mkk_bound_2007_vs_modern.png` | Moved to `results/figures/quark/exploratory/fig2_mkk_bound_2007_vs_modern.png` | Older comparison snapshot; not cited by the methodology note. |
| `results/figures/quark/fig2_mkk_bound_2007_vs_modern_pdg2024.pdf` | Moved to `results/figures/quark/exploratory/fig2_mkk_bound_2007_vs_modern_pdg2024.pdf` | Older comparison variant; not cited by the methodology note. |
| `results/figures/quark/fig2_mkk_bound_2007_vs_modern_pdg2024.png` | Moved to `results/figures/quark/exploratory/fig2_mkk_bound_2007_vs_modern_pdg2024.png` | Older comparison variant; not cited by the methodology note. |
| `results/figures/quark/fig2_mkk_bound_2sigma.pdf` | Moved to `results/figures/quark/exploratory/fig2_mkk_bound_2sigma.pdf` | Older bound-plot variant; not cited by the methodology note. |
| `results/figures/quark/fig2_mkk_bound_2sigma.png` | Moved to `results/figures/quark/exploratory/fig2_mkk_bound_2sigma.png` | Older bound-plot variant; not cited by the methodology note. |
| `results/figures/quark/rs_anarchy_cfw_comparison.png` | Moved to `results/figures/quark/exploratory/rs_anarchy_cfw_comparison.png` | PNG sibling of cited PDF; the methodology note cites the PDF only. |
| `results/figures/quark/rs_anarchy_gate_sensitivity.png` | Moved to `results/figures/quark/exploratory/rs_anarchy_gate_sensitivity.png` | PNG sibling of cited PDF; the methodology note cites the PDF only. |
| `results/figures/quark/rs_anarchy_max_ratio_vs_mkk.png` | Moved to `results/figures/quark/exploratory/rs_anarchy_max_ratio_vs_mkk.png` | PNG sibling of cited PDF; the methodology note cites the PDF only. |
| `results/figures/quark/rs_anarchy_mkk_min_hist.pdf` | Moved to `results/figures/quark/exploratory/rs_anarchy_mkk_min_hist.pdf` | Earlier perturbative/secondary-axis histogram; replaced in the note by `rs_anarchy_mkk_min_hist_gsstar.pdf`. |
| `results/figures/quark/rs_anarchy_mkk_min_hist.png` | Moved to `results/figures/quark/exploratory/rs_anarchy_mkk_min_hist.png` | Earlier perturbative/secondary-axis histogram; replaced in the note by `rs_anarchy_mkk_min_hist_gsstar.pdf`. |
| `results/figures/quark/rs_anarchy_mkk_min_hist_by_cvals.png` | Moved to `results/figures/quark/exploratory/rs_anarchy_mkk_min_hist_by_cvals.png` | PNG sibling of cited PDF; the methodology note cites the PDF only. |
| `results/figures/quark/rs_anarchy_mkk_min_hist_by_pdg_tightness.png` | Moved to `results/figures/quark/exploratory/rs_anarchy_mkk_min_hist_by_pdg_tightness.png` | PNG sibling of cited PDF; the methodology note cites the PDF only. |
| `results/figures/quark/rs_anarchy_mkk_min_hist_by_yprior.pdf` | Moved to `results/figures/quark/exploratory/rs_anarchy_mkk_min_hist_by_yprior.pdf` | Earlier perturbative-axis y-prior plot; replaced in the note by `rs_anarchy_mkk_min_hist_by_yprior_gsstar.pdf`. |
| `results/figures/quark/rs_anarchy_mkk_min_hist_by_yprior.png` | Moved to `results/figures/quark/exploratory/rs_anarchy_mkk_min_hist_by_yprior.png` | Earlier perturbative-axis y-prior plot; replaced in the note by `rs_anarchy_mkk_min_hist_by_yprior_gsstar.pdf`. |
| `results/figures/quark/rs_anarchy_mkk_min_hist_by_yprior_gsstar.png` | Moved to `results/figures/quark/exploratory/rs_anarchy_mkk_min_hist_by_yprior_gsstar.png` | PNG sibling of cited PDF; the methodology note cites the PDF only. |
| `results/figures/quark/rs_anarchy_mkk_min_hist_gsstar.png` | Moved to `results/figures/quark/exploratory/rs_anarchy_mkk_min_hist_gsstar.png` | PNG sibling of cited PDF; the methodology note cites the PDF only. |
| `results/figures/quark/rs_anarchy_pdg_pass_fraction.pdf` | Moved to `results/figures/quark/exploratory/rs_anarchy_pdg_pass_fraction.pdf` | Working summary plot; not cited by the methodology note. |
| `results/figures/quark/rs_anarchy_pdg_pass_fraction.png` | Moved to `results/figures/quark/exploratory/rs_anarchy_pdg_pass_fraction.png` | Working summary plot; not cited by the methodology note. |
| `results/figures/quark/rs_anarchy_per_system_pass_rate.png` | Moved to `results/figures/quark/exploratory/rs_anarchy_per_system_pass_rate.png` | PNG sibling of cited PDF; the methodology note cites the PDF only. |
| `results/figures/quark/yukawa_d_2sigma_vs_1sigma.pdf` | Moved to `results/figures/quark/exploratory/yukawa_d_2sigma_vs_1sigma.pdf` | Legacy Yukawa diagnostic; not cited by the methodology note. |
| `results/figures/quark/yukawa_d_2sigma_vs_1sigma.png` | Moved to `results/figures/quark/exploratory/yukawa_d_2sigma_vs_1sigma.png` | Legacy Yukawa diagnostic; not cited by the methodology note. |
| `results/figures/quark/yukawa_offdiag_spread.pdf` | Moved to `results/figures/quark/exploratory/yukawa_offdiag_spread.pdf` | Working diagnostic from `scripts/yukawa_per_element_anatomy.py`; the note cites `yukawa_size_envelope_vs_anarchic.pdf` instead. |
| `results/figures/quark/yukawa_offdiag_spread.png` | Moved to `results/figures/quark/exploratory/yukawa_offdiag_spread.png` | Working diagnostic from `scripts/yukawa_per_element_anatomy.py`; the note cites `yukawa_size_envelope_vs_anarchic.pdf` instead. |
| `results/figures/quark/yukawa_per_element_anatomy.pdf` | Moved to `results/figures/quark/exploratory/yukawa_per_element_anatomy.pdf` | Working diagnostic from `scripts/yukawa_per_element_anatomy.py`; not implicitly cited by the note. |
| `results/figures/quark/yukawa_per_element_anatomy.png` | Moved to `results/figures/quark/exploratory/yukawa_per_element_anatomy.png` | Working diagnostic from `scripts/yukawa_per_element_anatomy.py`; not implicitly cited by the note. |
| `results/figures/quark/yukawa_size_envelope_vs_anarchic.png` | Moved to `results/figures/quark/exploratory/yukawa_size_envelope_vs_anarchic.png` | PNG sibling of cited PDF; the methodology note cites the PDF only. |
| `results/figures/quark/yukawa_u_2sigma_vs_1sigma.pdf` | Moved to `results/figures/quark/exploratory/yukawa_u_2sigma_vs_1sigma.pdf` | Legacy Yukawa diagnostic; not cited by the methodology note. |
| `results/figures/quark/yukawa_u_2sigma_vs_1sigma.png` | Moved to `results/figures/quark/exploratory/yukawa_u_2sigma_vs_1sigma.png` | Legacy Yukawa diagnostic; not cited by the methodology note. |

### Moved Run A Scratch Exports

The unreferenced, ignored `results/figures/quark/runA/` scratch export directory was moved intact to `results/figures/quark/exploratory/runA/`:

- `rs_anarchy_gate_sensitivity_runA_8Mdraws.pdf`
- `rs_anarchy_gate_sensitivity_runA_8Mdraws.png`
- `rs_anarchy_max_ratio_vs_mkk_runA_8Mdraws.pdf`
- `rs_anarchy_max_ratio_vs_mkk_runA_8Mdraws.png`
- `rs_anarchy_mkk_min_hist_runA_8Mdraws.pdf`
- `rs_anarchy_mkk_min_hist_runA_8Mdraws.png`
- `rs_anarchy_pdg_pass_fraction_runA_8Mdraws.pdf`
- `rs_anarchy_pdg_pass_fraction_runA_8Mdraws.png`
- `rs_anarchy_per_system_pass_rate_runA_8Mdraws.pdf`
- `rs_anarchy_per_system_pass_rate_runA_8Mdraws.png`

## Status as of 2026-05-26 post-cleanup (C19)

Cleanup unit **C19** (Tier-5, final docs sweep) closes the R08-I3 INFO item
("Historical-snapshot figure dirs deliberately left in tree").

- `results/figures/quark_baseline_800k/` was already removed from the working
  tree by a downstream pass; no tracked files reference it (`git ls-files`
  returns no entries under that path).
- `results/figures/quark_pre_audit_constants/` remains in the working tree
  with 58 tracked files totaling ~5 MB. The original `## Summary` line above
  (item 36) noted these were "Considered acceptable forensic state for the
  rc1 window; revisit during repo-hygiene before tagging."

**C19 disposition:** the `quark_pre_audit_constants/` directory is **kept
in place** as a paper-archive companion. It is internally referenced by
provenance discussion of the FLAG-2024 / pre-audit bag-parameter constants
(see `docs/audits/bag_param_inventory.md` and the methodology-note rc1.1
section discussing the hole-#5 invalidation); removing it would lose the
visual record of the pre-audit exclusion boundaries. No methodology-note
figure (`\includegraphics`) target depends on this directory, so the LaTeX
build is unaffected. Tarballing-and-removing would be a follow-up housekeeping
gesture only and is **explicitly deferred** to a future repo-hygiene tag rather
than executed now, on the principle that the rc1 paper archive should keep
the pre-audit visual companion intact.

No `snapshots/` top-level directory exists on disk; the `.gitignore` line
`snapshots/` (line 86) refers to ad-hoc scratch staging under that name and
has no current footprint in the working tree.
