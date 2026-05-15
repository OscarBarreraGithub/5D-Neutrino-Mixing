# Zero-Pass Finite-Statistics Inventory

Date: 2026-05-15

Scope: Phase 2 hole #8, finite-statistics treatment of zero-PDG-pass scan
claims.  The inventory was seeded by:

```bash
grep -nE "zero PDG passes?|0/.*PDG|no.*pass|moreUV.*0|moreIR.*0|Run C.*0|factor-1\\.5.*0" docs/quark_scan_methodology_note.tex docs/audits/*.md docs/phase_logs/*.md
```

The Wilson-score upper limits use the repo helper
`quarkConstraints.finite_stats.wilson_upper_limit(k, n)`, with the
audit convention `z = 1.92`.

## Canonical Zero-Pass Ensembles

| Run | Run dir | N total | N pass observed | p_UL_95 | Implied physical bound |
|---|---|---:|---:|---:|---|
| Run 3 moreUV | `scan_outputs/rs_anarchy_run3_moreUV_20260515T085324` | 1,600,000 | 0 | 2.304e-6 | No more than 0.000230% of anarchic draws with the `c -> c + 0.05` pattern and default uniform Y-prior pass the factor-3/5 PDG gate, at 95% confidence. |
| Run 3 moreIR | `scan_outputs/rs_anarchy_run3_moreIR_20260515T085324` | 1,600,000 | 0 | 2.304e-6 | No more than 0.000230% of anarchic draws with the `c -> c - 0.05` pattern and default uniform Y-prior pass the factor-3/5 PDG gate, at 95% confidence. |
| Run C CFW-like gate | `scan_outputs/rs_anarchy_runC_20260515T085323` | 4,000,000 | 0 | 9.216e-7 | No more than 0.0000922% of anarchic draws with the baseline c-pattern and CFW-like wide/floored Y-prior pass the factor-1.5/2.5 PDG gate, at 95% confidence. |

## Provenance

| Run | Gate definition | c-pattern | Y-prior | Seed/tile provenance | Output hash |
|---|---|---|---|---|---|
| Run 3 moreUV | masses within factor 3, CKM within factor 3, J within factor 5 | `c_Q=(0.68,0.62,0.25)`, `c_u=(0.71,0.55,-0.45)`, `c_d=(0.71,0.66,0.60)` | `Re,Im Y_ij ~ U(-1.5,1.5)`, rejection sampler enforces `|Y_ij| >= 0.1` | `base_seed=20260506`; 8 M_KK tiles from 3 to 50 TeV; 200,000 draws/tile | `tile_summary.json` sha256 `3150aadc477183574d07d0be227edc3d646a08c81ce3e668fd5ef8daae3158cd` |
| Run 3 moreIR | masses within factor 3, CKM within factor 3, J within factor 5 | `c_Q=(0.58,0.52,0.15)`, `c_u=(0.61,0.45,-0.55)`, `c_d=(0.61,0.56,0.50)` | `Re,Im Y_ij ~ U(-1.5,1.5)`, rejection sampler enforces `|Y_ij| >= 0.1` | `base_seed=20260506`; 8 M_KK tiles from 3 to 50 TeV; 200,000 draws/tile | `tile_summary.json` sha256 `e6899f1af63821ab7afcd6389126d6522eff2f626ac154541c3e679890aa30d5` |
| Run C CFW-like gate | masses within factor 1.5, CKM within factor 1.5, J within factor 2.5 | `c_Q=(0.63,0.57,0.20)`, `c_u=(0.66,0.50,-0.50)`, `c_d=(0.66,0.61,0.55)` | `Re,Im Y_ij ~ U(-3,3)`, rejection sampler enforces `|Y_ij| >= 0.5` | `base_seed=20260506`; 8 M_KK tiles from 3 to 50 TeV; 500,000 draws/tile | `tile_summary.json` sha256 `92f57214f1e8164eac1dfffeec4c8e26154520d8f8543784f7ea5cead460cf3d` |

## Claim Locations

| File:line | Claim verbatim | Run dir | N total | N pass observed | Gate definition | p_UL_95 |
|---|---|---|---:|---:|---|---:|
| `docs/quark_scan_methodology_note.tex:736-738` | `moreUV ... returned zero PDG-passing draws over 1.6M attempts` | `scan_outputs/rs_anarchy_run3_moreUV_20260515T085324` | 1,600,000 | 0 | factor-3 masses/CKM, factor-5 J | 2.304e-6 |
| `docs/quark_scan_methodology_note.tex:736-738` | `moreIR ... returned zero PDG-passing draws over 1.6M attempts` | `scan_outputs/rs_anarchy_run3_moreIR_20260515T085324` | 1,600,000 | 0 | factor-3 masses/CKM, factor-5 J | 2.304e-6 |
| `docs/quark_scan_methodology_note.tex:752-754` | `In both cases the resulting 1.6M-draw scan returned exactly zero PDG-passing draws.` | `scan_outputs/rs_anarchy_run3_moreUV_20260515T085324` | 1,600,000 | 0 | factor-3 masses/CKM, factor-5 J | 2.304e-6 |
| `docs/quark_scan_methodology_note.tex:752-754` | `In both cases the resulting 1.6M-draw scan returned exactly zero PDG-passing draws.` | `scan_outputs/rs_anarchy_run3_moreIR_20260515T085324` | 1,600,000 | 0 | factor-3 masses/CKM, factor-5 J | 2.304e-6 |
| `docs/phase_logs/invalidation_gate_rerun.md:93` | `Run 3 moreUV | zero PDG passes | zero PDG passes | zero PDG passes | zero PDG passes` | `scan_outputs/rs_anarchy_run3_moreUV_20260515T085324` | 1,600,000 | 0 | factor-3 masses/CKM, factor-5 J | 2.304e-6 |
| `docs/phase_logs/invalidation_gate_rerun.md:94` | `Run 3 moreIR | zero PDG passes | zero PDG passes | zero PDG passes | zero PDG passes` | `scan_outputs/rs_anarchy_run3_moreIR_20260515T085324` | 1,600,000 | 0 | factor-3 masses/CKM, factor-5 J | 2.304e-6 |
| `docs/phase_logs/invalidation_gate_rerun.md:95` | `Run C | zero PDG passes | zero PDG passes | zero PDG passes | zero PDG passes` | `scan_outputs/rs_anarchy_runC_20260515T085323` | 4,000,000 | 0 | factor-1.5 masses/CKM, factor-2.5 J | 9.216e-7 |
| `docs/phase_logs/invalidation_gate_signoff.md:76-79` | `moreUV / moreIR / Run C still zero PDG passes ... all three return exactly 0` | `scan_outputs/rs_anarchy_run3_moreUV_20260515T085324` | 1,600,000 | 0 | factor-3 masses/CKM, factor-5 J | 2.304e-6 |
| `docs/phase_logs/invalidation_gate_signoff.md:76-79` | `moreUV / moreIR / Run C still zero PDG passes ... all three return exactly 0` | `scan_outputs/rs_anarchy_run3_moreIR_20260515T085324` | 1,600,000 | 0 | factor-3 masses/CKM, factor-5 J | 2.304e-6 |
| `docs/phase_logs/invalidation_gate_signoff.md:76-79` | `moreUV / moreIR / Run C still zero PDG passes ... all three return exactly 0` | `scan_outputs/rs_anarchy_runC_20260515T085323` | 4,000,000 | 0 | factor-1.5 masses/CKM, factor-2.5 J | 9.216e-7 |

## Grep Hits Not Treated As Canonical Zero-Pass Rows

| File:line | Matched text | Disposition |
|---|---|---|
| `docs/quark_scan_methodology_note.tex:487-493` | `purely anarchic draws never land within 2\sigma ... zero passing draws` | Absolute language is not tied to a canonical output with N/seed/hash in the current artifacts.  The methodology note is rewritten to avoid quoting it as a zero-pass measurement. |
| `docs/quark_scan_methodology_note.tex:506` | `nonetheless passes the PDG check` | False-positive grep hit, not a zero-pass claim. |
| `docs/phase_logs/invalidation_gate_signoff.md:52` | `rejected PDG-passing draws` | False-positive grep hit, not a zero-pass claim. |
| Other phase-log grep hits without zero-pass language | branch/test/status text | False-positive grep hits, not zero-pass claims. |
