# Phase 2 Hole #6 Implementation Log: Wilson-coefficient RG-running audit

Date: 2026-05-15  
Branch: `audit/wilson-rg`  
Base branch: `paper/quark-scan-2026q2`

## Summary

Audited `quarkConstraints/deltaf2.py::_evolve_wilsons` and the adjacent
`quarkConstraints/qcd_running.py` Wilson-running implementation against BMU
Delta F=2 leading-log conventions and the Phase 2 hole #5 bag-parameter inputs.

Result: numerical Wilson-RG discrepancies were found and fixed.  The invalidation
gate was already tripped by hole #5; this hole tightens it with an additional
greater-than-10% epsilon_K shift.

## Ten-item audit checklist findings

1. Basis identification  
   Answer: OK after documentation fix.  The code basis is BMU VLL/VRR plus
   scalar `O4_LR/O5_LR`, not literal BMU LR labels.  The map is
   `Q1_LR^BMU = 2 O5_LR`, `Q2_LR^BMU = O4_LR`.

2. Sign and factor-of-2 conventions  
   Answer: partially OK.  The RG factor-of-2 map is now explicit and tested.
   The tree-level Hamiltonian sign convention in `compute_delta_f2_wilsons`
   remains the repo-owned v1 convention; no matching-sign behavior was changed
   without broader approval.

3. Anomalous-dimension matrix  
   Answer: discrepancy found and fixed.  The paper branch used a non-BMU scalar
   LR matrix and inverse alpha ratio.  The audited code now uses
   `gamma_VLL = 4` and scalar `[C4,C5]` coefficient ADM
   `[[-16, 6], [0, 2]]`, mapped from BMU LR `[[2, 0], [12, -16]]`.

4. RG direction and endpoint  
   Answer: mixed.  Direction is fixed: `alpha_s(mu_high)/alpha_s(mu_low)` is
   used for Wilson evolution.  Endpoint remains globally `mu_had = 2 GeV`, which
   matches kaon `B_K` but not all FLAG BSM/B-system/D-system bag conventions.
   Endpoint migration was documented but not performed.

5. Threshold matching  
   Answer: discrepancy found and fixed.  Wilson RG now crosses `m_t`, `m_b`,
   and `m_c` as applicable; `3 TeV -> 2 GeV` runs through `n_f=6 -> 5 -> 4`.

6. Number of active flavours in alpha_s  
   Answer: OK after fix.  The Wilson wrapper now agrees with the shared
   threshold logic in `qcd.running.alpha_s` at one loop with continuous matching.

7. LO vs NLO  
   Answer: OK as LO; NLO not implemented.  Full BMU NLO/NDR evolution is left
   for an orchestrated upgrade.  The residual LR/scalar NLO effect is recorded
   as an order `10-30%` theory systematic.

8. Numerical sanity check  
   Answer: OK.  The standalone audit gives `C4_LR=1 -> 3.53816397486` and
   `C5_LR=1 -> (-0.894757448992, 0.853891627884)` for `3 TeV -> 2 GeV`, with
   max relative discrepancy `1.300e-16` against the closed-form LO reference.

9. Run direction for chirally flipped operators  
   Answer: discrepancy found and fixed.  High-scale `C5_LR` now feeds low-scale
   `C4_LR`; high-scale `C4_LR` does not generate `C5_LR` at LO in the scalar
   basis.

10. Operator mixing into matrix-element calculation  
    Answer: OK for the scalar `O4/O5` convention, with endpoint caveat.  The
    hadronic contractions pair `C4_LR` with `B4` and `C5_LR` with `B5` for
    kaons and generic mesons.

Full evidence and line references are in
`docs/audits/wilson_rg_inventory.md`; extracted BMU/FLAG reference values and
unit-vector numbers are in `docs/audits/wilson_rg_reference_values.md`.

## Code changes

Implementation and methodology commit: `7f71908`
(`physics(deltaf2): fix BMU LO Wilson RG running`).

- Fixed VLL/VRR Wilson evolution direction.
- Replaced the scalar LR anomalous-dimension matrix with the BMU-mapped
  `[C4_LR, C5_LR]` coefficient matrix.
- Added top-threshold handling to the Wilson-running alpha_s wrapper and
  threshold segmentation.
- Clarified `_evolve_wilsons` and matrix-element docstrings.
- Added `scripts/audit_wilson_rg.py`.
- Added `tests/test_wilson_rg_audit.py`.
- Updated affected Wilson-RG and Delta F=2 regression tests.

No bag parameters, CFW comparison files, or scan-output files were changed.

## Methodology-note update

Commit `7f71908` also added appendix subsection
`Wilson-coefficient RG running` to `docs/quark_scan_methodology_note.tex` and
rebuilt `docs/quark_scan_methodology_note.pdf`.

The note now states the basis, BMU-to-scalar LR map, LO anomalous dimensions,
threshold path, `3 TeV -> 2 GeV` reference factors, LO/NLO status, endpoint
caveat, and invalidation-gate impact.

## Verification

- `/n/home09/obarrera/.conda/envs/ising_bootstrap/bin/python scripts/audit_wilson_rg.py`
  -> max relative discrepancy `1.300e-16`.
- `/n/home09/obarrera/.conda/envs/ising_bootstrap/bin/python -m pytest -q tests/test_wilson_rg_audit.py tests/test_qcd_running.py tests/test_quark_deltaf2.py`
  -> `35 passed in 4.20s`.
- `pdflatex -interaction=nonstopmode quark_scan_methodology_note.tex` run twice
  from `docs/` -> exit 0 both times; `pdfinfo` reports `Pages: 17`.
- `git diff --check` -> clean.

## Invalidation gate

INVALIDATION_GATE_TRIPPED was already true from Phase 2 hole #5.

PHASE_2_H6_GATE_TIGHTENED: the Wilson-RG correction independently changes
epsilon_K by more than 10% at `M_KK = 3 TeV`.  For the standing `r=0.25`
benchmark, the post-hole-#5 pre-Wilson-audit epsilon_K ratio was
`0.5349556334878565`; after the Wilson-RG fix it is `2.367673073409882`, a
factor `4.425924179866831` increase.

Hole #6 ready for peer review.
