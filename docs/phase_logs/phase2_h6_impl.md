# Phase 2 Hole #6 Implementation Log: Wilson-coefficient RG-running audit

Date: 2026-05-15  
Branch: `audit/wilson-rg`  
Base branch: `paper/quark-scan-2026q2`

## Summary

Audited `quarkConstraints/deltaf2.py::_evolve_wilsons` and the adjacent
`quarkConstraints/qcd_running.py` Wilson-running implementation against BMU
Delta F=2 leading-log conventions and the Phase 2 hole #5 bag-parameter inputs.

Result: numerical Wilson-RG discrepancies were found and fixed.  A peer-review
BLOCKER then exposed a BMU-to-scalar LR sign bug in the first revision; this
revision harmonises the code to the conventional scalar-LR `O5` sign.
The invalidation gate was already tripped by hole #5; this hole still tightens
it with an additional greater-than-10% epsilon_K shift.

## Peer-review revision: Sign-convention chain audit

Diagnosis: **B, sign bug**.  The Wilson ADM used the positive
`Q1_LR^BMU = +2 O5_LR` map, while `_kaon_matrix_elements()` and
`_meson_matrix_elements()` used the conventional positive scalar `O5`
contraction with positive FLAG `B_5` inputs.  The fix keeps the conventional
matrix elements and bag inputs, flips the BMU-to-scalar map to
`Q1_LR^BMU = -2 O5_LR`, and changes the scalar-basis coefficient ADM from
`[[-16, 6], [0, 2]]` to `[[-16, -6], [0, 2]]`.

Action taken:

- Commit `7f71908` `physics(deltaf2): fix BMU LO Wilson RG running`
  established the Wilson-RG audit before peer review.
- Revision commit (this commit)
  `fix(deltaf2): correct BMU-to-scalar-LR sign convention in ADM and matrix elements`
  corrects the LR sign convention, updates fixtures, records the full
  sign-convention chain in `docs/audits/wilson_rg_inventory.md`, and refreshes
  stale broad-suite expectations exposed by the full test run.
- WARNING 2, endpoint mismatch for B/D systems and kaon BSM bags, is deferred
  to a follow-on hole.  This revision does not migrate per-system endpoints.
- WARNING 3 is left as-is because `docs/audits/*` has been added to the
  expected diff scope retroactively.
- NIT 4 is addressed in the methodology appendix by spelling out
  `m_t^{MSbar}(m_t) = 163.5 GeV`.

## Ten-item audit checklist findings

1. Basis identification  
   Answer: OK after documentation fix.  The code basis is BMU VLL/VRR plus
   conventional scalar `O4_LR/O5_LR`, not literal BMU LR labels.  The map is
   `Q1_LR^BMU = -2 O5_LR`, `Q2_LR^BMU = O4_LR`.

2. Sign and factor-of-2 conventions  
   Answer: sign bug found and fixed in the peer-review revision.  The RG
   factor-of-2 map and conventional O5 minus sign are now explicit and tested.
   The tree-level Hamiltonian sign convention in `compute_delta_f2_wilsons`
   remains the repo-owned v1 convention; no matching-sign behavior was changed
   without broader approval.

3. Anomalous-dimension matrix  
   Answer: discrepancy found and fixed.  The paper branch used a non-BMU scalar
   LR matrix and inverse alpha ratio.  The audited code now uses
   `gamma_VLL = 4` and scalar `[C4,C5]` coefficient ADM
   `[[-16, -6], [0, 2]]`, mapped from BMU LR `[[2, 0], [12, -16]]`.

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
   `C5_LR=1 -> (0.894757448992, 0.853891627884)` for `3 TeV -> 2 GeV`, with
   max relative discrepancy `1.300e-16` against the closed-form LO reference.

9. Run direction for chirally flipped operators  
   Answer: discrepancy found and fixed.  High-scale `C5_LR` now feeds low-scale
   `C4_LR`; high-scale `C4_LR` does not generate `C5_LR` at LO in the scalar
   basis.

10. Operator mixing into matrix-element calculation  
    Answer: OK for the scalar `O4/O5` convention, with endpoint caveat.  The
    hadronic contractions pair `C4_LR` with `B4` and conventional positive
    `C5_LR * <O5_LR>` with positive `B5` for kaons and generic mesons.

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

The note now states the conventional scalar basis, BMU-to-scalar LR map, LO
anomalous dimensions, threshold path, `3 TeV -> 2 GeV` reference factors,
LO/NLO status, endpoint caveat, and invalidation-gate impact.

## Verification

- `/n/home09/obarrera/.conda/envs/ising_bootstrap/bin/python scripts/audit_wilson_rg.py`
  -> `C4_LR=1 -> 3.53816397486`,
  `C5_LR=1 -> (0.894757448992, 0.853891627884)`, max relative discrepancy
  `1.300e-16`.
- `/n/home09/obarrera/.conda/envs/ising_bootstrap/bin/python -m pytest -q tests/test_wilson_rg_audit.py tests/test_qcd_running.py tests/test_quark_deltaf2.py`
  -> `35 passed in 2.86s` after fixture update.
- `/n/home09/obarrera/.conda/envs/ising_bootstrap/bin/python -m pytest -q`
  -> `538 passed, 1 skipped in 1040.47s`.
- `pdflatex -interaction=nonstopmode quark_scan_methodology_note.tex` run twice
  from `docs/` -> exit 0 both times; `pdfinfo` reports `Pages: 17` and
  `File size: 580512 bytes`.
- `git diff --check` -> clean.

## Invalidation gate

INVALIDATION_GATE_TRIPPED was already true from Phase 2 hole #5.

PHASE_2_H6_GATE_TIGHTENED: the revised Wilson-RG correction independently changes
epsilon_K by more than 10% at `M_KK = 3 TeV`.  For the standing `r=0.25`
benchmark, the post-hole-#5 pre-Wilson-audit epsilon_K ratio was
`0.5349556334878565`; after the revised Wilson-RG fix it is
`1.9286313761001348` (rounds to `1.93`, inside the peer-review
`1.93`--`2.37` bracket), a factor `3.6052174336881992` increase.
Multiplying by hole #5's `6.238805970149273` factor gives a cumulative
post-#5 plus post-revised-#6 invalidation-gate factor of
`22.492252048980177`.

The cumulative invalidation gate remains TRIPPED.
