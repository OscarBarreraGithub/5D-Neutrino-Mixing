# Wilson-Coefficient RG-Running Audit Inventory

Audit date: 2026-05-15  
Branch: `audit/wilson-rg`  
Code target: `quarkConstraints/deltaf2.py::_evolve_wilsons` plus
`quarkConstraints/qcd_running.py`

The audited code path computes tree-level Delta F=2 Wilson coefficients at
`M_KK`, evolves them with leading-log QCD, and contracts them with
hadronic matrix elements.  Exact reference numbers are in
`docs/audits/wilson_rg_reference_values.md`.

## Inventory

| Item | Finding | Evidence |
|---|---|---|
| Basis tag | Hybrid code basis: BMU VLL/VRR current-current operators plus scalar `O4_LR/O5_LR` paired with the code's `B4/B5` matrix elements. | `quarkConstraints/qcd_running.py:7`, `quarkConstraints/qcd_running.py:11`, `quarkConstraints/deltaf2.py:322`, `quarkConstraints/deltaf2.py:344`, `quarkConstraints/deltaf2.py:664` |
| Wilson matching | Tree-level `C1_VLL`, `C1_VRR`, `C4_LR`, and `C5_LR` are matched at `M_KK`; `C4_LR = -LR/M_KK^2`, `C5_LR = LR/(3 M_KK^2)`. | `quarkConstraints/deltaf2.py:313`, `quarkConstraints/deltaf2.py:341` |
| RG endpoint | Default endpoint is one global `mu_had = 2.0 GeV` for all systems. | `quarkConstraints/deltaf2.py:449`, `quarkConstraints/deltaf2.py:456`, `quarkConstraints/deltaf2.py:503`, `quarkConstraints/deltaf2.py:571`, `quarkConstraints/deltaf2.py:588` |
| LO/NLO flag | LO only: one-loop anomalous dimensions, one-loop alpha_s in `quarkConstraints/qcd_running.py`, continuous LO Wilson matching. | `quarkConstraints/qcd_running.py:3`, `quarkConstraints/qcd_running.py:20`, `quarkConstraints/qcd_running.py:44`, `quarkConstraints/qcd_running.py:85`, `quarkConstraints/qcd_running.py:226` |
| Active flavors | Fixed path includes top, bottom, and charm thresholds.  `3 TeV -> 2 GeV` crosses `n_f=6 -> 5 -> 4`; charm is below endpoint. | `quarkConstraints/qcd_running.py:36`, `quarkConstraints/qcd_running.py:65`, `quarkConstraints/qcd_running.py:144`, `quarkConstraints/qcd_running.py:276`, `tests/test_wilson_rg_audit.py:49` |
| Anomalous dimensions | Post-audit reference: `gamma_VLL = 4`; scalar LR coefficient ADM `[[-16, 6], [0, 2]]`, mapped from BMU LR `[ [2,0], [12,-16] ]`. | `quarkConstraints/qcd_running.py:44`, `quarkConstraints/qcd_running.py:54` |
| Evolution direction | Post-audit convention uses `alpha_s(mu_high)/alpha_s(mu_low)` for both VLL and LR. | `quarkConstraints/qcd_running.py:172`, `quarkConstraints/qcd_running.py:183`, `quarkConstraints/qcd_running.py:186`, `quarkConstraints/qcd_running.py:198` |
| Matrix elements | Kaon and generic meson contractions pair `C4_LR` with `B4` and `C5_LR` with `B5`. | `quarkConstraints/deltaf2.py:675`, `quarkConstraints/deltaf2.py:693`, `quarkConstraints/deltaf2.py:821`, `quarkConstraints/deltaf2.py:860` |
| alpha_s side feed | Matching couplings use the shared high-precision `qcd.alpha_s`; Wilson RG uses its own LO wrapper but was regression-tested against `qcd.running.alpha_s(..., n_loops=1, matching_loops=0)`. | `quarkConstraints/couplings.py:145`, `qcd/running.py:118`, `qcd/constants.py:28`, `tests/test_wilson_rg_audit.py:55` |

## Ten checklist findings

1. Basis identification  
   Answer: OK after documentation fix.  The code is not literal five-operator
   BMU notation; it uses BMU VLL/VRR and scalar `O4/O5` LR operators.  The
   mapping is now explicit: `Q1_LR^BMU = 2 O5_LR`,
   `Q2_LR^BMU = O4_LR`.  Evidence:
   `quarkConstraints/qcd_running.py:7`, `quarkConstraints/qcd_running.py:11`,
   `quarkConstraints/deltaf2.py:664`.  Recommended fix: no further code change
   for the basis label; keep the reference note.

2. Sign and factor-of-2 conventions  
   Answer: partially OK; residual convention risk documented.  The Wilson
   evolution now includes the factor-of-2 map between BMU `Q1_LR` and the
   scalar `O5_LR`.  The tree-level matching still sets `C4_LR = -LR/M_KK^2`
   and `C5_LR = LR/(3 M_KK^2)` in the repo-owned v1 convention, and the audit
   did not change Hermitian-conjugate or Hamiltonian sign conventions.  Evidence:
   `quarkConstraints/deltaf2.py:341`, `quarkConstraints/deltaf2.py:344`,
   `quarkConstraints/qcd_running.py:50`.  Recommended fix: if the orchestrator
   wants a full Hamiltonian-sign comparison to CFW/BBL, do it as a separate
   convention audit before changing matching signs.

3. Anomalous-dimension matrix  
   Answer: discrepancy found and fixed.  The pre-audit LR matrix on the paper
   branch was `[[8, -4], [-16/3, -28/3]]` and the evolution used the inverse
   alpha ratio.  The post-audit code uses `gamma_VLL = 4` and scalar
   `gamma_LR = [[-16, 6], [0, 2]]`, which is the BMU LO LR coefficient block
   conjugated into `[C4_LR, C5_LR]`.  Evidence:
   `quarkConstraints/qcd_running.py:44`, `quarkConstraints/qcd_running.py:54`,
   `quarkConstraints/qcd_running.py:179`, `quarkConstraints/qcd_running.py:193`.
   Recommended fix: completed; covered by `tests/test_wilson_rg_audit.py:28`.

4. RG direction and endpoint  
   Answer: mixed.  Direction is OK after the fix: Wilsons run from
   `matching_scale` down to `mu_had` using the BMU LO alpha ratio.  Endpoint
   remains a discrepancy for non-kaon systems: the global default is `2 GeV`,
   while hole #5 found `B_d/B_s` bags quoted at `mu = m_b` and D-system bags
   in comparable sources at `3 GeV`.  Kaon `B_K` matches `2 GeV`, but kaon
   BSM `B4/B5` are FLAG values quoted at `3 GeV`.  Evidence:
   `quarkConstraints/deltaf2.py:456`, `quarkConstraints/deltaf2.py:612`,
   `quarkConstraints/deltaf2.py:614`, `quarkConstraints/deltaf2.py:625`,
   `quarkConstraints/deltaf2.py:657`.  Recommended fix: do not change in this
   hole; a per-system endpoint migration would invalidate all scan outputs.

5. Threshold matching  
   Answer: discrepancy found and fixed for Wilson RG.  The paper-branch Wilson
   wrapper crossed only `m_b` and `m_c`, so `3 TeV -> 2 GeV` incorrectly used
   `n_f=5` above the top threshold.  The audited path now crosses `m_t`,
   `m_b`, and `m_c` as applicable; for kaons to `2 GeV`, the crossed sequence is
   `6 -> 5 -> 4`.  Evidence: `quarkConstraints/qcd_running.py:36`,
   `quarkConstraints/qcd_running.py:65`, `quarkConstraints/qcd_running.py:276`.
   Recommended fix: completed; test coverage at `tests/test_wilson_rg_audit.py:49`.

6. Number of active flavours in alpha_s  
   Answer: OK after Wilson wrapper fix.  At `m_t`, the wrapper switches
   `6 -> 5`; at `m_b`, `5 -> 4`; at `m_c`, `4 -> 3`.  The independent
   `qcd.running.alpha_s` infrastructure already has the same threshold list.
   Evidence: `quarkConstraints/qcd_running.py:65`, `qcd/constants.py:28`,
   `qcd/running.py:136`, `tests/test_wilson_rg_audit.py:55`.  Recommended fix:
   completed.

7. LO vs NLO  
   Answer: OK as an LO implementation; NLO not implemented.  The code uses
   one-loop ADM, one-loop alpha_s, and continuous LO threshold matching.
   BMU arXiv:hep-ph/0102316 provides NLO/NDR master formulae and emphasizes
   the largest RG effects in LR/scalar sectors.  A precise NLO comparison for
   `3 TeV -> 2 GeV` was not implemented in this hole because the full NLO
   matrices and scheme-dependent matching would be a larger upgrade.  Estimated
   NLO residual: order `10-30%` for LR/scalar coefficients at these scales.
   Evidence: `quarkConstraints/qcd_running.py:3`,
   `quarkConstraints/qcd_running.py:85`, `quarkConstraints/qcd_running.py:226`.
   Recommended fix: keep LO plus an explicit methodology note; decide later
   whether to add a full NLO WET evolution module.

8. Numerical sanity check  
   Answer: OK after the fix.  The standalone script gives, for
   `3 TeV -> 2 GeV`, `C4_LR=1 -> C4_LR(2 GeV)=3.53816397486` and
   `C5_LR=1 -> (C4_LR,C5_LR)=(-0.894757448992, 0.853891627884)`.  The maximum
   relative discrepancy vs the closed-form LO reference is `1.300e-16`.
   Evidence: `scripts/audit_wilson_rg.py:85`, `scripts/audit_wilson_rg.py:95`,
   `scripts/audit_wilson_rg.py:132`, `tests/test_wilson_rg_audit.py:28`.
   Recommended fix: completed.

9. Run direction for chirally flipped operators  
   Answer: discrepancy found and fixed.  The scalar LR block now has the large
   mapped BMU off-diagonal `gamma_45 = 6`, so high-scale `C5_LR` feeds
   low-scale `C4_LR`.  The previous paper-branch matrix also generated
   `C4 -> C5`, which is not the LO scalar-basis direction implied by the BMU
   map.  Evidence: `quarkConstraints/qcd_running.py:50`,
   `quarkConstraints/qcd_running.py:54`, `tests/test_qcd_running.py:243`.
   Recommended fix: completed.

10. Operator mixing into matrix-element calculation  
    Answer: OK for the scalar `O4/O5` convention, with endpoint caveat from
    item 4.  The hadronic contraction uses `C4_LR * O4_LR(B4)` and
    `C5_LR * O5_LR(B5)` for kaons and generic mesons.  Evidence:
    `quarkConstraints/deltaf2.py:675`, `quarkConstraints/deltaf2.py:677`,
    `quarkConstraints/deltaf2.py:692`, `quarkConstraints/deltaf2.py:860`.
    Recommended fix: completed docstring cleanup; do not alter bag constants.

## Invalidation-gate note

Hole #5 already tripped the gate.  Hole #6 tightens it: at `M_KK = 3 TeV`, the
audited `C4_LR` unit running factor changes from the paper-branch value
`1.96025649` (with an induced `C5_LR=-0.43912562`) to `3.53816397486` with no
`C4 -> C5` leakage.  For the existing `r=0.25` benchmark used in
`tests/test_quark_deltaf2.py`, the post-hole-5, pre-hole-6 epsilon_K ratio was
`0.5349556334878565`; after the Wilson-RG fix it is `2.367673073409882`, a
factor `4.425924179866831` increase.  This is an additional greater-than-10%
Wilson-RG shift on top of the bag-parameter invalidation.
