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
| Basis tag | Hybrid code basis: BMU VLL/VRR current-current operators plus conventional scalar `O4_LR/O5_LR` paired with the code's `B4/B5` matrix elements. | `quarkConstraints/qcd_running.py:7`, `quarkConstraints/qcd_running.py:11`, `quarkConstraints/deltaf2.py:322`, `quarkConstraints/deltaf2.py:344`, `quarkConstraints/deltaf2.py:670` |
| Wilson matching | Tree-level `C1_VLL`, `C1_VRR`, `C4_LR`, and `C5_LR` are matched at `M_KK`; `C4_LR = -LR/M_KK^2`, `C5_LR = LR/(3 M_KK^2)`. | `quarkConstraints/deltaf2.py:313`, `quarkConstraints/deltaf2.py:341` |
| RG endpoint | Default endpoint is one global `mu_had = 2.0 GeV` for all systems. | `quarkConstraints/deltaf2.py:449`, `quarkConstraints/deltaf2.py:456`, `quarkConstraints/deltaf2.py:503`, `quarkConstraints/deltaf2.py:571`, `quarkConstraints/deltaf2.py:588` |
| LO/NLO flag | LO only: one-loop anomalous dimensions, one-loop alpha_s in `quarkConstraints/qcd_running.py`, continuous LO Wilson matching. | `quarkConstraints/qcd_running.py:3`, `quarkConstraints/qcd_running.py:20`, `quarkConstraints/qcd_running.py:44`, `quarkConstraints/qcd_running.py:85`, `quarkConstraints/qcd_running.py:226` |
| Active flavors | Fixed path includes top, bottom, and charm thresholds.  `3 TeV -> 2 GeV` crosses `n_f=6 -> 5 -> 4`; charm is below endpoint. | `quarkConstraints/qcd_running.py:36`, `quarkConstraints/qcd_running.py:65`, `quarkConstraints/qcd_running.py:144`, `quarkConstraints/qcd_running.py:276`, `tests/test_wilson_rg_audit.py:49` |
| Anomalous dimensions | Revised post-review reference: `gamma_VLL = 4`; conventional scalar LR coefficient ADM `[[-16, -6], [0, 2]]`, mapped from BMU LR `[ [2,0], [12,-16] ]` with `Q1_LR^BMU = -2 O5_LR`. | `quarkConstraints/qcd_running.py:44`, `quarkConstraints/qcd_running.py:54` |
| Evolution direction | Post-audit convention uses `alpha_s(mu_high)/alpha_s(mu_low)` for both VLL and LR. | `quarkConstraints/qcd_running.py:172`, `quarkConstraints/qcd_running.py:183`, `quarkConstraints/qcd_running.py:186`, `quarkConstraints/qcd_running.py:198` |
| Matrix elements | Post-B3 M12-ready GGMS contractions pair `C4_LR` with `B4` and `C5_LR` with `B5`: `O1 = (1/3) f^2 m B1`, `O4 = (r_chi/4 + 1/24) f^2 m B4`, `O5 = (r_chi/12 + 1/8) f^2 m B5`. The earlier pre-B3 audit text used `2/3` for `O1` and had the LR coefficients swapped and doubled. | `quarkConstraints/deltaf2.py:1230`, `quarkConstraints/deltaf2.py:1231`, `quarkConstraints/deltaf2.py:1232`, `quarkConstraints/deltaf2.py:1451`, `quarkConstraints/deltaf2.py:1452`, `quarkConstraints/deltaf2.py:1453` |
| alpha_s side feed | Matching couplings use the shared high-precision `qcd.alpha_s`; Wilson RG uses its own LO wrapper but was regression-tested against `qcd.running.alpha_s(..., n_loops=1, matching_loops=0)`. | `quarkConstraints/couplings.py:145`, `qcd/running.py:118`, `qcd/constants.py:28`, `tests/test_wilson_rg_audit.py:55` |

## Sign-convention chain audit

Diagnosis: **B, sign bug**.  The initial hole #6 revision used the nonstandard
positive BMU-to-scalar map in the ADM, while the hadronic contraction and
FLAG `B_5` input were already being used with the conventional positive scalar
matrix-element sign.  This revision harmonises the code to the conventional
scalar-LR sign convention by flipping only the scalar ADM off-diagonal.

1. Operator definitions
   The original `deltaf2.py` implementation had only an implicit scalar
   `O4_LR/O5_LR` definition, deduced from the `B4/B5` contraction.  The current
   code now documents the explicit conventional scalar operators in
   `quarkConstraints/deltaf2.py:580` and `quarkConstraints/deltaf2.py:678`:
   `O4_LR = (bar h^alpha P_L q^alpha)(bar h^beta P_R q^beta)` and
   `O5_LR = (bar h^alpha P_L q^beta)(bar h^beta P_R q^alpha)`.  With BMU
   `Q1_LR` as a vector-LR operator, the Fierz map is
   `Q1_LR^BMU = -2 O5_LR`, `Q2_LR^BMU = O4_LR`; this is recorded in
   `quarkConstraints/qcd_running.py:7` and `quarkConstraints/qcd_running.py:50`.

2. Matrix-element contraction
   `_kaon_matrix_elements()` uses the post-B3 positive expressions
   `<O4_LR> = (r_chi / 4 + 1/24) f_K^2 m_K B_4_K` and
   `<O5_LR> = (r_chi / 12 + 1/8) f_K^2 m_K B_5_K` at
   `quarkConstraints/deltaf2.py:1231` and `quarkConstraints/deltaf2.py:1232`;
   generic B/D mesons use the same M12-ready form at
   `quarkConstraints/deltaf2.py:1452` and `quarkConstraints/deltaf2.py:1453`.
   The June-2026 B3 audit and 2026-07 full-repo audit corrected the stale
   pre-B3 text that put the large coefficient on `O5` and doubled both LR
   coefficients.  The contraction still uses the conventional positive scalar-B5
   sign, not a sign-flipped nonstandard O5 contraction.

3. Bag parameter `B_5`
   `B_5_K = 0.691` is the positive FLAG 2024 input at
   `quarkConstraints/deltaf2.py:620`.  Because the contraction is kept in the
   conventional scalar sign, the bag input is not sign-flipped.  A sign flip
   would only be required if the code chose the nonstandard
   `O5_nonstandard = -O5_conventional` convention, which it no longer does.

4. Wilson coefficient `C5_LR`
   Tree matching still produces `C5_LR = LR / (3 M_KK^2)` in the repo v1
   Hamiltonian convention at `quarkConstraints/deltaf2.py:344`.  RG running now
   maps conventional scalar coefficients through
   `C_BMU = (-C5_LR / 2, C4_LR)`, so the LO scalar coefficient ADM is
   `[[-16, -6], [0, 2]]` at `quarkConstraints/qcd_running.py:54`.  The evolved
   `C5_LR` then multiplies the conventional positive `<O5_LR>` contraction.

5. Propagation check
   The standalone audit now gives `C4_LR=1 -> C4_LR(2 GeV)=3.53816397486` and
   `C5_LR=1 -> (C4_LR,C5_LR)=(0.894757448992, 0.853891627884)`, with max
   relative discrepancy `1.300e-16`.  The `C4_LR=1` diagonal evolution is
   unchanged; the `C5_LR -> C4_LR` cross-feed has changed sign, as expected.

## Ten checklist findings

1. Basis identification  
   Answer: OK after documentation fix.  The code is not literal five-operator
   BMU notation; it uses BMU VLL/VRR and conventional scalar `O4/O5` LR
   operators.  The mapping is now explicit: `Q1_LR^BMU = -2 O5_LR`,
   `Q2_LR^BMU = O4_LR`.  Evidence:
   `quarkConstraints/qcd_running.py:7`, `quarkConstraints/qcd_running.py:11`,
   `quarkConstraints/deltaf2.py:670`.  Recommended fix: no further code change
   for the basis label; keep the reference note.

2. Sign and factor-of-2 conventions  
   Answer: sign bug found and fixed.  The Wilson evolution now includes the
   conventional minus sign and factor-of-2 map between BMU `Q1_LR` and scalar
   `O5_LR`.  The tree-level matching still sets `C4_LR = -LR/M_KK^2`
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
   `gamma_LR = [[-16, -6], [0, 2]]`, which is the BMU LO LR coefficient block
   conjugated into the conventional `[C4_LR, C5_LR]` basis.  Evidence:
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
   Answer: OK after the sign revision.  The standalone script gives, for
   `3 TeV -> 2 GeV`, `C4_LR=1 -> C4_LR(2 GeV)=3.53816397486` and
   `C5_LR=1 -> (C4_LR,C5_LR)=(0.894757448992, 0.853891627884)`.  The maximum
   relative discrepancy vs the closed-form LO reference is `1.300e-16`.
   Evidence: `scripts/audit_wilson_rg.py:85`, `scripts/audit_wilson_rg.py:95`,
   `scripts/audit_wilson_rg.py:132`, `tests/test_wilson_rg_audit.py:28`.
   Recommended fix: completed.

9. Run direction for chirally flipped operators  
   Answer: discrepancy found and fixed.  The scalar LR block now has the large
   mapped BMU off-diagonal `gamma_45 = -6`, so high-scale `C5_LR` feeds
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

## Status as of 2026-05-25 post-cleanup (C03)

The C03 cleanup unit (R04-I1..I4, Tier 1) closed the three Wilson-RG code /
documentation items assigned to it (`R04-I1`, `R04-I2`, `R04-I3`).  The four
items tagged `Recommended fix: defer / decide later` in the ten-item checklist
above remain explicitly deferred (R04-I4, reassigned to **C19** —
methodology-note polish — in `.orchestration/CLEANUP_PLAN.md` §C C19):

1. **Item 2 — Hamiltonian/CFW comparison.**  A full Hamiltonian-sign comparison
   to CFW/BBL was deferred to a separate convention audit; no change to the
   matching signs in this cleanup wave.
2. **Item 4 — Per-system RG endpoint migration.**  Non-kaon hadronic bag
   inputs (B_d, B_s at `mu = m_b`; D0 at `3 GeV`; kaon B4/B5 at `3 GeV`) are
   still evolved to the global default `mu_had = 2 GeV`.  Per-system endpoint
   alignment would invalidate all existing Delta F=2 scan outputs and is
   tracked as a separate orchestrated change.
3. **Item 7 — LO vs NLO.**  The QCD running stays LO with one-loop ADM,
   one-loop alpha_s, and continuous LO threshold matching.  A full NLO WET
   evolution module remains a larger upgrade; estimated NLO residual is
   `10-30%` for LR/scalar coefficients at `3 TeV -> 2 GeV`.
4. **Item 10 — Endpoint caveat for matrix-element contraction.**  The hadronic
   contraction is OK for the scalar `O4/O5` convention, with the per-system
   endpoint caveat from item 4 above; bag constants are not altered.

C03 explicitly verified, with no code change required, that the
`paper_0710_1869` LR-sign convention (R04-I3) remains consistent with the
post-C01 (R03-I1, FLAG 2024) canonical hadronic constants:

- Matching: `q4_lr = -(L*R) * 1/M_KK^2`, `q5_lr = +(L*R) / (3 M_KK^2)`
  (`quarkConstraints/paper_0710_1869/eft_deltaf2/matching_kkgluon.py:525-526`).
- Hadronic matrix elements: `<Q1>_GeV4 = +(2/3) f_K^2 m_K^2 B_K`,
  `<Q4>_GeV4 = +(1/2) m_K^2 f_K^2 R_chi B4`, and
  `<Q5>_GeV4 = +(1/6) m_K^2 f_K^2 R_chi B5`
  (`quarkConstraints/paper_0710_1869/eft_deltaf2/hadronic.py:765,2237,2248`).
- Contraction: `q4_lr * <Q4>_GeV4 + q5_lr * <Q5>_GeV4` with no extra sign
  and with the non-M12-ready GeV4 matrix elements divided by `2 m_K`
  (`quarkConstraints/paper_0710_1869/eft_deltaf2/observables.py:911-914`).
- BMU LR LO ADM stored as the Wilson-coefficient (upper-triangular)
  transpose `((2, 12), (0, -16))`
  (`quarkConstraints/paper_0710_1869/eft_deltaf2/rg.py:98-101`) is the
  transpose of the operator ADM `[[2, 0], [12, -16]]` cited in
  `quarkConstraints/qcd_running.py:50-57`, consistent with the Wilson-vs.
  operator ADM convention difference (see audit script
  `scripts/audit_wilson_rg.py:85-117`).

The R04-I3 closure adds no new sign-flip beyond the post-C01 cleanup; the
paper_0710_1869 slice and the repo-v1 `quarkConstraints/deltaf2.py` path
both use the conventional positive scalar `O4/O5` LR sign with positive bag
inputs.

The C03 code changes themselves are: (R04-I1) a defensive upper-triangular
guard added to `scripts/audit_wilson_rg.py::scalar_lr_segment_matrix`;
(R04-I2) the function-local `from .qcd_running import evolve_deltaf2_wilsons`
in `quarkConstraints/deltaf2.py::_evolve_wilsons` lifted to the module top
(unconditional per reviewer M-7; `qcd_running` does not import `deltaf2` at
module scope).  The audit script numeric output is unchanged
(`C4_LR=1 -> 3.53816397486`, `C5_LR=1 -> (0.894757448992, 0.853891627884)`,
max relative discrepancy `1.300e-16`).

## Invalidation-gate note

Hole #5 already tripped the gate.  Revised hole #6 tightens it: at `M_KK = 3 TeV`, the
audited `C4_LR` unit running factor changes from the paper-branch value
`1.96025649` (with an induced `C5_LR=-0.43912562`) to `3.53816397486` with no
`C4 -> C5` leakage.  For the existing `r=0.25` benchmark used in
`tests/test_quark_deltaf2.py`, the post-hole-5, pre-hole-6 epsilon_K ratio was
`0.5349556334878565`; after the revised Wilson-RG fix it is
`1.9286313761001348`, a factor `3.6052174336881992` increase.  Combined with
hole #5's `6.238805970149273` budget/input factor, the cumulative invalidation
gate factor is `22.492252048980177`.  This is an additional greater-than-10%
Wilson-RG shift on top of the bag-parameter invalidation.
