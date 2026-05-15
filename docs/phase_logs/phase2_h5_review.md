# Phase 2 Hole #5 Peer Review

Date: 2026-05-15
Reviewer: Phase 2 hole #5 peer reviewer

### Verdict

APPROVE-WITH-CONCERNS

### Item-by-item (1-12)

1. PASS - Branch state. `git status --short --branch` reports
   `paper/quark-scan-2026q2...origin/paper/quark-scan-2026q2`.
   `paper/quark-scan-2026q2`, `origin/paper/quark-scan-2026q2`,
   `audit/bag-inputs`, and `origin/audit/bag-inputs` all resolve to
   `695f35ea238c6c3dea7c173e554bdcf03fee654e`. `git merge-base --is-ancestor
   audit/bag-inputs paper/quark-scan-2026q2` exits 0, so the audit branch is
   FF-merged. `git merge-base --is-ancestor 1fec2c7 HEAD` also exits 0, so
   Phase 2 hole #4 remains in the chain.

2. PASS - Pre-audit constants. `git show 1fec2c7:quarkConstraints/deltaf2.py
   | grep -nE "B_1_K|B_4_K|B_5_K|EPSILON_K"` gives `B_1_K = 0.717`,
   `B_4_K = 0.78`, `B_5_K = 0.57`, and `EPSILON_K_SM = 1.81e-3`, matching
   the report's old-value claim.

3. PASS - Post-audit constants. Current `grep -nE
   "B_1_K|B_4_K|B_5_K|EPSILON_K" quarkConstraints/deltaf2.py` gives
   `B_1_K = 0.5503`, `B_4_K = 0.903`, `B_5_K = 0.691`, and
   `EPSILON_K_SM = 2.161e-3`.

4. PASS - Source verification, FLAG 2024. Web/source check:
   https://arxiv.org/abs/2411.04268 identifies the FLAG Review 2024 and its
   kaon `B_K` and BSM bag-parameter review. Independent PDF text extraction
   found `B_K^MS(2 GeV) = 0.5503(66)` and the Nf=2+1 BSM estimates
   `B_4 = 0.903(14)`, `B_5 = 0.691(14)`. The cited code values match the FLAG
   central values exactly, well within 1%.

5. PASS - Source verification, BGS 2020. Web/source check:
   https://arxiv.org/abs/1911.06822 reports the BGS SM prediction
   `epsilon_K = 2.16(6)(8)(15) x 10^-3`; the PDF text around Eq. (20) gives
   central value `2.161 x 10^-3` and combined form `2.16(18) x 10^-3`.
   The audit's `2.161e-3` central value is therefore source-consistent. No
   blocker on the central value.

6. PASS - 6x budget shrinkage arithmetic, with concern. Using the provided
   PDG/BGS central values, `|2.228e-3 - 2.161e-3| = 6.7e-5`. The old central
   budget was `2.228e-3 - 1.81e-3 = 4.18e-4`, and
   `4.18e-4 / 6.7e-5 = 6.23880597015`, matching the report's 6.24 factor.
   The audit acknowledges this as a central-value shrinkage, but it does not
   acknowledge or propagate an uncertainty on `|exp - SM|`. Important detail:
   using the BGS published total uncertainty `0.18e-3`, the uncertainty on the
   central gap is much larger than 13%; even PDG experimental error alone is
   `0.011e-3 / 0.067e-3 = 16%`. If the orchestrator intends a narrower SM-error
   convention, it needs to be documented explicitly.

7. FAIL - CKMfitter / UTfit cross-check documentation. I found no defense of
   choosing BGS as the paper's epsilon_K SM reference versus CKMfitter/UTfit or
   a literature-average budget. `rg -n "CKMfitter|UTfit|BGS|Brod|conservative"`
   finds the old CKMfitter comment/provenance and the BGS citation, but no
   explanation that BGS is intentionally conservative or why it is the right
   central value for final paper claims.

8. PASS - Methodology note update. `grep -nE
   "4\.18\\times|6\.7\\times|6\.7" docs/quark_scan_methodology_note.tex` finds
   the ratio definition using `6.7\times10^{-5}` and the later audit sentence
   saying the budget changed from `4.18\times10^{-4}` to `6.7\times10^{-5}`.
   That is consistent with the audited constants.

9. PASS - Yellow-row drift spot-check. I spot-checked `B_1_BD`: HPQCD 2019
   Table VI gives `B_d` operator-1 bag parameter `0.806(40)` in the
   scheme-matched `mu=m_b` table. FLAG 2024 quotes HPQCD 19A as the available
   Nf=2+1+1 B-mixing estimate, though FLAG's directly quoted average is the
   RGI `Bhat` quantity rather than the code's `B_1(m_b)` convention. The
   inventory's scheme-matched comparison is honest:
   `abs(0.87 - 0.806) / 0.806 = 7.94%`, i.e. Yellow.

10. PASS - Code scope. `git show --stat dc9c498` touches only
    `quarkConstraints/deltaf2.py` and `tests/test_quark_deltaf2.py`.
    `git show --stat 82a96f0` touches `docs/audits/bag_param_inventory.md`,
    `docs/quark_scan_methodology_note.tex`, and the rebuilt
    `docs/quark_scan_methodology_note.pdf`. No Wilson-coefficient RG, CFW, or
    scan-output drift was present.

11. PASS - Tests added and run. `tests/test_quark_deltaf2.py` contains
    `test_audited_deltaf2_hadronic_constants_match_selected_sources`.
    `pytest -q -k audited_deltaf2_hadronic` passed:
    `1 passed, 536 deselected in 157.32s`.

12. PASS - Invalidation-gate amplitude estimates. The `B_1_K` fixed-Wilson
    check gives `(0.5503 - 0.717) / 0.717 = -0.2324965`, so the 23.25%
    amplitude shift magnitude is consistent and the sign is a decrease in the
    one-operator NP amplitude. The report then says the epsilon_K central
    budget shrinks by 6.24x, increasing `epsilon_k_ratio` by 6.24x. That net
    direction is correct: the B1 amplitude decrease is outweighed by the much
    smaller central NP budget, so the net ratio increases by about
    `0.7675 * 6.2388 = 4.79` for that isolated O1 probe.

### Findings

1. WARNING - The BGS 2020 epsilon_K choice is not defended against CKMfitter,
   UTfit, or a literature-average central value. Recommended fix: add a short
   decision note explaining whether BGS is being used because it is the newest
   perturbative calculation, because it gives a conservative stronger bound, or
   because the paper intends to quote a central-value-only bound. Also state how
   claims change for a looser `1e-4` to `3e-4` NP budget.

2. WARNING - The epsilon_K NP budget is treated as an exact central-value gap.
   Recommended fix: propagate or bracket the uncertainty in
   `|epsilon_K^exp - epsilon_K^SM|` into the M_KK bound. For small relative
   budget errors, the mass-bound uncertainty is roughly half the budget
   uncertainty because the NP amplitude scales as `1/M_KK^2`; if the full BGS
   quoted uncertainty is used, the central gap is not a stable Gaussian budget
   and should be handled as a sensitivity/profiling choice rather than a single
   precise denominator.

### Physics decision flags

- BGS 2020 as epsilon_K^SM: the central value `2.161e-3` is source-correct and
  is defensible as a selected modern perturbative prediction. It is not yet
  documented as the appropriate choice for this paper's claims. Since BGS sits
  near the high-SM end, it leaves less NP room and produces stronger M_KK
  lower bounds than a CKMfitter/UTfit-style `1.9e-3` to `2.1e-3` SM prediction.
  That is conservative only for avoiding under-constraining NP; it is not
  conservative for quoting the weakest robust lower bound.

- Budget uncertainty: yes, this should be propagated or at least bracketed in
  the paper-facing M_KK bound. The audit does not currently acknowledge it. I
  also do not reproduce the prompt's `~13%` uncertainty from the PDG 2024
  experimental error plus the BGS published total error; the full BGS uncertainty
  is larger than the central gap. If a reduced SM-error convention is intended,
  it needs an explicit source and rationale.

### Final

Ready for Opus sign-off.

===PHASE_2_H5_REVIEW_END===
