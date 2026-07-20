# Phase 2 Hole #6 Peer Review

Date: 2026-05-15
Reviewer: Phase 2 hole #6 peer reviewer

### Verdict

REJECT-WITH-REVISIONS

### Item-by-item (1-13)

1. PASS - Branch state. `paper/quark-scan-2026q2`, `origin/paper/quark-scan-2026q2`, `audit/wilson-rg`, and `origin/audit/wilson-rg` all resolve to `561d848479ab058df9d384fe611e664721325ceb`; both ancestry checks exit 0.

2. FAIL - Diff scope. `git show --name-only 7f71908` also adds `docs/audits/wilson_rg_inventory.md` and `docs/audits/wilson_rg_reference_values.md`, beyond the stated whitelist; no bag, CFW, or scan-output files drifted.

3. PASS - Standalone audit script. `python scripts/audit_wilson_rg.py` prints `C4_LR: ... 3.53816397486` and `C5_LR: ... -0.894757448992, 0.853891627884`, with max relative discrepancy `1.300e-16`.

4. FAIL - BMU-to-scalar map. BMU defines `Q1_LR` as vector LR and `Q2_LR` as scalar LR, and maps `B1^LR -> B5`, `B2^LR -> B4`; it does not establish the positive-sign statement `Q1_LR^BMU = 2 O5_LR`. Standard scalar-basis Fierz evidence instead points to a minus sign for conventional `O5`.

5. FAIL - Scalar LR ADM. BMU gives operator ADM `[[2,12],[0,-16]]`; Wilson evolution uses its transpose, so the diagonal/eigenvalues are right, but the scalar off-diagonal is `+6` only if the failed positive `Q1=+2O5` map is adopted. With the conventional `Q1=-2O5` map it flips to `-6`.

6. PASS - Threshold handling. Code uses `_MT_DEFAULT = 163.5`, `_MB_DEFAULT = 4.18`, `_MC_DEFAULT = 1.27`, matching `qcd.constants.M_TOP_MS`, `M_BOTTOM`, `M_CHARM`; this is an MS-bar top threshold, not the top pole mass near 173 GeV.

7. PASS - C5-to-C4 direction. The code is triangular in `[C4,C5]`: high-scale `C5` feeds low-scale `C4`, while high-scale `C4` does not feed `C5`; the sign of the feed remains tied to items 4-5.

8. PASS - Independent numerical sanity check. Using the audit-script alpha values, `prod[(alpha_high/alpha_low)^(-16/(2 beta0))]` over `n_f=6,5,4` gives `1.41219*2.00505*1.24957 = 3.53816`, matching the code. The prompt's rough physical `alpha_s(2 GeV)=0.298` gives about `3.71`, a scheme/input difference rather than a code mismatch.

9. PASS - Gate-tightening factor. `0.535 -> 2.3677` is `4.4259x`; multiplied by hole #5's `6.2388x` gives `27.61x`. This is plausible in size for LR running plus the tighter epsilon_K budget, but the exact value is conditional on the LR sign convention.

10. PASS - Tests. `/n/home09/obarrera/.conda/envs/ising_bootstrap/bin/python -m pytest -q tests/test_wilson_rg_audit.py tests/test_qcd_running.py tests/test_quark_deltaf2.py` reports `35 passed in 4.91s`; collected test names directly cover unit vectors, thresholds, LR mixing, composition, and benchmark output.

11. PASS - Methodology appendix coverage. The new appendix states the basis, BMU-to-scalar map, LO ADM, thresholds, endpoints, LO/NLO status, unit factors, and endpoint caveat; the map/ADM content is blocked by items 4-5.

12. PASS - NLO theory uncertainty. BMU's own LO-vs-NLO tables show LR shifts at the 10-25% level between `m_t` and 2 GeV/1 GeV for representative inputs, so an order `10-30%` missing-NLO systematic is plausible.

13. PASS - LR endpoint caveat. The code does not correct the mismatch: all systems still run to `mu_had = 2 GeV`, while kaon BSM bags are quoted at 3 GeV, `B_{d,s}` at `m_b`, and D inputs around 3 GeV. This affects non-K bounds and should be tracked as a residual systematic or separate fix.

### Findings

1. BLOCKER - The BMU-to-scalar LR map is not sign-safe. BMU arXiv:hep-ph/0102316 defines `Q1_LR` and `Q2_LR` in Eq. (2.1), gives `B1^LR = B5`, `B2^LR = B4` in Eq. (7.7), and lists the one-loop LR ADM in Appendix A. Those facts support the B-index correspondence, but not the positive operator identity `Q1_LR^BMU = 2 O5_LR`. The standard scalar-basis Fierz relation for conventional `O5` carries the opposite sign, which would flip the scalar-basis coefficient ADM off-diagonal from `+6` to `-6` and reverse the sign of the induced `C5_high -> C4_low` contribution. Recommended fix: explicitly define `O4_LR` and `O5_LR` with color, chiral projectors, and sign; then redo the BMU conjugation and update code, tests, methodology, and reference values. If the current `+2` map is intentionally using a nonstandard `O5`, state that and make the matrix-element signs consistent.

2. WARNING - The global `2 GeV` endpoint is still a real physics mismatch. It is correct for `B_K`, but not for FLAG kaon BSM `B4/B5`, `B_{d,s}` bags at `m_b`, or D-system inputs near 3 GeV. Recommended fix: defer to a separate endpoint/bag-running hole, but do it before scan-output freeze because it changes all quoted bounds.

3. WARNING - The committed diff exceeds the requested file whitelist. The extra files are Wilson-audit documentation, not physics-output drift, but the review checklist explicitly said the diff should touch only the listed files. Recommended fix: either amend the whitelist in the orchestrator plan or split audit-reference docs into the sealing commit.

4. NIT - The top threshold language is easy to misread. The implementation uses `m_t(m_t) = 163.5 GeV`, consistent with the shared QCD module, while the prompt mentions a pole-like `m_t approx 173 GeV`. Recommended fix: call it `m_t^MSbar(m_t)` in the methodology note and audit logs.

### Physics-decision flags

- BMU mapping `Q1_LR^BMU = 2 O5_LR`, `Q2_LR^BMU = O4_LR`: not approved as written. `Q2 -> O4` and `B1^LR/B2^LR -> B5/B4` are supported; the positive sign in `Q1 -> +2 O5` is load-bearing and not established. Treat this as a blocker until the scalar operator convention is explicitly fixed.

- LO-only running: acceptable for an interim scan if the paper propagates a `10-30%` Wilson-RG systematic and states that the calculation is LO/NDR-aligned but not full BMU NLO. I do not require NLO before this hole signs off, but final paper claims should include the systematic band or a later NLO upgrade.

- B-system / D-system endpoint mismatch: defer out of hole #6 into a separate endpoint-alignment hole, because it changes all scan outputs. Do not freeze final B/D bounds until that hole either fixes the endpoints or documents a sensitivity bracket.

- Cumulative `~27x` epsilon_K ratio shift: arithmetically honest under the committed convention. The hole #5 budget shrink multiplies the hole #6 Wilson-RG change, so there is no generic cancellation term. However, the disputed `O5` sign changes the default benchmark ratio materially; a quick sign-flip check changes `2.37` to about `1.93`, so the exact cumulative factor is not final until the BMU/scalar map is resolved.

### Final

Send back to implementer because the BMU-to-scalar LR sign convention is unresolved and it controls the scalar ADM off-diagonal and induced `C5 -> C4` contribution.

===PHASE_2_H6_REVIEW_END===
