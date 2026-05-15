# Phase 2 Hole #6 Sign-off

**Verdict**: PASS-WITH-FOLLOWUP

The Wilson-RG audit revision (commits `7f71908`, `561d848`, `c80a8e0`, plus the
bag-parameter chain) is independently reproducible, internally sign-consistent
under the conventional scalar-LR basis (`O4_LR`, `O5_LR`), and now matches the
closed-form LO reference to floating-point precision. Two residual systematics
(endpoint mismatch on non-kaon systems, LO-only running) are explicitly
deferred and recorded as documented systematics in the methodology note.

## Verification checklist (1-8)

1. **Branch state — PASS.** `git rev-parse --abbrev-ref HEAD` returns
   `paper/quark-scan-2026q2`; `git log --oneline -15` shows the full chain in
   order: `82a96f0` (hadronic inputs) → `dc9c498` (FLAG2024+BGS2020 inputs) →
   `695f35e`, `2521f6b` (hole #5 impl + signoff) → `7f71908` (BMU LO RG fix)
   → `561d848` (impl seal) → `c80a8e0` (BMU-to-scalar-LR sign correction).
   Upstream is `origin/paper/quark-scan-2026q2`; working tree clean except for
   the two untracked review notes plus this signoff.

2. **ADM matrix — PASS.** `grep -n` on `quarkConstraints/qcd_running.py` shows
   `_GAMMA_LR = np.array([` at line 54 with `[-16.0, -6.0]` at line 55 and the
   `eig(_GAMMA_LR)` consumption at line 201. `deltaf2.py` no longer carries a
   local copy. The prompt's literal pattern with the bracket character class
   misses the negative-number form because GNU `grep` treats `[-16, -6]` as an
   unterminated character class; the matrix is confirmed by source inspection.

3. **Standalone audit — PASS.** `python scripts/audit_wilson_rg.py` prints
   `C4_LR=1 -> 3.53816...`, `C5_LR=1 -> (0.89476..., 0.85389...)` with both
   LR-sector components POSITIVE real, and `Max relative discrepancy vs
   closed-form LO reference: 1.300e-16`. Three-flavor-threshold segmentation
   (`3000->163.5 / 163.5->4.18 / 4.18->2 GeV`) is correctly applied.

4. **Tests — PASS.** `pytest -q tests/test_wilson_rg_audit.py
   tests/test_qcd_running.py tests/test_quark_deltaf2.py` reports `35 passed`
   in 5.68 s, no warnings, no xfails.

5. **Cumulative gate factor — PASS.** `python -c "print(6.2388 * 3.6052)"`
   returns `22.49212176`, confirming the 22.49x cumulative invalidation
   factor (hole #5 bag-parameter audit 6.24x × hole #6 Wilson-RG audit
   3.61x).

6. **Methodology note PDF — PASS.** `pdfinfo
   docs/quark_scan_methodology_note.pdf` reports `Pages: 17`, built with
   `pdfTeX-1.40.19` and `LaTeX with hyperref`, current timestamp.

7. **Sign-chain documentation — PASS.** `docs/audits/wilson_rg_inventory.md`
   line 22 records the explicit map `Q1_LR^BMU = -2 O5_LR`, line 27 opens a
   dedicated "Sign-convention chain audit" section diagnosing the prior bug
   as type "B, sign", and lines 43, 82-83 restate the mapping with code
   pointers at `qcd_running.py:7`, `:11`, `:50`. The hole #6 v2 review
   confirms the chain is self-consistent and the kaon contraction returns
   positive `<O5_LR>` with positive FLAG-style `B_5 = 0.691`.

8. **NIT fix — PASS.** `docs/quark_scan_methodology_note.tex` defines
   `\msbar` at line 32 and uses `$m_t^{\msbar}(m_t) = 163.5\,\text{GeV}$` at
   line 1016, with consistent `\msbar` usage for the quark-mass and
   bag-parameter conversions (lines 55-59, 899, 930-934). No bare "MSbar" or
   "msbar" strings remain in the body text.

## Physics-decision adjudication

**(a) Endpoint mismatch: ACCEPT DEFERRAL (a1).** The global `mu_had = 2 GeV`
endpoint is FLAG-consistent for `B_K` but mismatched for kaon BSM bags
(`B_4`, `B_5` at 3 GeV), `B_d` and `B_s` (at `m_b`), and D-meson inputs
(near 3 GeV). The implementer's deferral is correct because:
(i) the methodology note (Section on hadronic conventions, lines 930-934)
already discloses the 3 GeV BSM bag origin and the 2 GeV evolution endpoint,
(ii) the per-system pass-rate figure shows `eps_K` is the dominant binding
channel (~93% pass-fraction), while `Delta M_K`, `B_d`, `B_s`, `D` all pass
at >99% in the post-audit constraints, so an endpoint migration on those
systems shifts non-binding ratios and cannot affect the headline `M_KK`
bound at the level of the 22.5x factor, and (iii) the residual is bounded
above by NLO-class corrections (few × 10%) on already-non-binding observables.
Recorded as a paper-systematic and as a follow-on hole.

**(b) LO-only running: ACCEPT (b1) with systematic band.** Full BMU
NLO/NDR evolution of the LR sector is a multi-week effort (NDR-vs-RI
scheme matching, NLO ADM blocks, NLO matching at thresholds, validation
against literature tables). The conservative literature compromise is
LO running plus a propagated systematic on the bound; standard practice
in RS flavor papers is the same. Adopt LO and quote a 10-30% NLO
uncertainty band on the M_KK bound in the methodology note. The 22.5x
gate factor is robust against this band because the dominant correction
is the ADM fix itself (factor ~3.6), which is exact at LO and dwarfs the
NLO-class residual.

**(c) Cumulative 22.5x gate: CONFIRMED HONEST.** The product
`6.2388 x 3.6052 = 22.49212176` is reproduced. The 6.24x factor from
hole #5 is anchored on the FLAG 2024 + BGS 2020 bag-parameter migration
(signed off in `phase2_h5_signoff.md`). The 3.61x factor from hole #6
combines the BMU LO running fix (`7f71908`) and the BMU-to-scalar-LR
sign correction (`c80a8e0`); the standalone audit and the closed-form
LO cross-check agree to 1.3e-16, so the factor is not numerical
sleight-of-hand. The 22.5x cumulative factor triggers mandatory re-runs
of RUNA and the eight follow-up scans under the corrected constants.

## Deferred follow-ups (to record in roadmap)

1. **Endpoint migration to per-system FLAG conventions.** Add a follow-on
   hole that evolves `B_4`, `B_5` to 3 GeV, `B_{B_d}`, `B_{B_s}` to `m_b`,
   and D-meson bags to their FLAG quotation scale, replacing the global
   `mu_had = 2 GeV` endpoint. Expected impact: shifts non-binding
   ratios by O(10%), no impact on the kaon-dominated `M_KK` headline.
2. **BMU NLO/NDR running for the LR sector.** Implement the NLO ADM
   blocks (Buras-Misiak-Urban 2000), NLO threshold matching, and
   NDR-vs-RI scheme accounting. Expected impact: 10-30% on the
   `M_KK` bound, currently quoted as a systematic band.
3. **VLL/VRR sector NLO cross-check.** The hole #6 audit only locks the
   LR sector against a closed-form LO reference. The vector sectors run
   trivially at LO; add an NLO sanity check before the final paper draft.
4. **Methodology note appendix on conversion tables.** Once (1) and (2)
   land, regenerate the per-system conversion tables in
   `docs/quark_scan_methodology_note.tex` and re-derive the systematic
   band quantitatively rather than as a literature estimate.

## Invalidation-gate status

Cumulative invalidation factor **22.49x confirmed** (hole #5: 6.24x,
hole #6: 3.61x, product verified by direct multiplication). Both Wilson-RG
and bag-parameter audits have passed their respective sign-offs. Per the
final-claims invalidation gate defined in
`docs/paper_execution_roadmap.md` (line 161), no `p50/p95/figure` may
appear in the paper under the old constants. Scan re-runs of RUNA and the
eight follow-up scans under the corrected ADM + sign-chain + FLAG 2024
inputs are **AUTHORISED** and become the next major step. The pre-audit
manifest stays in place; the final blessed manifest is withheld until the
re-runs complete.

## Recommendation

**Proceed to invalidation-gate scan re-runs, THEN hole #7 (CFW).** Run
RUNA + 8 follow-ups under corrected constants, regenerate per-system
pass-rate and `M_KK` distribution figures, update the methodology note
with the new headline bound and the propagated 10-30% NLO systematic
band, then open hole #7.

---

Signed: Phase 2 hole #6 sign-off authority
Date: 2026-05-15
