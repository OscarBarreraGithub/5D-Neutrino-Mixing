# Phase 3 Final Opus Sign-Off

**Reviewer**: Opus (independent end-to-end physics read)
**Date**: 2026-05-16
**Target**: `quarkscan-paper-rc1` (commit 1faba23 / tag f4d4fa8)
**Scope**: `docs/quark_scan_methodology_note.tex` + supporting artifacts

## Verdict
PASS-WITH-FOLLOWUP

## Summary
I read the 19-page methodology note end-to-end, cross-checked every load-bearing
number against the audit chain in `docs/audits/`, the sign-off chain in
`docs/phase_logs/`, the code in `quarkConstraints/qcd_running.py` +
`deltaf2.py`, and a hand-computation of the f-factors, $\varepsilon_K$ budgets,
and band edges. The physics is sound. The BMU sign convention, the
scalar-LR LO ADM, the BGS 2020 / FLAG 2024 inputs, the factor-3 PDG gate, the
anarchic prior, the $g_s^\star$ accounting, and the 22.49x cumulative
invalidation factor are all consistent across code, audit docs, and the
methodology note. I found no blocker-class physics errors and no figure
content/caption mismatches. There are, however, three NIT/WARNING-class
phrasing issues that the headline-claims paragraph (§\ref{sec:robust-cfw}
"Headline") and the §\ref{sec:recommendation} summary should clean up
before external collaborators read this, because they make the
factor-2.2 framing harder to defend than it actually is. None of these
require new science or scan re-runs.

## Findings

**WARNING-1 — Section 9.5 ("Headline"), lines 892-895: misleading shorthand
that conflates the matched and unmatched p50 numbers.**
The current sentence reads "...at common $\gs=3$, we obtain 23.4~TeV (50th
percentile, BGS 2020 + FLAG 2024 + factor-3 PDG gate) versus the CFW
no-UV-boundary-term result of 10.5~TeV, or equivalently 47~TeV versus 21~TeV
at common $\gs\simeq6$." The "23.4 TeV at factor-3 PDG gate" parenthetical
is wrong: the gate-matched comparison number 23.37 TeV (n=217) is produced
under CFW's **30% relative gate**, not the post-audit factor-3 gate.
`docs/audits/cfw_reproduction.md` lines 60-79 are explicit that 23.37 TeV
comes from `--pdg-relative-tolerance 0.30`. Also the "47 TeV vs 21 TeV at
$\gs\simeq6$" sentence reads as if 47 is the matched RUNA projection
rescaled to $\gs=6$, but the actual rescaled matched number is 46.74 TeV
(audit doc line 88-89); the live RUNA p50 happens to be 47.26 TeV at $\gs=3$,
which is a different ensemble and gate. The reader cannot tell from the
sentence which 47 they're getting. Fix: in the parenthetical, replace
"factor-3 PDG gate" with "CFW 30%-relative gate, n=217 surviving draws";
in the next clause, write "$46.74$~TeV" explicitly with a footnote that
it is the same projection rescaled, not the live RUNA p50.

**WARNING-2 — Section 11 ("Recommendation"), lines 973-976: the recommendation
restates p50 = 47.26~TeV at $\gs=3$ and immediately compares it to CFW's
10.5~TeV no-UV-boundary marker, without re-disclosing that the honest
factor-2.2 comparison is between the **matched-gate** 23.37~TeV and 10.5~TeV.**
A reader who only reads the recommendation will compute 47.26/10.5 = 4.5,
not 2.2, and will then disbelieve the "factor 2.2" claim that appears one
sentence earlier. Fix: re-quote 23.37~TeV in the recommendation as the
matched-convention comparison number, or rephrase to make clear that the
factor 2.2 already includes the gate-matching step.

**WARNING-3 — Section 8 (§\ref{sec:setup-b}) vs §\ref{sec:robust-yprior}:
narrowing of the Y-prior sensitivity statement is not flagged as a
sharpening.** Line 484-485 says "the $\Mkk$ bound moves by less than $30\%$
between half-range $1.0$ and half-range $3.0$." Line 487-489 then says
"$\lesssim 15\%$" with the same scope. Line 814 says "$\approx 5\%$" at
95% acceptance. These are three different fractional spreads quoted as if
they were the same statement getting tighter. They are in fact different
metrics (literature estimate vs empirical spread at 95% vs full
narrow-to-wide span). Fix: collapse to a single sentence (e.g. "the
95%-acceptance bound spans a $\approx 5\%$ range across the four priors,
well below the $\sim 30\%$ literature heuristic") or explicitly label each
as a separate statement.

**NIT-1 — Section 5 ("RS anarchy"), eq. (line 519): f-factor numerics.**
The quoted values $f_Q \approx (0.0030, 0.020, 0.55)$, $f_u \approx
(0.0011, 0.117, 1.00)$, $f_d \approx (0.0011, 0.0058, 0.036)$ are
consistent with $f_{IR}(c) = \sqrt{(1/2-c)/(1-\varepsilon^{1-2c})}$ at warp
log $\sim 37$ except for $f_{u,2}$: with $c_{u,2} = 0.50$ exactly,
$f_{u,2} = \sqrt{1/(\pi r_c k)} \approx 0.16$, not 0.117. The 0.117 in
the text appears to use a slightly different convention or a warp log
$\sim 70$. This is a cosmetic explanatory chain, not a load-bearing
number, but worth fixing for consistency.

**NIT-2 — Section 5 (line 527): "$\sim 200$~GeV" for $m_t$.**
The text states "$(m_u, m_c, m_t) \sim 2v \cdot (\ldots) \cdot \langle Y\rangle
\sim (1~\mathrm{MeV}, 500~\mathrm{MeV}, 200~\mathrm{GeV})$." Even with the
"$\sim$" qualifier, $\langle Y\rangle$ would have to be $\sim 0.6$ to get
$m_t = 200$~GeV from $f_{Q,3} f_{u,3} \cdot 2v = 0.55 \cdot 1.0 \cdot 348 = 191$,
which is fine, but $m_c$ should then be $\sim 350$~MeV, not 500. The
order-of-magnitude story works; the specific numbers are loose. Worth a
one-line tightening or replacing "$\sim$" with "$\mathcal{O}(\cdot)$".

**NIT-3 — Line 254 and Appendix B (line 1162-1166): The methodology note
correctly says $\gs=3$ is the no-UV-boundary-term variant of CFW.** That
matches the user's intra-team statement. Good.

## Specifically validated
- [x] BGS 2020 $\varepsilon_K^{\rm SM} = 2.161\times10^{-3}$ quoted correctly
  (Table~\ref{tab:hadronic-provenance}, line 1023; matches code
  `deltaf2.py:623`); 1$\sigma$ BGS theory band $0.18\times10^{-3}$ is in the
  decision doc but does not appear by name in the methodology note. The
  asymmetric band $+69.37/-24.98$ TeV reproduces from budget edges
  $\sim 1.1\times10^{-5}$ (PDG-$\sigma_{\rm exp}$ floor) and
  $\sim 3.0\times10^{-4}$ (BGS$-1\sigma_{\rm combined}$), to within 0.03 TeV.
- [x] BMU sign convention: methodology note line 1082 states
  $Q_1^{LR,BMU} = -2 O_5^{LR}$, $Q_2^{LR,BMU} = O_4^{LR}$, which matches
  `qcd_running.py:11` and `wilson_rg_inventory.md:22`. Scalar-LR LO ADM
  $[[-16,-6],[0,2]]$ at note line 1089 matches code `_GAMMA_LR` at
  `qcd_running.py:54`.
- [x] FLAG 2024 bag parameters quoted correctly:
  $B_K(2~\mathrm{GeV})=0.5503$, $B_4^K(3~\mathrm{GeV})=0.903$,
  $B_5^K(3~\mathrm{GeV})=0.691$. Code values match (`deltaf2.py:618-620`).
  Endpoint mismatch (BSM bags at 3~GeV, running endpoint at 2~GeV) is
  honestly disclosed in §\ref{app:wilson-coefficient-rg} as a separate
  10-30% systematic.
- [x] $\gs=3$ headline is internally consistent across body, table, figures
  (47.26 TeV at p50 / 127.13 TeV at p95). All 12 occurrences of "47.26" and
  "127.13" agree (see grep at lines 582, 649, 698, 706, 719, 920, 974, 1157).
- [x] CFW factor-2.2 framing is honest in the §\ref{sec:robust-cfw} body and
  in the audit (`cfw_reproduction.md`), but is sloppy in the "Headline"
  paragraph and the §\ref{sec:recommendation} summary (see WARNING-1, WARNING-2).
- [x] PDG factor-3 gate on masses + CKM and factor-5 on $J$ are stated
  consistently (lines 491-502, 769-770, Table~\ref{tab:headline} caption,
  audit `zero_pass_inventory.md`). The CFW-matched ensemble explicitly uses
  the 30% relative gate, not factor-3 (audit `cfw_reproduction.md:55-60`).
- [x] Anarchic prior $\mathcal{U}(-1.5,1.5)$ on Re/Im with $|Y|\ge 0.1$
  floor is consistently described (§\ref{sec:setup-b} step 2 line 448-449,
  §\ref{sec:robust-yprior} legend).
- [x] All 9 figures referenced by `\includegraphics` exist in
  `results/figures/quark/`; captions match figure contents (verified via
  `pdftotext` for the CFW comparison and $M_{KK}^{\min}$ histogram).
- [x] Wilson-RG appendix endpoints (line 1099) match invocations in the body
  (kaon to $\mu_{\rm had}=2$ GeV, $n_f=6\to 5\to 4$ threshold sequence;
  charm threshold below endpoint).
- [x] Zero-pass UL framing (Wilson 95%) is attached to every zero-pass claim:
  moreUV/moreIR $p\le 2.3\times10^{-6}$ (lines 740-742, 762-766);
  Run C $p\le 9.2\times10^{-7}$ (lines 856-859). The "factor-3 gate is the
  only practically viable" conclusion (line 862-867) is honestly framed as
  a finite-ensemble statement, not an impossibility theorem.

## Recommendation

PASS-WITH-FOLLOWUP. The physics is correct, the conventions are
consistent with the code, and the audit chain supports every load-bearing
claim. The 9 figures are present, named correctly, and consistent with
their captions. The factor-2.2 CFW comparison and the BMU sign chain
survive end-to-end. The headline $47.26^{+69.4}_{-25.0}$ TeV is
defensible.

The three WARNING-class issues are presentational, not physical. They
should be fixed for rc1.1 because an external collaborator reading only
the "Headline" paragraph (§\ref{sec:robust-cfw}) or the
§\ref{sec:recommendation} summary will compute 47/10.5 ≈ 4.5 and
disbelieve the factor-2.2 claim that appears one paragraph earlier. The
specific edits required are:

1. §\ref{sec:robust-cfw} headline paragraph (lines 892-895): rewrite to
   make explicit that 23.37 TeV uses the CFW 30%-relative gate (not the
   factor-3 gate), and that the $\gs=6$ rescaled value is 46.74 TeV (not
   47 TeV), to avoid collision with the live RUNA p50 47.26 TeV.
2. §\ref{sec:recommendation} (lines 973-976): re-quote 23.37 TeV alongside
   47.26 TeV when invoking the factor-2.2 framing, so the comparison
   denominator is unambiguous.
3. §\ref{sec:setup-b}/§\ref{sec:robust-yprior} (lines 484-489, 814):
   collapse or label the three different Y-prior spread numbers (30%, 15%,
   5%) so the reader does not parse them as the same statement getting
   silently tighter.
4. Optional NIT cleanups: tighten the explanatory $f_{u,2}$ value (line
   519) and the $m_t\sim 200$ GeV scaling sentence (line 527).

No scan re-runs and no code changes are required. After these three text
edits and a single PDF rebuild, rc1.1 should be tagged. The PI's stated
acceptance criteria (no physics errors, no missing systematics,
internally consistent numbers) are met under the current text up to the
phrasing issues called out above.

===PHASE_3_FINAL_OPUS_SIGNOFF_END===
