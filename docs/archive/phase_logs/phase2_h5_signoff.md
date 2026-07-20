# Phase 2 Hole #5 Sign-off

Date: 2026-05-15
Authority: Claude Opus sign-off agent for Phase 2 hole #5
(hadronic bag-parameter audit), acting per the orchestrator delegation
in `docs/paper_execution_decisions.md`.
Branch: `audit/bag-inputs`, FF-merged into `paper/quark-scan-2026q2`.
Commits in chain: dc9c498, 82a96f0, 695f35e.

**Verdict**: PASS-WITH-FOLLOWUP

The audit's mechanical work (provenance inventory, code update, test,
methodology-note update, PDF rebuild) is correct and source-consistent.
The two outstanding peer-review WARNINGs are adjudicated below as
follow-up edits to the methodology note, not as blockers on hole #5
itself. The invalidation gate remains TRIPPED; scan re-runs are deferred
until hole #6 also signs off.

## Verification checklist

1. **`git log --oneline 1fec2c7..HEAD`** - PASS. Output:
   ```
   695f35e docs(paper): seal Phase 2 hole #5 (bag-param audit) impl log
   82a96f0 docs(paper): document hadronic input provenance
   dc9c498 physics(deltaf2): update kaon inputs to FLAG 2024 and BGS 2020
   ```
   Three commits, all expected for hole #5: code, docs, log.

2. **`pytest -q tests/test_quark_deltaf2.py`** - PASS. Output:
   `5 passed in 6.71s`. The new test
   `test_audited_deltaf2_hadronic_constants_match_selected_sources`
   is among the five.

3. **Methodology-note PDF rebuild via system pdflatex** - PASS.
   `/usr/bin/pdflatex -interaction=nonstopmode quark_scan_methodology_note.tex`
   from `docs/` exits 0; final line reports
   `Output written on quark_scan_methodology_note.pdf (16 pages,
   570549 bytes)`. Page count matches the impl log's claim of 16
   pages. No unresolved references. Aux files cleaned post-build.

4. **`docs/audits/bag_param_inventory.md` exists with >= 15 rows** -
   PASS. File exists (86 lines, 9103 bytes). Table contains 43 rows of
   `^|` lines, of which one is the header separator and one is the
   header itself; the data section contains 36 physical input rows
   plus a separate red-row amplitude table of three rows. Both
   table sections individually exceed the 15-row threshold; the
   inventory body has 36 >>= 15.

5. **`INVALIDATION_GATE_TRIPPED` honestly noted in
   `phase2_h5_impl.md`** - PASS. Line 62 reads
   `INVALIDATION_GATE_TRIPPED: B_1_K, B_4_K, B_5_K, and EPSILON_K_SM
   change epsilon_K amplitudes or ratios by more than 10% at
   M_KK = 3 TeV; affected final-claim scan is RUNA, because RUNA's
   final p50/p95/figure claims depend on pre-audit epsilon_k_ratio
   values.`. The impl log also lists the specific per-input
   amplitude shifts (23.25%, 15.77%, 21.23%, and 6.24x for the
   budget itself).

All five verification items pass. The audit deliverables are
internally consistent and source-verifiable.

## Physics-decision adjudication

### BGS 2020 epsilon_K^SM choice (WARNING 1)

**Decision: Option (a) - accept BGS 2020 (2.161e-3) as the central
value, and quote the M_KK epsilon_K-driven bound as a band whose
edges are set by the combined BGS theory uncertainty (0.18e-3), PDG
experimental uncertainty (0.011e-3), and the literature spread
across CKMfitter / UTfit / BGS (~0.15e-3).**

Justification: the orchestrator's audit philosophy in
`docs/paper_execution_decisions.md` ("prefer narrowing claims over
silently changing numbers") rules out option (b) (switching to a
softer UTfit central). BGS 2020 is not contradicted by literature
consensus; it sits at the high end of an internally-consistent band
of modern predictions. Replacing it with a softer value would be a
silent convention switch driven only by the desire for a looser
bound. Option (c) (band only, no central) is rejected as poor
presentation - flavor papers conventionally quote a central plus a
systematic band. (a) is the standard treatment: keep the most recent
modern perturbative input as the central, and disclose the
SM-prediction sensitivity as a band on the headline M_KK bound.
Full rationale and recommended prescription are recorded at
`docs/audits/epsilon_k_sm_decision.md`. (148 words.)

### Budget uncertainty propagation (WARNING 2)

**Decision: the M_KK lower bound must be quoted as a band, not a
single number, in the methodology note and in the paper text.**

Justification: the peer reviewer correctly notes that the BGS total
uncertainty (0.18e-3) is much larger than the central budget
(6.7e-5), so the budget is not a stable Gaussian denominator. The
correct treatment is therefore not to compute a one-sigma error on
M_KK from a Gaussian uncertainty, but to bracket the bound by
recomputing it at the budget edges. The mass-bound band roughly
covers M_KK_central * 0.4x (tight edge) to M_KK_central * 2.1x
(loose edge) under the
1 / sqrt(|Delta epsilon_K|) scaling, which is the headline number
the paper must disclose honestly. Concrete prescription is in
`docs/audits/epsilon_k_sm_decision.md`. (94 words.)

## Follow-up tasks for the orchestrator

These edits must complete before any final paper claim is frozen.
None of them are blockers on hole #5 itself, but they are blockers on
the final-claims invalidation gate clearing.

1. Edit `docs/quark_scan_methodology_note.tex` to quote the
   epsilon_K-driven M_KK lower bound as a band (central + asymmetric
   error bars) and cite `docs/audits/epsilon_k_sm_decision.md` for
   the sensitivity prescription. This edit is **safe to do before
   hole #6 signs off**, since it is purely a presentation change.

2. After hole #6 (Wilson-RG audit) signs off, re-run RUNA at three
   epsilon_K-budget edges: central (6.7e-5), low (~1e-5), high
   (~3e-4). Add a CLI flag to `scripts/run_rs_anarchy.py` such as
   `--epsilon-k-budget {central,low,high}` or a numeric override.
   Record each output file in `docs/artifact_manifest.md` with code
   SHA, command, sha256, and figure hashes.

3. Replace the single-number M_KK^min quotation in the methodology
   note and in any draft paper text with the band form: e.g.
   `M_KK^min (50%, g_s* = 3) = X^{+a}_{-b} TeV (epsilon_K budget,
   1 sigma)`.

4. Add a methodology-note paragraph (1-2 paragraphs of new text)
   describing the band construction. Cite BGS 2020 for the SM
   prediction, PDG 2024 for the experimental value, and FLAG 2024
   for the hadronic inputs. Cite the decision file.

5. Consider adding a sensitivity-band figure (M_KK^min vs
   epsilon_K-budget edge) to the figure set, especially if the
   1-sigma band straddles physically-distinct M_KK regimes (e.g.
   above vs below an EWPO floor).

## Invalidation gate status

The invalidation gate remains **TRIPPED**, exactly as
`docs/phase_logs/phase2_h5_impl.md` line 62 honestly records. The
following claims in the methodology note and any draft paper text
are gated:

- RUNA p50/p95 quoted M_KK lower bounds driven by epsilon_K.
- Any figure that uses the pre-audit `epsilon_k_ratio` denominator.
- The legacy single-number epsilon_K NP-budget figure of 4.18e-4.

Scan re-runs are **deferred** until hole #6 (Wilson-RG audit) also
signs off. Re-running RUNA twice (after #5 then again after #6)
would be wasteful, since hole #6 may change Wilson normalizations
and thresholds, forcing a third re-run. Per
`docs/paper_execution_decisions.md` stopping conditions, the
orchestrator does not stop for an audit returning "needs re-run";
the re-run feeds the final manifest, and the final manifest is
unsigned until both audits pass.

The gate clears when, in order: (1) hole #6 signs off; (2) RUNA is
re-run at three budget edges; (3) the methodology note quotes
M_KK^min as a band; (4) the final manifest records all re-runs.

## Recommendation

**Proceed to hole #6 (Wilson-RG audit).**

Hole #5 deliverables are correct, source-verified, and the two
open WARNINGs have been adjudicated into a concrete band-quote
prescription documented at
`docs/audits/epsilon_k_sm_decision.md`. The invalidation gate
remains tripped, but per the orchestrator's policy this is
expected: gate clearance is deferred until both hole #5 and hole #6
are signed off. No further work on hole #5 is required before
starting hole #6.

===PHASE_2_H5_SIGNOFF_END===
