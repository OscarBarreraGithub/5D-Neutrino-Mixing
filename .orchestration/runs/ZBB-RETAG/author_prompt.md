# Z to bb RE-TAG — make the minimal non-custodial constraint VETOING (codex author, gpt-5.x xhigh)
Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. AUTHOR task: first a SHORT plan, then implement. A codex reviewer AND an Opus reviewer will independently review; BOTH must APPROVE before commit.

## ESTABLISHED FACTS (triple-verified by proposal+review+Opus sign-off; see `.orchestration/runs/ZBB/*`)
- The MINIMAL non-custodial Z to bb prediction is ALREADY fully computed: the dominant gauge-KK profile shift to the b coupling is in z_delta_g_L/R_d[2,2] via the exact KK-tower overlap a(c) (rs_ew_spectrum.py ~789-917; rs_ew_couplings.py ~505-551), PLUS the small m_b^2 Casagrande admixture (rs_ew_couplings.py ~704-744, on in the scan). DO NOT add any new coupling/term and DO NOT change the computed R_b/A_b pulls. Adding the analytic m_t^2 term would DOUBLE-COUNT.
- T010 and T011 already compute real R_b/A_b pulls and set passes = (max abs_pull <= 1.0), i.e. they are veto-capable (T010.py ~512-555).
- They are held NON-vetoing ONLY because `_rs_zbb_matching_diagnostics` (T010.py ~452-469) UNCONDITIONALLY stamps `needs_human_physics` = custodial-variant status, which the scan's `tag_result` maps to tag="partial" -> advisory (run_full_catalog_scan.py ~951-954, _classify_results ~875-883).

## THE CHANGE (re-tag only; no physics-number change)
Gate the custodial `needs_human_physics` flag so it NO LONGER suppresses the minimal constraint. The MINIMAL piece (gauge tower + m_b^2 admixture) must be reported as a complete, VETOING constraint (so the scan tags it rigorous/minimal-complete and routes its HARD failures into the real veto accounting). The custodial / top-partner / BKT refinement must REMAIN explicitly deferred and recorded (a separate, clearly-labeled note/field that does NOT force the whole result to partial). In short: minimal Z to bb becomes live and vetoing; custodial stays a documented deferred refinement.
- Touch T010 and T011 (the diagnostics gate that emits needs_human) and, if needed, the scan's tag mapping so a constraint that is "minimal-complete with a deferred custodial refinement" is treated as rigorous/vetoing rather than partial/advisory. Keep the deferral HONEST and visible (e.g. a `custodial_variant_deferred=True` field), just not veto-suppressing.
- JUDGEMENT CALLS already decided (do not change): A_FB^{0,b} stays CONTEXT-ONLY (not a veto scalar) exactly as now; the Casagrande chiral cross-assignment (singlet->g_L, doublet->g_R) stays as currently documented.
- Do NOT change any other constraint's tagging. Do NOT alter the computed pulls, SM anchors, or coupling math.

## GATE TESTS
- T010/T011 now evaluate as vetoing-capable rigorous/minimal-complete (NOT partial) and a failing minimal point lands in the real HARD veto set; a passing point survives.
- The custodial refinement is still recorded (deferred flag present) but no longer forces non-vetoing.
- SM / universal-c limit: T010/T011 recover SM (no spurious veto).
- Determinism; finite; no change to the numeric R_b/A_b pulls vs before (assert the predicted/ratio values are byte-identical to pre-change for a fixed point — only the tag/veto-routing changes).
- `python -m pytest tests/ -q` stays green (update any test that asserted the old partial/advisory tag, enumerating the change).

## OUTPUT (<=16 lines)
Short plan; exact files/lines changed; how the minimal-vs-custodial split is now represented; proof the numeric pulls are unchanged (same predicted/ratio) and only the tag/veto-routing changed; the gate-test results; pytest counts. End with: ZBB-RETAG-AUTHOR-DONE.
