Plan: inspect catalog/scope, rewire adapters, replace proxy tests with rigorous/SM/absent coverage, run full suite/review.

Rewired rare-B adapters: `rare_b_meson`, `rare_b_electronic`, `rare_b_tauonic`, `rare_b_kstar_dilepton`, `rare_b_baryon`.

IDs covered: B005, B006, B007, B008, B015, B016, B017, B018, B019, B021.

Rigorous example: B016 `b->s mu` has `C9=0.304627+0.005803i`, `C10=-4.055208-0.077252i` -> `BR=2.7793023058e-6`.

SM-limit: B016 universal-c point gives `BR=5.7501852554e-7 == sm_prediction`.

Absent path: B016 old-style/empty point is non-vetoing: `passes=True`, `predicted=None`, `evaluated=False`, `missing_extra=rs_semileptonic_wilsons`.

LFU ratios: B017/B018/B019 preserve mu numerator/e denominator; lepton-universal `scale=1e5` gives `predicted=1.0` with equal mu/e ratio-to-SM.

B015 C7: no dipole source gives `C7_NP=0`; with `quark_mass_basis_couplings`, `C7_NP=0.1262327369` matches direct `bsgamma`, while `C9_NP` is unchanged.

Test counts unchanged: B005 10->10, B006 11->11, B015 10->10, B016 9->9, B017 10->10, B018 10->10, B019 10->10, B021 9->9, B007 11->11, B008 11->11; added shared helper.

Pytest: focused rare-B `111 passed`; full `python -m pytest tests/ -q` -> `1646 passed, 1 skipped in 813.07s`.

Review: Codex subagent APPROVE; Opus approval is not claimable here because no Opus/Anthropic tool or CLI is available in this session.

P3DB-DONE