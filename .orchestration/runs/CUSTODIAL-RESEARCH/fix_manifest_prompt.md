# Custodial PR1 — small fix: EW001 allowlist-extras manifest (codex, gpt-5.x xhigh). SMALL focused task.
Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. One test fails:
`tests/test_full_catalog_scan_harness.py::test_quark_only_allowlist_matches_candidate_verification_and_drops_deferred_leptons`.
Cause: the custodial PR1 change made `flavor_catalog_constraints/primary/top_higgs_ew/EW001.py` read an extra `rs_ew_couplings` (via `_resolve_ew_model` reading `point.get_extra("rs_ew_couplings").metadata["ew_model"]`). The AST-based allowlist test compares each constraint's actual `get_extra(...)` reads to the manifest `QUARK_ONLY_ALLOWLIST_EXTRAS` in `scripts/run_full_catalog_scan.py`; EW001's manifest entry is now missing `rs_ew_couplings` (test: "Left contains one more item: 'rs_ew_couplings'").

FIX: update the `QUARK_ONLY_ALLOWLIST_EXTRAS["EW001"]` tuple in `scripts/run_full_catalog_scan.py` to include `rs_ew_couplings`, in the correct position to match the AST read order the test computes. CONFIRM `rs_ew_couplings` is NOT in `QUARK_ONLY_FORBIDDEN_EXTRAS` (it is a quark-side extra the quark-only point already provides, NOT a swept-lepton extra) so no false lepton-dependence is introduced. Do NOT change EW001 physics or any other manifest entry. If the `_resolve_ew_model` read order vs the existing reads matters for the tuple ordering, match exactly what `_candidate_get_extra_usage` extracts.

VERIFY: `python -m pytest tests/test_full_catalog_scan_harness.py -q` GREEN, and `python -m pytest tests/ -q` GREEN (env: source ~/.bashrc && conda activate ising_bootstrap && export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH").

OUTPUT (<=8 lines): the exact manifest line change; confirm rs_ew_couplings is quark-side/not-forbidden; pytest counts (harness file + full). End with: MANIFEST-FIX-DONE.
