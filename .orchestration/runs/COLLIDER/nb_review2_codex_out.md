Rendered: 11 figure cells, all 11 have embedded `image/png`; 0 error outputs, 0 tracebacks.
New run: `.py/.ipynb` have 0 refs to `wq_quarkonly_1M_full`; `SCAN_ROOT` points to `wq_quarkonly_1M_20128400`.
Cache/run: new-run cache is 1,000,000 rows; 500 JSONL summaries sum to 1,000,000; allowlist is 46.
Schema: T010/T011 are `rigorous`/`HARD`; lepton-free CR proxy columns present through CR013, with CR005/006/009/011/014 deferred.
Narrative: matches new data; no stale “rigorous reach ~1-3 TeV” or “Zbb advisory” claim found. `~1-3 TeV` appears only for subdominant ΔF=2.
Caveats present: non-custodial Zb_L, custodial relaxation deferred, and ~25-30 TeV pending literature cross-check.
Plots/code: Zbb is in veto-eligible rigorous panels, not advisory; CR* are proxy carvers; M_KK >= 5 TeV frame with 4 TeV axis floor; C3 veto is per-r, not pooled.
Spot-check: T010 veto fraction 15/20/30 TeV = 1.000 / 0.8918 / 0.000; rigorous/veto survival = 0 through 15, 0.108 at 20, 1.000 at 30.
Scope/tests: changed files are notebook, exported notebook `.py`, and regenerated figs only; no production code changes. `pytest tests/ -q`: 1716 passed, 1 skipped.
NB-REVIEW2: APPROVE