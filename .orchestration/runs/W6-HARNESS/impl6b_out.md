Short plan executed: read W6 plan/W6a hook, implemented harness, added tests, ran fresh 10k smoke, ran full pytest, completed dual review.

Driver: `scripts/run_full_catalog_scan.py`; harness tests: `tests/test_full_catalog_scan_harness.py`.

Reviews: Codex APPROVE and Opus APPROVE, recorded in `.orchestration/runs/W6-HARNESS/rev6b_*_out.md`.

Smoke output: `.orchestration/runs/W6-HARNESS/smoke_w6b_fresh/` with 10,000 rows, 533 evaluated, 9,467 skipped.

Post-cache timing: 0.385822 s/draw; 7.23868 s/evaluated point.

1e8 extrapolation: 1.0717e4 core-hours per draw basis; 2.0107e5 core-hours per evaluated-point basis.

Constraints: 98 evaluated and 98 active per evaluated point; exception rate 0.0485437.

Top rigorous HARD vetoes: B022 533, K004 463, L001 181, B003/B004/B023/K001 138.

Top proxy HARD vetoes: CR009 533, CR006 391, B016/CR001/CR005/CR012/CR013/EW001 265, K010 256, B015 205.

Universal-c SM sanity: passed; rigorous=[], proxy=[].

Pytest: `python -m pytest tests/ -q` -> 1706 passed, 1 skipped.

W6B-DONE.