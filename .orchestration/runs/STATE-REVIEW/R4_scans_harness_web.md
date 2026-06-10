Wrote the R4 section to `.orchestration/runs/STATE-REVIEW/R4.md` and verified it has exactly the four required buckets and ends with `SECTION-READY`.

Verification run:
- `LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH" pytest -q tests/test_full_catalog_scan_harness.py` -> 18 passed
- `LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH" pytest -q tests/test_wq_quarkonly_comparison.py` -> 2 passed

Note: `.orchestration/runs/STATE-REVIEW/` is currently untracked, including both `FORMAT.md` and the new `R4.md`.