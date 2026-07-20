# Phase 2 Hole #4 Sign-off

**Verdict**: PASS

**Verification checklist**:
1. PASS - Current branch is `paper/quark-scan-2026q2`; `git status -sb` shows tracking `origin/paper/quark-scan-2026q2` with no ahead/behind markers (HEAD = `7e679ea`, matches `origin/paper/quark-scan-2026q2`).
2. PASS - `pytest -q tests/test_quark_fit.py` reports `14 passed in 5.17s` with no xfailed / no failed entries.
3. PASS - `grep "xfail_strict" pyproject.toml setup.cfg pytest.ini` returns `pyproject.toml:xfail_strict = true` (the `setup.cfg`/`pytest.ini` exit 2 is simply because those files do not exist; the setting is correctly placed in `pyproject.toml`).
4. PASS - `git log --oneline b0675cd..HEAD` shows exactly the five expected commits in the expected order: `16e6b36` (fixture pre), `2871fce` (seed fix), `559b851` (strict-xfail config), `fd83d96` (fixture post), `7e679ea` (report seal).
5. PASS - `git show 2871fce -- quarkConstraints/benchmarks.py` shows only the `default_spurion_seed()` numeric values being re-derived (singular values, two `RotationParameters` triples, two zeroed right rotations). `overall_scale=2.8` is preserved. Nothing surprising.
6. PASS - `git show --stat 2871fce` touches exactly two files: `quarkConstraints/benchmarks.py` (+30/-6) and `tests/test_quark_fit.py` (-20, the two xfail decorators). No physics constants outside the spurion seed are touched; no `deltaf2.py`, `mass_running.py`, `pdg_quark_masses.py`, or scan/figure files are affected.
7. PASS - `tests/baselines/pre-fix-xfail-output.txt` (6053 bytes) holds the `--runxfail` `12 passed, 2 failed` baseline with the quotient-equivalence and orientation-false failure traces; `tests/baselines/post-fix-test-output.txt` (439 bytes) holds the post-fix `14 passed` capture. Both are non-empty and contain meaningful before/after evidence.

**Outstanding nits**:
- Cosmetic: the impl log (`phase2_h4_impl.md`, Fix Path paragraph) states that, with the repaired seed, "both the base and shifted runs terminate after one evaluation." The reviewer's independent shell reproduction showed `nfev = 26/26` with `fit_orientation=True` and `nfev = 3/3` with `fit_orientation=False`. Residual / canonical-vector invariances all hold to the documented precision (max mass diff ~1e-14, max CKM diff ~1e-13, max canonical-vector diff 0.0), so the science claim is intact and the test passes are real - only the "one evaluation" sentence is slightly inaccurate. Flag for a one-line correction during the Phase 3 polish pass; not blocking.

**Reasoning**: All seven verification items pass with concrete evidence. The seed re-derivation is scoped strictly to `default_spurion_seed()`, the two now-stale xfail decorators are removed (consistent with the new strict policy), `xfail_strict = true` is in place to prevent regression, and the pre/post baseline fixtures provide reproducible evidence. The lone reviewer WARNING is a doc-prose nit about optimizer `nfev` wording, not a science issue.

**Recommendation**: Proceed to hole #5 (bag params). Note the cosmetic `nfev` wording correction for the Phase 3 polish backlog.
