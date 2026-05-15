# Phase 2 Hole #4 Implementation Log: xfail repair

## Diagnosis

The pre-fix `--runxfail` fixture is `tests/baselines/pre-fix-xfail-output.txt`. It captured `12 passed, 2 failed` for `tests/test_quark_fit.py` when the two xfails were forced to run.

The first failure was in `test_fit_quark_sector_is_invariant_under_quotient_directions`, at the final canonical-vector equality assertion:

```python
np.testing.assert_allclose(
    encode_quark_fit_canonical_vector(base_solution.seed),
    encode_quark_fit_canonical_vector(shifted_solution.seed),
    rtol=0.0,
    atol=5e-6,
)
```

The two fitted canonical vectors differed in 6 of 14 entries, with max absolute difference `0.015689995721`. The earlier residual comparisons did not fail in this baseline run: max up/down mass-residual difference was `8.6259e-10`, and max CKM-residual difference was `7.6440e-09`, both below the test tolerances.

The second failure was in `test_fit_orientation_false_remains_deterministic_and_restricted`, at the first mass equality assertion:

```python
np.testing.assert_allclose(
    base_solution.result.masses_up,
    shifted_solution.result.masses_up,
    rtol=0.0,
    atol=2e-4,
)
```

The orientation-disabled base and shifted fits landed in distinct minima. The up-sector masses were `[1.200731e-03, 6.044365e-01, 1.557261e+02]` vs `[1.071363e-03, 5.907030e-01, 2.490268e+02]`, with top-mass difference `93.30065292` GeV. Their scores were `0.5449157281` and `0.9413364162`.

The PDG target machinery is now `default_quark_targets()` -> `pdg_up_down_arrays_at_scale(DEFAULT_QUARK_FIT_SCALE_GEV)` -> `run_msbar_mass(...)`, giving PDG-2024 MS-bar masses at `mu = m_t(m_t) = 163.5 GeV`:

```text
up   = [0.001189543183065, 0.6022859509993, 162.4227471603]
down = [0.002588357852040, 0.05149179982249, 2.731983679692]
```

The old `default_spurion_seed()` was still the pre-PDG seed. Under these PDG targets it started at score `6.305944155`, so tiny quotient-chart differences between the base seed and the shifted seed were enough for least-squares to settle into nearby but different reported minima. The invariant itself is still physically correct: common overall scale, right rotations, and `2pi` angle shifts are quotient redundancies, so quotient-equivalent seeds should give the same fitted observables and reported canonical seed when the benchmark seed is near a unique PDG-target minimum. This is not an obsolete pre-PDG test; the stored seed was obsolete.

## Fix Path

Chose path A: re-derive the seed.

I ran a deterministic least-squares fit under `default_quark_targets()` from the legacy benchmark seed. The raw optimum terminated by `gtol` after 53 evaluations with score `1.2093216839e-10`. I stored that PDG-target optimum in `default_spurion_seed()` while preserving `overall_scale = 2.8` by dividing the fitted physical spectra by 2.8.

With the repaired seed, both quotient-equivalent starts are already at the same PDG-target minimum. In both `fit_orientation=True` and `fit_orientation=False`, the base and shifted runs terminate after one evaluation, with max mass-residual difference `1.11e-14`, max CKM-residual difference `1.75e-13`, and max reported canonical-vector difference `0.0`. The two xfail decorators were removed so the original assertions now run normally.

## Commits

- `16e6b36` `test(fixtures): capture pre-fix xfail behavior for spurion-seed tests`
- `2871fce` `fix(benchmarks): re-derive spurion seed under PDG targets`
- `559b851` `chore(pytest): enable strict xfail policy`
- `fd83d96` `test(fixtures): capture post-fix passing test output`

## Pytest Counts

- Pre-fix forced xfail run: `12 passed, 2 failed` from `pytest -ra --runxfail tests/test_quark_fit.py`.
- Post-fix acceptance: `14 passed` from `pytest -q tests/test_quark_fit.py`; zero xfails and zero failures.
- Post-fix captured run: `14 passed` from `pytest -ra tests/test_quark_fit.py`.
- Targeted quark regression suite: `77 passed` from the requested nine-file pytest command.

## Strict Xfail Policy

Before this task, `pyproject.toml` had no `[tool.pytest.ini_options]` block and no global `xfail_strict` setting. After this task, `pyproject.toml` sets:

```toml
[tool.pytest.ini_options]
xfail_strict = true
```

## Deviations

No plan deviations. The methodology note, scan outputs, and figures were not touched.

Hole #4 ready for peer review.
