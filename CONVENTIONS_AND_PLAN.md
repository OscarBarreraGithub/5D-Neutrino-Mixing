# Conventions and Consistency Plan

Date: 2026-02-04

This document records the **canonical physics conventions**, the current **mismatches**, and the **fix plan** so the codebase stays consistent.

## 0) Decisions (Explicit and Final)

These are the explicit choices to keep the code **PDG-aligned** while staying consistent with the Perez–Randall setup already used in this repo.

- **Bulk mass parameter**: `c = M5 / k` (no `+ 1/2` shift).
  - Rationale: this is the convention implicitly used by the overlap formulas and the KK solver (`alpha = |c + 1/2|`).
- **Default curvature scale**: `k = M_Pl = 1.2209e19 GeV` (not reduced Planck).
  - Rationale: matches existing `warpConfig/baseParams.py` and the current benchmark results against Perez–Randall.
- **PMNS convention**: PDG parameterization, with Majorana phases on the right:
  - `U_PMNS = U(θ12, θ23, θ13, δ) * diag(1, e^{i α/2}, e^{i β/2})`
  - Use the PDG central values for angles and δ; keep δ values split by ordering as in `neutrinos/neutrinoValues.py`.
- **Mass/Yukawa formulas**: keep the factor-of-2 normalization exactly as used in the code:
  - `m_Ei = 2 v k f_L Y_Ei f_Ei`
  - `m_νi = (2 k^2 v^2 f_L^2 f_N^2)/((f_N^UV)^2 M_N) * Y_Ni^2`
  - `Ybar = 2k * Y`

## 1) Canonical Conventions (Approved Source of Truth)

These are the conventions used by the current implementation and the standard RS formulas used throughout this repo:

- **Bulk mass parameter**: `c = M5 / k` (dimensionless)
- **Bessel order**: `alpha = |c + 1/2|`
- **Warp factor**: `epsilon = Lambda_IR / k`
- **IR overlap**:
  - `f_IR^2(c) = (1/2 - c) / (1 - epsilon^(1-2c))`
- **UV overlap**:
  - `f_UV^2(c) = (1/2 - c) / (epsilon^(2c-1) - 1)`
- **Limit at c -> 1/2**:
  - `f_IR^2 = f_UV^2 = -1 / (2 ln epsilon)`
- **Charged lepton masses (IR Higgs)**:
  - `m_Ei = 2 v k f_L Y_Ei f_Ei`
- **Neutrino seesaw (universal limit)**:
  - `m_nu_i = (2 k^2 v^2 f_L^2 f_N^2) / ((f_N^UV)^2 M_N) * Y_Ni^2`
- **Rescaled Yukawa**:
  - `Ybar = 2k * Y`
- **PMNS convention**:
  - In the charged-lepton mass basis: `Y_N = V_PMNS * diag(Y_Ni)`

## 2) Mismatch Audit (Resolved in Worktree)

Status: All items below have been corrected in this worktree. Kept here for traceability.

1. **Bulk-mass definition mismatch**
   - `derivations/conventions.tex` and `PROJECT_STATUS.md` define `c = M5/k + 1/2`,
     but the **code** and overlap formulas use `c = M5/k`.
   - This is inconsistent with the overlap formulas as written.

2. **c -> 1/2 limit missing factor 1/2**
   - `derivations/conventions.tex` states `f_IR^2 -> -1/ln epsilon`.
   - Correct limit is `-1/(2 ln epsilon)` (matches code).

3. **Y_N_bar example range mismatch**
   - `yukawa/neutrino.py` docstring says `Y_N_bar ~ O(0.01–0.1)`.
   - Running the benchmark gives `Y_N_bar ~ [0.20, 0.43, 1.02]`.

4. **Import path bug**
   - `neutrinos/massConstraints.py` uses `from neutrinoValues import ...`, which
     fails when imported as a package from repo root.

## 3) Fix Plan (Detailed)

1. **Standardize the definition of c**
   - Update docs to: `c = M5/k`.
   - Update `derivations/conventions.tex`, `PROJECT_STATUS.md`, and `CLAUDE.md`.

2. **Correct the c -> 1/2 limit**
   - Update `derivations/conventions.tex` to `-1/(2 ln epsilon)`.

3. **Align documentation examples with actual outputs**
   - Update `yukawa/neutrino.py` example range.
   - Re-check `PROJECT_STATUS.md` / `README.md` examples for consistency.

4. **Packaging cleanup**
   - Add `__init__.py` for `neutrinos/`, `warpConfig/`, etc.
   - Fix absolute/relative imports (notably `neutrinos/massConstraints.py`).
   - Remove `sys.path` hack in `yukawa/compute_yukawas.py` once packaging is stable.

5. **Verification script / test**
   - Add a small `scripts/benchmark_perez_randall.py` (or a `tests/` check)
     that validates the benchmark point:
     - `f_L ~ 0.016`
     - `f_N ~ 0.48`
     - `f_N^UV ~ 1.2e-4`
     - `Y_E_bar ~ [2.9, 4.4, 5.4]`
     - `Y_N_bar ~ [0.20, 0.43, 1.02]`

6. **Dependencies**
   - Add `requirements.txt` or `pyproject.toml` with `numpy`, `scipy`.

## 4) Verification Checklist

- Benchmark reproduces Perez–Randall point within ~1–2%.
- `c -> 1/2` overlap limit matches `-1/(2 ln epsilon)` numerically.
- `Y_N Y_N^dagger` uses the same PMNS convention as `neutrinoValues.get_pmns()`.
- All modules import cleanly from repo root.

## 5) Open Questions

None at this time. If any new ambiguity appears, add it here before changing conventions.
