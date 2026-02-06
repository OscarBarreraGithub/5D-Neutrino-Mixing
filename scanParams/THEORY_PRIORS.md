# scanParams Theory Priors and Scan Scope

Date: 2026-02-05
Branch: `codex/scanner-theory-priors`

This note defines scanner priors before large cluster runs. The objective is to keep model flexibility while preserving the logic used in Perez-Randall style RS lepton scans.

## 1) Baseline Conventions

- RS convention in this repo: `c = M5/k`.
- `k` is an explicit config input in GeV and is recorded in every output row.
- v1 production runs fix `k = 1.2209e19 GeV` (unreduced Planck), but plumbing remains fully explicit.
- `Lambda_IR` is scanned in `3000` to `7000` GeV.
- Neutrino ordering is fixed to `normal` in v1.

## 2) UV Scale Convention (Important)

To avoid ambiguity when changing `k`, UV parameters are represented dimensionlessly:

- scan `MN_over_k = M_N / k`,
- derive `M_N = MN_over_k * k` internally.

Scanner policy:

- all UV scales are stored as ratios to `k`,
- scanner internals do not assume or auto-convert reduced vs unreduced Planck conventions.

## 3) KK Scale Mapping Convention

`Lambda_IR` and `M_KK` are connected explicitly in metadata.

Use:

- `M_KK = xi_KK * Lambda_IR`,
- v1 default `xi_KK = 1.0`.

Because LFV scales as `(3 TeV / M_KK)^4`, changing `xi_KK` by O(1) factors materially changes pass/fail rates.

Systematics note:

- if later matching to specific first-KK conventions (gauge KK, fermion KK, etc.), typical mappings can imply `xi_KK ~ O(2-3)` and must be tracked as a named systematic run.

## 4) LFV Convention and Basis Definition

The scanner treats LFV input as the experimental bound on `BR(mu->e gamma)`, then derives the internal NDA coefficient `C`.

Code anchor in this repo:

- `flavorConstraints/muToEGamma.py` defines `PREFAC_BR` and
  `coefficient_from_br_limit(br_limit, prefactor=PREFAC_BR)`.

Definition used in this pipeline:

- `Ybar_N = 2k * Y_N_matrix`,
- `Y_N_matrix = U_PMNS * diag(Y_N_i)`,
- LFV is evaluated in the charged-lepton mass basis.

Current module relation:

`BR(mu->e gamma) ~= PREFAC_BR * |(Ybar_N Ybar_N^dagger)12|^2 * (3 TeV / M_KK)^4`,

which implies:

`C = sqrt(BR_limit / PREFAC_BR)`.

Reference values with current `PREFAC_BR = 4e-8`:

- Perez-Randall paper-era limit: `BR < 1.2e-11` -> `C ~= 0.0173` (often rounded to `0.02`).
- MEG II published limit (2025): `BR < 1.5e-13` -> `C ~= 0.00194`.

Required LFV metadata for reproducibility:

- `lfv_model` (example: `mu_to_e_gamma_nda_v1`),
- `br_limit`,
- `prefac_br`,
- `C_derived`,
- `reference_scale_GeV` (numeric, v1 default `3000.0`),
- `xi_KK`.

## 5) Oscillation Inputs and Acceptance Policy

v1 policy (explicit):

- use central oscillation inputs from `neutrinos/neutrinoValues.py`,
- do not scan oscillation parameters in v1,
- no oscillation chi2 in v1; acceptance is based on the model constraints already implemented in the scan pipeline.

Future extension modes can add sampled NuFIT windows or full oscillation chi2, but they must be a separately named run mode.

## 6) LMFV-Compatible c-Parameterization

The key restriction is not "no c scans"; it is "no arbitrary independent c scans".

Perez-Randall style logic implies a strong degeneracy prior for lepton doublets:

- use one shared `c_L0` in v1 (generation-universal),
- enforce a primary wavefunction prior `max_fL_ratio` (v1 default `1.1`) where
  pairwise LH overlap ratios must satisfy `max_ij[f_Li/f_Lj] <= max_fL_ratio`,
- derive `delta_cL_max` from the exact `f_IR(c, epsilon)` implementation for
  each `(k, Lambda_IR, c_L0)` and record it in metadata,
- if splittings are enabled later, apply the derived `delta_cL_max` rather than
  a hard-coded universal percentage on `c_L`.

For right-handed neutrinos in v1:

- use one shared `c_N0` (generation-universal),
- optional extension later: percent-level splittings.

For charged singlets:

- scan `c_E1, c_E2, c_E3` in moderate windows,
- implement ordering as a labeling convention by sorting sampled values
  descending (`c_E1 >= c_E2 >= c_E3`), not by reject-until-ordered sampling.

## 7) v1 Scan/Fixed Decision

| Parameter | v1 decision | Prior / value | Notes |
|---|---|---|---|
| `Lambda_IR` | scan | `3000` to `7000` GeV | outer geometry loop |
| `k` | fixed but explicit | `1.2209e19 GeV` | keep configurable in schema |
| `xi_KK` | fixed but explicit | `1.0` | maps `Lambda_IR` to `M_KK` |
| `c_L0` | scan | uniform linear `0.52` to `0.70` | universal in v1; `delta_cL_max` derived from `max_fL_ratio` |
| `c_N0` | scan | uniform linear `0.15` to `0.45` | universal in v1; avoids most pathological IR-heavy corner |
| `c_E1,c_E2,c_E3` | scan | each uniform linear `0.45` to `0.90` | sample then sort descending |
| `MN_mode` | fixed in baseline | `fixed_ratio` + preset `PR_benchmark` | means `MN_over_k = 0.1` |
| `MN_over_k` | default fixed | `0.1` | optional exploratory log scan |
| `lightest_nu_mass` | scan | log-uniform `1e-5` to `3e-2` eV | NO broad mode |
| `ordering` | fixed | `normal` | v1 consistency choice; global preference is dataset-dependent |
| `majorana_alpha,beta` | fixed | `0, 0` | see Sec. 9 scope |

Exploratory `MN_mode = scan_ratio` (optional): `MN_over_k` log-uniform over `1e-6` to `1`.

## 8) Lightest Neutrino Mass Prior

Using this repo's `compute_masses` for normal ordering:

- minimum sum at `m_lightest = 0`: `sum m_nu = 0.05898 eV`,
- if `sum m_nu < 0.12 eV`: max `m_lightest = 0.03008 eV`,
- if `sum m_nu < 0.072 eV`: max `m_lightest = 0.00868 eV`.

Scanner recommendation:

- broad mode: log-uniform `m_lightest in [1e-5, 3e-2] eV`,
- strict mode: log-uniform `m_lightest in [1e-5, 8.7e-3] eV`,
- treat exactly `m_lightest = 0` as a dedicated benchmark point, not part of the log prior.

Note: the strict `0.072 eV` option is scenario-dependent (dataset and cosmology-model dependent), so it should be a named run mode rather than a hard global default.

## 9) Majorana Phase Scope

In the current pipeline, LFV uses `Y_N Y_N^dagger` with
`Y_N = U_PMNS diag(y_i)` and Majorana phases on the right of `U_PMNS`.
For this observable, Majorana phases cancel.

So fixing phases to zero is safe for v1 acceptance under current filters.

Scope caveat: this does not imply phases are irrelevant for observables such as
`0nu beta beta`, leptogenesis-related quantities, or other phase-sensitive probes.

## 10) Anarchic Yukawa Prior and Score (Scanner Mode)

For this model, "anarchic" means random O(1) matrix structure with a separate overall neutrino scale.

v1 sampling definition:

- sample raw complex matrices `Ytilde_E`, `Ytilde_N` entrywise i.i.d.,
- magnitudes: log-uniform in `[1/3, 3]`,
- phases: uniform in `[0, 2*pi)`,
- no enforced matrix structure (not Hermitian, not triangular),
- set `Y_N = yN_overall * Ytilde_N`,
- scan `yN_overall` in a modest band (example `0.01` to `0.2`).

Minimal v1 anarchic score (documented default):

- `S = -w_band * P_band - w_cond * log(kappa(Ytilde_N))^2 - w_fit * chi2_total`,
- `P_band` penalizes entries outside `[1/3, 3]`,
- `kappa` is matrix condition number,
- `chi2_total` is optional in v1 (set `w_fit = 0` if no chi2 mode is active).

Default score weights must be stored in metadata (`w_band`, `w_cond`, `w_fit`).

Implementation note for this branch:

- scanner acceptance applies the anarchic score to the solved Yukawa point by
  decomposing `Ybar_N = yN_overall * Ytilde_N` with geometric-mean `yN_overall`,
- no additional random-matrix draw is used in the acceptance filter.

## 11) Cluster Reproducibility Metadata Contract

Every run should store at minimum:

- code identity: `git_commit` (full hash), `dirty_tree`,
- randomness: `rng_seed_global`, and deterministic per-sample seed (`rng_seed_sample` or deterministic derivation rule),
- scan priors and mode labels: all knobs in Sec. 7,
- c_L degeneracy metadata: `max_fL_ratio`, derived `delta_cL_max_*` fields,
- LFV metadata: all fields in Sec. 4,
- Yukawa sampling metadata: magnitude prior, phase prior, i.i.d. policy, any rescaling policy,
- score metadata: score formula id/version and weight values,
- acceptance policy metadata: oscillation input mode from Sec. 5.

## 12) Open Choices Before Production Cluster Run

1. Use broad or strict cosmology mode (`0.12` vs `0.072` eV sum scenario).
2. Keep `c_N0` lower bound at `0.15` (efficient default) or extend to `0.0` in exploratory runs.
3. Keep `MN_mode = fixed_ratio` with `PR_benchmark` (`MN_over_k = 0.1`) for first production, or run exploratory `scan_ratio` mode.
4. Keep `xi_KK = 1.0` fixed initially, or add a controlled systematic variation run.

## 13) References

- Perez, G. and Randall, L. *Natural Neutrino Masses and Mixings from Warped Geometry* (2008): https://arxiv.org/abs/0805.4652
- Csaki, C. et al. *Neutrino masses in models with warped extra dimensions* (2000): https://www.sciencedirect.com/science/article/pii/S0370269300012925
- Chen, M. C. and Yu, J. H. *Minimal Flavor Violation in the Lepton Sector of the Randall-Sundrum Model* (2008): https://www.sciencedirect.com/science/article/pii/S0370269308004745
- MEG II Collaboration, published limit (2025): https://arxiv.org/abs/2504.15711
- NuFIT 6.0 global analysis (2024): https://arxiv.org/abs/2410.05380
- Planck 2018 cosmology (A&A 2020): https://www.aanda.org/articles/aa/full_html/2020/09/aa33910-18/aa33910-18.html
- DESI BAO + CMB + SN mass-sum update (2024): https://www.osti.gov/pages/biblio/2371418
