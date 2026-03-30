# Quark Constraints Implementation Plan

Reference paper:
- [arXiv:0710.1869](https://arxiv.org/abs/0710.1869)
- [PDF](https://arxiv.org/pdf/0710.1869)
- Local copy: [`0710.1869v1.pdf`](0710.1869v1.pdf)

## Converged first step

Implement a narrow quark-sector MFV module that can reproduce the paper's core
numerical claim:

- RS bulk masses are built from two anarchic Yukawas, not from arbitrary
  independent quark-sector flavor matrices.
- The down sector becomes parametrically protected when the dial `r` is small.
- The model remains viable if quark masses and CKM mixing can be fit while the
  leading down-sector flavor-violation proxy stays suppressed.

The first deliverable should be a deterministic fit-and-diagnostic module, not
a full hadronic-physics package.

## Implement now

1. A quark-sector parameterization matching the paper's spurion structure:
   - `C_u,d ~ Y_u,d^\dagger Y_u,d`
   - `C_Q ~ r Y_u Y_u^\dagger + Y_d Y_d^\dagger`
2. A profile/eigenvalue helper that converts quark bulk-mass eigenvalues into
   overlap factors using the existing `warpConfig.wavefuncs` formulas.
3. A quark fit routine that takes target quark masses and CKM inputs and returns:
   - `c_Q`, `c_u`, `c_d`
   - `F_Q`, `F_u`, `F_d`
   - fitted Yukawa singular values / mixing angles
   - CKM residuals and mass residuals
4. A cheap flavor-violation proxy based on the paper's scaling estimate:
   - down-sector `Delta F = 2` suppression summary
   - the `r`-dependent hierarchy indicator
   - an up-sector `D - Dbar` proxy if the same fit machinery already exposes it

## Later

- Full kaon phenomenology, including precise `epsilon_K` and meson-mixing
  predictions.
- Full EDM calculations.
- Operator-basis and RG-evolution machinery for flavor observables.
- Combined quark+lepton global scans before the quark fit layer is validated.
- Precise exclusion claims tied to current experiments without sourcing
  external EFT and hadronic inputs.

## Paper-supported vs external-source-required

Paper-supported:

- the MFV spurion structure for `C_Q`, `C_u`, `C_d`
- the `r -> 0` down-sector protection limit
- the quark-mass relation `m_u,d ~= 2 v F_Q Y_u,d F_u,d`
- the CKM scaling logic `V_CKM ~ f_Qi / f_Qj`
- the `h_RS`-style scaling proxy
- the benchmark guidance `|r| ~ 0.1-0.4`, `y ~ 3`, and suppressed down-sector
  flavor violation

External-source-required:

- precise `Delta F = 2` Wilson coefficients and operator basis
- hadronic matrix elements and RG running for kaon, `B`, and `D` observables
- full EDM formulas and current bound interpretation
- any modern claim that the model is allowed or excluded at a given KK scale

## Required formulas, inputs, outputs

Use the same overlap formula already used elsewhere in the repo:

- `f^2(c) = (1/2 - c) / (1 - epsilon^(1 - 2c))`

Use the paper's quark mass and mixing relations:

- `m_u,d ~= 2 v F_Q Y_u,d F_u,d`
- `V_CKM ~ f_Qi / f_Qj`

Use the paper's leading `Delta F = 2` estimate as a fast proxy:

- `h_RS ~ 0.5 * (3 TeV / m_KK)^2 * (f_Q3 / 0.3)^4`

Use the paper's reported numerical guidance as validation targets:

- preferred `|r| ~ 0.1-0.4`
- representative Yukawa size `y ~ 3`
- down-sector suppression of order `O(0.25)` in the favored window

Inputs should be:

- matching-scale quark masses
- CKM parameters
- `r`
- overall Yukawa scale `y`
- optionally random seed and benchmark choice

Outputs should be:

- fit residuals for masses and CKM
- overlap eigenvalues and bulk-mass eigenvalues
- a compact down-sector suppression summary
- a boolean or score for whether the point is consistent with the paper's
  flavor-protection story

## Why this matters

This is the smallest numerically useful step because it tests the paper's
central claim in a way this repo can support now:

- the flavor problem is being moved from "generic RS excluded" to a controlled
  low-dimensional `r`-dependent fit problem,
- the existing `warpConfig` and `scanParams` infrastructure can be reused
  instead of re-invented,
- the output will tell us whether the paper's claimed `r` window is robust in
  the repo's conventions,
- and it creates a clean bridge from the current lepton-sector scans to a
  future unified quark + lepton flavor module.

## Recommended repo integration

- Put the numerical code in a future `quarkConstraints/fit.py` and
  `quarkConstraints/proxies.py`.
- Reuse `warpConfig.wavefuncs` for profiles and `diagonalization` for SVD
  / basis rotation.
- Add a future scan wrapper only after the fit routine is validated on the
  paper's benchmark window.
