# μ→eγ Constraint (NDA Dipole)

This module implements the **brane-localized dipole NDA estimate** used in
Perez & Randall (arXiv:0805.4652) for the μ→eγ bound.

## Definitions

In the charged-lepton mass basis:

- Neutrino Yukawa matrix:
  \[
  Y_N = U_{\rm PMNS} \;\mathrm{diag}(Y_{N_1},Y_{N_2},Y_{N_3})
  \]

- Rescaled Yukawa:
  \[
  \bar Y_N \equiv 2k\,Y_N
  \]

- Product entering LFV:
  \[
  \bar Y_N \bar Y_N^\dagger = (2k)^2\,U\,\mathrm{diag}(Y_{N_i}^2)\,U^\dagger
  \]

We take the off-diagonal element
\(|(\bar Y_N \bar Y_N^\dagger)_{12}|\) as the relevant LFV spurion.

## NDA Dipole Estimate (Perez–Randall)

Perez & Randall estimate the μ→eγ branching ratio from a **brane-localized
(TeV-brane) dipole operator** as

\[
\mathrm{BR}(\mu\to e\gamma) \simeq 4\times10^{-8}
\left|(\bar Y_N \bar Y_N^\dagger)_{12}\right|^2
\left(\frac{3~\mathrm{TeV}}{M_{KK}}\right)^4.
\]

Using the (paper-era) experimental limit
\(\mathrm{BR}(\mu\to e\gamma) < 1.2\times10^{-11}\),
they quote the bound

\[
\left|(\bar Y_N \bar Y_N^\dagger)_{12}\right| \le C\left(\frac{M_{KK}}{3~\mathrm{TeV}}\right)^2,
\qquad C=0.02.
\]

This **C = 0.02** is the default in `check_mu_to_e_gamma()`.

### Updated limit (MEG II 2024)

The MEG II 2024 result gives
\(\mathrm{BR}(\mu\to e\gamma) < 7.5\times10^{-13}\) (90% CL).
Using the same NDA formula above:

\[
C_{\rm MEGII} = \sqrt{\frac{7.5\times10^{-13}}{4\times10^{-8}}} \approx 4.33\times10^{-3}.
\]

You can pass this as `C` explicitly or compute it via
`coefficient_from_br_limit()`.

## API

- `check_mu_to_e_gamma(yukawa_result, C=C_PAPER, reference_scale=3000)`
- `check_mu_to_e_gamma_raw(Y_N_bar, pmns, M_KK, C=C_PAPER, reference_scale=3000)`
- `coefficient_from_br_limit(br_limit, prefactor=4e-8)`

## Notes

- The **bulk-Higgs** estimate in the paper is much weaker
  (they quote \(C \simeq 0.3\)). If you want that case, pass a larger `C`.
- The check is purely the NDA dipole bound; it does **not** include full loop
  factors or flavor model details.
