**RS EDM Estimator Spec**

**Core Physics**
Tree-level KK-gluon exchange does **not** generate a quark EDM. In the repo’s mass-basis object, `couplings.left_down` and `right_down` are Hermitian, so diagonal entries are real and
\[
\mathrm{Im}\left[(G_L^d)_{11}(G_R^d)_{11}\right]=0.
\]
The leading RS result is instead the one-loop Higgs/KK-fermion dipole. Agashe-Perez-Soni show the KK-gluon loop is aligned with the down mass matrix and has vanishing imaginary part at leading order, while the Higgs/KK-fermion spurion is not aligned:
\[
C_{d_n}\propto
\left[F_Q\left(Y_uY_u^\dagger+Y_dY_d^\dagger\right)Y_dF_d\right]_{11}.
\]
See their Eqs. 48-51 and discussion of the RS CP problem. They estimate \(d_n/e\sim 10^{-24}\,\mathrm{cm}(3\,\mathrm{TeV}/M_{KK})^2\) for anarchic phases. [Agashe-Perez-Soni](https://escholarship.org/content/qt50c8s3km/qt50c8s3km.pdf)

Define the point-level down EDM invariant:
\[
A_d=U_{Ld}^\dagger F_QY_dF_dU_{Rd},\quad
B_d=U_{Ld}^\dagger F_Q(Y_uY_u^\dagger+Y_dY_d^\dagger)Y_dF_dU_{Rd},
\]
\[
S_d=B_{d,11}/A_{d,11},\qquad I_d=\mathrm{Im}(S_d).
\]
This is the dominant CP diagnostic. Add the up analogue \(S_u\) for the neutron formula.

Optional diagnostic from the already stored KK-gluon matrices:
\[
I_d^{G}=\sum_{j=s,b}\frac{m_j}{m_d}\,
\mathrm{Im}\left[(G_L^d)_{1j}(G_R^d)_{j1}\right].
\]
This is rephasing invariant, but should be reported as `kk_gluon_subleading`; central RS matching should default it to zero unless a deliberate EFT threshold estimate is requested.

**Inputs**
Existing:
`couplings.left_down/right_down`, `couplings.M_KK`, `couplings.alpha_s`; `f_Q`, `f_d`, `f_u`; `U_L_d`, `U_R_d`, optionally `U_L_u`, `U_R_u`; `masses_down`, `masses_up`.

Available in fitted path:
`QuarkFitResult.bulk_state.F_Q/F_u/F_d`, `Y_u_bulk_basis`, `Y_d_bulk_basis`, rotations and masses.

Available in Bauer forward generator but not stored in current parquet:
local `Y_u`, `Y_d`, `f_Q`, `f_d`, rotations, masses. The instrumented parquet stores only `G_L/G_R` off-diagonal summaries, so it cannot compute the dominant Higgs/KK-fermion EDM unless `S_d`, `S_u`, or raw Yukawas/rotations are added.

Missing for production accuracy:
KK-fermion spectrum, loop functions, wrong-chirality Yukawas \(\tilde Y_{u,d}\), Higgs localization prescription, brane counterterm estimate, and dipole RG helper. Koenig-Neubert-Straub emphasize wrong-chirality Yukawas and dipole constraints in composite/warped models. [KNS](https://arxiv.org/pdf/1403.2756)

**Matching Convention**
Use KNS dipole operators:
\[
Q_{qq\gamma}=\frac{em_q}{16\pi^2}\bar q\sigma^{\mu\nu}P_RqF_{\mu\nu},\quad
Q_{qqg}=\frac{g_sm_q}{16\pi^2}\bar q\sigma^{\mu\nu}T^aP_RqG^a_{\mu\nu}.
\]
Then
\[
d_q=\frac{em_q}{8\pi^2}\mathrm{Im}\,C_{qq\gamma},\quad
\tilde d_q=\frac{g_sm_q}{8\pi^2}\mathrm{Im}\,C_{qqg}.
\]
[KNS Eqs. 11-16](https://arxiv.org/pdf/1403.2756)

Central estimator:
\[
C_{d\gamma}(M_{KK})=\frac{\kappa_\gamma}{2M_F^2}S_d,\quad
C_{dg}(M_{KK})=\frac{\kappa_g}{2M_F^2}S_d,
\]
with defaults `M_F=M_KK`, `kappa_gamma=kappa_chromo=1.0`, uncertainty flag `O(3)`.

Run \((C_\gamma,C_g)\) to \(\mu_l=1\) or \(2\) GeV with LL mixing:
\[
C_\gamma(\mu_l)=\eta_{\gamma\gamma}C_\gamma(M)+\eta_{\gamma g}C_g(M),\quad
C_g(\mu_l)=\eta_{gg}C_g(M).
\]

**Hadronic Translation**
Use QCD sum-rule estimate:
\[
d_n=(1\pm0.5)\left[1.4(d_d-0.25d_u)+1.1e(\tilde d_d+0.5\tilde d_u)\right].
\]
This is the standard Pospelov-Ritz formula under PQ/removed \(\bar\theta\). [Pospelov-Ritz](https://arxiv.org/abs/hep-ph/0010037)

Current neutron bound: PDG 2026 lists
\[
|d_n|<1.8\times10^{-26}\ e\,\mathrm{cm}\quad(90\%\ \mathrm{CL}),
\]
from Abel 2020, with measured value
\[
d_n=(0.0\pm1.1_{\rm stat}\pm0.2_{\rm sys})\times10^{-26}\ e\,\mathrm{cm}.
\]
[PDG Live](https://pdgprod.lbl.gov/pdgprod/pdgLive/DataBlock.action?node=S017EDM), [Abel 2020](https://arxiv.org/abs/2001.11966)

**Function Signature**
```python
def edm_from_point(
    *,
    couplings,
    Y_u=None, Y_d=None,
    f_Q=None, f_u=None, f_d=None,
    U_L_u=None, U_R_u=None, U_L_d=None, U_R_d=None,
    masses_up=None, masses_down=None,
    M_F=None,
    mu_low=1.0,
    include_up=True,
    include_kk_gluon_diagnostic=True,
    kappa_gamma=1.0,
    kappa_chromo=1.0,
    uncertainty="O3",
) -> dict:
    ...
```

Return keys:
`d_d_e_cm`, `d_u_e_cm`, `tilde_d_d_cm`, `tilde_d_u_cm`, `d_n_e_cm`, `d_n_chromo_e_cm`, `S_d`, `S_u`, `I_d_G`, `bound_ratio`, `uncertainty_flags`, `missing_inputs`.

**Algorithm**
1. Validate shapes and use `M_KK = couplings.M_KK`; set `M_F=M_KK` unless supplied.
2. If Yukawas/rotations/f-factors exist, compute `S_d`; if `include_up`, compute `S_u`.
3. Match `C_gamma`, `C_g` at `M_F` using the NDA equations above.
4. Evolve dipoles to `mu_low`.
5. Convert to `d_q/e` and `tilde_d_q`; multiply GeV\(^{-1}\) by `1.973269804e-14` to get cm.
6. Apply neutron formula; return central and `hadronic_scale_factor=(0.5,1.5)`.
7. If only `G_L/G_R` are available, return no central Higgs EDM and only `I_d_G` diagnostic.

**Sanity Checks**
1. Rephase \(d_i\to e^{i\alpha_i}d_i\): `S_d`, `S_u`, and `I_d_G` must not change.
2. Set all Yukawas real: EDM invariants should vanish up to numerical noise.
3. An anarchic \(Y_*\sim3\), \(M_{KK}=3\) TeV point with \(\mathrm{Im}S_d=O(1-10)\) should give \(d_n\sim10^{-26}-10^{-24}\ e\,\mathrm{cm}\), with explicit `O3` uncertainty.