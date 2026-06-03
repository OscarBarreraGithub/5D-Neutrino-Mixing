# Independent design: RS electroweak-coupling sector

Independence note: this design was produced from the requested source files and
derivation files only. I did not search for or read any
`rs_ew_sector_design*.md` file.

## Design target

Replace the catalog's current electroweak new-physics proxies with a single
RS electroweak sector builder that produces mass-basis Z, W, W'/Z'/gamma'
couplings, the Wilson coefficients derived from them, and explicit partial
status for genuinely loop-hard pieces. The validated SM-side and observable
arithmetic in the existing `physics_adapters` should remain the owner of
branching-ratio, pull, and budget calculations.

The contract should be:

- tree-level neutral-current and charged-current RS matching is computed once
  per point, stored on `ParameterPoint.extras`, and reused by every constraint;
- adapters consume rigorous couplings or Wilson bundles, not overlap proxies;
- unknown extras still fail loudly in `make_point`;
- complex amplitudes and Wilson coefficients stay in diagnostics or typed
  extras, never in `ConstraintResult.predicted`, `ratio`, or `budget`.

## 1. Physics inventory

### Conventions and existing machinery

Use the repo conventions:

- Geometry: `epsilon = Lambda_IR / k`, `z_h = 1/k`, `z_v = 1/Lambda_IR`,
  `L = warp_log = log(k/Lambda_IR)` from `warpConfig/baseParams.py`.
- Bulk mass convention: `c = M_5/k`.
- Zero-mode endpoint overlaps:
  `F_IR(c) = f_IR(c, epsilon)` and `F_UV(c) = f_UV(c, epsilon)` from
  `warpConfig/wavefuncs.py`.
- Gauge KK roots: `solve_kk(species="gauge", bc="NN", geometry, exact=True)`,
  with `m_n = x_n Lambda_IR`, from `solvers/bessel.py`.
- PMNS and lepton Yukawas: `compute_all_yukawas()` and `get_pmns()`.

Derivation references:

- `derivations/conventions.tex`: geometry, c convention, f factors, gauge NN
  quantization, charged-lepton masses, seesaw, PMNS.
- `derivations/zero_modes/README.md`: `f_IR`/`f_UV` are done; full zero-mode
  `g(z)` profiles are still marked NEEDED.
- `derivations/kk_modes/README.md`: gauge NN masses are specified; gauge/fermion
  profile normalization is still marked NEEDED.
- `derivations/kk_modes/bessel_equation.tex`: Bessel roots, brane Yukawa
  normalization, and zero-mode mass formula.
- `derivations/flavor/README.md`: `mu_to_e_gamma.tex` is still NEEDED, so
  loop-level lepton dipoles cannot be promoted to rigorous RS yet.

Important gap: the repo does not currently contain a completed derivation file
for RS gauge-KK couplings, Z-coupling shifts, W shifts, or dipole loop
matching. The design below gives the implementation formula, but the first
implementation PR should add a derivation note for these formulas before
shipping production matching.

### KK electroweak spectrum

Quantity needed:

- `m_Wprime_n_gev`, `m_Zprime_n_gev`, `m_gammprime_n_gev` for the heavy EW
  vectors.
- `m_W_gev`, `m_Z_gev` for the light eigenstates after EWSB.

Formula:

1. Compute geometry:
   `geom = get_warp_params(k=k, Lambda_IR=Lambda_IR)`.
2. Compute gauge NN roots:
   `m_gauge_n = solve_kk("gauge", "NN", geom, n_roots=N)[0]`.
3. Build the charged vector mass matrix in basis
   `(W_zero, W_KK1, ..., W_KKN)`:

   ```text
   M^2_charged[m,n] = m_KK[m]^2 delta_mn
                      + (v^2/4) g_2^2 chi_m(z_v) chi_n(z_v)
   ```

   with `m_KK[0] = 0` and `chi_0(z_v) = 1` in the zero-mode normalization.

4. Build the neutral vector mass matrix in basis
   `(W3_zero, B_zero, W3_KK1, B_KK1, ...)`:

   ```text
   M^2_neutral[(a,m),(b,n)] =
       m_a,m^2 delta_ab delta_mn
       + (v^2/4) q_ab chi_m(z_v) chi_n(z_v)

   q = [[g_2^2, -g_2 g_Y],
        [-g_2 g_Y, g_Y^2]]
   ```

5. Diagonalize the matrices. Identify photon, W, Z by smallest eigenvalues and
   SM-like eigenvectors; the first mostly-KK charged and neutral states define
   `kk_ew_mass_gev`, `m_wprime_1_gev`, `m_zprime_1_gev`.

The exact `chi_n(z)` normalization is a small missing helper over the existing
Bessel solver. It should be implemented from the normalization formula already
listed in `derivations/kk_modes/README.md`, then cached by
`(epsilon, n_roots)`.

### Gauge-KK couplings to zero-mode fermions

Quantity needed:

- Mode-by-mode interaction-basis overlap factors `Omega_n(c)` and mass-basis
  coupling matrices for all quark and lepton chiralities.

Exact formula:

```text
Omega_n(c) = integral_{z_h}^{z_v} dz w_0(c,z) chi_n(z)

w_0(c,z) = normalized zero-mode probability density,
           integral dz w_0(c,z) = 1
```

The endpoint normalization must be consistent with `f_IR(c, epsilon)`. The
profile helper should be derived so that evaluating the normalized zero mode at
`z_v` reproduces `f_IR` and at `z_h` reproduces `f_UV` in the repo convention.

For a gauge group `X` and a chiral species `f_A`:

```text
G_X,fA,n^mass = g_X q_X(f_A) U_fA^\dagger diag(Omega_n(c_fA_i)) U_fA
```

where `q_X` is `T3` for SU(2)_L doublets and `Y` for hypercharge. Singlets have
`T3 = 0`. For left quarks, use separate rotations into up and down mass bases:

```text
G_X,uL,n = g_X q_X(Q_L) U_L_u^\dagger diag(Omega_n(c_Q)) U_L_u
G_X,dL,n = g_X q_X(Q_L) U_L_d^\dagger diag(Omega_n(c_Q)) U_L_d
```

and analogously for right up, right down, charged leptons, and neutrinos. This
is the rigorous replacement for using KK-gluon overlap matrices as electroweak
proxies.

Fast-scan implementation: precompute `Omega_n(c)` on the actual c grid, store a
small interpolation table, and keep only the first few modes or a validated KK
sum. The 100M-point path must not numerically integrate profiles per point.

### Z-fermion coupling shifts

Quantity needed:

- Dimensionless matrices in the repo Z-pole convention:

  ```text
  L_Z = g_Z Z_mu fbar gamma^mu (g_L P_L + g_R P_R) f
  delta_g_Z_fA[species][i,j] = g_Z^RS_fA[i,j] / g_Z - g_Z^SM_fA delta_ij
  ```

  Required species:
  `u_L`, `u_R`, `d_L`, `d_R`, `e_L`, `e_R`, `nu_L`.

Rigorous computing formula:

1. Build and diagonalize the neutral vector mass matrix above.
2. For the light Z eigenvector `Z = sum_{X,n} V_Z[X,n] X_n`, compute:

   ```text
   g_Z^RS(f_A) = sum_{X,n} V_Z[X,n] G_X,fA,n^mass
   delta_g_Z_fA = g_Z^RS(f_A)/g_Z - g_Z^SM(f_A) I
   ```

3. Use `g_Z^SM = T3 - Q sin^2(theta_W)` for left doublets and
   `g_Z^SM = -Q sin^2(theta_W)` for singlets, matching
   `quarkConstraints.zpole`.

This covers:

- diagonal Z-pole shifts for `Zbb`, `Zcc`, and lepton Z-pole observables;
- off-diagonal down FCNC-Z for `Z -> q_i q_j`;
- off-diagonal charged-lepton Z couplings for `Z -> e mu`, `e tau`, `mu tau`;
- the tree-level Z-penguin component of `s -> d l l`, `b -> q l l`,
  `c -> u l l`, and `b -> s nu nubar`.

Classic RS Zbb:

- Include both gauge-vector mixing and fermion-KK mixing contributions.
- Gauge-vector mixing is covered by the vector diagonalization formula.
- Fermion-KK mixing is the additional tree-level effect from diagonalizing the
  EWSB mass matrices including zero and fermion KK modes:

  ```text
  M_f(v) =
      [[m_00, m_0K],
       [m_K0, M_KK + m_KK-Higgs]]

  g_Z^mass = U_L^\dagger g_Z^interaction U_L   and   U_R^\dagger g_Z^interaction U_R
  delta_g_Zbb = g_Z^mass[zero,zero]/g_Z - g_Z^SM
  ```

  This is needed for the top-Yukawa enhanced `Z b_L b_L` shift. It requires
  fermion KK profiles and brane Yukawa overlap matrices. Since
  `derivations/kk_modes/README.md` marks fermion normalization as NEEDED, the
  first rigorous implementation should either add this derivation or explicitly
  label `Zbb_fermion_mixing_included=False`. Do not silently call a gauge-only
  Zbb result "classic RS Zbb".

### Flavor off-diagonal Z couplings

Quantity needed:

- `delta_g_Z_dL`, `delta_g_Z_dR`, `delta_g_Z_uL`, `delta_g_Z_uR` for FCNC.
- `delta_g_Z_eL`, `delta_g_Z_eR` for charged LFV.

Formula:

```text
delta_g_Z_fA[i,j] = (1/g_Z) sum_{X,n} V_Z[X,n] g_X q_X(f_A)
                    (U_fA^\dagger diag(Omega_n(c_fA)) U_fA)[i,j]
                    - g_Z^SM(f_A) delta_ij
```

For off-diagonal entries the SM term is zero. These matrices are the direct
replacement for:

- `zpole_down_fcnc_coupling_proxy`;
- `z_lfv_coupling_proxy`;
- the Z-like pieces in rare `C9/C10`, `Y`, and `X` Wilson proxies.

### Heavy neutral-vector four-fermion coefficients

Rare decays should not use only the shifted light-Z coupling. Heavy `Z'`,
`gamma'`, and mixed neutral KK eigenstates also generate direct contact terms.
After diagonalizing the neutral vector matrix, for each heavy neutral eigenstate
`V_a`:

```text
C_AB^{q_i q_j l_a l_b} =
    (g_Z^2 / m_Z^2) delta_g_Z_qA[i,j] g_Z_lB[a,b]
    + sum_{V_heavy} g_V_qA[i,j] g_V_lB[a,b] / M_V^2
```

Units: `GeV^-2`. `A,B in {L,R}` for charged leptons; neutrinos only have `B=L`.
This bundle is the rigorous input for semileptonic and LFV rare decays.

### C9/C10 matching

For `q_j -> q_i l_a l_b`, define `lambda_t` as the adapter already does. From
the four-fermion coefficients above:

```text
C9_NP  = -pi/(sqrt(2) G_F alpha lambda_t)
         * (C_LL + C_LR)
C10_NP = -pi/(sqrt(2) G_F alpha lambda_t)
         * (C_LR - C_LL)
C9p_NP  = -pi/(sqrt(2) G_F alpha lambda_t)
          * (C_RL + C_RR)
C10p_NP = -pi/(sqrt(2) G_F alpha lambda_t)
          * (C_RR - C_RL)
```

This replaces the current "one Z-like boson at M_KK" proxy in:

- `rare_b_dilepton` and exclusive/inclusive `b -> s l l` adapters;
- `rare_kaon_dilepton` and LFV rare-kaon dilepton adapters;
- `rare_charm_dilepton` and LFV rare-charm dilepton adapters.

Scalar/pseudoscalar Wilsons remain zero unless a separate Higgs/radion sector
is implemented.

### b -> s nu nubar matching

Use the same four-fermion coefficient, with neutrino left-handed couplings:

```text
X_NP_L = C_LL^{s b nu nu} / g_SM^2
X_NP_R = C_RL^{s b nu nu} / g_SM^2
g_SM^2 = 4 G_F^2 M_W^2 / (2 pi^2)
```

This replaces `RARE_B_NUNU_RS_MATCHING_ASSUMPTION_V1` and feeds the existing
`epsilon`, `eta`, `R_K`, and `R_Kstar` response.

### Charged-current and W shifts

Quantity needed:

- `delta_g_W_ud_L`, `delta_g_W_ud_R`: 3x3 complex matrices in the convention
  `L_W = g_2/sqrt(2) W^+ ubar gamma^mu [(V + delta_g_L) P_L + delta_g_R P_R] d`.
- `delta_g_W_lnu_L`, `delta_g_W_lnu_R`: lepton charged-current matrices.
- `delta_G_F_over_G_F`: real scalar.

Computing formula:

1. Diagonalize the charged vector mass matrix.
2. For the light W eigenvector:

   ```text
   g_W^RS(u_i d_j)_A = sum_n V_W[n] G_W,fA,n^mass[i,j]
   delta_g_W_ud_L = sqrt(2)/g_2 * g_W^RS(u_L d_L) - V_CKM
   delta_g_W_ud_R = sqrt(2)/g_2 * g_W^RS(u_R d_R)
   ```

3. The measured muon-decay normalization is:

   ```text
   G_F^meas = G_F^0 * (1 + delta_G_F_over_G_F)
   ```

   where the shift is computed from light-W coupling shifts plus direct W'
   exchange in `mu -> e nu nu`.

4. Semileptonic extraction shifts:

   ```text
   A(d_j -> u_i l_a nu_a) / A_SM =
       1
       + delta_g_W_ud_L[i,j] / V_ij
       + delta_g_W_lnu_L[a,a]
       + C_Wprime_LL[i,j,a,a] / C_SM
       - delta_G_F_over_G_F
       + right-handed/scalar terms
   ```

This makes EW002, EW003, and K018 point-dependent. EW003 still has
inclusive/exclusive covariance and form-factor interpretation risk, but the RS
charged-current amplitude itself becomes rigorous.

### Lepton and neutrino profiles

Quantity needed:

- `c_L[3]`, `c_E[3]`, `c_N[3]`.
- `f_L_IR[3]`, `f_E_IR[3]`, `f_N_IR[3]`, `f_N_UV[3]`.
- `Y_E_bar[3,3]`, `Y_N_bar[3,3]`, PMNS, charged-lepton rotations.

Current repo path:

- `compute_all_yukawas()` supports universal `c_L` and `c_N`, diagonal
  `c_E`, diagonal charged leptons, and `Y_N = PMNS diag(Y_N_i)`.

Required extension:

- Accept vector `c_L` and `c_N`, or explicitly mark universal mode.
- Store rotations even when they are identity, so the same builder works for
  future non-universal lepton scans.
- Build `Y_N_bar Y_N_bar^\dagger` only as a diagnostic/dipole spurion, not as
  a rigorous loop prediction.

### Dipole coefficients

Needed by:

- `mu -> e gamma`, `tau -> mu gamma`, `tau -> e gamma`;
- `b -> s gamma` and exclusive radiative B modes;
- charged-lepton EDMs and quark EDM/cEDM inputs.

Honest status:

- Not fully rigorous from current derivations. `derivations/flavor/README.md`
  explicitly marks `mu_to_e_gamma.tex` as NEEDED, and `derivations/integrals`
  marks triple overlaps/loop integrals as future.
- The builder can compute the spurions and profile inputs required by a loop
  matcher:

  ```text
  lepton_dipole_spurion_LL = Y_N_bar Y_N_bar^\dagger
  quark_dipole_profile_data = c_Q,c_u,c_d,U_L/R,Y_u,Y_d,fermion_KK_masses
  ```

- A rigorous dipole object must contain low-scale Wilsons only after a real
  one-loop matcher is implemented:

  ```text
  C7_bsgamma(M_KK), C8_bsgamma(M_KK), C7(mu_b), C8(mu_b)
  A_L/R(l_i -> l_j gamma)
  d_l/e = Im(C_l^dipole)
  qEDM, qCEDM, Weinberg operator
  ```

Until then, these constraints remain partial or NEEDS-HUMAN.

### Higgs LFV

Needed by T018-T020:

```text
Y_h^mass = U_L^\dagger (d M_e(v) / d v) U_R
BR(h -> l_i l_j) = m_h (|Y_ij|^2 + |Y_ji|^2)/(8 pi Gamma_h)
```

Honest status:

- With only diagonal zero-mode charged leptons, `Y_h^mass` is diagonal.
- Off-diagonal Higgs LFV requires full fermion-KK mixing and Higgs-profile
  assumptions. This is tree-level but not available from the current builder.
- Mark T018-T020 partial until the fermion-KK mass-matrix block is added.

### Mu-e conversion and LFV three-body decays

Tree-level vector pieces become rigorous from the neutral-current bundle:

```text
G_AB(l_i -> 3 l_j) =
    [ (g_Z^2/m_Z^2) delta_g_Z_eA[j,i] g_Z_eB[j,j]
      + sum_V g_V_eA[j,i] g_V_eB[j,j]/M_V^2
      + box_AB ] / (4 G_F/sqrt(2))
```

For coherent mu-e conversion:

```text
g_LV^q = C_LV^q / (4 G_F/sqrt(2))
g_RV^q = C_RV^q / (4 G_F/sqrt(2))

g_LV^p = 2 g_LV^u + g_LV^d
g_LV^n = g_LV^u + 2 g_LV^d
```

and similarly for `R`. Scalar coefficients need Higgs/radion matching and stay
partial. Dipole-contact interference remains phase-sensitive until dipoles are
rigorous.

### EDMs

Honest status:

- Charged-lepton EDMs need CP-odd loop dipole matching.
- Neutron, mercury, atomic EDMs need both RS quark EDM/cEDM/Weinberg matching
  and hadronic/nuclear matrix elements.
- These should remain NEEDS-HUMAN after the RS-EW tree-level sector. The new
  schema can carry explicit low-energy Wilsons if supplied, but the builder
  must not invent them.

## 2. ParameterPoint schema

Keep `ParameterPoint` frozen. Add keys only through `KNOWN_EXTRA_KEYS`.

### Existing keys to keep

| Extra key | Type | Units | Status |
|---|---|---:|---|
| `quark_mass_basis_couplings` | `QuarkMassBasisCouplings` | mixed | Keep for KK-gluon and backward compatibility. |
| `kk_gluon_mass_gev` | `float` | GeV | Keep. |
| `kk_ew_mass_gev` | `float` | GeV | Fill with first physical EW gauge KK mass, not `Lambda_IR`. |
| `lepton_mass_basis_couplings` | `RSLeptonMassBasisCouplings` | mixed | Reuse existing placeholder and fill it. |

### New keys

| Extra key | Type | Units | Real/complex policy |
|---|---|---:|---|
| `rs_ew_spectrum` | `RSEWSpectrum` | GeV | Real arrays only. |
| `rs_ew_couplings` | `RSEWMassBasisCouplings` | dimensionless and GeV^-2 | Complex matrices allowed in typed extra. |
| `rs_semileptonic_wilsons` | `RSSemileptonicWilsonBundle` | dimensionless Wilsons | Complex Wilsons allowed in typed extra; adapters put copies in diagnostics. |
| `rs_charged_current` | `RSChargedCurrentCouplings` | dimensionless | Complex matrices allowed. |
| `rs_dipole_wilsons` | `RSDipoleWilsonBundle` | dimensionless or GeV^-1 | Optional/partial; absent unless loop matcher ran. |
| `rs_higgs_yukawas` | `RSHiggsYukawaCouplings` | dimensionless | Optional/partial; complex matrix. |

### Dataclass sketches

`RSEWSpectrum`

```python
@dataclass(frozen=True)
class RSEWSpectrum:
    model_label: str
    k_gev: float
    lambda_ir_gev: float
    epsilon: float
    warp_log: float
    n_gauge_modes: int
    gauge_roots_x: np.ndarray        # shape (N,), real
    gauge_masses_gev: np.ndarray     # shape (N,), real
    m_w_gev: float
    m_z_gev: float
    m_wprime_gev: np.ndarray         # shape (N,), real
    m_zprime_gev: np.ndarray         # shape (N,), real
    m_gammaprime_gev: np.ndarray     # shape (N,), real
    exact_bessel: bool
    profile_normalization: str
```

`RSLeptonMassBasisCouplings`

```python
@dataclass(frozen=True)
class RSLeptonMassBasisCouplings:
    model_label: str
    c_L: np.ndarray                  # shape (3,), real
    c_E: np.ndarray                  # shape (3,), real
    c_N: np.ndarray                  # shape (3,), real
    f_L_IR: np.ndarray               # shape (3,), real
    f_E_IR: np.ndarray               # shape (3,), real
    f_N_IR: np.ndarray               # shape (3,), real
    f_N_UV: np.ndarray               # shape (3,), real
    Y_E_bar: np.ndarray              # shape (3,3), complex
    Y_N_bar: np.ndarray              # shape (3,3), complex
    U_e_L: np.ndarray                # shape (3,3), complex
    U_e_R: np.ndarray                # shape (3,3), complex
    U_nu_L: np.ndarray               # shape (3,3), complex
    pmns: np.ndarray                 # shape (3,3), complex
    lfv_dipole_spurion: np.ndarray   # Y_N_bar Y_N_bar^dagger, complex
    m_kk_ew_gev: float
    matching_status: Mapping[str, str]
```

`RSEWMassBasisCouplings`

```python
@dataclass(frozen=True)
class RSEWMassBasisCouplings:
    model_label: str
    spectrum: RSEWSpectrum
    z_delta_g_L_u: np.ndarray        # 3x3 complex, dimensionless
    z_delta_g_R_u: np.ndarray
    z_delta_g_L_d: np.ndarray
    z_delta_g_R_d: np.ndarray
    z_delta_g_L_e: np.ndarray
    z_delta_g_R_e: np.ndarray
    z_delta_g_L_nu: np.ndarray
    z_total_g_L_u: np.ndarray        # SM + delta, dimensionless
    z_total_g_R_u: np.ndarray
    z_total_g_L_d: np.ndarray
    z_total_g_R_d: np.ndarray
    z_total_g_L_e: np.ndarray
    z_total_g_R_e: np.ndarray
    z_total_g_L_nu: np.ndarray
    neutral_contact_LL: Mapping[str, np.ndarray]  # GeV^-2 matrices by transition
    neutral_contact_LR: Mapping[str, np.ndarray]
    neutral_contact_RL: Mapping[str, np.ndarray]
    neutral_contact_RR: Mapping[str, np.ndarray]
    includes_heavy_neutral_exchange: bool
    includes_light_z_shift: bool
    includes_fermion_kk_mixing: bool
    diagnostics: Mapping[str, Any]
```

`RSChargedCurrentCouplings`

```python
@dataclass(frozen=True)
class RSChargedCurrentCouplings:
    model_label: str
    delta_g_W_ud_L: np.ndarray       # 3x3 complex, dimensionless
    delta_g_W_ud_R: np.ndarray       # 3x3 complex, dimensionless
    delta_g_W_lnu_L: np.ndarray      # 3x3 complex, dimensionless
    delta_g_W_lnu_R: np.ndarray      # 3x3 complex, dimensionless or zero
    delta_G_F_over_G_F: float        # real
    charged_contact_LL: Mapping[str, complex]  # GeV^-2 by flavor channel
    diagnostics: Mapping[str, Any]
```

`RSSemileptonicWilsonBundle`

```python
@dataclass(frozen=True)
class RSSemileptonicWilsonBundle:
    model_label: str
    b_to_s_ll: Mapping[str, complex]     # C9,C10,C9p,C10p by lepton pair
    b_to_d_ll: Mapping[str, complex]
    s_to_d_ll: Mapping[str, complex]     # Y_NP_L/R and y7V/y7A
    c_to_u_ll: Mapping[str, complex]
    b_to_s_nunu: Mapping[str, complex]   # X_NP_L/R
    lfv_llqq: Mapping[str, complex]
    matching_scale_gev: float
    diagnostics: Mapping[str, Any]
```

`RSDipoleWilsonBundle` should be absent by default. When present, it must
state `matching_status="loop-matched"` and contain actual low-scale Wilsons,
not NDA spurions.

## 3. Point builder logic

Add a new construction path:

```python
def build_from_rs_ew_inputs(
    *,
    quark_fit_result: QuarkFitResult | None = None,
    quark_couplings: QuarkMassBasisCouplings | None = None,
    lepton_inputs: RSLeptonSweepInputs,
    Lambda_IR: float,
    k: float = MPL,
    n_gauge_modes: int = 1,
    ew_model: Literal["minimal"] = "minimal",
    include_fermion_kk_mixing: bool = False,
    raw: Any = None,
) -> ParameterPoint:
    ...
```

Implementation steps:

1. Compute geometry with `get_warp_params`.
2. Compute EW gauge roots with `solve_kk("gauge", "NN", geometry, n_roots)`.
3. Fill `kk_ew_mass_gev` with the first physical EW gauge KK mass:
   `gauge_masses_gev[0]`, not `Lambda_IR`.
4. Build or accept quark mass-basis data:
   - If `quark_fit_result` is supplied, use `bulk_state.c_Q,c_u,c_d,F_Q,F_u,F_d`
     and rotations `U_L_u,U_L_d,U_R_u,U_R_d`.
   - If only `QuarkMassBasisCouplings` is supplied, keep old quark extra but
     mark exact EW overlap reconstruction unavailable unless the needed c
     vectors and rotations are attached.
5. Build lepton data:
   - Current compatible mode: call `compute_all_yukawas()` for universal
     `c_L`, `c_N`, diagonal `c_E`, `M_N`, lightest mass, ordering, Majorana
     phases.
   - Store `Y_E_bar` as a diagonal 3x3 matrix and `Y_N_bar = PMNS diag(Y_N_bar_i)`.
   - Store identity charged-lepton rotations unless a future non-diagonal
     charged-lepton Yukawa input is supplied.
6. Precompute or look up `Omega_n(c)` for each c value.
7. Build charged and neutral vector mass matrices and diagonalize them.
8. Compute `RSEWMassBasisCouplings`, `RSChargedCurrentCouplings`, and
   `RSSemileptonicWilsonBundle`.
9. Return `make_point(raw=raw, quark_mass_basis_couplings=..., kk_gluon_mass_gev=...,
   kk_ew_mass_gev=..., lepton_mass_basis_couplings=..., rs_ew_spectrum=...,
   rs_ew_couplings=..., rs_charged_current=..., rs_semileptonic_wilsons=...)`.

### New sweep inputs required

Current quark-only `build_from_quark_couplings` is insufficient. The scan row
must add:

| Input | Type | Reason |
|---|---|---|
| `Lambda_IR`, `k` | floats, GeV | geometry and physical EW KK mass. |
| `n_gauge_modes` | int | KK truncation and validation. |
| `c_L` or `c_L[3]` | real | lepton doublet overlaps and LFV Z/W couplings. |
| `c_E[3]` | real | charged-lepton singlet overlaps. |
| `c_N` or `c_N[3]` | real | neutrino profiles and seesaw. |
| `M_N` or `M_N[3]` | GeV | UV Majorana mass. |
| `lightest_nu_mass`, `ordering`, `majorana_alpha`, `majorana_beta` | scalar | PMNS/seesaw. |
| `Y_E_bar[3,3]` optional | complex | needed for non-diagonal charged-lepton scans. |
| `Y_N_bar[3,3]` optional | complex | lets future scans bypass universal seesaw inversion. |
| quark `c_Q,c_u,c_d` and rotations | arrays | exact EW overlaps; current coupling-only object is not enough. |
| `g_2`, `g_Y`, `sin2_theta_w`, `v` policy | real | coupling normalization. |
| `include_fermion_kk_mixing` | bool | whether Zbb/Higgs LFV include fermion-KK mixing. |

For a 100M-point scan, all grid-valued c inputs should be integer-coded or
cached so overlap table lookups are O(1).

## 4. Matching map

The existing adapters should add a rigorous branch:

```text
if point has rs_semileptonic_wilsons / rs_ew_couplings:
    use rigorous input
else:
    keep current proxy or return missing-extra note
```

| Constraint family | Current proxy | New rigorous input | Formula | Status after design |
|---|---|---|---|---|
| Z-pole T010-T012 | `zbb_coupling_shift_proxy` from quark KK-gluon overlaps | `rs_ew_couplings.z_delta_g_L/R_d`, `z_delta_g_L/R_u` | `R_b,A_b,A_FB` from existing `zpole_shifted_couplings(SM, delta_g)` | Fully rigorous only if `includes_fermion_kk_mixing=True`; gauge-only is partial for classic Zbb. |
| Z-LFV T015-T017 | caller lepton overlap proxy mapped to `(m_Z/M_KK)^2 overlap` | `z_delta_g_L/R_e[i,j]` | Existing `z_lfv_branching_fraction_from_couplings(delta_g_L, delta_g_R)` | Fully rigorous for tree-level LFV Z once lepton rotations/profiles are present. |
| FCNC-Z T014 | down-overlap proxy | `z_delta_g_L/R_d[i,j]` | Existing FCNC Z width using actual off-diagonal light-Z couplings | Fully rigorous. |
| Higgs-LFV T018-T020 | explicit off-diagonal Higgs Yukawa proxy | `rs_higgs_yukawas.Y_h_mass[i,j]` | `BR = m_h(|Y_ij|^2+|Y_ji|^2)/(8 pi Gamma_h)` | Partial until fermion-KK Higgs mass-matrix block exists. |
| `b -> s l l`, `b -> d l l` | one Z-like boson with SM lepton charges | `rs_semileptonic_wilsons.b_to_s_ll`, `.b_to_d_ll` | C9/C10 formula above, then existing rare-B observable cores | Fully rigorous for vector/axial tree EW; scalar, nonlocal charm, and full angular/global-fit effects remain partial. |
| `s -> d l l` | one Z-like `Y_NP` proxy | `rs_semileptonic_wilsons.s_to_d_ll` | `Y_NP` or y7V/y7A from neutral contact coefficients | Fully rigorous for direct short-distance vector/axial; ChPT interference/CPC long-distance remains partial. |
| `c -> u l l` | one Z-like C9/C10 proxy | `rs_semileptonic_wilsons.c_to_u_ll` | C9/C10 formula with `lambda_b` | Fully rigorous for short-distance vector/axial; resonance windows remain partial. |
| `b -> s nu nu` | one Z-like `X_NP` proxy | `rs_semileptonic_wilsons.b_to_s_nunu` | `X_NP_L/R = C_L/R / g_SM^2` | Fully rigorous for tree EW. |
| Lepton dipoles L001/L007/L008 | `Y_N_bar Y_N_bar^dagger` NDA proxy | `rs_dipole_wilsons.lepton_A_L/R` | Existing BR formula from explicit chiral dipole amplitudes | Partial until one-loop RS dipole matching exists. |
| `b -> s gamma` and exclusive radiative B | overlap C7/C8 proxy plus LL running | `rs_dipole_wilsons.b_to_s_gamma` | Existing C7/C8 running and rate, fed by real loop Wilsons | Partial until one-loop RS C7/C8 matching exists. |
| Mu-e conversion L003-L005 | explicit low-energy scalar/vector/dipole proxies | vector coefficients from `rs_ew_couplings`, dipoles from `rs_dipole_wilsons`, scalar from Higgs/radion block | KKO formula already implemented | Partial: vector tree piece rigorous; dipole/scalar/interference remain partial. |
| LFV three-body L002/L009/L010/L023 | lepton-overlap Z proxy plus explicit boxes | `z_delta_g_e` and heavy neutral contacts; optional boxes | `G_AB` formula above | Partial: Z/contact tree piece rigorous; boxes and dipole phases remain partial. |
| EDMs E001/E002/E004/E006-E009 | explicit low-energy CP-odd proxy or stub | `rs_dipole_wilsons.edm`, qEDM/qCEDM/Weinberg if loop matched | Observable converters already exist | Still NEEDS-HUMAN without loop and hadronic/nuclear matrix elements. |
| Charged-current EW002 | no RS shift | `rs_charged_current`, `delta_G_F_over_G_F` | corrected first-row `Vij_eff` sum | Fully rigorous for tree charged-current shifts. |
| Charged-current EW003 | no RS shift, data tension only | `rs_charged_current` | corrected inclusive/exclusive `Vcb,Vub` amplitudes | Partial because inclusive/exclusive theory covariance is not an RS-coupling issue. |
| Charged-current K018 | no RS shift | `rs_charged_current`, `delta_G_F_over_G_F` | `|V_us|_eff = |V_us| * |1 + delta_CC(K_l3)|` | Fully rigorous for tree charged-current shifts. |

### Named-scope count

Counting the process IDs explicitly in the request's affected families, not
the full catalog:

- Fully rigorous RS tree matching after the neutral-current and charged-current
  blocks: 26 process IDs.
- Partial after this design because they also need dipoles, Higgs/radion,
  boxes, long-distance windows, or inclusive/exclusive covariance: 17 process
  IDs.
- Still NEEDS-HUMAN after this design because loop plus hadronic/nuclear
  matrix elements are missing: 7 EDM process IDs.

This count is intentionally about RS new-physics matching. Some constraints
with fully rigorous short-distance matching still carry honest
observable-side long-distance diagnostics.

## 5. Phasing and risk

### Recommended implementation order

1. Add derivation notes for gauge profiles, gauge-KK overlaps, vector
   diagonalization, and Z/W coupling normalization.
2. Implement `RSEWSpectrum`, exact gauge roots, gauge profile normalization,
   and `Omega_n(c)` table caching. Validate KK root spacing first.
3. Implement light-Z and heavy-neutral contact couplings for quarks only.
   Rewire T010-T014 and same-flavor rare B/K/charm C9/C10 inputs.
4. Fill `lepton_mass_basis_couplings` and charged-lepton/neutrino neutral
   currents. Rewire Z-LFV, LFV rare K/charm, and b->s nu nu.
5. Implement charged-current W/W' and `G_F` shifts. Rewire EW002, K018, then
   EW003 with a conservative covariance note.
6. Add fermion-KK mass-matrix mixing for classic Zbb and Higgs LFV. Promote
   Zbb from gauge-only partial to full tree-level status.
7. Add real loop matchers for lepton dipoles, b->s gamma, and EDM Wilsons.

### What remains NEEDS-HUMAN

- Lepton and quark dipole loop functions until derived and validated.
- Neutron/atomic/mercury EDM hadronic and nuclear matrix-element choices.
- Long-distance charm and kaon interference pieces, resonance windows, and
  full rare-B angular/global-fit likelihoods.
- Model choices outside minimal RS: custodial embeddings, brane kinetic terms,
  extended gauge groups, and radion/Higgs localization.

### Top correctness risks

1. Coupling normalization and signs: the existing Z-pole and rare-decay cores
   use precise conventions. Centralizing vector diagonalization and storing
   dimensionless `delta_g` is the main mitigation.
2. Double counting: rare decays must include light-Z shifts plus direct heavy
   vector exchange once, and must not also apply the old proxy.
3. Basis rotations: lepton PMNS, charged-lepton rotations, and quark up/down
   rotations must be explicit in the builder. Universal `c_L` hides LFV; a
   non-universal scan must not accidentally rotate in the wrong basis.
4. `M_KK` convention: `kk_ew_mass_gev` must be physical root times
   `Lambda_IR`; old code sometimes used `M_KK = Lambda_IR` as bookkeeping.
5. 100M-point cost: profile integrals and diagonalization must be cached or
   reduced to validated lookup tables.

## 6. Validation strategy

### Unit-level physics checks

- KK roots: verify `solve_kk("gauge","NN")` returns first root near the known
  NN gauge value used in the repo's LFV helper, and that masses scale linearly
  with `Lambda_IR`.
- Profile normalization: numerically integrate each `w_0(c,z)` to 1 and check
  endpoint values reproduce `f_IR` and `f_UV`.
- Gauge universality: for universal c values and identity rotations,
  off-diagonal neutral currents vanish to numerical precision.
- Decoupling: all `delta_g` and contact Wilsons scale as `v^2/M_KK^2` or
  `1/M_KK^2` and go to zero at large `Lambda_IR`.

### Known RS checks

- Zbb sign and magnitude: compare gauge-only and gauge-plus-fermion-mixing
  `delta_g_L^b`, `delta_g_R^b` against standard minimal-RS benchmark results
  at `M_KK = 3 TeV`, using the same `c_Q3`, `c_b`, `c_t` and Yukawa input.
- KK mass spacing: compare first few gauge masses with Bessel-root expectations
  and with the `GAUGE_KK_ROOT_NN` convention currently documented in
  `flavorConstraints/muToEGamma.py`.
- Lepton LFV null test: in the current universal `c_L` and diagonal
  charged-lepton scan, tree-level `Z e mu` must vanish unless non-universal
  lepton inputs or fermion-KK mixing are enabled.

### Catalog integration checks

- SM sides unchanged: with RS extras absent or all new couplings set to zero,
  every affected adapter reproduces its current SM-limit prediction.
- Proxy replacement tests: create one synthetic point where old proxy and new
  contact coefficient are numerically set equal; the adapter result should be
  identical. Then switch to a real RS point and assert the diagnostics no
  longer contain `NEEDS-HUMAN-PHYSICS` for tree-level matching.
- Complex discipline: Wilson coefficients appear in typed extras and
  diagnostics; `ConstraintResult` scalar fields remain real.
- Regression sweep: run a small grid over `Lambda_IR` and representative c
  values and check monotonic decoupling of Z-pole, FCNC-Z, and semileptonic
  ratios.
