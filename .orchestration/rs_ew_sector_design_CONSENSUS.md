# RS Electroweak-Coupling Sector - Consensus Design

Status: consensus design only. No production code changes are implied here.

Inputs reconciled:

- Design A: `.orchestration/rs_ew_sector_design.md`
- Design B: `.orchestration/runs/G1-DESIGN/codex_design.md`
- Repo derivations and code contracts: `derivations/conventions.tex`,
  `derivations/zero_modes/README.md`, `derivations/kk_modes/README.md`,
  `derivations/kk_modes/bessel_equation.tex`,
  `derivations/flavor/README.md`,
  `flavor_catalog_constraints/base.py`,
  `flavor_catalog_constraints/point_builder.py`,
  `quarkConstraints/zpole.py`, and the rare-decay/dipole adapter cores.

## Consensus Position

The two designs agree on the central physics object: RS electroweak proxy
matching is currently standing in for non-universal gauge overlap effects after
electroweak symmetry breaking. The leading compact form is

```text
delta g_A^f = s_Z * g_A^SM(f) * (m_Z^2 / M_KK^2)
              * [U_A^dagger diag(a(c_Ai, epsilon) - a_ref) U_A]_{ij}
              + delta g_A,fermion-KK^f
```

in the existing Z-pole convention

```text
L_Z = g_Z Z_mu fbar gamma^mu (g_L P_L + g_R P_R) f,
delta g_A^f = g_A^RS - g_A^SM       # dimensionless, additive
g_L^SM = T3 - Q sin^2 theta_W
g_R^SM = -Q sin^2 theta_W
g_Z = g / cos theta_W
```

where `s_Z` is fixed by the vector mass-matrix diagonalization and by the Zbb
benchmark, not by hand. The sign check is: for an IR-composite `b_R` relative
to light down singlets, the implemented minimal-RS `delta g_R^b` must be
negative and of order `10^-3` at `M_KK ~ 3 TeV`.

Where A and B differ in form, the consensus picks B's mode-by-mode
spectrum/coupling construction as the production formula. The first
implementation of `a(c)-a_ref` is the numerical overlap/vector-diagonalization
limit, cached for scans. Reason: `derivations/conventions.tex` and
`derivations/kk_modes/bessel_equation.tex` already pin the geometry, gauge NN
quantization, Bessel roots, and normalization structure, while
`derivations/zero_modes/README.md` and `derivations/kk_modes/README.md` still
mark full zero-mode profiles and KK profile normalization as needed. The repo
does not yet contain a completed `gauge_kk_coupling.tex` derivation for the
closed-form prefactor in A. Therefore closed-form `a(c)` prefactors are
deferred until a zero-mode/gauge-profile normalization derivation is added to
`derivations/` and checked against the numerical overlap builder.

## 1. Physics Inventory

### Geometry and Zero-Mode Inputs

Pinned by `derivations/conventions.tex` and `warpConfig/baseParams.py`:

```text
epsilon = Lambda_IR / k = z_h / z_v
z_h = 1/k
z_v = 1/Lambda_IR
L = log(k / Lambda_IR)
c = M_5 / k
```

Pinned by `derivations/conventions.tex`,
`derivations/zero_modes/README.md`, and `warpConfig/wavefuncs.py`:

```text
f_IR(c)^2 = (1/2 - c) / (1 - epsilon^(1 - 2c))
f_UV(c)^2 = (1/2 - c) / (epsilon^(2c - 1) - 1)
```

`f_IR` and `f_UV` are implemented. Full normalized `g_0(c,z)` profiles are not
yet implemented and must be added before numerical overlap matching is called
production.

### Electroweak KK Spectrum

Use the existing gauge NN Bessel problem:

```text
m_n = x_n Lambda_IR
J_0(x_n) Y_0(epsilon x_n) - J_0(epsilon x_n) Y_0(x_n) = 0
```

implemented as:

```python
masses, extras = solve_kk("gauge", "NN", geometry, n_roots=N, exact=True)
```

`kk_ew_mass_gev` is the first physical electroweak gauge KK mass
`masses[0]`, not the bookkeeping `Lambda_IR`. In the small-epsilon limit,
`x_1 ~ 2.405` and the first tower ratios follow the usual NN gauge Bessel
spacing. This is grounded in `derivations/conventions.tex`,
`derivations/kk_modes/bessel_equation.tex`, and `solvers/bessel.py`.

Production builder:

```text
M_charged^2[m,n] = m_KK,m^2 delta_mn
                   + (v^2/4) g_2^2 chi_m(z_v) chi_n(z_v)

M_neutral^2[(X,m),(Y,n)] =
    m_X,m^2 delta_XY delta_mn
    + (v^2/4) [[g_2^2, -g_2 g_Y], [-g_2 g_Y, g_Y^2]]_XY
      chi_m(z_v) chi_n(z_v)
```

Diagonalizing these matrices identifies photon, light `W`, light `Z`, and the
heavy `W'`, `Z'`, `gamma'` modes. This is more general than treating all EW
KK vectors as one mass and is the agreed production path.

### Gauge-KK Overlap Form Factors

Numerical production definition:

```text
Omega_n(c) = integral_{z_h}^{z_v} dz w_0(c,z) chi_n(z)
integral dz w_0(c,z) = 1
```

where `chi_n(z)` is the normalized gauge profile and `w_0(c,z)` is the
normalized fermion zero-mode probability density. Endpoint checks must
reproduce `f_IR(c)` and `f_UV(c)`.

Cached leading form, obtained from the numerical overlap builder first:

```text
a(c, epsilon) = KK-sum / first-mode-reduced overlap form factor
delta g_A^f = s_Z * g_A^SM(f) * (m_Z^2/M_KK^2) * [a(c_A)-a_ref]
```

`a_ref` is the universal piece absorbed into measured `g_Z`, `m_Z`, and
`G_F`; it must be subtracted before storing any Z-pole shift. A closed-form
`a(c)` is not accepted for production until it is derived from the same
normalized overlap integral in `derivations/` and checked against Tier-B
numerical quadrature. The exact closed-form prefactor is not pinned in current
`derivations/`.

### Interaction-Basis and Mass-Basis Z Couplings

For each chiral species:

```text
G_X,fA,n^mass = g_X q_X(f_A) U_A^dagger diag(Omega_n(c_Ai)) U_A
```

The light-Z coupling in the existing `zpole` convention is:

```text
g_Z^RS(f_A) = (1/g_Z) * sum_{X,n} V_Z[X,n] G_X,fA,n^mass
delta g_Z,fA = g_Z^RS(f_A) - g_A^SM(f_A) I
```

Required mass-basis matrices:

```text
z_delta_g_L_u, z_delta_g_R_u
z_delta_g_L_d, z_delta_g_R_d
z_delta_g_L_e, z_delta_g_R_e
z_delta_g_L_nu
```

Left quark doublets use `c_Q` but different rotations:

```text
u_L: U_L_u^dagger diag(...) U_L_u
d_L: U_L_d^dagger diag(...) U_L_d
```

Charged leptons use `U_e_L`, `U_e_R`. Active neutrino Z couplings come from
the lepton doublet profile `c_L` and the PMNS/active-neutrino basis. The
right-handed neutrino `c_N` is an electroweak singlet and does not get an SM
Z charge; it enters seesaw and dipole spurions only. This resolves A's
`a_neutrino = a(c_N)` proposal in favor of B's gauge-charge construction.

### Heavy Neutral Contact Terms

Rare decays should not be represented only by shifted light-Z exchange once
the heavy spectrum is built. Store the physical contact coefficients:

```text
C_AB(light Z)^{q_i q_j l_a l_b} =
    (g_Z^2 / m_Z^2)
    * [ (g_qA^SM delta_ij + delta g_Z,qA[i,j])
        (g_lB^SM delta_ab + delta g_Z,lB[a,b])
        - g_qA^SM g_lB^SM delta_ij delta_ab ]

C_AB(total) = C_AB(light Z)
              + sum_{V heavy} g_V,qA[i,j] g_V,lB[a,b] / M_V^2
```

Equivalently, the light-Z piece contains all three NP products:
`g_qA^SM delta_ij * delta g_lB[a,b]`,
`delta g_qA[i,j] * g_lB^SM delta_ab`, and
`delta g_qA[i,j] * delta g_lB[a,b]`. This full product is required for LFV
channels such as `mu -> e` conversion, `mu -> 3e`, and `tau -> 3l`, where the
quark current can be SM-diagonal while the lepton current is LFV. Units are
`GeV^-2`. The heavy sum is required for true `Z'/gamma'` contact matching and
collider/contact consistency.

### Semileptonic Wilsons

For `q_j -> q_i l_a l_b`, with `C_AB` in `GeV^-2`:

```text
C9_NP   = -pi/(sqrt(2) G_F alpha lambda_t) * (C_LL + C_LR)
C10_NP  = -pi/(sqrt(2) G_F alpha lambda_t) * (C_LR - C_LL)
C9p_NP  = -pi/(sqrt(2) G_F alpha lambda_t) * (C_RL + C_RR)
C10p_NP = -pi/(sqrt(2) G_F alpha lambda_t) * (C_RR - C_RL)
```

Use the adapter's own CKM factor and SM input bundle (`lambda_t` for
`b -> q`, `lambda_b` for `c -> u`, the kaon convention for `s -> d`). This
replaces the current "one Z-like boson at `M_KK`" proxy.

The new rare-decay path must not reuse the existing `_wilson_prefactor` or any
helper that already inserts `1/M_KK^2`, `v^2/M_KK^2`, or a single proxy-boson
scale. The `C_AB` contacts are physical `GeV^-2` coefficients and already
carry their own `1/m_Z^2` or `1/M_V^2`. Adapters should compute dimensionless
`C9_NP`, `C10_NP`, `C9p_NP`, and `C10p_NP` with the formulas above, then pass
those additive NP Wilsons into the existing `rare_*` cores at the same point
where proxy Wilsons are currently consumed. There is no second
`1/M_KK^2` multiplication in the rare core.

For `b -> s nu nubar` and `s -> d nu nubar`:

```text
X_NP_L = C_LL^{q_i q_j nu nu} / g_SM^2
X_NP_R = C_RL^{q_i q_j nu nu} / g_SM^2
g_SM^2 = 4 G_F^2 M_W^2 / (2 pi^2)
```

This matches the existing `rare_b_nunu.py` and `rare_kaon_snd.py`
normalization. As above, `X_NP_L/R` is an additive dimensionless short-distance
function passed to the existing `rare_*_nunu` cores without `_wilson_prefactor`
or any additional `1/M_KK^2`.

### Charged Current and G_F Shifts

Charged-vector diagonalization gives:

```text
L_W = g_2/sqrt(2) W^+ ubar gamma^mu [(V_CKM + delta g_W_ud_L) P_L
                                      + delta g_W_ud_R P_R] d

delta g_W_ud_L = sqrt(2)/g_2 * g_W^RS(u_L d_L) - V_CKM
delta g_W_ud_R = sqrt(2)/g_2 * g_W^RS(u_R d_R)
```

Minimal SU(2)_L has no right-handed W coupling before fermion-KK/custodial
extensions, so `delta_g_W_ud_R` is zero unless such a block is explicitly
enabled. Muon decay defines:

```text
G_F^meas = G_F^0 * (1 + delta_G_F_over_G_F)
```

where the shift includes light-W coupling shifts and direct `W'` exchange.
Semileptonic extraction shifts use:

```text
A(d_j -> u_i l_a nu_a) / A_SM =
    1
    + delta_g_W_ud_L[i,j]/V_ij
    + delta_g_W_lnu_L[a,a]
    + C_Wprime_LL[i,j,a,a] / C_SM
    - delta_G_F_over_G_F
    + optional right-handed/scalar terms
```

### Lepton/Neutrino Yukawa and Profile Inputs

Current repo path: `yukawa.compute_all_yukawas()` accepts scalar `c_L`,
scalar `c_N`, vector `c_E[3]`, and `M_N`; returns `Y_E`, `Y_E_bar`, `Y_N`,
`Y_N_bar`, unbarred `Y_N_matrix`, `f_L`, `f_E`, `f_N`, `f_N_UV`, and PMNS-built
neutrino Yukawas. Contract for the EW builder: the returned `Y_N_bar` is the
legacy vector output, while `Y_N_matrix` is the unbarred `(3,3)` PMNS-basis
matrix. The builder must define and store
`Y_N_bar_matrix = 2*k*Y_N_matrix`, store the PMNS matrix used for that
construction, and store identity `U_e_L`, `U_e_R` in the current
charged-lepton mass-basis mode. This is grounded in
`derivations/conventions.tex` and `derivations/flavor/README.md`.

Consensus builder must store identity charged-lepton rotations in the current
diagonal mode and allow future non-diagonal `Y_E_bar` inputs to provide
nontrivial `U_e_L`, `U_e_R`.

### Dipoles

The tree-level EW sector can provide profile and spurion inputs:

```text
lepton_dipole_spurion = Y_N_bar_matrix Y_N_bar_matrix^dagger
quark_dipole_profile_data = c_Q,c_u,c_d,U_L/R,Y_u,Y_d,KK masses
```

It must not fabricate rigorous `C7`, `C8`, lepton dipole amplitudes, EDMs, or
Weinberg-operator coefficients. `derivations/flavor/README.md` marks
`mu_to_e_gamma.tex` as needed, and the finite KK loop coefficients are not
present. Therefore dipoles are partial until a loop matcher exists.

### Higgs LFV

Tree formula once a fermion-KK/Higgs block exists:

```text
Y_h^mass = U_L^dagger (d M_e(v) / d v) U_R
BR(h -> l_i l_j) = m_h (|Y_ij|^2 + |Y_ji|^2) / (8 pi Gamma_h)
```

With the current diagonal zero-mode charged-lepton builder, this is diagonal.
Off-diagonal Higgs LFV needs fermion-KK mixing and Higgs localization
assumptions, so it is partial.

### Oblique Parameters

Current proxy:

```text
Delta S ~ c_S v^2/M_KK^2
Delta T ~ pi L/(2 cos^2 theta_W) v^2/M_KK^2
```

Consensus replacement is the leading KK gauge-sum result using
`rs_ew_spectrum` and vector profiles. The leading `v^2/M_KK^2` piece becomes
RS-derived, but custodial embeddings and brane kinetic terms remain
model-dependent diagnostics. The default is documented `minimal-RS` with an
explicit diagnostic flag; custodial protection, extended gauge groups, and
brane kinetic terms are human model choices and must not be silently enabled.

## 2. ParameterPoint Schema

`ParameterPoint` stays frozen. `make_point` remains fail-loud via
`KNOWN_EXTRA_KEYS`. Complex matrices and Wilsons may live in typed extras and
in `ConstraintResult.diagnostics`; scalar fields on `ConstraintResult` remain
real only.

### Agreed Extra Keys

| Key | Type | Units | Policy |
|---|---|---:|---|
| `quark_mass_basis_couplings` | existing `QuarkMassBasisCouplings` | mixed | Keep for KK gluon and backward compatibility. |
| `kk_gluon_mass_gev` | `float` | GeV | Keep. |
| `kk_ew_mass_gev` | `float` | GeV | Reuse existing key; fill with first physical EW gauge KK root. |
| `lepton_mass_basis_couplings` | `RSLeptonMassBasisCouplings` | mixed | Reuse existing placeholder; folds in A's proposed `lepton_overlaps`. |
| `rs_ew_spectrum` | `RSEWSpectrum` | GeV | New; real spectrum and mixing metadata. |
| `rs_ew_couplings` | `RSEWMassBasisCouplings` | dimensionless, GeV^-2 | New; complex mass-basis Z and neutral-contact matrices. |
| `rs_semileptonic_wilsons` | `RSSemileptonicWilsonBundle` | dimensionless | New; complex Wilsons by transition. |
| `rs_charged_current` | `RSChargedCurrentCouplings` | dimensionless, GeV^-2 | New; complex W shifts and contacts. |
| `rs_dipole_wilsons` | `RSDipoleWilsonBundle` | dimensionless or GeV^-1 | Optional; absent unless a real loop matcher ran. |
| `rs_higgs_yukawas` | `RSHiggsYukawaCouplings` | dimensionless | Optional; absent or partial until fermion-KK/Higgs block exists. |

The new required extras `rs_ew_spectrum`, `rs_ew_couplings`,
`rs_semileptonic_wilsons`, and `rs_charged_current`, plus optional extras
`rs_dipole_wilsons` and `rs_higgs_yukawas` whenever they are emitted, MUST be
added to `KNOWN_EXTRA_KEYS`. Otherwise `make_point` correctly raises `KeyError`
and the new builder must fail loud.

No separate `lepton_overlaps` extra is added. It would duplicate
`lepton_mass_basis_couplings` and increase fail-loud schema churn.

### Dataclass Sketches

Typed keys are required so tests can assert full family coverage instead of
checking ad hoc string names:

```python
FamilyIndex = int  # must be 0, 1, or 2
QuarkSector = Literal["u", "d"]
LeptonSector = Literal["e", "nu"]
ChiralityKey = Literal["LL", "LR", "RL", "RR"]  # quark chirality first
QuarkTransitionKey = tuple[QuarkSector, FamilyIndex, FamilyIndex]
# ("d", 1, 2) means b -> s; ("u", 0, 1) means c -> u.
LeptonPairKey = tuple[LeptonSector, FamilyIndex, FamilyIndex]
# ("e", a, b) means l_b -> l_a current; ("nu", alpha, beta) is active nu.
NeutralContactKey = tuple[QuarkTransitionKey, LeptonPairKey, ChiralityKey]
RareWilsonName = Literal["C9", "C10", "C9p", "C10p", "X_L", "X_R"]
RareWilsonKey = tuple[QuarkTransitionKey, LeptonPairKey, RareWilsonName]
LFVLLQQKey = tuple[QuarkTransitionKey, LeptonPairKey, ChiralityKey]
```

For neutrino final states, only lepton-left entries are physical in
minimal-RS: `LL` and `RL` contacts, or `X_L` and `X_R` Wilsons. Missing
right-handed-neutrino contact keys must be absent by schema, not silently set
to zero. Family-coverage tests should require all `3x3` charged-lepton pairs,
all `3x3` active-neutrino pairs where relevant, and all supported quark
transitions for every adapter family.

```python
@dataclass(frozen=True)
class RSEWSpectrum:
    model_label: str
    k_gev: float
    lambda_ir_gev: float
    epsilon: float
    warp_log: float
    n_gauge_modes: int
    gauge_roots_x: np.ndarray        # (N,), real
    gauge_masses_gev: np.ndarray     # (N,), real
    m_w_gev: float
    m_z_gev: float
    m_wprime_gev: np.ndarray         # (N,), real
    m_zprime_gev: np.ndarray         # (N,), real
    m_gammaprime_gev: np.ndarray     # (N,), real
    charged_vector_mixing: np.ndarray
    neutral_vector_mixing: np.ndarray
    exact_bessel: bool
    profile_normalization: str
```

```python
@dataclass(frozen=True)
class RSLeptonMassBasisCouplings:
    model_label: str
    c_L: np.ndarray                  # (3,), real; scalar inputs broadcast
    c_E: np.ndarray                  # (3,), real
    c_N: np.ndarray                  # (3,), real; scalar inputs broadcast
    f_L_IR: np.ndarray               # (3,), real
    f_E_IR: np.ndarray               # (3,), real
    f_N_IR: np.ndarray               # (3,), real
    f_N_UV: np.ndarray               # (3,), real
    Y_E_bar: np.ndarray              # (3,3), complex
    Y_N_bar: np.ndarray              # (3,), complex, legacy vector output
    Y_N_matrix: np.ndarray           # (3,3), complex, unbarred PMNS basis
    Y_N_bar_matrix: np.ndarray       # (3,3), complex = 2*k*Y_N_matrix
    U_e_L: np.ndarray                # (3,3), complex
    U_e_R: np.ndarray                # (3,3), complex
    U_nu_L: np.ndarray               # (3,3), complex
    pmns: np.ndarray                 # (3,3), complex
    lfv_dipole_spurion: np.ndarray   # Y_N_bar_matrix Y_N_bar_matrix^dagger
    m_kk_ew_gev: float
    matching_status: Mapping[str, str]
```

```python
@dataclass(frozen=True)
class RSEWMassBasisCouplings:
    model_label: str
    spectrum: RSEWSpectrum
    convention: str                  # "zpole_dimensionless_additive"
    a_ref: float
    form_factors: Mapping[str, np.ndarray]  # real arrays by species
    z_delta_g_L_u: np.ndarray        # (3,3), complex, dimensionless
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
    neutral_contacts: Mapping[NeutralContactKey, complex]  # GeV^-2
    includes_heavy_neutral_exchange: bool
    includes_light_z_shift: bool
    includes_fermion_kk_mixing: bool
    diagnostics: Mapping[str, Any]
```

```python
@dataclass(frozen=True)
class RSChargedCurrentCouplings:
    model_label: str
    delta_g_W_ud_L: np.ndarray       # (3,3), complex, dimensionless
    delta_g_W_ud_R: np.ndarray
    delta_g_W_lnu_L: np.ndarray
    delta_g_W_lnu_R: np.ndarray
    delta_G_F_over_G_F: float
    charged_contact_LL: Mapping[tuple[FamilyIndex, FamilyIndex, FamilyIndex, FamilyIndex], complex]
    # key is (u_i, d_j, charged_l_a, nu_alpha), GeV^-2
    diagnostics: Mapping[str, Any]
```

```python
@dataclass(frozen=True)
class RSSemileptonicWilsonBundle:
    model_label: str
    b_to_s_ll: Mapping[RareWilsonKey, complex]
    b_to_d_ll: Mapping[RareWilsonKey, complex]
    s_to_d_ll: Mapping[RareWilsonKey, complex]
    c_to_u_ll: Mapping[RareWilsonKey, complex]
    b_to_s_nunu: Mapping[RareWilsonKey, complex]
    s_to_d_nunu: Mapping[RareWilsonKey, complex]
    lfv_llqq: Mapping[LFVLLQQKey, complex]
    matching_scale_gev: float
    diagnostics: Mapping[str, Any]
```

`RSDipoleWilsonBundle` is absent by default. If present, it must state
`matching_status="loop-matched"` and contain actual high- and low-scale
Wilson coefficients, not NDA spurions. `RSHiggsYukawaCouplings` is likewise
optional and must state whether fermion-KK mixing and Higgs-profile effects
were included.

## 3. Point Builder Logic and Sweep Inputs

Add a new construction path in `flavor_catalog_constraints.point_builder`:

```python
def build_from_rs_ew_inputs(
    *,
    quark_fit_result: QuarkFitResult | None = None,
    quark_couplings: QuarkMassBasisCouplings | None = None,
    lepton_inputs: RSLeptonSweepInputs,
    Lambda_IR: float,
    k: float = MPL,
    n_gauge_modes: int = 1,
    ew_model: Literal["minimal_rs", "custodial_protected", "custom"] = "minimal_rs",
    include_fermion_kk_mixing: bool = False,
    include_loop_dipoles: bool = False,
    raw: Any = None,
) -> ParameterPoint:
    ...
```

`ew_model="minimal_rs"` is the documented default and must emit diagnostics
stating that no custodial protection or brane kinetic terms were included.
`custodial_protected` and nonzero brane kinetic terms are explicit human model
choices, not automatic fallbacks.

Build path:

1. Compute geometry with `get_warp_params(k, Lambda_IR)`.
2. Compute gauge NN roots with `solve_kk("gauge", "NN", geometry, n_roots)`.
3. Fill `kk_ew_mass_gev = gauge_masses_gev[0]`.
4. Build `RSEWSpectrum`, including charged/neutral vector mass matrices and
   mixing metadata when profile normalization is available.
5. Build or accept quark data. Exact EW overlaps need `c_Q`, `c_u`, `c_d` and
   `U_L_u`, `U_L_d`, `U_R_u`, `U_R_d`; if only the old
   `QuarkMassBasisCouplings` object is supplied, keep it for compatibility but
   mark exact EW reconstruction unavailable unless these values are attached.
6. Build lepton data with `compute_all_yukawas()` in current-compatible mode:
   scalar `c_L`, vector `c_E[3]`, scalar `c_N`, and scalar `M_N`. Broadcast
   scalar `c_L`, `c_N`, `M_N` to shape `(3,)` in the stored dataclass. Store
   returned vector `Y_N_bar`, unbarred `Y_N_matrix`, derived
   `Y_N_bar_matrix = 2*k*Y_N_matrix`, PMNS, and identity `U_e_L`, `U_e_R`
   unless non-diagonal charged-lepton Yukawas are supplied.
7. Compute `Omega_n(c)` and/or numerical-overlap `a(c)` for every needed c
   value from cached profile tables. Do not use a closed-form `a(c)` prefactor
   until a derivation file pins it.
8. Build `RSEWMassBasisCouplings`, including light-Z shifts and full
   product-minus-SM neutral contacts plus heavy neutral exchange.
9. Build `RSSemileptonicWilsonBundle` from the contact coefficients, using
   each adapter family's normalization and no `_wilson_prefactor`.
10. Build `RSChargedCurrentCouplings`.
11. Add `rs_dipole_wilsons` and `rs_higgs_yukawas` only if the corresponding
    matchers actually ran.
12. Return `make_point(...)`; every emitted extra must be registered in
    `KNOWN_EXTRA_KEYS`, and unknown extras still raise `KeyError`.

New sweep inputs the current quark-only builder lacks:

| Input | Type | Why |
|---|---|---|
| `Lambda_IR`, `k` | floats, GeV | Geometry and physical EW KK mass. |
| `n_gauge_modes` | int | KK truncation and validation. |
| `c_L` or `c_L[3]` | real | Lepton doublet Z/W profiles and active neutrinos. |
| `c_E[3]` | real | Charged-lepton singlet profiles. |
| `c_N` or `c_N[3]` | real | Sterile-neutrino profiles for seesaw/dipole spurions. |
| `M_N` or `M_N[3]` | GeV | UV Majorana mass scale(s). |
| `lightest_nu_mass`, `ordering` | scalar/string | Neutrino spectrum for seesaw inversion. |
| `majorana_alpha`, `majorana_beta` | radians | PMNS/seesaw phases. |
| optional `Y_E_bar[3,3]` | complex | Future non-diagonal charged-lepton scans. |
| optional `Y_N_matrix[3,3]` or `Y_N_bar_matrix[3,3]` | complex | Future direct neutrino-Yukawa scans; builder must not confuse these with vector `Y_N_bar`. |
| quark `c_Q,c_u,c_d` plus rotations | arrays | Exact EW overlaps; old KK-gluon matrices alone are insufficient. |
| `g_2`, `g_Y`, `sin2_theta_w`, `v` policy | real | Coupling normalization and vector diagonalization. |
| `ew_model` | enum | Defaults to documented `minimal_rs`; custodial/BKT variants require explicit human choice. |
| `include_fermion_kk_mixing` | bool | Zbb/Higgs LFV promotion. |
| `include_loop_dipoles` | bool | Whether rigorous dipole Wilsons are present. |

For high-volume scans, c-grid overlap values and vector diagonalizations must
be cached; no per-point quadrature in a 100M-point path.

## 4. Matching Map

Status labels:

- FULL: rigorous RS tree-level matching in the agreed scope.
- PARTIAL: rigorous tree/profile piece exists but another loop, scalar,
  box, long-distance, covariance, or model-choice ingredient remains.
- HUMAN: still needs human theory inputs after this EW sector.

| Family / adapters | Current proxy | New input | Formula | Status |
|---|---|---|---|---|
| Z-pole `T010-T012` (`zpole`, `zpole_charm`) | `zbb/zcc` overlap proxy and `kk_ew_mass_gev` fallback | `rs_ew_couplings.z_delta_g_*`, `kk_ew_mass_gev` | Vector-diagonalized `shifted_couplings(sm, delta_g)` path; old proxy helpers bypassed | FULL only after fermion-KK Zbb block; gauge-only is PARTIAL for classic Zbb |
| FCNC-Z `T014` | `ZPOLE_DOWN_FCNC_PROXY_V1` | `z_delta_g_L/R_d[i,j]` | Off-diagonal mass-basis light-Z couplings | FULL |
| Z-LFV `T015-T017` | caller lepton `delta_g` proxy | `z_delta_g_L/R_e[i,j]` | Existing `z_lfv_branching_fraction_from_couplings` | FULL for tree LFV Z |
| Higgs-LFV `T018-T020` | off-diagonal Higgs Yukawa proxy | `rs_higgs_yukawas.Y_h_mass` | `m_h(|Y_ij|^2+|Y_ji|^2)/(8 pi Gamma_h)` | PARTIAL until fermion-KK/Higgs block |
| `b -> s ll`, `b -> d ll` (`B005-B008`, `B015-B019`, `B021`) | one Z-like boson and KK-gluon overlap | `rs_semileptonic_wilsons.b_to_s_ll`, `.b_to_d_ll` | C9/C10 from full neutral contacts; no `_wilson_prefactor` | FULL for vector/axial tree EW; exclusive form factors, C7, scalar, nonlocal charm remain PARTIAL diagnostics |
| `b -> s nu nubar` (`B022-B023`) | one Z-like `X_NP` proxy | `rs_semileptonic_wilsons.b_to_s_nunu` | `X_NP_L/R = C_L/R / g_SM^2` from full contacts; no extra `1/M_KK^2` | FULL for tree EW |
| Charged-current beauty (`B009`, `B025-B026`) | charged-current stress/LFU proxy | `rs_charged_current` | corrected W/W' amplitude and `delta_G_F` formula | `B009`, `B025` FULL for tree amplitude; `B026` PARTIAL due inclusive/exclusive covariance |
| `s -> d nu nubar` (`K004-K005`) | one Z-like `X_NP` proxy | `rs_semileptonic_wilsons.s_to_d_nunu` | `X_NP_L/R = C_L/R / g_SM^2` from full contacts; no extra `1/M_KK^2` | FULL for tree EW |
| `s -> d ll` (`K006`, `K008-K010`, `K012-K013`) | one Z-like `Y_NP` or `y7V/y7A` proxy | `rs_semileptonic_wilsons.s_to_d_ll` | Buras/ISU kaon mapping from full neutral contacts; no `_wilson_prefactor` | FULL for direct short-distance vector/axial; long-distance/interference diagnostics remain PARTIAL where present |
| kaon LFV (`K019-K021`) | quark proxy plus lepton LFV proxy | `lfv_llqq` Wilsons from neutral contacts | Full light-Z product plus heavy neutral contacts into LFV `Y_NP`/form-factor mapping | FULL for tree short-distance; form-factor/q2 simplifications remain diagnostics |
| charged-current kaon (`K017-K018`) | `kk_ew_mass_gev` charged-current proxy | `rs_charged_current`, `delta_G_F` | corrected `K_l2/K_l3` amplitudes | FULL for tree amplitude |
| rare charm (`C004-C008`) | one Z-like c-u and LFV proxy | `rs_semileptonic_wilsons.c_to_u_ll`, `lfv_llqq` | C9/C10 from full contacts with rare-charm CKM convention; no `_wilson_prefactor` | FULL for short-distance vector/axial; resonance/long-distance windows PARTIAL |
| lepton dipoles (`L001`, `L007-L008`) | `Y_N_bar Y_N_bar^dagger` NDA proxy | `rs_dipole_wilsons` only if loop matched; otherwise lepton spurion diagnostics | chiral dipole BR formulas | PARTIAL until one-loop RS dipole matching |
| LFV three-body (`L002`, `L009-L010`, `L023`) | Z-penguin/box proxy | `z_delta_g_e` plus heavy neutral contacts | `G_AB` from full light-Z product plus heavy vectors | PARTIAL: tree vector/contact rigorous; boxes/dipole phases missing |
| mu-e conversion (`L003-L005`) | explicit vector/scalar/dipole proxy | vector from `rs_ew_couplings`, optional dipoles/scalars | nucleon vector combinations from full `llqq` contacts, including `g_q^SM * delta g_l^LFV` | PARTIAL: vector tree rigorous; scalar/dipole/interference missing |
| muonium conversion (`L006`) | lepton four-point proxy | neutral lepton contacts | full light-Z product plus heavy neutral lepton-current contact | PARTIAL until box/contact convention is fully pinned |
| `b -> s gamma`, radiative B (`B011-B014`) | KK-gluon overlap as C7/C8 proxy | `rs_dipole_wilsons.b_to_s_gamma` only if loop matched | existing C7/C8 running and rate | PARTIAL until one-loop C7/C8 matcher |
| top FCNC EW/dipole/Higgs (`T001-T008`) | Z/photon/gluon/Higgs proxy masses and overlaps | Z part from `rs_ew_couplings`; dipole/Higgs optional bundles | top FCNC adapters' effective-coupling formulas | PARTIAL overall; Z tree part can be rigorous, dipole/Higgs need extra blocks |
| oblique `EW001` | `Delta S,T ~ v^2/M_KK^2` proxy | `rs_ew_spectrum` gauge-sum | leading KK gauge contribution | PARTIAL: leading minimal-RS piece only |
| `EW002`, `EW003` | missing/CKM stress proxy | `rs_charged_current`, `delta_G_F` | CKM extraction shifts | `EW002` FULL for tree shifts; `EW003` PARTIAL due covariance/extraction theory |
| collider EW masses/contact (`CR005-CR007`, `CR009`, `CR012-CR014`, `kk_graviton_resonance`) | `kk_ew_mass_gev` or `M_KK` fallback | `kk_ew_mass_gev`, `rs_ew_spectrum`, neutral contacts where applicable | first physical EW gauge root; contact scale from heavy neutral vectors | FULL for mass resolution; contact/diboson couplings PARTIAL if model-specific |
| EDMs (`E001-E002`, `E004`, `E006-E009`) | explicit CP-odd proxy/stub | optional `rs_dipole_wilsons.edm`, qEDM/qCEDM/Weinberg only after loop matching | EDM observable converters | HUMAN |

### Agreed Count

Firm count for the named affected families after arbitration:

```text
FULL after the agreed tree-level scope, including fermion-KK Zbb promotion: 26
PARTIAL after the same agreed tree-level scope: 17
STILL-HUMAN EDM/matrix-element cases: 7
```

This reconciles the design-A range and design-B count as follows: before the
fermion-KK Zbb/Higgs subphase, the full count is `24` and the two Zbb-sensitive
entries remain partial. After that subphase, the agreed full count is `26`.
The partial count is fixed at `17` because B's explicit process-ID count is the
better-defined scope; A's `13-16` range was a family-level estimate.
Custodial/BKT variants, loop-normalized dipoles, and EDM promotions are not in
these FULL counts until the human-input items below are resolved.

## 5. Implementation Phasing

Each phase is independently testable.

1. Derivation pins: add derivation notes for normalized gauge profiles,
   zero-mode `g_0(c,z)`, gauge-KK overlaps, vector diagonalization, and
   Z/W coupling normalization. This must precede production matching.
2. Spectrum and overlap kernel: implement `RSEWSpectrum`,
   `kk_ew_mass_gev = x_1 Lambda_IR`, profile normalization, `Omega_n(c)`,
   and cached numerical-overlap `a(c)`.
3. Quark neutral currents: build quark mass-basis Z matrices and heavy neutral
   contacts. Rewire Z-pole/FCNC-Z and same-flavor rare B/K/charm vector
   Wilsons through the vector-diagonalized path, bypassing old proxy helpers.
4. Lepton neutral currents: fill `lepton_mass_basis_couplings`, charged-lepton
   and active-neutrino Z matrices, LFV neutral contacts, Z-LFV, LFV rare
   K/charm, `mu -> 3e`, `tau -> 3l`, `mu-e` conversion, and
   `b -> s nu nubar` with the full light-Z product.
5. Charged-current sector: build W/W' shifts, `delta_G_F_over_G_F`, and
   charged contacts. Rewire `EW002`, `K017/K018`, `B009`, `B025`, then `EW003`
   with a covariance note.
6. Fermion-KK tree mixing and Higgs block: add classic Zbb fermion-mixing
   contribution and Higgs LFV mass-matrix matching. This is what promotes the
   full count from 24 to 26.
7. Loop matchers: add lepton dipoles, `b -> s gamma` C7/C8, top dipoles, and
   EDM Wilsons. Only after this may `rs_dipole_wilsons` drive full dipole/EDM
   predictions.

## 6. Normalization and Convention Pin

This is the highest-risk item in both designs.

### Storage Convention

`rs_ew_couplings.z_delta_g_*` stores dimensionless additive shifts in the
`quarkConstraints.zpole` convention:

```text
L_Z = g_Z Z_mu fbar gamma^mu (g_L P_L + g_R P_R) f
zpole.shifted_couplings(sm, delta_g_left=stored_delta)
```

Do not multiply these stored `delta_g` matrices by `g_Z` before passing them
to `zpole`.

For rare decays, use the Wilson bundle or physical contacts. If constructing a
light-Z contact from stored `delta_g`, use the full product minus the SM
product:

```text
C_AB(light Z) =
    (g_Z^2 / m_Z^2)
    * [ (g_qA^SM delta_ij + delta_g_qA[i,j])
        (g_lB^SM delta_ab + delta_g_lB[a,b])
        - g_qA^SM g_lB^SM delta_ij delta_ab ]
```

where `g_qA^SM`, `g_lB^SM`, and stored `delta_g` values are all dimensionless
couplings in the same parentheses convention. This retains LFV
`g_q^SM * delta_g_l` terms when the quark current is flavor diagonal. Do not
feed new `delta_g` values through old proxy helpers that multiply by `g_Z/2`
again, and do not route contact-derived Wilsons through `_wilson_prefactor`.

### Required Universal-Profile Tests

There are two distinct tests:

1. Full SM-subtraction test. Force every relevant chiral representation to
   have `a(c)=a_ref`, or use a truly universal overlap across all chiral reps
   with the same universal subtraction, and use identity rotations. Then:

```text
all z_delta_g_* diagonal shifts = 0
all z_delta_g_* off-diagonal entries = 0
delta_g_W_ud_L = 0 relative to CKM
delta_G_F_over_G_F = 0 after universal subtraction
all affected Z-pole/FCNC/semileptonic adapters reproduce their existing
SM-limit or zero-NP predictions
```

2. Per-species family-universal test. Set `c_Q`, `c_u`, `c_d`, `c_L`, and
   `c_E` family-universal within each gauge-charge species, but do not require
   every species to equal `a_ref`. Then all FCNC and LFV entries vanish, while
   diagonal Z-pole and charged-current shifts may remain as universal
   per-representation shifts.

Together these tests catch double-counting of the universal overlap piece,
wrong placement of `g_Z`, and false assumptions that family universality alone
removes all diagonal Z-pole shifts.

### Required Zbb Sign/Magnitude Cross-Check

Use the new vector-diagonalized path, not the old `zbb/zcc` proxy helpers, on
a minimal non-custodial benchmark with light down singlets UV-localized and
`b_R` more IR-localized (`c_bR < c_lightR`). The implementation must return:

```text
a(c_bR) > a(c_lightR)
delta g_R^b < 0
|delta g_R^b| ~ 10^-3 for M_KK ~ 3 TeV and O(1) IR compositeness
all shifts scale as 1/M_KK^2
```

If vector diagonalization gives the opposite sign, the implementation has a
sign or convention mismatch relative to the existing validated `zpole` core.
If the magnitude is off by an order of magnitude, check whether `M_KK` used
`Lambda_IR` instead of the physical root `x_1 Lambda_IR`, whether `g_Z` was
double-counted, and whether `a_ref` was omitted.

## 7. Firm Scope After Arbitration

1. Exact closed-form `a(c)` prefactor: A proposes a closed-form expression
   proportional to `f_IR(c)^2` and a log/universal subtraction. The current
   `derivations/` tree does not pin the full `g_0(c,z)` and gauge profile
   normalization needed to verify that prefactor. Production uses the
   numerical overlap derivation first. Closed form is deferred until a new
   derivation file validates it.
2. Classic Zbb completeness: A treats the gauge-induced Z shift as sufficient
   for full Zbb rigor. B requires fermion-KK mixing and custodial/embedding
   flags for the classic `Z b_L b_L` and `Z b_R b_R` story. The consensus
   adopts B for production status: gauge-only is a rigorous gauge-sector
   contribution but not a complete classic-Zbb prediction.
3. Model-dependent custodial and brane-kinetic choices: the minimal builder can
   compute a documented `minimal_rs` answer and emit diagnostics, but
   custodial protection, extended gauge groups, and brane kinetic terms are
   not determined by current inputs. The custodial/BKT variant is a human
   decision and must not be silently selected.
4. Dipole absolute normalization: both designs agree the profile/spurion
   structure can be computed, but the finite one-loop KK coefficient is not in
   `derivations/`. Dipoles remain PARTIAL/absent and must not be promoted
   without a loop matcher.
5. Observable-side long-distance/covariance treatment: the RS short-distance
   Wilsons can be rigorous while rare-B angular fits, charm resonance windows,
   kaon interference, and inclusive/exclusive CKM covariance remain outside
   this sector. These should stay explicit diagnostics rather than being
   hidden in the EW builder.

## 8. HUMAN-INPUT ITEMS surfaced by design review

1. Custodial/BKT model choice: whether to stay with diagnostic `minimal_rs` or
   provide explicit custodial embeddings, extended gauge group, and brane
   kinetic terms.
2. Dipole-loop normalization: finite KK-loop matcher and matching convention
   for lepton dipoles, `C7/C8`, top dipoles, and related CP-odd operators.
3. EDM items: which EDM operator basis, CP-phase inputs, hadronic/nuclear
   matrix-element conventions, and promotion criteria to use once loop
   matching exists.
