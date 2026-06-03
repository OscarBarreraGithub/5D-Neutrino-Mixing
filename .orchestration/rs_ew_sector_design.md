# RS Electroweak-Coupling Sector — Implementation Design

**Status:** DESIGN ONLY (no code written). Author: RS flavor-physics architect pass.
**Goal:** Add a rigorous "RS electroweak (EW) coupling sector" so the flavor-catalog
constraints that currently consume NEW-PHYSICS *proxies* (the "G1/G2 gap") can be
re-wired to RS-derived inputs computed from the existing zero-mode/KK machinery
(`warpConfig/wavefuncs.py`, `warpConfig/baseParams.py`, `solvers/bessel.py`).

This doc is grounded in:
- `CLAUDE.md` conventions: `c = M_5/k`, `alpha = |c + 1/2|`, `f_IR^2 = (1/2-c)/(1-eps^{1-2c})`,
  `f_UV^2 = (1/2-c)/(eps^{2c-1}-1)`, seesaw + PMNS in the charged-lepton mass basis.
- `derivations/kk_modes/bessel_equation.tex`: the unified Sturm-Liouville KK problem,
  gauge BC `s=2, M_phi^2=0, alpha=nu=0`, the normalization weight `w(y)=e^{(2-s)sigma}`,
  the cylinder-function `N_n^2`, and the gauge mass formula
  `m_n ~ (n + alpha/2 - 1/4) pi k e^{-pi k R}` (for the gauge tower `x_1 ~ 2.405` exact).
- `flavor_catalog_constraints/base.py` (`ParameterPoint` frozen container, real-only
  `ConstraintResult` numeric fields, complex -> `diagnostics`), `point_builder.py`
  (`KNOWN_EXTRA_KEYS`, fail-loud `make_point`, `build_from_quark_couplings`).
- The proxy sites in `flavor_catalog_constraints/physics_adapters/` and their
  `quarkConstraints` cores. Every proxy shares one structural form (see §0).

---

## 0. The single structural pattern behind every proxy

Reading the `*_RS_*_PROXY_V1` / `*_RS_MATCHING_ASSUMPTION_V1` strings across
`quarkConstraints/` shows that **almost all** EW proxies stand in for the *same*
missing object: the RS **Z-coupling shift** induced by gauge-KK mixing. The proxies
all take the form

```
delta g_f          ~  proxy_strength * (m_Z / M_KK)^2 * (overlap_f - overlap_light)    [zpole Zbb]
delta g_ij (FCNC)  ~  proxy_strength * (m_Z / M_KK)^2 * overlap_ij                     [zpole down-FCNC, zpole_lfv]
C9/C10, C_L^nu     ~  "single Z-like boson, Delta_q = Delta_lepton = g/(2 c_W)"        [rare_b/k/charm dilepton, b->s nu nu, rare kaon]
C7/C8 (dipole)     ~  KK-gluon overlap used as flavor-overlap proxy, LL-run           [bsgamma, lepton mu->e gamma]
Delta S, Delta T   ~  c_S v^2/M_KK^2 ; pi L/(2 cW^2) v^2/M_KK^2                        [oblique_stu]
```

The rigorous replacement is a *single* well-defined computation: the **non-universal
Z-coupling shift from integrating out the KK gauge tower**, which in RS is the
classic result (Agashe-Contino-Pomarol; Casagrande et al. 0807.4937; Gherghetta TASI):

```
delta g_f / g_Z^SM  =  (m_Z^2 / M_KK^2) * [ a_f(c_f) - <a> ]                              (rigorous, leading)
```

where `a_f(c_f)` is the **gauge-fermion overlap form factor** — the warped integral of the
flavor-universal KK gauge profile against the fermion zero-mode squared — and `<a>` is its
universal (flavor-blind) piece reabsorbed into `g_Z`, `m_Z`. Crucially `a_f` is a pure
function of `c_f` and `epsilon`, computable **today** from `f_IR/f_UV` plus one extra
warped integral. This single object feeds Z-pole, FCNC-Z, b->s l l, s->d l l, c->u l l,
b->s nu nu, Z-LFV, and (with a different chirality projector) the dipole pieces.

So the design is: **build `a_f(c)` once, rotate it to the mass basis, and expose the
resulting diagonal + off-diagonal `delta g_L/R` matrices** for quarks, charged leptons,
and neutrinos. Everything else is re-using existing SM-side machinery already validated
in the constraint files.

---

## 1. Physics inventory (RS-EW quantities the constraints need)

### (a) KK gauge masses (W', Z' = first KK W/Z)
Compute with the existing solver:
```python
masses, _ = solve_kk(species="gauge", bc="NN", geometry=get_warp_params(k, Lambda_IR), n_roots=3)
M_KK_gauge = masses[0]            # = x_1 * Lambda_IR, x_1 ~ 2.4048 (J_0 first zero) in eps->0 limit
```
- `M_W1 ~ M_Z1 ~ x_1 * Lambda_IR` at leading order (the W/Z KK masses are degenerate up
  to `O(m_Z^2/M_KK^2)` electroweak-symmetry-breaking splittings — record both but treat
  the splitting as a diagnostic, not a separate solve, in G1).
- Derivation: `bessel_equation.tex` §"KK masses" with `alpha=0` (gauge, `s=2`).
- This replaces the *mass proxy* `kk_ew_mass_gev := M_KK` fallback in
  `collider_resonance.resolve_kk_ew_mass_gev` and `kk_graviton_resonance` with an actual
  first-gauge-KK root, and supplies the `M_KK` that normalizes every `(m_Z/M_KK)^2`.

### (b) Z-fermion coupling shifts delta g_L^f / delta g_R^f (the core object)
The KK gauge zero-mode is flat (`c_gauge -> 1/2` analog, `nu=0`); the first KK mode is
IR-peaked. Its overlap with a fermion zero mode of bulk mass `c` defines the form factor
`a(c)`. With the IR-overlap normalization already in the repo, the leading non-universal
shift is

```
a(c)         = warped overlap integral of the (normalized) 1st-KK gauge profile against |f0(c,y)|^2
delta g_L^i  = g_Z^SM * (m_Z^2 / M_KK^2) * ( a(c_{L_i}) - a_ref )    (left doublets: quarks & leptons)
delta g_R^i  = g_Z^SM * (m_Z^2 / M_KK^2) * ( a(c_{E_i}) - a_ref )    (right charged leptons)
```

Two computation tiers:
- **Tier-A (closed form, G1):** use the standard RS analytic form-factor. For an
  IR-localized fermion the leading overlap reduces to a simple function of `f_IR^2(c)`
  and `log(1/eps)`; for a UV-localized fermion it is exponentially suppressed. Concretely
  `a(c) ≈ f_IR^2(c) * [ (1/2)/(3 - 2c)  -  (universal log piece) ]` (the exact prefactor
  comes from the `∫ z (J_0 + b Y_0)^2 |f0|^2` integral in `bessel_equation.tex`). This is
  a one-page integral over the *existing* normalized profiles — implementable now.
- **Tier-B (numeric, G2):** evaluate the same overlap by quadrature using the
  cylinder-function gauge profile from `solve_kk(..., return_extras=True)['b']` and the
  fermion `f_n` from the same Sturm-Liouville solution. Use Tier-B as the validation
  oracle for Tier-A (see §6).

Sign convention: for `c_b < c_light` (b more IR-composite) `a(c_b) > a_light`, giving the
classic **negative** `delta g_R^b` / shift to `Z b bbar` — the sign the proxy was hand-
tuned to mimic. This MUST reproduce (§6).

Derivation references: `bessel_equation.tex` (profiles, normalization integral); add a
new `derivations/flavor/gauge_kk_coupling.tex` capturing the overlap integral (currently a
README stub).

### (c) Flavor off-diagonal Z couplings (FCNC-Z)
The shifts in (b) are diagonal in the *gauge/interaction* basis but `c`-dependent, hence
**non-universal**. Rotating to the fermion mass basis with the same unitary matrices used
by the Yukawa fit produces off-diagonal entries:
```
[delta g_L]^mass = U_L^dagger  diag( delta g_L^i )  U_L         (and likewise R with U_R)
delta g_{ij} = off-diagonal entry, i != j   (drives FCNC-Z)
```
- Quarks: `U_L_d, U_R_d, U_L_u, U_R_u` already exist in `QuarkMassBasisCouplings`
  provenance (`fit_result.U_*`); the EW form factor `a(c)` replaces the KK-gluon overlap
  `g_s^2 f^2` used there. **Same rotation, different coupling.**
- Leptons: requires the charged-lepton mass-basis rotations (new — see §3). In the
  charged-lepton mass basis (CLAUDE.md convention), neutrino mixing is the PMNS matrix.
- This replaces `ZPOLE_DOWN_FCNC_PROXY_V1`, `ZPOLE_LFV_PROXY_V1`, the e-mu Z-like
  couplings in `rare_charm_lfv_*`, and feeds the Z-penguin pieces of b->s l l / s->d l l.

### (d) Charged-current / W shifts (G_F, |V_ij| matching)
The analogous KK-W overlap shifts the effective 4-fermion charged-current strength and
introduces non-universal corrections to extracted CKM/PMNS elements and `G_F`:
```
delta G_F / G_F   ~  (m_W^2 / M_KK^2) * (a(c_{L,i}) + a(c_{L,j}) - 2 a_ref)
delta |V_ij|/|V_ij| ~ same form-factor combination
```
- Feeds EW002/EW003 (G_F / W mass / oblique) and K018-type unitarity/charged-current
  stress tests, replacing the charged-current *stress proxy* in `semileptonic_lfu`.
- Tractable in Tier-A from the same `a(c)`; flagged PARTIAL because the full custodial /
  brane-kinetic-term structure that protects `Z b_L b_L` is model-dependent (see §5).

### (e) Lepton & neutrino bulk profiles (c_L, c_E, c_N -> f's)
Needed for LFV dipoles, mu-e conversion, EDMs, Z-LFV. Already computable:
```python
eps = get_warp_params(k, Lambda_IR)['epsilon']
fL  = f_IR(c_L, eps);  fE = f_IR(c_E, eps);  fN = f_IR(c_N, eps);  fN_uv = f_UV(c_N, eps)
```
plus the charged-lepton mass-basis rotations `U_L^e, U_R^e` and the PMNS matrix from the
`yukawa` module (`compute_all_yukawas` already returns `Y_E_bar`, `Y_N_bar`, `Y_N_matrix`,
and the PMNS used to build it).

### (f) Dipole coefficients (C7, lepton dipole) — partially tractable
The RS one-loop dipole from KK fermion + KK gauge/Higgs exchange. Two pieces:
- **Tractable (G2):** the *chirality-flipping overlap structure* `(f_{L_i} f_{E_j})`
  weighted by the KK mass and the relevant Yukawa insertion gives the leading
  flavor/scale dependence of `C7` (quark `b->s gamma`) and the lepton dipole
  (`mu->e gamma`, `tau->l gamma`). The repo's `flavorConstraints.muToEGamma` already
  encodes the NDA structure `BR ~ |(Y Y^dagger)_{ij}|^2 (ref/M_KK)^4`; rigorize it by
  feeding the *actual* RS overlaps + mass-basis rotation instead of caller-supplied
  spurions, and by computing the loop function prefactor.
- **HONESTLY loop-hard (stays NEEDS-HUMAN):** the finite `O(1)` loop coefficient
  (KK-fermion sum, Goldstone/`A_5`, custodial cancellations) is genuinely loop-level and
  scheme-dependent. G1/G2 fixes the *flavor structure, scale, and sign*; the absolute
  `O(1)` normalization remains a documented theory factor. EDMs (CP-odd dipoles) are in
  the same boat and additionally need the complex phase structure of the 5D Yukawas.

### Inventory headline list
1. KK gauge masses `M_W1, M_Z1` (bessel solver).
2. Gauge-fermion overlap form factor `a(c, eps)` (closed-form Tier-A + numeric Tier-B).
3. Diagonal Z-coupling shifts `delta g_{L,R}^{q,e,nu}` (interaction basis).
4. Mass-basis Z-coupling matrices (diagonal + off-diagonal FCNC-Z) for quarks, charged
   leptons, neutrinos.
5. Charged-current / W form-factor shift (`G_F`, `|V_ij|`).
6. Lepton/neutrino zero-mode `f`'s + charged-lepton `U_L^e, U_R^e`, PMNS.
7. Dipole flavor-overlap structures for `C7` and lepton dipole (flavor+scale rigorous;
   `O(1)` loop factor flagged).
8. Oblique `S, T` from the KK gauge sum (rigorous leading `v^2/M_KK^2`, custodial flagged).

---

## 2. ParameterPoint schema additions

Keep the frozen-container + fail-loud contract. Add one rich object per sector plus a few
scalars, mirroring how `quark_mass_basis_couplings` works. Reuse the two existing
placeholders (`kk_ew_mass_gev`, `lepton_mass_basis_couplings`).

New / re-used `KNOWN_EXTRA_KEYS`:

| key | type | units | real/complex | notes |
|---|---|---|---|---|
| `kk_ew_mass_gev` *(exists)* | `float` | GeV | real | now FILLED = first KK gauge root `M_Z1` |
| `rs_ew_couplings` *(new)* | `RSEWCouplings` dataclass | — | matrices complex | the §1(b–d) object; **the** new container |
| `lepton_mass_basis_couplings` *(exists)* | `LeptonMassBasisCouplings` dataclass | — | matrices complex | now FILLED (charged-lepton + neutrino Z/dipole couplings) |
| `lepton_overlaps` *(new)* | `LeptonOverlaps` dataclass | — | real `f`'s, complex `U` | `fL, fE[3], fN, fN_uv, U_L_e, U_R_e, pmns` |

New frozen dataclasses (mirroring `QuarkMassBasisCouplings`; all matrices `np.ndarray`
shape `(3,3)`, complex; constraints read scalars into `ConstraintResult` and put matrices
into `diagnostics`):

```
@dataclass(frozen=True)
class RSEWCouplings:
    M_KK_gauge_gev: float          # first KK gauge root (Z')
    m_z_gev: float                 # SM Z mass (for the (m_Z/M_KK)^2 normalization)
    g_z_sm: float                  # g/cos theta_W
    # interaction-basis diagonal form factors
    a_left_quark:  np.ndarray (3,) # a(c_{Q_i})
    a_right_up:    np.ndarray (3,)
    a_right_down:  np.ndarray (3,)
    a_left_lepton: np.ndarray (3,) # a(c_{L_i})
    a_right_lepton:np.ndarray (3,) # a(c_{E_i})
    a_neutrino:    np.ndarray (3,) # a(c_{N_i}) (UV-localized -> tiny)
    a_ref: float                   # universal piece reabsorbed
    # mass-basis Z-coupling shift matrices (diagonal + off-diagonal)
    dgL_down:  np.ndarray (3,3); dgR_down: ...; dgL_up: ...; dgR_up: ...
    dgL_lepton: np.ndarray (3,3); dgR_lepton: ...
    dgL_neutrino: np.ndarray (3,3)
    # charged-current
    delta_gf_over_gf: float
    method: str                    # "tierA_closed_form" | "tierB_numeric"
    provenance: Mapping[str, Any]  # c-values, eps, solver roots, sign-check flags
```

`LeptonMassBasisCouplings` carries the lepton analog of `QuarkMassBasisCouplings`
(`U_L_e, U_R_e`, the mass-basis dipole overlap matrices, `M_KK`, `pmns`) so the lepton
dipole / mu-e-conversion / Z-LFV adapters can consume a real object instead of
caller-supplied spurions.

Contract: `make_point` still rejects unknown keys; every new key gets a one-line comment
in `KNOWN_EXTRA_KEYS`; the contract test that pins "every key a constraint reads is
declared" continues to hold.

---

## 3. point_builder logic

Add `build_from_lepton_sector(...)` and extend the unified builder. New REQUIRED sweep
inputs the current quark-only builder lacks: **`c_L, c_E[3], c_N, M_N`** (lepton bulk
masses + Majorana scale), plus the existing `k, Lambda_IR`. These already exist as the
signature of `yukawa.compute_all_yukawas`, so the builder wraps it.

Flow:
```
geometry = get_warp_params(k, Lambda_IR)
eps      = geometry['epsilon']

# (i) KK gauge mass — fills kk_ew_mass_gev
M_KK_gauge = solve_kk("gauge","NN", geometry, n_roots=1)[0][0]

# (ii) form factors a(c) — Tier-A closed form from f_IR/f_UV (+ Tier-B numeric oracle)
a_left_lepton = a_of_c(c_L, eps);  a_right_lepton = a_of_c(c_E, eps); ...

# (iii) Yukawa + rotations (reuse existing module)
yuk = compute_all_yukawas(Lambda_IR, c_L, c_E, c_N, M_N, ...)
U_L_e, U_R_e, pmns = yuk.U_L_e, yuk.U_R_e, yuk.pmns      # (expose if not already returned)

# (iv) diagonal interaction-basis shifts, then rotate to mass basis
dgL_lepton = U_L_e^dagger @ diag(g_z*(m_z/M_KK)^2*(a_left_lepton - a_ref)) @ U_L_e
# ... same for R, neutrino, and (reusing fit_result.U_*) the quark sectors

point = make_point(
    raw=scan_row,
    quark_mass_basis_couplings=qc,           # existing
    kk_ew_mass_gev=float(M_KK_gauge),        # now real, not a proxy
    rs_ew_couplings=RSEWCouplings(...),      # new
    lepton_mass_basis_couplings=LeptonMassBasisCouplings(...),
    lepton_overlaps=LeptonOverlaps(...),
)
```
`build_from_quark_couplings` stays; the EW path is additive. The quark form factors reuse
`fit_result.U_L_d/U_R_d/U_L_u/U_R_u` already threaded through the quark couplings object,
so no quark-fit change is needed — only the *coupling kernel* `g_s^2 f^2 -> g_z (m_z/M_KK)^2 a(c)`.

---

## 4. Matching map (current proxy -> rigorous input)

`R` = becomes FULLY rigorous under G1/G2; `P` = still PARTIAL (documented residual).

### Z-pole (T010–T013) and FCNC-Z (T014)
| adapter / core | current proxy | new rigorous input | formula |
|---|---|---|---|
| `zpole.zbb_coupling_shift_proxy` (`ZPOLE_RS_ZBB_PROXY_V1`) | `dg_b = strength*(m_Z/M_KK)^2*(ovl_b-ovl_light)` | `rs_ew_couplings.dgR_down[2,2], dgL_down[2,2]` | `dg = g_z (m_z/M_KK)^2 (a(c_b)-a_ref)` **R** |
| `zpole` R_q/A_q/A_FB (T010–T013) | same proxy fed into validated pseudo-obs | diagonal `dgL/R` for all quarks+leptons | reuse existing `shifted_couplings`/`asymmetry` **R** |
| `zpole` down-FCNC (`ZPOLE_DOWN_FCNC_PROXY_V1`, T014) | `dg_ij = strength*(m_Z/M_KK)^2*ovl_ij` | `rs_ew_couplings.dgL_down[i,j]` (mass-basis rotation) | off-diag of rotated matrix **R** |

### Z-LFV (T015–T017)
| `zpole_lfv*` (`ZPOLE_LFV_PROXY_V1`) | caller spurion `dg_emu` | `rs_ew_couplings.dgL_lepton[i,j], dgR_lepton[i,j]` | rotated lepton Z matrix **R** |

### Higgs-LFV (T018–T020)
| `higgs_lfv` (`HIGGS_LFV_RS_PROXY_V1`) | caller Higgs-Yukawa spurions `Y_ij` | off-diagonal RS Higgs Yukawa from `f_{L_i} f_{E_j}` misalignment + `U^e` rotation | `Y_ij^h ~ (v/M_KK) * misalignment(f_L,f_E)` **P** (overall width/`O(1)` flagged) |

### b->s l l / s->d l l / c->u l l (rare B/K/C dilepton)
| `rare_b_dilepton` (C9/C10, `..._RS_MATCHING_ASSUMPTION_V1`) | "single Z-like boson, Delta_q=Delta_l=g/(2 c_W)" | Z-penguin from `dgL_down[2,1]` + leptonic SM-Z charges (now derived, not assumed) | `C9,C10 ∝ dg_{sb} * (Z-l vertex)` **R** for Z-penguin piece; box/photon-penguin **P** |
| `rare_kaon_dilepton`/`rare_kaon_snd` | same Z-like assumption (s->d) | `dgL_down[1,0]` | **R** (Z-penguin), long-distance **P** |
| `rare_charm_dilepton` | Z-like c->u | `dgL_up[1,0], dgR_up[1,0]` | **R** (Z-penguin); GIM/long-distance **P** |
| `rare_charm_lfv_*` (`..._PROXY_V1`) | caller e-mu Z-like | `dgL_lepton/dgR_lepton` off-diag + up-sector FCNC | **R** short-distance |

### b->s nu nu
| `rare_b_nunu` (`..._RS_MATCHING_ASSUMPTION_V1`) | "Z-like, Delta_q=Delta_nu=g/(2 c_W)" | `dgL_down[2,1]` + neutrino Z charge (`dgL_neutrino` ~ 0, UV) | `C_L^nu ∝ dg_{sb}` **R** (Z-penguin); box **P** |

### Dipoles (mu->e gamma / tau->l gamma; b->s gamma)
| `lepton` (mu->e gamma, `LEPTON_DIPOLE_PROXY_ASSUMPTION_V1`) | caller `Y_N_bar`,PMNS,M_KK -> NDA `BR~|(YY†)_{ij}|^2(ref/M_KK)^4` | `lepton_mass_basis_couplings`: real RS overlaps + `U^e` rotation | flavor+scale **R**; `O(1)` loop factor **P** |
| `bsgamma` (C7, `BSGAMMA_RS_MATCHING_ASSUMPTION_V1`) | KK-gluon overlap as b-s flavor proxy, LL-run | mass-basis chirality-flip overlap `(f_{L_b} f_{E_s})` weighted | C7 flavor/scale **R**; finite C7/C8 `O(1)` **P** |

### mu-e conversion + LFV 3-body
| `mu_e_conversion` (`MU_E_CONVERSION_PROXY_V1`) | caller nucleon `g_LV/RV/LS/RS` | dipole from `lepton_mass_basis_couplings` + vector from `dgL/R_lepton` (Z-exchange to quarks) | dipole+vector-Z **R**; scalar/tensor **P** |
| `lfv_three_body*` (`LFV_THREE_BODY_PROXY_V1`) | caller Z-penguin + box amplitudes | `dgL/R_lepton` off-diag for the Z-penguin `G_AB^Z` | Z-penguin **R**; genuine box **P** |

### EDMs
| `edm`/`atomic`/`neutron`/`mercury` (`EDM_RS_MATCHING_GAP_V1`) | caller CP-odd dipole coeff | complex phases of RS lepton/quark dipole overlaps | **P** — stays NEEDS-HUMAN (loop-level + complex 5D Yukawa phases); G1 only supplies the *real* overlap magnitude as a diagnostic envelope |

### Charged-current (EW002/EW003/K018) + oblique (EW001)
| `semileptonic_lfu` (`..._RS_MATCHING_ASSUMPTION_V1`) | KK-gluon charged-current stress proxy | `rs_ew_couplings.delta_gf_over_gf` + per-generation `a(c_L)` | LFU ratios **R** (leading); form-factor **P** |
| `oblique_stu` (`OBLIQUE_STU_RS_PROXY_V1`) | `dS=c_S v^2/M_KK^2`, `dT=pi L/(2cW^2) v^2/M_KK^2` | KK gauge-sum `S,T` from actual `M_KK_gauge` + profiles | leading `v^2/M_KK^2` **R**; custodial/brane-kinetic `O(1)` **P** |

### Count
- **FULLY rigorous (R) under G1/G2:** Z-pole (4: T010–T013), FCNC-Z (1: T014), Z-LFV (3:
  T015–T017), b->s l l Z-penguin (B family ~3), s->d l l Z-penguin (K family ~3), c->u l l
  (C family ~2), b->s nu nu (1–2), charged-current LFU leading (2–3), KK gauge mass
  consumers (CR/graviton, ~2). **≈ 24–26 constraints become fully rigorous** for their
  short-distance/Z-mediated piece.
- **PARTIAL (R short-distance, P residual):** all dipoles (mu->e gamma, tau->l gamma,
  b->s gamma; ~5), Higgs-LFV (3), mu-e conversion (1–3), LFV 3-body (3), oblique (1).
  **≈ 13–16 constraints become partially rigorous** (flavor/scale/sign rigorous, `O(1)`
  loop or hadronic factor flagged).
- **Stays NEEDS-HUMAN even after G1/G2:** all EDMs (electron/neutron/mercury/atomic; ~4)
  — genuine loop + complex 5D-Yukawa-phase items; the finite `O(1)` dipole loop
  coefficients; lattice/long-distance hadronic inputs (already SM-side, out of scope).

---

## 5. Phasing & risk

### Recommended implementation order (unlocks most first)
1. **Block 1 — KK gauge mass + form factor `a(c)` (Tier-A).** Pure function of
   `c, eps`; no rotations. Unlocks the *normalization* (`M_KK_gauge`, `(m_Z/M_KK)^2`) and
   the diagonal Z-pole shifts. Highest leverage: Z-pole (T010–T013) + Zbb sign check.
2. **Block 2 — quark mass-basis Z matrices.** Reuse existing `U_L/R` rotations; swap the
   coupling kernel. Unlocks FCNC-Z (T014), all b/k/charm dilepton Z-penguins, b->s nu nu.
   These reuse already-validated SM cores -> large constraint count for small new code.
3. **Block 3 — lepton sector (`c_L,c_E,c_N`, `U^e`, PMNS).** New sweep inputs + builder
   path. Unlocks Z-LFV, LFV 3-body Z-penguin, lepton dipole flavor structure, mu-e
   conversion vector piece.
4. **Block 4 — charged-current + oblique leading order.** Unlocks EW002/EW003/K018 leading,
   EW001 leading.
5. **Block 5 — Tier-B numeric overlap** as validation oracle + dipole `O(1)` study (mostly
   for confidence, not new constraints).

### Stays NEEDS-HUMAN even after G1/G2
- Charged-lepton & quark **EDMs** (complex 5D Yukawa phases + loop).
- The **finite `O(1)` loop coefficient** of every dipole (C7/C8 absolute normalization,
  lepton dipole prefactor) — flavor/scale rigorous, magnitude is a documented theory factor.
- **Custodial / brane-kinetic-term** structure protecting `Z b_L b_L` and shifting `T` —
  genuinely model-dependent; G1 uses the minimal (non-custodial) value and flags it.
- Hadronic long-distance pieces (charm loops in b->s l l, `K_L->pi0 e e` interference) —
  already SM-side, out of scope.

### Top 3 correctness risks
1. **Coupling normalization / convention consistency** with the *existing* validated
   `zpole`/`rare_*` cores. Those cores expect `delta_g_left/right` defined relative to a
   specific `g_Z`, `sin^2 theta_W`, and charge convention (`ZPoleSMInputs`,
   `shifted_couplings`). The new `a(c)` -> `delta g` map MUST land in the *same*
   normalization, or every ratio is silently off by `g_Z`/`cos theta_W` factors. Pin with
   a unit test that the SM limit (`a(c)=a_ref`) gives exactly zero shift and reproduces the
   existing SM pseudo-observables.
2. **Basis / PMNS rotation for leptons (double-rotation hazard).** Work strictly in the
   charged-lepton mass basis (CLAUDE.md). The lepton Z matrices use `U_L^e/U_R^e`;
   neutrino mixing is PMNS. Mixing these up (rotating the lepton Z shift by PMNS, or
   forgetting to rotate at all) produces spurious or vanishing FCNC. Cross-check against
   the quark path which already does the analogous rotation correctly.
3. **Double-counting the universal piece.** `a_ref` (the flavor-blind overlap) is absorbed
   into the measured `g_Z, m_Z, G_F`. If the diagonal shift uses `a(c)` instead of
   `a(c)-a_ref`, a large *universal* shift contaminates flavor-conserving observables (and
   the oblique `S,T`). The subtraction must be explicit and tested (flat/universal `c`
   -> all `delta g = 0`).

---

## 6. Validation strategy

1. **Reproduce the classic RS `Z b bbar` shift.** With `c_b < c_light` and a custodial-off
   minimal setup, the rigorous `delta g_R^b` must be **negative** and of order
   `(m_Z/M_KK)^2 * O(f_IR^2(c_b))`; magnitude `~10^-3` for `M_KK ~ few TeV`. Compare sign
   and order to the literature (Agashe-Contino-Pomarol; Casagrande et al.) — the proxy was
   hand-built to mimic this, so it is the anchor test.
2. **KK mass spacing.** `solve_kk("gauge","NN")` first roots must follow the gauge tower
   `x_n ≈ 2.40, 5.52, 8.65` (`J_0` zeros) in the `eps->0` limit, i.e. `M_n/M_1 ≈ 2.3, 3.6`.
   Already exercised by the bessel solver; add a regression pin.
3. **SM limit / universal-c test.** Set all `c` equal -> all `delta g_{ij} = 0` (diagonal
   and off-diagonal). Every Z-pole / FCNC / dilepton prediction must collapse to its
   existing validated SM value (the constraint files already pin the SM side).
4. **Tier-A vs Tier-B agreement.** The closed-form `a(c)` must match the numeric quadrature
   overlap (Block 5) to the `O(m_Z^2/M_KK^2)` truncation accuracy across `c in [0.3,0.9]`.
5. **Decoupling.** All shifts scale as `1/M_KK^2`; verify `delta g`, FCNC BRs -> 0 as
   `Lambda_IR -> infinity`, and that the constraints' `ratio` saturations track it.
6. **Cross-sector consistency.** The same `a(c)` used in quark Z-pole and in b->s l l must
   give a `dg_{sb}` consistent with the `dg_b` diagonal under one rotation — a single
   internal consistency assertion catches normalization drift between adapters.

---

## Appendix: files touched (when implemented)
- `flavor_catalog_constraints/point_builder.py` — add keys + `build_from_lepton_sector`.
- `flavor_catalog_constraints/physics_adapters/` — re-wire the listed adapters to read
  `rs_ew_couplings` / `lepton_mass_basis_couplings` instead of caller proxies.
- NEW `warpConfig/` helper or `quarkConstraints/rs_ew.py` — `a_of_c(c, eps)` (Tier-A) +
  numeric Tier-B overlap; `RSEWCouplings`, `LeptonMassBasisCouplings`, `LeptonOverlaps`.
- `yukawa/compute_yukawas.py` — expose `U_L_e, U_R_e, pmns` on `YukawaResult` if absent.
- `derivations/flavor/gauge_kk_coupling.tex` — promote README stub to the overlap-integral
  derivation backing `a(c)` and `delta g`.
