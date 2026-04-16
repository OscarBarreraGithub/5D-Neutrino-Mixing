# Quark-Sector Scan: Complete Assumptions and Pipeline

This document traces every assumption, formula, and convention that goes into
producing the publication figures `fig1_exclusion_boundaries` and
`fig2_mkk_bound_2007_vs_modern`. It follows the computation in the order it
actually executes, from geometry definition through to the final pixel.

---

## Step 0: Define the Warped Geometry

The spacetime is a slice of AdS_5 (Randall-Sundrum I) compactified on an
S_1/Z_2 orbifold, bounded by two 3-branes:

- **UV brane** at the orbifold fixed point y = 0 (Planck-scale physics)
- **IR brane** at y = pi r_c (TeV-scale physics)

The 5D metric in conformal coordinates (mostly-plus signature) is:

```
ds^2 = (k z)^{-2} (eta_{mu nu} dx^mu dx^nu + dz^2)
```

where eta = diag(-1, +1, +1, +1) and z runs from z_h = 1/k (UV brane) to
z_v = 1/Lambda_IR (IR brane). In proper-distance coordinates:

```
ds^2 = e^{-2 k |y|} eta_{mu nu} dx^mu dx^nu + dy^2
```

Two parameters fix the entire geometry:

| Parameter | Symbol | Default | Meaning |
|-----------|--------|---------|---------|
| AdS curvature | k | 1.2209 x 10^19 GeV (= M_Pl) | Sets the UV scale |
| IR scale | Lambda_IR | scanned (500 -- 20000 GeV) | Sets the TeV-scale brane |

From these, everything else is derived:

- **Warp factor**: epsilon = Lambda_IR / k (typically ~ 10^{-15} for
  Lambda_IR ~ few TeV)
- **Orbifold radius**: r_c = ln(k / Lambda_IR) / (pi k). The warp factor
  and radius are related by epsilon = e^{-pi k r_c}.
- **UV brane position** (conformal): z_h = 1/k ~ 10^{-19} GeV^{-1}
- **IR brane position** (conformal): z_v = 1/Lambda_IR ~ 10^{-4} GeV^{-1}

**Fixed external parameter**:
- **Electroweak VEV**: v = 174 GeV. This is the Higgs VEV, not the reduced
  VEV v/sqrt(2) = 246/sqrt(2). The code uses 174 GeV throughout.

**Assumption 0a**: The Higgs is localized on (or very near) the IR brane. This
is why the Yukawa interaction is proportional to the fermion zero-mode value
at the IR brane, and why TeV-scale new physics is connected to the IR scale.

**Assumption 0b**: k is fixed at M_Pl = 1.2209 x 10^19 GeV. The scan does
not vary k. This means the warp factor epsilon and the hierarchy are controlled
entirely by Lambda_IR.

**Assumption 0c**: Backreaction of bulk fields on the geometry is neglected.
The metric is taken to be the pure AdS_5 slice throughout.

---

## Step 1: Zero-Mode Profiles (f-factors)

Each 5D fermion field has a bulk Dirac mass M_5. The dimensionless bulk mass
parameter is c = M_5 / k. The zero-mode wavefunction in the extra dimension
is determined by c: it satisfies the 5D Dirac equation with the AdS background,
and the boundary conditions select a single normalizable zero mode per chirality.

The physical quantity that enters all 4D couplings is the value of this
normalized zero-mode wavefunction at the IR brane. For left-handed zero modes
(the ones relevant for Higgs overlap), this "f-factor" is:

```
f_IR^2(c, epsilon) = (1/2 - c) / (1 - epsilon^{1 - 2c})
f_IR(c, epsilon) = sqrt(max(f_IR^2, 0))
```

The formula arises from integrating the square of the zero-mode profile
z^{2-c} over the extra dimension to normalize it, then evaluating at z = z_v.
The max(..., 0) is a numerical safeguard against round-off.

**Physical behavior**:
- c > 1/2: The exponent (1-2c) is negative, so epsilon^{1-2c} is huge
  (recall epsilon ~ 10^{-15}). The denominator is dominated by -epsilon^{1-2c},
  and f_IR^2 ~ (c - 1/2) * epsilon^{2c-1}, which is exponentially small.
  *These fermions are UV-localized and have tiny overlap with the IR-brane Higgs.*
- c < 1/2: The exponent (1-2c) is positive, so epsilon^{1-2c} ~ 0. Then
  f_IR^2 ~ (1/2 - c), which is O(1).
  *These fermions are IR-localized and have large Yukawa overlap.*
- c = 1/2: The formula is 0/0. The limit is f_IR^2 = 1/(-2 ln epsilon).
  For epsilon ~ 10^{-15}, this gives f_IR^2 ~ 1/69 ~ 0.014.

**Numerical examples** (at epsilon = 3000 / 1.2209e19 ~ 2.46 x 10^{-16}):
- c = 0.72 (light quark): f_IR ~ 4.4 x 10^{-8}
- c = 0.51 (near flat): f_IR ~ 0.11
- c = 0.30 (heavy quark): f_IR ~ 0.45

This enormous dynamic range — 10 orders of magnitude between c = 0.72 and
c = 0.30 — is the core of the RS flavor mechanism. Mass hierarchies come from
geometry, not from hierarchical Yukawa couplings.

**Convention**: c = M_5 / k, so positive c means the bulk mass and the AdS
curvature have the same sign. This is the standard convention in the RS flavor
literature (Agashe, Contino, Perez 2005; Csaki et al. 2009).

There is also a UV-brane overlap f_UV, used for neutrino Majorana masses in
the lepton sector, but it does not enter the quark scan.

---

## Step 2: Choose the Scan Grid

The dense scan (the one used for the publication figures) sweeps a 3D grid:

| Parameter | Grid | Range | Points |
|-----------|------|-------|--------|
| r | log-uniform | 0.02 -- 2.0 | 100 |
| overall_scale | linear | 1.5 -- 6.0 | 20 |
| Lambda_IR | log-uniform | 500 -- 20000 GeV | 50 |

Total: 100 x 20 x 50 = 100,000 scan points.

The grid is enumerated deterministically: Lambda_IR (outer) -> overall_scale
(middle) -> r (inner). Each point is fitted independently using
default_spurion_seed() as the starting point (there is no warm-starting from
the previous r). Points are distributed across SLURM shards for parallelism;
the result is independent of shard count or execution order.

**Fixed parameters at scan time**:
- xi_KK = 1.0. This means M_KK = Lambda_IR in the scan's bookkeeping
  convention. The Wilson coefficients use M_KK = Lambda_IR as the propagator
  mass. (See Step 16b for the post-hoc conversion to a physical KK-gluon mass.)
- g_s* = 3.0. The enhanced KK-gluon coupling, not the perturbative g_s. This
  value is the one-loop estimate from Gedalia, Grossman, Nir, Perez (2009).
  It is passed directly into the coupling computation at scan time, so all
  stored Wilson coefficients and ratios already include this enhancement.
- k = 1.2209 x 10^19 GeV (= M_Pl).
- v = 174 GeV (electroweak VEV).

**The initial fit seed** is the deterministic default from the benchmarks
module:
```
up_singular_values   = [0.003,  0.12, 1.0]
down_singular_values = [0.006, 0.035, 0.45]
overall_scale        = 2.8
up_left   = (theta12=0.02,  theta13=0.002,  theta23=0.03,  delta=0.8)
down_left = (theta12=0.23,  theta13=0.0035, theta23=0.041, delta=1.2)
```
These values were chosen by hand to be in the right ballpark for SM-like quark
masses. The raw seed also includes nonzero right rotations (e.g.,
up_right.theta12 = 0.04), but these are set to identity during the
canonicalization step at the start of the fit (Step 10). For a new scan point,
the seed is either this default (at the start of a new Lambda_IR /
overall_scale slice) or the converged result from the previous r value.

---

## Step 3: Construct the Yukawa Spurions

Y_u and Y_d are complex 3x3 matrices. They are the 5D Yukawa coupling matrices
that appear in the action on the IR brane as L = Y_{ij} Q_i H u_j (schematic).
In the scan, they are parametrized via a singular value decomposition:

```
Y_u = overall_scale * U_left^u . diag(sigma_u) . (U_right^u)^dagger
Y_d = overall_scale * U_left^d . diag(sigma_d) . (U_right^d)^dagger
```

where:
- sigma_u = (sigma_u1, sigma_u2, sigma_u3) are 3 positive singular values,
  ordered light-to-heavy. These control the "hierarchy" within each Yukawa.
- sigma_d = (sigma_d1, sigma_d2, sigma_d3) similarly for the down sector.
- U_left^{u,d} are CKM-like 3x3 unitary matrices, each with 4 real
  parameters: three Euler angles (theta_12, theta_13, theta_23) and one
  CP-violating phase (delta). These control the "orientation" of the Yukawa
  in flavor space.
- U_right^{u,d} are fixed to the 3x3 identity matrix (not fitted). This is
  a restriction on the parameter space — see below.
- overall_scale is a common real positive prefactor multiplying both Y_u and
  Y_d. It controls the overall magnitude of the Yukawa entries.

The unitary matrices use the standard CKM (PDG) parametrization:
```
U = [[c12 c13,                           s12 c13,                        s13 e^{-i delta}],
     [-s12 c23 - c12 s23 s13 e^{i delta}, c12 c23 - s12 s23 s13 e^{i delta}, s23 c13       ],
     [ s12 s23 - c12 c23 s13 e^{i delta},-c12 s23 - s12 c23 s13 e^{i delta}, c23 c13       ]]
```
where c_ij = cos(theta_ij), s_ij = sin(theta_ij). All angles are wrapped
into [-pi, pi).

**Canonicalization**: Before the fit begins, the seed is "canonicalized":
overall_scale is absorbed into the singular values (multiplied in), and
overall_scale is set to 1.0. So the fitter works with "physical singular
values" = overall_scale * sigma. The right rotations are set to identity.
This removes a redundancy (rescaling all singular values by a common factor
and compensating with overall_scale would give the same Y).

After canonicalization, the actual objects being constructed are:
```
Y_u = U_left^u . diag(physical_sigma_u) . I = U_left^u . diag(physical_sigma_u)
Y_d = U_left^d . diag(physical_sigma_d) . I = U_left^d . diag(physical_sigma_d)
```

**Consequence of U_right = I**: Since the right rotations are identity:
- Y_u^dagger Y_u = diag(sigma_u)^2 (diagonal)
- Y_d^dagger Y_d = diag(sigma_d)^2 (diagonal)
- Y_u Y_u^dagger = U_left^u . diag(sigma_u)^2 . (U_left^u)^dagger (generally dense)
- Y_d Y_d^dagger = U_left^d . diag(sigma_d)^2 . (U_left^d)^dagger (generally dense)

So C_u and C_d (defined next) are automatically diagonal, while C_Q is
generally dense. This means all off-diagonal flavor structure in the left-handed
quark sector — and therefore all tree-level FCNCs — come from the combination
of the two left rotations U_left^u and U_left^d through C_Q.

**Assumption 3a**: Right rotations are frozen to identity. This restricts the
Yukawa space: a general complex 3x3 matrix has 18 real parameters, but this
parametrization has only 7 per sector (3 singular values + 4 angles), plus 1
shared overall_scale = 15 total. The full space would have 36 + 1 = 37.

---

## Step 4: Build the Bulk Mass Matrices

From the Yukawa spurions, three Hermitian matrices are constructed:

```
C_u = Y_u^dagger . Y_u           [controls RH up-type quark localization]
C_d = Y_d^dagger . Y_d           [controls RH down-type quark localization]
C_Q = r . (Y_u . Y_u^dagger) + (Y_d . Y_d^dagger)   [controls LH doublet localization]
```

**Physical motivation for C_Q**: In the RS NMFV framework, the 5D bulk mass
matrix for the left-handed quark doublets is assumed to be a function of the
Yukawa spurions (since the Yukawas are the only sources of flavor breaking).
At leading order in the spurion expansion:

```
C_Q = alpha_Q . I + beta_Q . Y_d Y_d^dagger + gamma_Q . Y_u Y_u^dagger + ...
```

The parameter r = gamma_Q / beta_Q measures the relative weight of the up-type
vs down-type Yukawa contribution. The code uses the simplified form with
beta_Q = 1 (normalized out) and no flavor-universal alpha_Q term:

```
C_Q = r . (Y_u Y_u^dagger) + (Y_d Y_d^dagger)
```

C_Q is Hermitian and positive semi-definite (it is a sum of positive
semi-definite matrices with non-negative coefficients).

**The role of r**: When r = 0, C_Q = Y_d Y_d^dagger. If we work in the basis
where Y_d is diagonal (which we can always choose), then C_Q is also diagonal.
In that basis, the left-handed overlap matrix F_Q is diagonal, and there are
no tree-level FCNCs in the down sector. This is the "alignment" limit.

When r > 0, the Y_u Y_u^dagger term introduces off-diagonal entries (since Y_u
contains CKM-like rotations relative to Y_d). These off-diagonal entries in C_Q
misalign the left-handed quark basis from the down-Yukawa eigenbasis, generating
tree-level FCNCs.

**Key difference from the standard NMFV parametrization**: The literature
(Csaki et al. 0907.0474) writes C_Q = alpha_Q I + beta_Q Y_d Y_d^dagger +
gamma_Q Y_u Y_u^dagger. The alpha_Q I term provides a flavor-universal baseline
bulk mass (typically alpha_Q ~ 0.5-0.7 to localize quarks near the UV brane).
This scan omits alpha_Q I. Its role is instead played by the BulkMassMap
(Step 5), which maps the eigenvalues of C_Q onto a bounded physical range
[0.30, 0.72]. The BulkMassMap gives c_Q = 0.72 when the eigenvalue is zero
(equivalent to alpha_Q = 0.72), but for nonzero eigenvalues the relationship
is nonlinear (sigmoid), not linear as it would be with explicit alpha_Q.

**Note on C_u and C_d**: Because U_right = I (Step 3), both C_u and C_d are
diagonal: C_u = diag(sigma_u^2), C_d = diag(sigma_d^2). Only C_Q carries
off-diagonal structure.

---

## Step 5: Diagonalize and Map to Physical Bulk Masses

### 5a: Hermitian eigendecomposition

Each of C_Q, C_u, C_d is diagonalized using numpy.linalg.eigh (which is
specialized for Hermitian matrices and guaranteed to return real eigenvalues):

```
C_X = rotation_X . diag(lambda_X) . rotation_X^dagger
```

The eigenvalues lambda_X are sorted in ascending order (lightest generation
first, heaviest third). The rotation_X are the corresponding unitary eigenvector
matrices, with columns ordered to match.

For C_u and C_d, which are already diagonal (see Step 4), the eigenvalues are
just the diagonal entries (sigma^2) and the rotations are (up to column
permutation) the identity.

For C_Q, the eigenvalues and rotations are nontrivial. The rotation_Q matrix
is the key object: it defines the basis transformation from the "spurion basis"
(where Y_u, Y_d were constructed) to the "bulk-mass eigenbasis" (where the
bulk masses are diagonal and the f-factors can be applied).

### 5b: Nonlinear map from eigenvalues to bulk masses

The eigenvalues lambda_X (non-negative reals) are mapped to physical bulk mass
parameters c_X through a bounded monotone sigmoid:

```
c_X = c_uv - (c_uv - c_ir) . lambda / (lambda + eigen_scale)
```

With the defaults c_uv = 0.72, c_ir = 0.30, eigen_scale = 1.0, this becomes:

```
c_X = 0.72 - 0.42 . lambda / (lambda + 1)
```

This map has the following properties:
- lambda = 0  =>  c = 0.72  (UV-localized, light fermion)
- lambda = 1  =>  c = 0.51  (intermediate)
- lambda -> infinity  =>  c -> 0.30  (IR-localized, heavy fermion)
- The derivative is dc/dlambda = -0.42 / (lambda + 1)^2, which is always
  negative (larger eigenvalue = smaller c = more IR-localized).
- The map is bounded: c is always in [0.30, 0.72], regardless of the
  eigenvalue.

**Why this map is needed**: Without it, there is no guaranteed connection
between the spurion eigenvalues and a physically meaningful bulk mass. The
zero-mode profile formula f_IR(c, epsilon) requires c to be in a specific
range to produce sensible results (c < 0 gives imaginary f-factors; c > ~0.8
gives essentially zero overlap). The BulkMassMap enforces this constraint
by construction.

**Assumption 5a**: The window [0.30, 0.72] is a design choice. c_uv = 0.72
means that a fermion with zero spurion eigenvalue (i.e., decoupled from the
Yukawa structure) sits at c = 0.72, which is UV-localized and very light
(f_IR ~ 4 x 10^{-8}). c_ir = 0.30 means that the most IR-localized fermion
has c = 0.30 (f_IR ~ 0.45). The top quark, which needs f_IR ~ O(1), would
need c ~ 0.3 or below; the BulkMassMap cannot produce c < 0.30. This forces
the top Yukawa singular value to be ~ 2.4 to reproduce m_t = 172 GeV (see
Step 8), which is perturbative but larger than the "anarchic O(1)" expectation.

**Assumption 5b**: The sigmoid shape (lambda/(lambda+1)) is a convenience
choice. In the standard parametrization with explicit (alpha_Q, beta_Q), the
map would be linear: c_{Q,i} = alpha_Q + beta_Q * lambda_i. The sigmoid
compresses large eigenvalues, preventing any c from falling below c_ir.

---

## Step 6: Compute Zero-Mode Overlaps

With the physical bulk masses c_Q, c_u, c_d (each a length-3 vector, one per
generation), the zero-mode IR-brane overlaps are computed element-wise using
the formula from Step 1:

```
F_Q[i] = f_IR(c_Q[i], epsilon)     for i = 1, 2, 3  (generations)
F_u[i] = f_IR(c_u[i], epsilon)
F_d[i] = f_IR(c_d[i], epsilon)
```

These are real, positive numbers stored as 1D arrays of length 3. They will
be used as diagonal matrices in the mass formula.

Because of the BulkMassMap range [0.30, 0.72], all c values are > 1/2 except
possibly the third generation. In practice:
- First-generation quarks (c ~ 0.7): F ~ 10^{-7} to 10^{-8}
- Second-generation quarks (c ~ 0.55-0.65): F ~ 10^{-3} to 10^{-5}
- Third-generation quarks (c ~ 0.30-0.45): F ~ 0.1 to 0.45

---

## Step 7: Rotate Yukawas to the Bulk-Mass Basis

The diagonalization of C_Q, C_u, C_d in Step 5 defines three new bases. The
Yukawas must be rotated into the bulk-mass eigenbasis before the mass matrices
can be built, because the overlap matrices F_Q, F_u, F_d are diagonal only in
that basis.

The rotation is:
```
Y_u^{bulk} = rotation_Q^dagger . Y_u . rotation_u
Y_d^{bulk} = rotation_Q^dagger . Y_d . rotation_d
```

where rotation_Q, rotation_u, rotation_d are the unitary eigenvector matrices
from Step 5.

**Important**: This rotation does NOT diagonalize the Yukawas. It
diagonalizes the bulk mass matrices. The rotated Yukawas Y_u^{bulk} and
Y_d^{bulk} are generally dense 3x3 complex matrices. Off-diagonal entries in
these rotated Yukawas, combined with the non-degenerate diagonal overlap
factors, are what ultimately produce off-diagonal entries in the mass matrices
and drive the CKM mixing.

**Why rotation_u and rotation_d matter here**: Even though C_u and C_d are
diagonal (so rotation_u and rotation_d are trivially identity up to column
ordering), the column ordering from the eigenvalue sort in Step 5a can
permute the right-hand side. In practice, since the singular values are
ordered light-to-heavy, and the eigenvalues of C_u = diag(sigma_u^2) are
also ordered by the sort, these rotations reduce to the identity or at
most a column permutation.

---

## Step 8: Build 4D Mass Matrices

The effective 4D quark mass matrices in the bulk-mass eigenbasis are:

```
M_u = 2v . diag(F_Q) . Y_u^{bulk} . diag(F_u)
M_d = 2v . diag(F_Q) . Y_d^{bulk} . diag(F_d)
```

where v = 174 GeV is the electroweak VEV.

**Origin of the factor 2v**: With v = 174 GeV, the prefactor 2v = 348 GeV.
Note that 2 x 174 = 348 ~ sqrt(2) x 246, where 246 GeV is the standard
electroweak VEV (v_SM = 246 GeV, with the Higgs field normalization
<H> = v_SM / sqrt(2) = 174 GeV). So the mass formula is equivalently:

```
M = sqrt(2) . v_SM . diag(F_Q) . Y^{bulk} . diag(F)
```

The factor of sqrt(2) arises from the normalization of the Higgs doublet
VEV relative to the mass term in the Lagrangian. The 5D Yukawa couplings Y
in this code are dimensionless (they absorb the factor of k that appears
in the lepton-sector formula m_E = 2vk f_L Y_E f_E used elsewhere in this
repo). What matters is that the target quark masses (Step 10) are defined
consistently with this convention.

**Structure of M_u and M_d**: The left overlap diag(F_Q) is the same for both
sectors (since Q is the SU(2) doublet containing both u_L and d_L). The right
overlaps diag(F_u) and diag(F_d) differ. The Yukawas Y_u^{bulk} and
Y_d^{bulk} are generally dense.

The physical quark mass for generation i is roughly:
```
m_i ~ 2v . F_Q[i] . Y^{bulk}_{ii} . F_{u or d}[i]
```
for the approximately diagonal case. The exponential hierarchy in the F-factors
produces the quark mass hierarchy without requiring hierarchical Yukawa entries.

---

## Step 9: Extract Physical Masses and CKM

Each mass matrix is decomposed via SVD (biunitary diagonalization) using the
convention:

```
M = U_L . diag(singular values) . U_R^dagger
```

so that U_L^dagger . M . U_R = diag(m_1, m_2, m_3) with m_i >= 0. The code
uses numpy's SVD: M = U . diag(s) . V^dagger, then identifies U_L = U,
U_R = V.

The singular values are the physical quark masses. They are sorted in ascending
order (m_1 < m_2 < m_3 corresponds to u < c < t or d < s < b).

The CKM matrix is constructed from the left-handed rotation matrices:

```
V_CKM = (U_L^u)^dagger . U_L^d
```

This is the standard definition: V_CKM rotates from the down mass eigenbasis
to the up mass eigenbasis for left-handed quarks.

**CKM observables extracted** (4 real numbers):
1. |V_us| = |V_CKM[0,1]| (Cabibbo angle)
2. |V_cb| = |V_CKM[1,2]|
3. |V_ub| = |V_CKM[0,2]|
4. J = Im(V_CKM[0,1] . V_CKM[1,2] . V_CKM[0,2]* . V_CKM[1,1]*)
   — the Jarlskog invariant, a rephasing-invariant measure of CP violation

---

## Step 10: The Fit (Optimization Loop)

Steps 3 through 9 are wrapped inside a Trust Region Reflective (TRF) optimizer
(scipy.optimize.least_squares). The fit adjusts the Yukawa parameters until the
predicted quark masses and CKM match the observed values.

### 10a: Free parameters (14 total)

The optimizer works in a "canonical chart" with 14 real coordinates:

1-3. ln(physical_sigma_u[1]), ln(physical_sigma_u[2]), ln(physical_sigma_u[3])
     — log of the up-sector singular values (the log ensures positivity)
4-6. ln(physical_sigma_d[1]), ln(physical_sigma_d[2]), ln(physical_sigma_d[3])
     — log of the down-sector singular values
7-10. theta_12^u, theta_13^u, theta_23^u, delta^u
     — four parameters of the up-sector left rotation
11-14. theta_12^d, theta_13^d, theta_23^d, delta^d
     — four parameters of the down-sector left rotation

All angles are in [-pi, pi) after canonicalization.

**Fixed during fit**: r, overall_scale (absorbed into singular values as 1.0),
Lambda_IR, k, v, BulkMassMap parameters, right rotations (= I).

### 10b: Target observables (10 total)

| Observable | Target value | Source |
|------------|-------------|--------|
| m_u | 0.0013 GeV | Repo-owned target set |
| m_c | 0.62 GeV | Repo-owned target set |
| m_t | 172.0 GeV | Repo-owned target set |
| m_d | 0.0028 GeV | Repo-owned target set |
| m_s | 0.057 GeV | Repo-owned target set |
| m_b | 2.86 GeV | Repo-owned target set |
| \|V_us\| | 0.225 | CKM parametrization |
| \|V_cb\| | 0.0415 | CKM parametrization |
| \|V_ub\| | 0.00361 | CKM parametrization |
| J | ~3.08 x 10^{-5} | From (theta12=0.2274, theta13=0.00368, theta23=0.0415, delta=1.196) |

The target CKM is built from a CKM-like unitary with angles
(theta_12=0.2274, theta_13=0.00368, theta_23=0.0415, delta=1.196). The four
CKM observables (|V_us|, |V_cb|, |V_ub|, J) are extracted from this matrix.

**Note on mass targets**: These are repo-owned "effective" quark masses, not
MS-bar masses at any single consistent renormalization scale. Some values
(m_u = 1.3 MeV, m_d = 2.8 MeV) are close to MS-bar at mu = 2 GeV; others
(m_c = 0.62 GeV) are closer to MS-bar at mu ~ 3 TeV; m_b = 2.86 GeV is
close to MS-bar at mu = m_b. The values were chosen for scan stability and
to give physically reasonable bulk mass parameters, not for precise
consistency with a single renormalization scheme. Since the fit determines
the bulk masses (and thus the f-factors and FCNC couplings) from these
targets, any systematic shift in the targets would propagate into the
constraint ratios. This is a source of unquantified systematic uncertainty.

### 10c: Residual vector (10 components)

The optimizer minimizes a vector of 10 residuals:

```
residuals = [
  ln(m_u^pred / 0.0013),    ln(m_c^pred / 0.62),     ln(m_t^pred / 172.0),
  ln(m_d^pred / 0.0028),    ln(m_s^pred / 0.057),     ln(m_b^pred / 2.86),
  (|V_us|^pred - 0.225) / 0.225,
  (|V_cb|^pred - 0.0415) / 0.0415,
  (|V_ub|^pred - 0.00361) / 0.00361,
  (J^pred - J^target) / |J^target|
]
```

Mass residuals are in log-space (dimensionless), CKM residuals are
fractional deviations normalized by the target value.

**Fit score** = sqrt(mean(residuals^2)). This is the RMS of all 10 residuals.
A score of 0.1 means ~ 10% average deviation across all observables.

### 10d: Optimizer settings

- **Method**: Trust Region Reflective (TRF) (scipy least_squares default, "trf" method)
- **max_nfev**: 120 function evaluations in the scan config (overriding the
  fit_quark_sector function default of 200). Each evaluation runs Steps 3-9
  once.
- **Error handling**: If any step inside the residual vector computation
  raises a ValueError or LinAlgError (e.g., negative eigenvalues from an
  extreme parameter combination), the residual returns a vector of 1e6
  for all 10 components. This steers the optimizer away from unphysical
  regions without crashing.

### 10e: Convergence and quality

After the optimizer finishes (either converging or exhausting max_nfev):

- The best-fit parameters are decoded back into a seed.
- The best-fit point is re-evaluated to get final masses, CKM, and score.
- If fit_score > 0.1, the point is flagged as poorly converged and will be
  discarded at plotting time (Step 16a).

---

## Step 11: Mass-Basis KK-Gluon Couplings

After the fit converges, we have: (1) the diagonal overlap profiles F_Q, F_u,
F_d, and (2) the unitary rotation matrices U_L^u, U_R^u, U_L^d, U_R^d from the
SVD of the mass matrices. These are combined to build the KK-gluon coupling
matrices in the physical quark mass basis.

The first KK-gluon excitation couples to zero-mode quarks proportionally to
their IR-brane wavefunction values squared. In the bulk-mass eigenbasis, this
coupling is diagonal: g_s* . diag(F^2). To get the coupling in the mass basis,
we rotate:

```
(overlap)_{ij}^{mass basis} = (U^dagger . diag(F^2) . U)_{ij}
```

where U is the appropriate left or right rotation matrix from the SVD.

Four coupling matrices are built:

```
g_L^{down}[i,j] = g_s* . (U_L^d)^dagger . diag(F_Q^2) . U_L^d    [i,j]
g_R^{down}[i,j] = g_s* . (U_R^d)^dagger . diag(F_d^2) . U_R^d    [i,j]
g_L^{up}[i,j]   = g_s* . (U_L^u)^dagger . diag(F_Q^2) . U_L^u    [i,j]
g_R^{up}[i,j]   = g_s* . (U_R^u)^dagger . diag(F_u^2) . U_R^u    [i,j]
```

Each is a Hermitian 3x3 matrix (enforced by a numerical projection
(M + M^dagger)/2 to suppress round-off). The diagonal entries [i,i] are the
flavor-conserving couplings; the off-diagonal entries [i,j] with i != j drive
flavor-changing neutral currents.

Note: F_Q appears in both the LH-down and LH-up couplings because both u_L
and d_L live in the same SU(2) doublet Q, which has the same bulk mass and
therefore the same overlap profile. The difference between g_L^{down} and
g_L^{up} comes entirely from the different rotation matrices U_L^d vs U_L^u.

**On g_s***: The code uses g_s* = 3.0 (passed at scan time). This replaces
the perturbative QCD coupling g_s = sqrt(4 pi alpha_s(M_KK)). In RS, the
first KK-gluon has a wavefunction peaked near the IR brane, enhancing its
coupling to IR-localized fermions by a factor sqrt(2 k pi r_c) ~ 6-8 at tree
level. The one-loop corrected value is g_s* ~ 3 (Gedalia et al. 2009).

Since g_s* enters the Wilson coefficients as g_s*^2, the difference between
g_s* = 1 (perturbative) and g_s* = 3 is a factor of 9 in all constraint
ratios. This is a major assumption.

---

## Step 12: Tree-Level Wilson Coefficients

The KK-gluon is integrated out at the matching scale mu = M_KK, producing
four Delta F = 2 four-quark operators in the BMU basis (Buras, Misiak, Urban,
hep-ph/0005183). For a meson system involving quark flavors i and j:

```
prefactor = 1 / M_KK^2

C1_VLL =  (g_L[i,j])^2         . prefactor / 6
C1_VRR =  (g_R[i,j])^2         . prefactor / 6
C4_LR  = -(g_L[i,j] . g_R[i,j]) . prefactor
C5_LR  =  (g_L[i,j] . g_R[i,j]) . prefactor / 3
```

where g_L[i,j] and g_R[i,j] are the (complex) off-diagonal entries of the
coupling matrices from Step 11.

**The Wilson coefficients are complex**: Since g_L[i,j] and g_R[i,j] are
off-diagonal entries of Hermitian matrices, they are generally complex (not
real). The products g_L^2, g_R^2, and g_L*g_R are complex numbers. The
imaginary parts carry CP-violating information and are essential for the
epsilon_K constraint.

**Color factor origin**: These come from the color algebra of tree-level
single-KK-gluon exchange in the t-channel. Each vertex carries a color
generator T^a, giving the product:
  (T^a)_{alpha beta} (T^a)_{gamma delta} = 1/2 delta_{alpha delta}
  delta_{gamma beta} - 1/(2 N_c) delta_{alpha beta} delta_{gamma delta}
After Fierz rearrangement into the BMU basis:
- The VLL/VRR operators pick up 1/(2 N_c) = 1/6 for N_c = 3.
- The LR operators arise from Fierz-rearranging the (V-A)x(V+A) structure
  into scalar/pseudoscalar operator pairs. The specific coefficients -1 and
  +1/3 for C4_LR and C5_LR follow from the Fierz identity for the LR color
  contractions (see e.g., Buras, Weak Hamiltonian review, Eq. 6.6).

**Why only 4 operators**: Tree-level single vector boson exchange produces
only vector current structures at the vertices (gamma^mu P_L or P_R). This
gives VLL, VRR, and LR operators but NOT scalar LL/RR operators (SLL, SRR),
which would require scalar current insertions. The full BMU basis has 8
operators, but SLL/SRR are absent at tree level. They could appear at one
loop but are not included in this analysis.

**Flavor indices for each meson system**:

| System | Sector | (i,j) | Coupling matrices used |
|--------|--------|-------|------------------------|
| epsilon_K | down | (0,1) = d,s | g_L^{down}, g_R^{down} |
| Delta M_K | down | (0,1) = d,s | g_L^{down}, g_R^{down} |
| Delta M_{B_d} | down | (0,2) = d,b | g_L^{down}, g_R^{down} |
| Delta M_{B_s} | down | (1,2) = s,b | g_L^{down}, g_R^{down} |
| Delta M_{D^0} | up | (0,1) = u,c | g_L^{up}, g_R^{up} |

epsilon_K and Delta M_K use the same Wilson coefficients (same quark pair);
they differ only in which observable is computed from M_12 (Step 15).

**Assumption 12a**: Tree-level matching only. No one-loop corrections at the
matching scale (no threshold corrections, no scheme conversion beyond tree
level). This is the standard approximation in the RS FCNC literature.

**Assumption 12b**: Only the first KK-gluon contributes. Higher KK modes are
not summed. For the leading exclusion bound this is typically a ~10-20%
correction.

**Assumption 12c**: M_KK = Lambda_IR is used as both the propagator mass and
the matching scale. The physical first KK-gluon mass is m_g^{(1)} = 2.45 x
Lambda_IR (Step 16b corrects for this at plotting time).

---

## Step 13: QCD Running to the Hadronic Scale

The Wilson coefficients are evaluated at the matching scale M_KK (multi-TeV).
The hadronic matrix elements (Step 14) are evaluated at mu_had = 2 GeV. The
Wilson coefficients must be evolved from M_KK down to mu_had using the QCD
renormalization group.

### 13a: VLL and VRR operators (multiplicative running)

The VLL and VRR operators do not mix with any other operators under QCD. Their
anomalous dimension is gamma_VLL = 6(N_c - 1)/N_c = 4 for N_c = 3.

For a single segment with n_f active flavors between scales mu_high and mu_low:

```
C_VLL(mu_low) = C_VLL(mu_high) . [alpha_s(mu_low) / alpha_s(mu_high)]^{gamma_VLL / (2 beta_0)}
```

where beta_0 = (33 - 2 n_f) / 3 is the one-loop beta function coefficient.

Since gamma_VLL > 0 and alpha_s(mu_low) > alpha_s(mu_high), the ratio > 1
and the VLL coefficient is ENHANCED at lower scales.

### 13b: LR operators (2x2 matrix mixing)

The C4_LR and C5_LR operators mix under QCD running. The anomalous dimension
matrix is:

```
gamma_LR = [[ 8,      -4     ],
            [-16/3,  -28/3   ]]
```

For each n_f segment, the evolution matrix is:

```
U_LR = M . diag([alpha_s_low/alpha_s_high]^{d_i/(2 beta_0)}) . M^{-1}
```

where d_i are the eigenvalues of gamma_LR and M is the eigenvector matrix.
The dominant eigenvalue of gamma_LR is large and positive, giving a significant
enhancement (~3-5x) of the LR operators at low scales. This is the main reason
why epsilon_K is such a powerful constraint in RS models: the LR operators,
which receive the chirality-enhanced hadronic matrix elements (proportional to
(m_K/(m_s + m_d))^2 ~ 25), are further amplified by QCD running.

### 13c: Multi-segment evolution with flavor thresholds

The running from M_KK to mu_had = 2 GeV crosses two quark thresholds:

```
M_KK  -->  m_b = 4.18 GeV   (n_f = 5)
m_b   -->  m_c = 1.27 GeV   (n_f = 4)
m_c   -->  mu_had = 2 GeV   (n_f = 3)   [only if mu_had < m_c]
```

Wait — since mu_had = 2 GeV > m_c = 1.27 GeV, only the m_b threshold is
crossed. The segments are:

```
M_KK  -->  m_b = 4.18 GeV   (n_f = 5)
m_b   -->  2.0 GeV          (n_f = 4)
```

At each threshold, the Wilson coefficients are matched continuously (no finite
matching corrections at leading order).

The alpha_s running within each segment uses a one-loop formula:

```
alpha_s(mu) = alpha_s(mu_ref) / [1 + beta_0 alpha_s(mu_ref) / (2 pi) . ln(mu/mu_ref)]
```

with alpha_s(M_Z = 91.19 GeV) = 0.1179 as the reference value, and threshold
matching at m_b and m_c.

The cumulative evolution factors for VLL/VRR are multiplied across segments;
the 2x2 LR matrices are multiplied in order.

**Assumption 13a**: Leading-log anomalous dimensions only. The operator mixing
is LO (one-loop gamma). This is adequate for the current precision.

**Assumption 13b**: The alpha_s running within the QCD RG module
(qcd_running.py) uses a one-loop formula with alpha_s(M_Z) = 0.1179 as
reference. The separate qcd/ module (used for the coupling computation in
Step 11) uses 4-loop running with alpha_s(M_Z) = 0.1180. This introduces
an O(few-percent) inconsistency between the coupling value at M_KK and the
evolution from M_KK. In practice this is subdominant compared to the g_s*
uncertainty.

**Assumption 13c**: Continuous threshold matching (no finite corrections at
m_b or m_c). This is the standard LO prescription.

---

## Step 14: Hadronic Matrix Elements

The RG-evolved Wilson coefficients at mu_had = 2 GeV are contracted with
hadronic matrix elements to obtain M_12^{NP}, the new-physics contribution to
the meson mixing amplitude.

The general formula for a pseudoscalar meson P with valence quarks q_1, q_2 is:

```
M_12^NP = C1_VLL . ME_VLL + C1_VRR . ME_VLL + C4_LR . ME_LR4 + C5_LR . ME_LR5
```

Note: VRR has the same matrix element as VLL by parity symmetry of QCD. This
means <P-bar|O_VRR|P> = <P-bar|O_VLL|P>.

The matrix elements are:

```
ME_VLL = (2/3) . f_P^2 . m_P . B_1

r_chi = (m_P / (m_q1 + m_q2))^2

ME_LR4 = (r_chi/6 + 1/4) . f_P^2 . m_P . B_4
ME_LR5 = (r_chi/2 + 1/12) . f_P^2 . m_P . B_5
```

The factor r_chi is the chiral enhancement factor. For kaons:
r_chi = (0.49761 / (0.0934 + 0.00467))^2 ~ 25.2.
This is why the LR operators dominate epsilon_K: their matrix elements are
enhanced by ~25 relative to the VLL matrix element.

### Hadronic inputs by meson system

**Kaon (K^0)** — valence quarks d, s:

| Input | Value | Source |
|-------|-------|--------|
| f_K | 155.7 MeV | PDG 2024 |
| m_K | 497.61 MeV | PDG 2024 |
| m_s(2 GeV) | 93.4 MeV | FLAG 2024 |
| m_d(2 GeV) | 4.67 MeV | FLAG 2024 |
| B_1^K | 0.717 | FLAG/lattice average |
| B_4^K | 0.78 | ETM 2013 |
| B_5^K | 0.57 | ETM 2013 |

**B_d** — valence quarks d, b:

| Input | Value | Source |
|-------|-------|--------|
| f_{B_d} | 190.0 MeV | FLAG 2024 |
| m_{B_d} | 5279.72 MeV | PDG 2024 |
| m_b(m_b) | 4.18 GeV | PDG (MS-bar) |
| m_d(2 GeV) | 4.67 MeV | FLAG 2024 |
| B_1^{B_d} | 0.87 | FLAG 2024 |
| B_4^{B_d} | 1.02 | FLAG 2024 |
| B_5^{B_d} | 0.96 | FLAG 2024 |

**B_s** — valence quarks s, b:

| Input | Value | Source |
|-------|-------|--------|
| f_{B_s} | 230.3 MeV | FLAG 2024 |
| m_{B_s} | 5366.92 MeV | PDG 2024 |
| m_s(2 GeV) | 93.4 MeV | FLAG 2024 |
| B_1^{B_s} | 0.87 | FLAG 2024 |
| B_4^{B_s} | 1.02 | FLAG 2024 |
| B_5^{B_s} | 0.96 | FLAG 2024 |

**D^0** — valence quarks u, c:

| Input | Value | Source |
|-------|-------|--------|
| f_D | 212.0 MeV | FLAG 2024 |
| m_{D^0} | 1864.84 MeV | PDG 2024 |
| m_c(m_c) | 1.27 GeV | PDG (MS-bar) |
| m_u(2 GeV) | 2.16 MeV | FLAG 2024 |
| B_1^D | 0.75 | Lattice (less precise) |
| B_4^D | 1.0 | Estimated |
| B_5^D | 1.0 | Estimated |

**Assumption 14a**: Bag parameters B_4^D and B_5^D for D^0 mixing are set to
1.0. No reliable lattice determination exists. This introduces O(1) uncertainty
in the D^0 constraint.

**Assumption 14b**: The quark masses entering the chiral enhancement factor
r_chi are MS-bar masses at their respective self-scales, NOT all at mu = 2 GeV.
For B mesons, m_b(m_b) = 4.18 GeV is used rather than m_b(2 GeV) ~ 5.4 GeV;
for D^0, m_c(m_c) = 1.27 GeV. This is inconsistent with evaluating the matrix
elements at mu_had = 2 GeV. Using m_b(2 GeV) instead would change r_chi for
B_d from ~1.6 to ~1.0, an O(40%) shift in the LR matrix element. For the B
and D systems the LR contribution is subdominant (r_chi ~ 1-2, vs ~25 for
kaons), so this inconsistency has limited impact on the final bounds, but it
is a known systematic.

**Assumption 14c**: The matrix element formulae above are the "vacuum
insertion approximation" corrected by the bag parameters B_i. The B_i
encode deviations from the naive factorization <P|O|P> = <P|qq|0><0|qq|P>.
Lattice QCD determines the B_i.

**Normalization convention**: M_12^NP as computed here is the FULL off-diagonal
Hamiltonian matrix element <P-bar|H_eff|P>, with dimensions of GeV. It does
NOT include a 1/(2 m_P) normalization that some references fold into the
definition of M_12. Consistently, the mass difference is Delta m = 2 |M_12|,
so the constraint |M_12^NP| <= Delta m / 2 is equivalent to |M_12^NP| <=
|M_12^exp|.

---

## Step 15: Observable Evaluation and Acceptance

### 15a: epsilon_K (CP violation in kaon mixing)

The NP contribution to epsilon_K is:

```
epsilon_K^NP = |kappa_epsilon / (sqrt(2) . Delta m_K)| . |Im(M_12^NP)|
```

with:
- kappa_epsilon = 0.94 (accounts for long-distance and higher-order
  corrections, Buras et al.)
- Delta m_K = 3.484 x 10^{-15} GeV (K_L - K_S mass difference, PDG)

The NP budget is the gap between the experimental value and the SM prediction:

```
budget = |epsilon_K^exp - epsilon_K^SM| = |2.228e-3 - 1.81e-3| = 4.18e-4
```

Ratio to bound: epsilon_K^NP / budget. **Passes if ratio <= 1.0.**

**Assumption 15a-i**: The SM prediction epsilon_K^SM = 1.81 x 10^{-3} is taken
from the CKMfitter central value. The experimental value is 2.228 x 10^{-3}
(PDG). The difference of 4.18 x 10^{-4} is the allowed "room" for NP. This is
a generous budget — some analyses use a tighter SM prediction.

**Assumption 15a-ii**: The formula takes |Im(M_12^NP)|, discarding the sign
of the NP CP-violating phase. This means we do not account for constructive vs
destructive interference between NP and SM. If the NP phase happens to
partially cancel the SM contribution, the actual constraint would be weaker;
if it enhances it, stronger. Using the absolute value is a standard
simplification for model-independent bounds.

### 15b: Delta M_K (kaon mass splitting)

The NP contribution to the kaon mass difference is:

```
|M_12^NP|_K / (Delta m_K / 2) = ratio
```

Budget = Delta m_K / 2 = 1.742 x 10^{-15} GeV.

The constraint is that |M_12^NP| should not exceed |M_12^exp| = Delta m_K / 2.
This is conservative because the SM contribution to Delta m_K is dominated by
poorly-known long-distance effects, so the entire experimental value is treated
as a budget for NP (rather than requiring |M_12^NP| < |M_12^exp - M_12^SM|,
which would need a reliable SM prediction). **Passes if ratio <= 1.0.**

### 15c: B_d, B_s, D^0 mixing (mass differences)

For each system, M_12^NP is computed using the system-specific hadronic
parameters. The budget definitions are:

| System | Budget formula | Numerical value |
|--------|----------------|-----------------|
| B_d | max(Delta m^exp / 2, \|Delta m^exp - Delta m^SM\| / 2) | max(1.667e-13, 1.33e-14) = **1.667 x 10^{-13} GeV** |
| B_s | max(Delta m^exp / 2, \|Delta m^exp - Delta m^SM\| / 2) | max(5.844e-12, 6.0e-15) = **5.844 x 10^{-12} GeV** |
| D^0 | Delta m^exp / 2 (SM long-distance dominated) | **3.125 x 10^{-15} GeV** |

The max(...) formula is conservative: it uses whichever is larger — half the
experimental value, or half the experiment-SM difference. For B_d and B_s, the
experimental value dominates (the SM predictions are close to experiment).

Ratio to bound: |M_12^NP| / budget. **Passes if ratio <= 1.0.**

### 15d: Experimental inputs for budgets

| System | Delta m^exp | Delta m^SM |
|--------|-------------|------------|
| K | 3.484 x 10^{-15} GeV | (enters only via epsilon_K) |
| B_d | 3.334 x 10^{-13} GeV | 3.6 x 10^{-13} GeV (CKMfitter) |
| B_s | 1.1688 x 10^{-11} GeV | 1.17 x 10^{-11} GeV (CKMfitter) |
| D^0 | 6.25 x 10^{-15} GeV | Long-distance dominated (HFLAV) |

### 15e: Overall acceptance

A scan point is **accepted** if the ratio to bound is <= 1.0 for ALL FIVE
systems simultaneously: epsilon_K, Delta M_K, Delta M_{B_d}, Delta M_{B_s},
Delta M_{D^0}.

Additionally, the fit must have converged with fit_score <= 0.1 (meaning the
quark masses and CKM are reproduced well enough that the Wilson coefficients
are meaningful).

**Note on the code structure**: The legacy deltaf2 module has 4 input entries
(epsilon_K, B_d, B_s, D^0). Delta M_K was added as a 5th system in the modern
scan pipeline, which evaluates it separately via evaluate_delta_mk_with_running.
The JSONL output from the modern pipeline includes all 5 systems in the
ratio_to_bound_by_system dict, which is what the plotting script reads.

---

## Step 16: Post-Processing for Figures

The raw scan stores results in JSONL format (one JSON object per scan point).
Before plotting, several transformations are applied.

### 16a: Quality filter

Points with fit_score > 0.1 are discarded entirely. A poor fit means the
quark masses and/or CKM were not reproduced, so the zero-mode profiles,
rotation matrices, and Wilson coefficients are unreliable.

Points where the fit did not converge (fit_success = False or
fit_converged = False) are also discarded.

### 16b: Mass-axis convention (Lambda_IR -> m_g^{(1)})

The scan uses M_KK = Lambda_IR (xi_KK = 1) as the propagator mass in all
Wilson coefficient formulas (Step 12). The physical first KK-gluon mass is
larger:

```
m_g^{(1)} = xi_kk . Lambda_IR
```

where xi_kk = 2.448687 is the first root of the gauge-sector KK mass
eigenvalue equation with Neumann-Neumann boundary conditions on a
warped interval.

To correctly relabel the axis, we must account for the fact that the scan used
Lambda_IR (not m_g^{(1)}) in the 1/M_KK^2 propagator. Since the Wilson
coefficients scale as 1/M_KK^2, a point at Lambda_IR with ratio R should have
ratio R / xi_kk^2 when interpreted as being at m_g^{(1)} = xi_kk . Lambda_IR:

```
axis_scale = xi_kk / xi_scan       (= 2.448687 / 1.0 = 2.448687)
M_KK_plot  = Lambda_IR . axis_scale
ratio_plot = ratio_stored / axis_scale^2
```

The net effect: a point that was marginal (ratio = 1) in the raw scan
becomes well within bounds (ratio = 1/xi_kk^2 ~ 1/6) after the correction,
because the physical propagator mass is 2.45x larger than what the scan used.
The exclusion boundary shifts to a lower Lambda_IR (where the raw ratio was
~xi_kk^2 ~ 6). In terms of m_g^{(1)}, this lower Lambda_IR multiplied by
xi_kk partially compensates, so the boundary m_g^{(1)} is similar to (but
not exactly equal to) the raw-scan Lambda_IR boundary.

### 16c: Grid aggregation (common best point)

At each unique (r, M_KK) grid point, there are up to 20 different
overall_scale values (from the scan grid). The plotter selects the single
overall_scale that minimizes:

```
max(ratio_epsilon_K, ratio_K, ratio_{B_d}, ratio_{B_s}, ratio_{D^0})
```

across all five systems simultaneously. All five system ratios are then read
from that same point.

**Why this matters**: An earlier version of the code picked the best
overall_scale independently for each system. That was overly optimistic because
the five contour lines came from five different Yukawa structures, not one
consistent MFV point. With the common best point, all contours are physically
consistent.

**What overall_scale does here**: A larger overall_scale means larger Yukawa
entries, which pushes the bulk masses c_Q toward c_ir = 0.30 (more IR-localized)
and increases the overlap F_Q. This generally increases BOTH the quark masses
AND the flavor-changing couplings. The fit compensates by adjusting the singular
values. The net effect of sweeping overall_scale at fixed (r, Lambda_IR) is to
explore different tradeoffs between the diagonal masses and the off-diagonal
FCNCs. Picking the best overall_scale finds the tradeoff that minimizes FCNCs
while still reproducing the quark masses.

---

## Step 17: Figure 1 — System-by-System Exclusion Boundaries

**What it shows**: Contour lines in the (r, m_g^{(1)}) plane where each meson
system's ratio_to_bound = 1, with green/red shading for the combined
allowed/excluded regions.

**Construction**:
1. Take the aggregated data from Step 16c (one best point per (r, M_KK)).
2. Convert M_KK from GeV to TeV.
3. Work in log-log space: compute log10(r) and log10(M_KK/TeV).
4. For each of the five systems and the combined max-ratio, interpolate
   log10(ratio) onto a regular 300 x 300 grid in (log10 r, log10 M_KK)
   using scipy.interpolate.griddata with **cubic** interpolation. Values
   below 10^{-10} are clipped to 10^{-10} before taking log10.
5. For each system, draw a contour at log10(ratio) = 0 (ratio = 1) with the
   system-specific color and linestyle.
6. Shade the combined allowed/excluded regions using the max-ratio
   interpolation: green where log10(max_ratio) < 0 (all systems pass), red
   where log10(max_ratio) > 0 (at least one system fails).
7. Axis ticks are formatted to show physical values (not log10 values).

**Contour line styles**:
| System | Color | Style |
|--------|-------|-------|
| epsilon_K | Purple (#7B2D8E) | Solid |
| Delta M_K | Blue (#2166AC) | Dashed |
| Delta M_{B_d} | Orange (#D95F02) | Dash-dot |
| Delta M_{B_s} | Red (#C0392B) | Dotted |
| Delta M_{D^0} | Green (#1B7837) | Dash-dot-dot |

**Interpretation**: At a given r, any m_g^{(1)} above all contour lines is
allowed by all five meson constraints. The contour that sits highest at each r
is the binding constraint. Typically epsilon_K dominates for r ~ 0.1-0.5, and
Delta M_K takes over at large r.

**Important caveat**: This is a best-case envelope. At each grid point, the
fitter found the single most favorable Yukawa structure, and the plotter
picked the most favorable overall_scale. The contour shows the lowest M_KK at
which any fitted point survives — it is an existence bound, not a statement
about typical NMFV behavior.

---

## Step 18: Figure 2 — M_KK Lower Bound: 2007 vs Modern

**What it shows**: For each r value, the minimum m_g^{(1)} among all accepted
points, comparing modern (2024+) constraints to rescaled 2007-era constraints.

**Construction**:
1. For each unique r value in the scan grid, collect all (M_KK, ratios) pairs
   across all Lambda_IR and overall_scale values.
2. For the **modern** bound: check each pair's max(all 5 ratios) <= 1.
   Record the smallest M_KK among passing points. If no point passes at any
   M_KK, mark this r as fully excluded.
3. For the **2007** bound: multiply each system's ratio by the corresponding
   BOUND_RATIO factor (which makes the constraint looser, since modern bounds
   are tighter). Then apply the same acceptance check and find the minimum
   passing M_KK.
4. Plot both curves as a function of r.
5. Shade the region between the two curves to visualize the improvement since
   2007.

**Note on consistency with Figure 1**: Figure 1 uses _aggregate_best_scale to
pick a common best overall_scale at each (r, M_KK) grid point before drawing
contours. Figure 2 uses _compute_min_mkk_by_r, which iterates over individual
scan rows and checks max(all 5 ratios) <= 1 per row — so it does NOT mix
different overall_scale values across systems. The difference from Figure 1 is
that Figure 2 does not first collapse each (r, M_KK) grid point to a single
common best-scale representative; it considers every stored row individually
when searching for the minimum passing M_KK.

**2007 rescaling factors** (BOUND_RATIOS):

These factors represent the ratio of the modern experimental bound to the
2007-era bound. Multiplying the modern ratio by this factor gives the
2007-equivalent ratio:

| System | Factor | Origin | Interpretation |
|--------|--------|--------|----------------|
| epsilon_K | 0.70 | Tighter lattice + CKM inputs | Modern 30% tighter |
| Delta M_K | 0.70 | Same improvement as epsilon_K | Modern 30% tighter |
| Delta M_{B_d} | 0.67 | Better lattice + exp | Modern 33% tighter |
| Delta M_{B_s} | 0.37 | Huge improvement in B_s measurements | Modern 63% tighter |
| Delta M_{D^0} | 0.106 | D^0 mixing barely measured in 2007 | Modern ~90% tighter |

Example: if a point has modern B_s ratio = 0.9 (barely passing), its 2007
equivalent is 0.9 x 0.37 = 0.33 (easily passing). So the modern constraint
is much more stringent for B_s mixing.

**Interpretation**: The purple shaded band between the two curves shows how
much constraints have tightened since 2007. At small r (where epsilon_K
dominates), the improvement is modest (~30%). At large r (where B_s and D^0
dominate), the improvement is dramatic.

---

## Summary of Key Assumptions

1. **Geometry**: Slice of AdS_5 with k = M_Pl = 1.2209 x 10^19 GeV,
   Higgs on the IR brane, no backreaction.

2. **Yukawa parametrization**: SVD with right rotations fixed to identity.
   This forces C_u, C_d to be diagonal. Only C_Q has off-diagonal structure.
   The full Yukawa space is not explored.

3. **Bulk mass construction**: C_Q = r Y_u Y_u^dagger + Y_d Y_d^dagger with
   no explicit alpha_Q identity term. The BulkMassMap (sigmoid mapping
   eigenvalues to c in [0.30, 0.72]) substitutes for alpha_Q nonlinearly.

4. **Fit procedure**: 14-parameter Trust Region Reflective (TRF) optimization with
   max 120 evaluations. Finds the best Yukawa structure at each scan point,
   not a representative sample. Each point starts from the same deterministic
   seed.

5. **KK-gluon coupling**: g_s* = 3.0 (one-loop enhanced, Gedalia et al. 2009).
   Enters as g_s*^2 in all Wilson coefficients. The tree-level estimate is
   g_s* ~ 6; the perturbative QCD value is g_s ~ 1. This single parameter
   changes all bounds by a factor of up to 36.

6. **Matching**: Tree-level only, at scale M_KK = Lambda_IR. Only the first
   KK-gluon mode. No loop corrections.

7. **QCD running**: LO anomalous dimensions for Delta F=2 operator mixing,
   one-loop alpha_s running with threshold matching at m_b and m_c.
   Hadronic evaluation scale mu_had = 2 GeV.

8. **Hadronic inputs**: FLAG 2024 / PDG 2024 for bag parameters and decay
   constants. D^0 bag parameters B_4, B_5 estimated at 1.0. Kaon B_4, B_5
   from ETM 2013.

9. **Acceptance budgets**: epsilon_K uses the exp-SM gap (4.18 x 10^{-4}).
   Delta M_K uses half the experimental value. B_d, B_s use
   max(exp/2, |exp-SM|/2). D^0 uses exp/2.

10. **Axis convention**: Publication figures use m_g^{(1)} = 2.4487 Lambda_IR
    (first gauge KK eigenvalue with NN BCs, from Bessel root) with a
    corresponding 1/xi_kk^2 ratio correction.

11. **Best-case envelope**: Both figures show the most favorable point at each
    (r, M_KK). This is an existence bound ("NMFV can allow M_KK this low"),
    not a characterization of typical NMFV behavior. A sampling/histogram
    analysis is needed for the latter.

12. **Target quark masses**: The fit targets (m_u = 1.3 MeV, m_c = 620 MeV,
    m_t = 172 GeV, m_d = 2.8 MeV, m_s = 57 MeV, m_b = 2.86 GeV) are
    repo-owned values, not standard PDG pole or MS-bar masses at a canonical
    scale. They are chosen for scan stability.
