**Design Doc: `\mu \to e\gamma` In This RS Framework**

No code changes were made. I read the repo implementation first: [flavorConstraints/muToEGamma.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavorConstraints/muToEGamma.py), [flavor_catalog_constraints/physics_adapters/lepton.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/physics_adapters/lepton.py), and [flavor_catalog_constraints/primary/charged_lepton/L001.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/charged_lepton/L001.py).

One correction up front: the canonical Agashe-Blechman-Petriello paper is `hep-ph/0606021`, not `hep-ph/0606293`.

**1. Current Repo State**

The repo currently implements the Perez-Randall lepton-MFV NDA proxy, not a generic anarchic RS one-loop dipole.

The core formula in `flavorConstraints/muToEGamma.py` is:

```text
BR(mu -> e gamma) ~= 4e-8 * |(Ybar_N Ybar_N^\dagger)_{e mu}|^2 * (3 TeV / M_KK)^4
```

with the paper-era bound encoded as:

```text
|(Ybar_N Ybar_N^\dagger)_{e mu}| < C * (M_KK / 3 TeV)^2,   C = 0.02
```

This matches the Perez-Randall NDA structure. The adapter and `L001` catalog constraint already flag that the current `ParameterPoint` does not yet carry the full mass-basis lepton coupling object needed for a true one-loop charged-lepton dipole. The existing implementation needs only `Y_N`, PMNS, and an `M_KK` convention; it does not use a full anarchic charged-lepton Yukawa matrix `Y_E`, separate `c_{L_i}`, `c_{E_i}`, diagonalization matrices, KK fermion spectra, Higgs localization, or custodial representation data.

**2. Literature Survey**

**Agashe, Blechman, Petriello, “Probing the Randall-Sundrum geometric origin of flavor with lepton flavor violation”, hep-ph/0606021**  
Source: https://arxiv.org/abs/hep-ph/0606021

ABP is the early canonical RS lepton LFV analysis. They treat tree-level LFV and loop-induced dipoles. The effective dipole is written as

```text
L_eff ⊃ C_R/(2 m_mu) * \bar e_R sigma^{mu nu} F_{mu nu} mu_L
      + C_L/(2 m_mu) * \bar e_L sigma^{mu nu} F_{mu nu} mu_R .
```

Their important conclusion is model-dependent calculability: for a Higgs exactly localized on the IR brane, the dipole is UV sensitive/cutoff dependent; for a bulk Higgs, the Higgs/KK-fermion loop is calculable. The relevant loop has internal KK fermions and Higgs fields. Its flavor structure is not aligned with the charged-lepton mass matrix, so it generates LFV after rotation to the mass basis.

For a bulk Higgs, ABP find a schematic structure of the form

```text
C_L ~ e m_mu /(32 pi^2 M_KK^2) * [Delta_R Delta_2 Delta_L]_{e mu},
```

where the `Delta` matrices encode overlap and Yukawa insertions. In anarchic RS language, this is the charged-lepton-Yukawa-driven Higgs/KK-fermion dipole. ABP also emphasize the classic RS tension: tree-level LFV falls with increasing anarchic Yukawa size because the zero-mode wavefunctions become more UV localized, while dipoles grow with Yukawa size.

**Csaki, Falkowski, Weiler, “The Flavor of the Composite Pseudo-Goldstone Higgs”, 0804.1954**  
Source: https://arxiv.org/abs/0804.1954

CFW is mostly a quark-sector/anarchic-RS flavor paper, but it is canonical for the RS-GIM logic that also controls leptons. They show how off-diagonal gauge couplings scale with fermion profile factors:

```text
g_{ij}^{KK} ~ f_i f_j
```

after rotation to the mass basis. They conclude that fully anarchic warped flavor is highly constrained unless there are flavor symmetries or moderate tuning. For this project, CFW is useful less as a `\mu -> e gamma` calculator and more as the reference for how overlap factors and anarchic rotations suppress flavor violation.

**Blanke, Buras, Duling, Gori, Weiler, “Delta F = 2 observables and fine-tuning in a warped extra dimension with custodial protection”, 0809.1073**  
Source: https://arxiv.org/abs/0809.1073

Blanke et al are the canonical full-scan custodial RS flavor reference. They focus on quark `\Delta F=2` and rare processes, but the computational lesson is directly relevant: once realistic rotations, custodial structure, and full gauge sectors are included, naive estimates can move by order-one to order-of-magnitude factors, and accidental cancellations/tunings matter. For `\mu -> e gamma`, this supports designing an explicit mass-basis interface rather than treating the NDA number as a universal RS bound.

**Perez & Randall, “Natural Neutrino Masses and Mixings from Warped Geometry”, 0805.4652**  
Source: https://arxiv.org/abs/0805.4652

This is exactly the formula currently coded. Perez-Randall impose a lepton-MFV-like structure in which charged-lepton flavor violation vanishes without the neutrino Yukawa spurion. In the charged-lepton mass basis,

```text
Y_N -> V_PMNS diag(Y_N)
S_N = Ybar_N Ybar_N^\dagger .
```

They estimate

```text
BR(mu -> e gamma)_IR
  ~ 4e-8 * |(Ybar_N Ybar_N^\dagger)_{e mu}|^2 * (3 TeV / M_KK)^4 .
```

The dominant chirality has the structure

```text
C_L ∝ m_mu * (Y_N Y_N^\dagger)_{e mu} / M_KK^2,
C_R ~ (m_e / m_mu) C_L .
```

The zero-mode overlap factors combine with the charged-lepton Yukawa insertion to reproduce `m_mu`; the flavor violation is entirely in `Y_N Y_N^\dagger`. This is not the generic charged-lepton-anarchy dipole.

**Chen & Yu, “Minimal Flavor Violation in the Lepton Sector of the Randall-Sundrum Model”, 0804.2503**  
Source: https://arxiv.org/abs/0804.2503

Chen-Yu study lepton MFV in RS. The key point is aligned with Perez-Randall: if the charged-lepton sector is flavor aligned and neutrino Yukawas are the only LFV spurion, then dangerous tree-level charged LFV can be heavily suppressed, and `\mu -> e gamma` is controlled by neutrino spurions. This is a different theory assumption from generic anarchic charged-lepton Yukawas.

**Csaki, Grossman, Tanedo, Tsai, “Warped penguin diagrams”, 1004.2037**  
Source: https://arxiv.org/abs/1004.2037

This paper computes RS dipole/penguin diagrams using a mixed position/momentum formalism for an IR-brane Higgs. The leading dipole operator is written schematically as

```text
R'^2 * e/(16 pi^2) * v/sqrt(2)
* f_{L_i} [ a_{k l} Y_{i k} Y^\dagger_{k l} Y_{l j}
            + b_{ij} Y_{ij} ]
* f_{E_j}
* \bar L_i sigma^{mu nu} E_j F_{mu nu}.
```

This cleanly separates the two generic structures:

```text
Higgs / KK-fermion:   f_L [Y_E Y_E^\dagger Y_E] f_E
Gauge / aligned:      f_L [Y_E] f_E
```

The single-Yukawa piece is aligned with the mass matrix in the universal limit and generates LFV only through misalignment from non-universal profiles and full rotations. The three-Yukawa Higgs term is the dominant anarchic dipole for moderate or large `Y_*`.

**Agashe et al, “Warped Dipole Completed, with a Tower of Higgs Bosons”, 1412.6468**  
Source: https://arxiv.org/abs/1412.6468

Agashe-Azatov-Cui-Randall-Son show that in bulk/narrow-bulk Higgs RS, KK Higgs modes are finite, unsuppressed, and comparable to the SM-Higgs contribution. This matters because a complete one-loop RS dipole is not just “zero-mode Higgs plus first KK fermion.” The KK Higgs tower can shift the coefficient by order one and is required for a 5D-covariant calculation. Their conclusion reinforces that a true explicit calculation needs a Higgs-localization prescription.

**Beneke, Moch, Rohrwild, “Lepton flavour violation in RS models with a brane- or nearly brane-localized Higgs”, 1508.01705**  
Source: https://arxiv.org/abs/1508.01705

This is the most complete practical reference for charged-lepton dipoles in minimal and custodially protected RS. It gives the dimension-six operators

```text
L_eff ⊃ 1/T^2 [
  a^B_{ij} (\bar L_i sigma^{mu nu} E_j) Phi B_{mu nu}
+ a^W_{ij} (\bar L_i tau^A sigma^{mu nu} E_j) Phi W^A_{mu nu}
] + h.c.
```

and after EWSB/mass rotation,

```text
alpha^A = U_L^\dagger (c_W a^B - s_W a^W) U_R .
```

The low-energy amplitude is

```text
L = m_mu A_R \bar e sigma^{mu nu} F_{mu nu} P_R mu
  + m_mu A_L \bar e sigma^{mu nu} F_{mu nu} P_L mu + h.c.

BR(mu -> e gamma)
  = m_mu^5 / (4 pi Gamma_mu) * (|A_L|^2 + |A_R|^2).
```

Beneke et al identify three important contributions: Higgs exchange, gauge-boson exchange, and loop-induced effects from other dimension-six operators. The dominant anarchic contribution for moderate/large `Y_*` is the Higgs/KK-fermion term,

```text
a^H_{ij} ~ e/(192 pi^2) * f_{L_i} [Y_E Y_E^\dagger Y_E]_{ij} f_{E_j},
```

up to model-dependent minimal/custodial and Higgs-localization coefficients. The gauge contribution is closer to

```text
a^g_{ij} ~ f_{L_i} Y_{E,ij} f_{E_j} * A_{ij},
```

and would be aligned if `A_{ij}` were universal. Its LFV part is suppressed but provides a floor at small `Y_*`.

Their useful estimates are:

```text
BR(mu -> e gamma)_Higgs ~ 5e-9 * Y_*^4 * (1 TeV / T)^4

BR(mu -> e gamma)_gauge ~ 0.5e-11 * (1 TeV / T)^4   minimal RS
                         ~ 5e-11   * (1 TeV / T)^4   custodial RS
```

These are charged-lepton-anarchy estimates, not Perez-Randall neutrino-spurion estimates.

**Custodial RS References**

The foundational custodial symmetry references are Agashe-Delgado-May-Sundrum, https://arxiv.org/abs/hep-ph/0308036, and Agashe-Contino-Da Rold-Pomarol, https://arxiv.org/abs/hep-ph/0605341. For lepton LFV specifically, custodial embeddings and their consequences are treated in Csaki et al 1004.2037 and Beneke et al 1508.01705. Custodial protection changes the gauge-sector contribution and can add additional Yukawa structures, but it does not eliminate the need to compute the dipole in the charged-lepton mass basis.

**3. Dipole Structure**

Use the mass-basis electromagnetic dipole convention:

```text
L_eff =
m_mu A_R \bar e sigma^{mu nu} F_{mu nu} P_R mu
+ m_mu A_L \bar e sigma^{mu nu} F_{mu nu} P_L mu + h.c.

BR(mu -> e gamma)
= m_mu^5 / (4 pi Gamma_mu) * (|A_L|^2 + |A_R|^2).
```

At the RS matching scale, the SMEFT-like form is:

```text
L_eff ⊃ 1/T^2 [
  a^B_{ij} (\bar L_i sigma^{mu nu} E_j) Phi B_{mu nu}
+ a^W_{ij} (\bar L_i tau^A sigma^{mu nu} E_j) Phi W^A_{mu nu}
] + h.c.
```

After EWSB and rotations,

```text
alpha^A = U_L^\dagger (c_W a^B - s_W a^W) U_R .
```

For `\mu -> e gamma`, the relevant entries are approximately:

```text
A_R ∝ alpha^A_{e mu} * v / (sqrt(2) T^2 m_mu),
A_L ∝ alpha^{A *}_{mu e} * v / (sqrt(2) T^2 m_mu).
```

The flavor spurion depends on the lepton model.

**Generic anarchic charged-lepton RS**

The dominant Higgs/KK-fermion structure is:

```text
a^H_{ij} ~ e/(16 pi^2) * f_{L_i}
          [Y_E Y_E^\dagger Y_E]_{ij}
          f_{E_j}
          * O(1 model coefficient).
```

The zero-mode charged-lepton mass matrix is schematically:

```text
m^E_{ij} ~ v * f_{L_i} Y^E_{ij} f_{E_j}.
```

Thus the dipole has the same external profile factors as the mass matrix, but a different internal Yukawa structure. After diagonalizing the mass matrix, the three-Yukawa term is misaligned and generates LFV. In anarchic estimates this gives the familiar scaling:

```text
A(mu -> e gamma)_Higgs ∝ Y_*^2 / T^2,
BR ∝ Y_*^4 / T^4.
```

Gauge/KK-vector loops instead have a mostly single-Yukawa structure:

```text
a^g_{ij} ~ e/(16 pi^2) * f_{L_i} Y^E_{ij} f_{E_j} * A_{ij}.
```

If `A_{ij}` is universal, this is aligned with the mass matrix and gives no LFV after rotation. LFV comes from non-universal profile-dependent pieces, so it is typically smaller than the Higgs term but nonzero.

**Perez-Randall / lepton-MFV RS**

In the Perez-Randall setup, the charged-lepton sector is aligned and the LFV spurion is neutrino-Yukawa driven:

```text
S_N = Ybar_N Ybar_N^\dagger,
Ybar_N = 2 k Y_N.
```

In the charged-lepton mass basis,

```text
Y_N = V_PMNS diag(Y_N1, Y_N2, Y_N3),
S_N = V_PMNS diag(Ybar_N1^2, Ybar_N2^2, Ybar_N3^2) V_PMNS^\dagger.
```

The dominant chirality is:

```text
C_L ∝ m_mu * (S_N)_{e mu} / M_KK^2,
C_R ~ (m_e / m_mu) C_L.
```

The overlap factors enter as

```text
f_L * Y_N Y_N^\dagger * Y_E * f_E,
```

but `Y_E f_L f_E` is fixed by `m_mu`, so the final dominant chirality is controlled by `m_mu S_N`, not by an additional `sqrt(m_e/m_mu)` suppression.

**When Each Spurion Dominates**

`Y_N Y_N^\dagger` dominates only if the charged-lepton sector is aligned or protected so that charged-lepton-Yukawa LFV is absent/subleading. That is the Perez-Randall / lepton-MFV assumption and is what the repo currently implements.

`Y_E Y_E^\dagger Y_E` dominates in generic anarchic charged-lepton RS with an IR or narrow-bulk Higgs, especially for `Y_* ~ 1-3`. This is the standard explicit one-loop result in the later dipole literature.

If both sectors are present and not symmetry-separated, the physical amplitude is a coherent sum:

```text
A_{e mu} = A^{Higgs/Y_E}_{e mu}
         + A^{gauge/Y_E}_{e mu}
         + A^{neutrino/Y_N}_{e mu}
         + ...
```

so the code should report individual pieces and the total only after a physicist fixes phases, basis, and model assumptions.

**4. Theoretical Comparison To Perez-Randall NDA**

The repo’s NDA estimate agrees with Perez-Randall in three important ways:

```text
BR ∝ |S_N|^2,
S_N = (Ybar_N Ybar_N^\dagger)_{e mu},
BR ∝ M_KK^{-4}.
```

The paper-era coefficient follows from

```text
BR ~= 4e-8 * |S_N|^2 * (3 TeV / M_KK)^4.
```

Using the old `1.2e-11` bound gives

```text
|S_N| < sqrt(1.2e-11 / 4e-8) * (M_KK / 3 TeV)^2
      ~= 0.017 * (M_KK / 3 TeV)^2,
```

which is the repo’s rounded `C = 0.02`.

Where NDA and explicit calculations differ:

```text
Perez-Randall NDA:
  S_N = Ybar_N Ybar_N^\dagger.
  Assumes lepton-MFV/aligned charged leptons.
  Unknown O(1) brane coefficient absorbed into 4e-8.

Generic anarchic one-loop:
  Dominant spurion is usually Y_E Y_E^\dagger Y_E.
  Gauge loops add a smaller but irreducible floor.
  Exact coefficient depends on Higgs localization, KK Higgs tower,
  minimal vs custodial representation, and full mass-basis rotations.
```

Therefore the Perez-Randall NDA is not uniformly conservative or aggressive. It is model-specific.

For example, with `M_KK = 3 TeV` and `|S_N| = 0.01`, Perez-Randall gives roughly

```text
BR ~ 4e-12.
```

Beneke’s generic charged-lepton-anarchy Higgs estimate gives

```text
BR_H ~ 5e-9 * Y_*^4 * (1 TeV / T)^4.
```

For `T = 3 TeV`, this is roughly

```text
BR_H ~ 6e-11 * Y_*^4.
```

So for `Y_* ~ 1`, the generic charged-lepton-anarchy Higgs contribution can already exceed the Perez-Randall example by an order of magnitude; for `Y_* ~ 3`, it can be orders of magnitude larger. Conversely, if the model is truly lepton-MFV with small `Y_N`, the generic charged-lepton contribution is absent by assumption and the Perez-Randall NDA is the right leading proxy.

**5. Ranked Computational Options**

**Option A: Keep NDA As-Is**

What it captures:

```text
Perez-Randall lepton-MFV neutrino-spurion dipole.
BR = 4e-8 |(Ybar_N Ybar_N^\dagger)_{e mu}|^2 (3 TeV / M_KK)^4.
```

Inputs needed:

```text
Y_N eigenvalues or Y_N matrix,
PMNS,
M_KK convention.
```

Missing:

```text
Generic charged-lepton anarchy,
Y_E Y_E^\dagger Y_E Higgs/KK-fermion loop,
gauge-loop floor,
mass-basis charged-lepton rotations,
Higgs-localization dependence,
custodial model dependence.
```

Expected accuracy:

```text
Order-of-magnitude only, and only inside Perez-Randall/lepton-MFV assumptions.
```

Implementation cost:

```text
Already done.
```

Verdict:

```text
Useful baseline, not the best final physics model if the lepton sector is generic anarchic RS.
```

**Option B: Perez-Randall NDA With Explicit Mass-Basis Spurion**

What it captures:

```text
Same Perez-Randall normalization, but computes the correct mass-basis
S_N = U_{eL}^\dagger Ybar_N Ybar_N^\dagger U_{eL}
instead of assuming the repo's simplified PMNS/eigenvalue structure.
```

Inputs needed from `ParameterPoint`:

```text
Full Y_N or Ybar_N,
charged-lepton left rotation U_{eL},
M_KK/T convention,
possibly PMNS only if the model fixes Y_N = V_PMNS diag(Y_N).
```

Missing:

```text
All generic Y_E-driven anarchic dipoles.
```

Expected accuracy:

```text
Best low-cost implementation for the Perez-Randall model.
Still NDA-normalized and O(1)-uncertain.
```

Implementation cost:

```text
Low. Requires a well-defined lepton mass-basis extra on ParameterPoint.
```

Verdict:

```text
Best immediate cleanup of the current repo logic.
```

**Option C: Beneke-Calibrated Mass-Basis Spurion Estimator**

What it captures:

```text
Dominant charged-lepton-anarchy Higgs/KK-fermion dipole:
  f_L [Y_E Y_E^\dagger Y_E] f_E

Gauge-loop floor:
  mostly aligned single-Yukawa term with profile-dependent misalignment

Optional Perez-Randall neutrino-spurion term:
  Y_N Y_N^\dagger
```

Inputs needed from `ParameterPoint`:

```text
Full complex 3x3 Y_E,
full complex 3x3 Y_N if the neutrino term is included,
c_{L_i}, c_{E_i},
charged-lepton mass diagonalization matrices U_L, U_R,
M_KK or T convention,
Higgs localization: exact brane / narrow bulk / bulk,
minimal vs custodial RS flag,
Y_* normalization convention.
```

Missing:

```text
Exact 5D propagator integrals,
exact KK Higgs tower summation,
full gauge-sector model dependence,
Barr-Zee and all dimension-six mixing pieces except as optional estimates.
```

Expected accuracy:

```text
Correct dominant spurions and parametric behavior.
Likely factor-few to order-of-magnitude depending on Higgs/custodial assumptions.
Much more physically faithful than Perez-Randall NDA for generic anarchic leptons.
```

Implementation cost:

```text
Medium. Requires real lepton-sector data in ParameterPoint.
```

Verdict:

```text
Best practical approach for this repo if the physicist wants generic anarchic RS LFV.
```

**Option D: Full Explicit One-Loop RS Dipole**

What it captures:

```text
Full Beneke-style one-loop matching:
a^B, a^W, gauge exchange, Higgs exchange, KK Higgs tower,
custodial representations, profile integrals, rotations, and operator mixing.
```

Inputs needed:

```text
Everything in Option C,
plus complete model representation data,
5D gauge couplings,
Higgs profile beta or exact-brane prescription,
KK spectrum or 5D propagator machinery,
wrong-chirality Yukawa prescription,
custodial partner content.
```

Missing:

```text
Only higher-order corrections and UV-completion-dependent brane counterterms.
```

Expected accuracy:

```text
Best available theory calculation.
Needed for publication-grade numerical predictions in a generic RS lepton model.
```

Implementation cost:

```text
High to very high.
Not supportable from the current ParameterPoint/YukawaResult alone.
```

Verdict:

```text
Best ultimate calculation, but not the best next implementation step unless the full lepton model is finalized.
```

**6. Recommended Design**

The best computational path is a two-tier design.

First, define a proper mass-basis lepton-dipole input object on `ParameterPoint`. It should contain:

```text
Y_E_bar       full 3x3 charged-lepton 5D Yukawa matrix
Y_N_bar       full 3x3 neutrino 5D Yukawa matrix, if present
c_L[3], c_E[3], c_N[3] or model-specific neutrino localization data
F_L[3], F_E[3] evaluated consistently with the repo's normalization
U_eL, U_eR    charged-lepton diagonalization matrices
m_e, m_mu, m_tau reproduction diagnostics
M_KK or T, with explicit convention
Higgs localization flag
minimal/custodial model flag
```

Then implement two reported dipole estimates from the same interface:

```text
1. PR-LMFV term:
   BR_PR = 4e-8 |(U_eL^\dagger Ybar_N Ybar_N^\dagger U_eL)_{e mu}|^2
           (3 TeV / M_KK)^4

2. Generic charged-lepton-anarchy term:
   alpha_H ~ coefficient * U_eL^\dagger
             [F_L Ybar_E Ybar_E^\dagger Ybar_E F_E]
             U_eR
```

The code should label them separately:

```text
mu_to_e_gamma_pr_lmfv_proxy
mu_to_e_gamma_charged_yukawa_spurion_proxy
mu_to_e_gamma_gauge_floor_proxy
```

The total should not be advertised as a unique physical prediction until the physicist chooses whether the model is Perez-Randall lepton-MFV or generic anarchic charged leptons. These are different theories, not just different numerical approximations.

**7. Open Physics Decisions For The Physicist**

1. Is the intended lepton model Perez-Randall/lepton-MFV, generic anarchic charged leptons, or a hybrid with both `Y_E` and `Y_N` LFV?

2. Are `Y_E` and `Y_N` stored in the charged-lepton mass basis, the 5D interaction basis, or a neutrino/seesaw basis? The code must not guess this.

3. Should the LFV neutrino spurion be

```text
Y_N Y_N^\dagger
```

in left-handed lepton space, as Perez-Randall assumes, or something modified by right-handed neutrino localization/seesaw structure?

4. Should the charged-lepton dipole include the generic anarchic spurion

```text
Y_E Y_E^\dagger Y_E
```

or is that forbidden/aligned by the model’s flavor symmetry?

5. What is the KK scale convention: `Lambda_IR = T`, first gauge KK mass `m_1 ~= 2.45 T`, or a model-specific fermion KK mass?

6. Is the Higgs exactly brane-localized, narrow-bulk, or true bulk? Exact brane Higgs dipoles have prescription/UV sensitivity; narrow-bulk/bulk Higgs requires KK Higgs treatment.

7. Is the model minimal RS or custodially protected RS? Custodial leptons alter gauge contributions and can add partner-Yukawa structures.

8. Are `c_{L_i}` universal, fixed by neutrino mixing, or scanned independently? This strongly affects whether PMNS-like mixing comes from `Y_N` or charged-lepton rotations.

9. Are wrong-chirality Yukawas equal to the ordinary 5D Yukawas, independent anarchic matrices, or absent by construction?

10. Should the constraint use only central estimates, or should it scan unknown O(1) loop coefficients/phases to expose accidental cancellations?

**Executive Recommendation**

1. The current repo correctly implements the Perez-Randall lepton-MFV NDA proxy, but that is not the generic anarchic RS `\mu -> e\gamma` calculation.  
2. The best practical next approach is Option C: a Beneke-calibrated mass-basis spurion estimator, with the current Perez-Randall term retained as a separately labeled baseline.  
3. It needs a real `ParameterPoint` lepton-sector object: full `Y_E`, full `Y_N`, `c_L/c_E`, profile factors, `U_eL/U_eR`, KK-scale convention, Higgs localization, and minimal/custodial flag.  
4. Compared to Perez-Randall NDA, the explicit anarchic result usually replaces `Y_N Y_N^\dagger` with the dominant charged-lepton spurion `Y_E Y_E^\dagger Y_E`, plus a smaller gauge floor.  
5. If the physicist confirms Perez-Randall lepton-MFV, use Option B; if they confirm generic anarchic leptons, use Option C now and reserve full Beneke-style one-loop matching for final publication-grade validation.  

MUEG-RESEARCH-DONE