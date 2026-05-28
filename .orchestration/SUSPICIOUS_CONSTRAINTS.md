# Suspicious Constraints

## K003 Re(eps'/eps)

K003 is registered as a PRIMARY kaon constraint, but its evaluator is an
explicit non-vetoing placeholder. The repository currently has Delta F = 2
kaon-mixing machinery for K001/K002 only; the requested audit found no
existing K -> pi pi, Delta S = 1, eps'/eps, or QCD-penguin implementation in
`quarkConstraints/`, `scripts/`, or `flavor_catalog_constraints/physics_adapters/`.

Full evaluation would need:

- Delta S = 1 Wilson-coefficient matching for the RS NP sources.
- RG evolution for the relevant QCD penguin, electroweak penguin, and
  chromomagnetic operators.
- K -> pi pi hadronic matrix elements, including the lattice inputs needed
  to connect the Wilsons to Re(eps'/eps).
- A convention audit so the result is not mixed with the existing Delta F = 2
  epsilon_K / Delta m_K machinery.

Constraints with the same broad issue are likely K005, K006, K012, and K013:
rare or nonleptonic kaon modes require Delta S = 1 amplitudes, form factors
or matrix elements, and usually electroweak-penguin or long-distance treatment
that is not present in the current physics core. K017 is suspicious for a
different reason: it is a charged-current LFU ratio, not a rare Delta S = 1
mode.

## K004 K+ -> pi+ nu nubar

K004 is registered as a PRIMARY kaon constraint, but its evaluator is an
explicit non-vetoing placeholder. This is a Delta S = 1 rare decay, not a
Delta F = 2 kaon-mixing observable and not the K003 eps'/eps placeholder.

Full evaluation would need:

- K -> pi nu nu Wilson matching, including the relevant BMU operator mapping.
- CKM elements and a documented SM-input convention.
- SM--NP interference treatment for the charged branching ratio.
- A Standard Model reference prediction from Buras, Schwab, and Uhlig,
  hep-ph/0508165, using the SM-NNLO benchmark
  BR_SM(K+ -> pi+ nu nubar) ~ 8.4e-11.

## K005 K_L -> pi^0 nu nubar

K005 is registered as a PRIMARY kaon constraint, but its evaluator is an
explicit non-vetoing placeholder. This is a Delta S = 1 rare decay and the
neutral-kaon mode is a Standard-Model CP-violation-pure observable, with no
charm contribution to lowest order.

Full evaluation would need:

- Delta S = 1 K -> pi nu nubar Wilson matching and operator-basis mapping.
- CKM inputs, SM--NP interference conventions, and the CP-odd K_L projection.
- A neutral-mode branching-ratio evaluator, separate from the charged K004
  formula.
- A Standard Model reference prediction from Buras, Gorbahn, Haisch, and
  Nierste, hep-ph/0603079, using BR_SM(K_L -> pi^0 nu nubar) ~ 3e-11.
- The direct KOTO 2024 upper limit BR(K_L -> pi^0 nu nubar) < 2.2e-9 at
  90% CL.

## K006 K_L -> mu mu

K006 is registered as a PRIMARY kaon constraint, but its evaluator is an
explicit non-vetoing placeholder. K_L -> mu mu is long-distance dominated
(~99% from 2-photon intermediate state), short-distance SM ~9e-10 from
Buras et al; full evaluation deferred.

## K008 K_L -> pi0 e+ e-

K008 is registered as a PRIMARY kaon constraint, but its evaluator is an
explicit non-vetoing placeholder. This is a Delta S = 1 rare semileptonic
mode with direct/indirect CP interference and long-distance inputs from
K_S -> pi0 e+ e- and K_L -> pi0 gamma gamma; it is not covered by the
Delta F = 2 kaon-mixing machinery used for K001/K002.

Full evaluation would need:

- s -> d e+ e- Wilson matching and operator-basis mapping.
- CKM inputs plus SM--NP interference conventions.
- Direct/indirect CP separation using the K_S -> pi0 e+ e- input.
- CP-conserving long-distance treatment tied to K_L -> pi0 gamma gamma.
- A branching-ratio evaluator for K_L -> pi0 charged-lepton modes.

## K009 K_L -> pi0 mu+ mu-

K009 is registered as a PRIMARY kaon constraint, but its evaluator is an explicit non-vetoing placeholder because rare semileptonic Delta S = 1 K_L -> pi0 mu+ mu- evaluation is deferred pending s -> d mu+ mu- matching, CP-interference, and long-distance inputs.

## K010 K_S -> pi0 e+ e-

K010 is registered as a PRIMARY kaon constraint, but its evaluator is an explicit non-vetoing placeholder because rare semileptonic Delta S = 1 K_S -> pi0 e+ e- evaluation is deferred pending s -> d e+ e- matching, SM-NP interference conventions, form-factor treatment, phase-space extrapolation, and branching-ratio machinery.

## K012 K_S -> mu mu

K012 is registered as a PRIMARY kaon constraint, but its evaluator is an explicit non-vetoing placeholder because rare dimuon Delta S = 1 K_S -> mu+ mu- evaluation is deferred pending s -> d mu+ mu- matching, SM-NP interference conventions, K_S-specific long-distance two-photon treatment, time-dependent neutral-kaon interference conventions, and branching-ratio machinery.

## K013 K_L -> pi0 gamma gamma

K013 is registered as a PRIMARY kaon constraint, but its evaluator is an explicit non-vetoing placeholder because radiative ChPT evaluation is deferred pending chiral perturbation theory O(p^4)+O(p^6) amplitudes and photon matrix elements; Delta S = 1 Wilson matching alone is not enough.

## K017 R_K charged-kaon LFU ratio

K017 is registered as a PRIMARY kaon constraint, but its evaluator is an
explicit non-vetoing placeholder. R_K is a tree-level charged-kaon
lepton-universality ratio, not a Delta S = 1 rare decay. Full evaluation is
deferred pending modified W-lepton coupling treatment in 5D RS, radiative QED
corrections matched to the inclusive K_l2 convention, and any new tree-level
four-fermion contact-operator matching.

## K018 K_l3 |V_us| f_+(0)

K018 is registered as a PRIMARY kaon constraint, but its evaluator is an
explicit non-vetoing placeholder. K_l3 is a tree-level charged-current
semileptonic kaon source observable for |V_us| f_+(0), not a Delta S = 1 rare
decay. Full evaluation is deferred pending source-level |V_us| extraction,
5D RS effects on the W-quark-quark coupling, radiative/isospin corrections,
lattice f_+(0) provenance, and any new four-fermion contact-operator matching.

## EDM family (E001-E009) — all deferred

The current physics core has no CP-violating dipole machinery. All 7 EDM
constraints are non-vetoing placeholders:

- **E001** Electron EDM |d_e| < 4.1e-30 e·cm (JILA HfF+ 2023)
- **E002** Muon EDM |d_μ| (BNL 2008, PSI muEDM aiming 1e-21)
- **E004** Neutron EDM |d_n| < 1.8e-26 e·cm (nEDM PSI 2020)
- **E006** Hg-199 atomic EDM < 7.4e-30 e·cm (Graner 2016)
- **E007** Ra-225 + Xe-129 atomic EDMs (Schiff moment enhancement)
- **E008** Quark chromo-EDM ~d_q from hadronic EDMs
- **E009** Weinberg three-gluon CP-violating operator

Full evaluation requires: CP-violating dipole Wilson matching, 5D RS
one-loop contributions to imaginary parts of dipole operators, hadronic
matrix elements (Schiff moments for atomic, QCD sum rules for neutron),
and theta-bar / Peccei-Quinn treatment for hadronic species. None of
this is in the physics core yet.

## Charged lepton family (L001-L023) — all deferred

The current physics core has μ→eγ machinery (used in scanParams) but lacks
LFV 4-fermion operators, μ→e conversion form factors, μ-onium / antimuonium
treatment, and ν-trident amplitudes. All 11 charged-lepton constraints are
non-vetoing placeholders.

Deferred subgroups:
- L001/L007/L008: LFV dipole decays (l → l' γ); MEG II 2025 BR(μ→eγ) < 1.5e-13.
- L002/L009/L010: LFV 4-lepton decays (μ→3e, τ→3μ, τ→3e).
- L003/L004/L005: μ → e conversion on Al / Au / Ti nuclei (SINDRUM-II Au < 7e-13).
- L006: Muonium-antimuonium oscillation (P_MMbar).
- L023: Neutrino trident production (CCFR; Altmannshofer-Gori-Pospelov 1406.2332).

Full evaluation requires: LFV dipole + 4-fermion Wilson matching at the
EW scale, RG to low energies, nuclear-overlap form factors for μ→e
conversion (Kitano-Koike-Okada), Z' tridents for L023, and 5D-RS
contributions to all of the above (KK gauge / Higgs / radion loops).

## Charm family (C001-C008) — all deferred

The current physics core has ΔC=2 D-mixing machinery via deltaf2.py
(adapter could plug in for C001/C002), but ΔC=1 rare-decay machinery
(D → ll, D → π ll, LFV charm) is absent. All 8 charm constraints are
non-vetoing placeholders for now.

- **C001/C002** D-mixing (ΔC=2): could be promoted to full implementation
  using the existing deltaf2.delta_m_D evaluator. Tracked for follow-up.
- **C003** Direct CPV: ΔC=1, requires NP penguins.
- **C004/C005** D⁰ → ll: ΔC=1 rare leptonic.
- **C006/C008** D LFV: requires LFV operators.
- **C007** D → π ll: rare semileptonic.

## Beauty family (B001-B034) — all deferred (B001-B004 promotable)

22 beauty placeholders. B001-B004 are ΔB=2 mixing observables and could be
promoted to full impl using deltaf2.py (Δm_d, Δm_s exist; CP phases would
need Im-part adapter). The other 18 are ΔB=1 (rare leptonic, semileptonic,
radiative, LFU, baryonic, invisible, penguin) and need new machinery.

Promotion candidates: B001 (Δm_d), B002 (sin2β), B003 (Δm_s), B004 (φ_s).

## Top/Higgs/EW family (EW001-EW003 + T001-T020) — all deferred

20 placeholders. EW oblique (S, T, U) parameters, CKM unitarity tensions,
top FCNC (t→qV, t→Hq), Z-pole asymmetries (bb, cc), LFV Z + LFV Higgs
decays. All require new EW-precision and FCNC machinery (radion + KK
gauge mixing for STU; top FCNC dipole + chromo Wilson coefficients;
LFV-boson 5D RS contributions).

## Collider-RS family (CR001-CR014) — all deferred

14 placeholders. Direct-collider bounds on RS-sector states: KK gluon
(dijet/ttbar), KK Z'/W' (dilepton/W'->lν), KK graviton (diphoton/dilepton),
custodial top partners (T_{5/3}, T, B), singlet/doublet VLQs, DY EFT,
VBS, diboson/diphoton resonances, four-top. All require LHC direct-search
recasting + 5D-RS mass-coupling map.
