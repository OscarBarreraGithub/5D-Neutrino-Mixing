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

Constraints with the same broad issue are likely K005, K006, K012, K013, and
K017: rare or nonleptonic kaon modes require Delta S = 1 amplitudes, form
factors or matrix elements, and usually electroweak-penguin or long-distance
treatment that is not present in the current physics core.

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
