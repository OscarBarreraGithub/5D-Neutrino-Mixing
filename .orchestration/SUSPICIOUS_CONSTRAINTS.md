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
