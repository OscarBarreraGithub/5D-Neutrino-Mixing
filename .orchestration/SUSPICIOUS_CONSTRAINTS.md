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
