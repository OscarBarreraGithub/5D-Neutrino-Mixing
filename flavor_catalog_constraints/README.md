# flavor_catalog_constraints

A plug-in constraint system for the RS flavor scan. Each catalogued
flavor-physics observable becomes **one self-contained file** that wraps
existing physics, loads its experimental anchor from the catalog, and returns
a standardized pass/fail result (the registry currently discovers 103 = 95
PRIMARY + 8 SECONDARY). Adding a constraint is dropping in a file; deleting one
is removing a file. No file imports another constraint, and one broken
constraint cannot break the rest.

Beyond the ΔF=2 adapter (`physics_adapters/deltaf2.py`), the package wires
RS-electroweak observables (oblique S,T,U via `quarkConstraints/oblique_stu.py`,
Z→bb T010/T011, off-diagonal Z FCNC T014, EW001), top/Higgs and collider
recasts, and semileptonic Wilson coefficients, via `rs_ew_builder.py` and
`point_builder.py`. Runtime tag policy (`rigorous` / `proxy` / `partial` /
`stub`, plus HARD/SOFT/INFO severity) is decided by the scan harness
(`scripts/run_full_catalog_scan.py:tag_result`), not by static metadata.

## Layout

```
flavor_catalog_constraints/
├── base.py              # contract types: ParameterPoint, ConstraintResult,
│                        #   Severity, ConstraintLevel, ConstraintProtocol
├── anchors.py           # schema-flex YAML loader -> typed Anchor (loud on miss)
├── registry.py          # @register decorator, auto-discovery, evaluate_all
├── point_builder.py     # ParameterPoint construction + KNOWN_EXTRA_KEYS
├── physics_adapters/    # the ONLY modules that import physics code
│   └── deltaf2.py       #   wraps quarkConstraints/deltaf2.py
├── primary/<family>/<ID>.py     # PRIMARY-tier constraints
├── secondary/<family>/<ID>.py   # SECONDARY-tier constraints
├── TEMPLATE.py          # copy this to create a new constraint
└── README.md
```

Families: `beauty`, `charged_lepton`, `charm`, `collider_rs`,
`edm_neutrino`, `kaon`, `top_higgs_ew`.

Each constraint file has a sidecar in
`flavor_catalog/processes/<family>/<ID>.yaml` (primary) or
`.../secondary/<family>/<ID>.yaml` (secondary). The registry refuses to
register a constraint whose sidecar is missing.

## The contract

Every constraint is a plain class named `Constraint` decorated with
`@register`. It declares three attributes and one method:

```python
@register
class Constraint:
    process_id = "K001"          # ^[A-Z]+[0-9]+$, matches the sidecar
    severity   = Severity.HARD   # HARD / SOFT / INFO
    observable = "epsilon_K"

    def evaluate(self, point: ParameterPoint) -> ConstraintResult: ...
```

`level` (PRIMARY/SECONDARY) and `family` are **derived from the file's
path** by the decorator — never declared by the author, so a file is
portable between tiers/families just by moving it.

`ConstraintResult` numeric fields (`predicted`, `sm_prediction`,
`experimental`, `ratio`, `budget`) are **real floats or None**. Complex
amplitudes, Wilson coefficients, etc. go in `diagnostics`.

**Severity semantics:** HARD = observed bound, a failure vetoes the
point. SOFT = SM-vs-exp tension / projection, advisory. INFO = never
vetoes.

## Adding a constraint

1. `cp TEMPLATE.py primary/<family>/<ID>.py`
2. Set `process_id`, `severity`, `observable`.
3. Set `_ANCHOR_CANDIDATES` to the key(s) your sidecar's
   `pdg_or_equivalent` block actually uses (the loader tries them in
   order and raises `AnchorError` if none match — never a silent
   default).
4. Reach physics through a `physics_adapters.*` wrapper (add one if
   needed — adapters are append-only). Never import `quarkConstraints` /
   `flavorConstraints` / `qcd` directly.
5. `cp ../../tests/constraints/test_constraint_template.py
   ../../tests/constraints/test_<ID>.py` and fill in.

That's it — discovery finds the file automatically.

## Removing a constraint

Delete the one file (and its per-constraint test). Nothing references it
by name; the registry simply stops discovering it.

## Isolation guarantees

- **Plug-in / atomic:** a constraint is one file with no cross-constraint
  imports. Cross-constraint coordination is limited to the shared extras
  keys in `point_builder.KNOWN_EXTRA_KEYS`.
- **Import isolation:** `discover()` imports each leaf module in its own
  `try/except`; a failure is recorded in `import_failures()`, never
  raised, so one broken file does not break discovery.
- **Evaluation isolation:** `evaluate_all()` wraps every `evaluate`
  call; a raising constraint yields a failing `ConstraintResult` with
  the exception in `diagnostics`, and the rest still run.
- **Adapter moat:** all physics imports live in `physics_adapters/`, so
  an upstream signature change is absorbed in one place.

## ParameterPoint construction

```python
from flavor_catalog_constraints import make_point, empty_point, build_from_quark_couplings

empty_point()                                   # no extras (smoke tests)
make_point(kk_gluon_mass_gev=3000.0)            # hand-build; unknown keys raise
build_from_quark_couplings(couplings)           # real quark-sector path
```

`make_point` rejects any keyword not in `KNOWN_EXTRA_KEYS`, so a typo'd
input fails loudly at construction instead of surfacing as a silent
`None` inside a constraint. Constraints read inputs via
`point.get_extra("key")`.

## Running the tests

```
python -m pytest tests/constraints/
```

The global contract test (`test_contract.py`) passes on an empty
constraint set and scales as real constraints land: it checks discovery,
protocol conformance, metadata invariants, sidecar existence, exception
isolation, and the loud-failure behavior of the anchor loader.
