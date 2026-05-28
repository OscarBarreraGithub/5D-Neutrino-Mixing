"""Thin import-boundary wrappers over shared physics modules.

The directive (``.orchestration/PHASE3_CONSTRAINTS_DIRECTIVE.md``)
requires that modifying one constraint cannot break another. Shared
physics modules are the risk vector: 20 different kaon and B-mixing
constraints all sit on top of ``quarkConstraints.deltaf2``, so a
signature change there could ripple across the whole zoo.

The adapter layer here is the isolation moat. Conventions:

- One adapter module per shared physics module. Constraints import
  from these adapters; they never import ``quarkConstraints.*`` or
  ``flavorConstraints.*`` directly.
- Adapters return the existing physics-module dataclasses
  **unchanged** (R1 CR-3). The adapter's only job is the import
  boundary — it does not repackage signatures.
- Adapters are **append-only** in the codex loop (R1 M-6). A
  constraint may add a new wrapper function; it may not modify the
  signature of an existing one. Modifying an existing wrapper is a
  cross-constraint change that requires a separate sign-off.

Today this package contains:

- :mod:`flavor_catalog_constraints.physics_adapters.deltaf2` — wraps
  ``quarkConstraints/deltaf2.py`` (used by K001).

Future adapters will be added on demand:

- ``meg`` → ``flavorConstraints/muToEGamma.py``
- ``modern_phenomenology`` → ``quarkConstraints/modern/phenomenology.py``
- ``edm`` → (TBD)
- ``collider`` → (TBD)
"""
