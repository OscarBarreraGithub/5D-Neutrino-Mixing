"""Physics-adapter boundary: the ONLY place that imports physics modules.

Why this layer exists
---------------------
Many constraints sit on top of a handful of shared physics modules
(``quarkConstraints.deltaf2``, ``flavorConstraints``, ``qcd``, ...). If
each constraint imported those directly, a signature change in one
physics function would ripple across every constraint that touches it —
violating the "modifying one constraint cannot break another" rule.

Conventions
-----------
- One adapter module per shared physics module.
- Constraints import ONLY from ``flavor_catalog_constraints.physics_adapters.*``;
  never from ``quarkConstraints.*`` / ``flavorConstraints.*`` / ``qcd`` directly.
- Adapters are **append-only**: add new wrapper functions freely, but do
  not change the signature of an existing wrapper (that is a
  cross-constraint change). Absorbing an upstream signature change
  happens *inside* the existing wrapper body, leaving its signature
  stable for callers.
- Adapters re-export the upstream result dataclasses unchanged so the
  result shape is owned by the physics module, not duplicated here.
"""

from __future__ import annotations

__all__: list[str] = []
