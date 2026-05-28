"""PRIMARY-tier constraints.

Each subpackage is a family (``beauty``, ``charged_lepton``, ``charm``,
``collider_rs``, ``edm_neutrino``, ``kaon``, ``top_higgs_ew``); each
leaf module is a single constraint registered via
``@register_constraint``.

Family ``__init__.py`` files are intentionally empty so that removing
a constraint is a single ``rm`` of the leaf module — no dangling
imports.
"""
