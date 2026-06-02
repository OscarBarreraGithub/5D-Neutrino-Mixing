Plan executed: studied scaffold/CR010/CR012/K001/K004/CR014.yaml, appended adapter helpers, added CR014, added tests, ran suite.

Physics: HARD mass-vs-limit proxy for CMS-B2G-25-005 four-top two-lepton top-philic vector `Z'`, `Gamma/m=50%`; `m_proxy = kk_ew_mass_gev or M_KK`.

Reused/built: reused `collider_resonance` mass comparison; appended CR014-specific adapter helpers only.

Budget: `850 GeV` observed 95% CL from `CR014.yaml`, converted to `0.85 TeV`.

Validation: independent core recomputation checks `1.0 TeV -> ratio 0.85`; safe/excluded tests pass at `0.90/0.80 TeV`.

SM/reference diagnostics: YAML-loaded ATLAS/CMS four-top measurements are retained for context, not veto.

NEEDS-HUMAN-PHYSICS: sigma*BR, width dependence, top-philic couplings, four-top acceptance, and SM-four-top background recast flagged in docstring + diagnostics.

Files changed: `physics_adapters/collider_resonance.py`, `primary/collider_rs/CR014.py`, `tests/constraints/primary/collider_rs/test_CR014.py`.

Pytest: `python -m pytest tests/constraints/ -q` -> `1041 passed in 33.54s`.

Unrelated worktree changes were present and left untouched.