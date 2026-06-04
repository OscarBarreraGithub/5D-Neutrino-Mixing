Updated [.orchestration/runs/W6-HARNESS/plan.md](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/.orchestration/runs/W6-HARNESS/plan.md).

Recorded fixes:
- Per-tile precomputed `a(c)` spline plus `Omega_n(c)` vector lookup; no reliance on `_omega_cache`/`_a_cache`.
- Required future builder hook: `spectrum=` plus `a_of_c=`/`omega_of_c=` injection.
- Honest veto tagging: `rigorous | proxy | partial | stub`, dual `survives_all_HARD`, separate rigorous/proxy exclusions, `hard_not_evaluated`, `coverage_complete`.
- Real `QuarkFitResult`; new driver reusing only `run_rs_anarchy.py` mp/tiling/seeding style.
- Smoke gate must measure real post-cache timing.

Realistic target recorded: `1.2-3 s/point`; `1e8` is about `4-8` core-years and only proceeds if smoke confirms. No production code changed.