1. NIT: ΔF=2 amplitude/QCD-running checks are N/A for CR006; the verdict path is a collider mass lower-bound comparison, not `M12` or RG evolution. Active code uses `kk_charged_current_prediction_from_m_kk_gev` plus `evaluate_collider_resonance_limit`: `CR006.py:415-420`, `collider_resonance.py:166-168`.

2. NIT: Shared resolver docstring still says “neutral electroweak KK”/CR005, but CR006 uses it for charged `W_KK`; behavior is OK as a common EW mass proxy, but wording should be generalized. `physics_adapters/collider_resonance.py:280-290`, used at `CR006.py:365-368`.

3. NIT: Active budget is consistent: `6.0 TeV` PDG/ATLAS SSM `W' -> e nu` 95% CL lower limit, ratio `m_limit/m_WKK`; snapshots also match `5.6`, `5.7`, `4.3`, `3.9 TeV` diagnostics. `CR006.yaml:80-167`, `pdg2025_wprime_searches.txt:14-19`, `CR006.py:467-480`.

4. NIT: `sigma*BR` recast is correctly not used for the HARD verdict and is marked `NEEDS-HUMAN-PHYSICS`; correct physics needs production, BR, width/interference, acceptance, and mass-dependent limits. `CR006.py:17-23`, `CR006.py:459-464`.

PHYSICS-OK