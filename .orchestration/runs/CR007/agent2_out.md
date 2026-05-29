1. SHOULD-FIX: HARD verdict uses `m_GKK = kk_ew_mass_gev, kk_gluon_mass_gev, or M_KK` (`kk_graviton_resonance.py:54`, `:57`, `:60`; `CR007.py:379`, `:437`), but CR007.yaml says the limit constrains the first spin-2 graviton mass and “does not by itself constrain the KK gluon … without a model spectrum” (`CR007.yaml:288`). Correct physics needs a dedicated `m_GKK` or explicit spectrum map; repo `M_KK` can be bookkeeping `Lambda_IR` (`scales.py:43-45`), while first gauge root is `2.4487` (`scales.py:21`) and first graviton root is ~`3.83`, so pass/fail can shift by O(1.6-3.8).

2. NIT: ΔF=2 amplitude checks are N/A, not wrong: CR007 has no `M12`, Δm, or CP observable; it calls the collider mass-limit path (`CR007.py:437-442`), not real/imaginary mixing amplitudes.

3. NIT: QCD running check is N/A: no Wilson coefficients or `deltaf2` evaluator are involved; the active verdict is the collider ratio `m_limit / m_pred` (`quarkConstraints/collider_resonance.py:166-168`), so no non-running ΔF=2 path is being used.

4. NIT: Anchor numbers match the YAML/snapshots: active CMS bulk-RS `1.4 TeV` at `k/Mbar_Pl=0.5` (`CR007.yaml:129-150`; PDG snapshot `:82-86`); diagnostics retain generic `4.8`, `4.78`, and ATLAS bulk `2.3 TeV` (`CR007.py:462-480`).

5. NIT: Units/severity/σ×BR handling are internally consistent: result mass fields are TeV, HARD ratio is `limit/mass`, and the missing full `sigma*BR`, width, BR, acceptance, and limit-curve recast is explicitly flagged (`CR007.py:485-486`).

PHYSICS-NEEDS-FIXES