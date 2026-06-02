1. NIT: ΔF2 amplitude check is N/A for CR014; no \(M_{12}\), Re/Im, or CP observable is evaluated. Verdict is a collider mass-threshold comparison only: `CR014.py:690-696`, `quarkConstraints/collider_resonance.py:166-168`.

2. NIT: QCD running check is N/A; CR014 does not use Wilson coefficients or any DeltaF2 path. Its verdict calls the collider resonance evaluator, not a non-running flavor evaluator: `CR014.py:690-696`, `physics_adapters/collider_resonance.py:619-646`.

3. NIT: Budget is correct for the documented proxy: strongest applicable vector, Γ/m=0.50, observed mass-exclusion row is 850 GeV = 0.85 TeV, with ratio \(m_\text{limit}/m_\text{proxy}\). See `CR014.yaml:88-108`, `CR014.py:504-545`, `CR014.py:581-587`.

4. NIT: Width/sigmaBR/acceptance/SM-background limitations are correctly flagged as NEEDS-HUMAN-PHYSICS in docstring and diagnostics: `CR014.py:16-23`, `CR014.py:743-750`.

5. NIT: Anchor numbers match local snapshots/YAML: CMS 850 GeV observed, 1000 GeV expected at 50% width; PDG four-top 22.5+6.6/-5.5 fb and 17.7+4.4/-4.0 fb; ATLAS 21-119 fb range. See `CR014.yaml:88-181`.

6. NIT: Severity and units are consistent: HARD observed 95% CL bound; GeV input converted to TeV. Spot check: 800 GeV fails ratio 1.0625, 850 GeV passes ratio 1.0, 900 GeV passes ratio 0.9444. Focused tests: 8 passed.

PHYSICS-OK