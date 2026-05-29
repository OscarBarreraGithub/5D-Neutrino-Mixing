1. NIT — ΔF2 amplitude check is N/A: CR001 computes no `M_12`; it compares `m_gKK` to a direct ttbar resonance limit, so no magnitude-vs-imaginary-part issue applies. `flavor_catalog_constraints/primary/collider_rs/CR001.py:377-411`, `quarkConstraints/collider_resonance.py:166-203`.

2. NIT — QCD running check is N/A: no Wilson coefficients or `mu_had=2 GeV` evolution belong to this collider mass proxy; no deltaF2 non-running path is used for the verdict. Running effect: not defined for this observable. `quarkConstraints/collider_resonance.py:1-15`.

3. NIT — Budget is physics-defensible for the documented proxy: YAML active limit is CMS2026 5.5 TeV at 95% CL from excluded interval `0.5 < m(g_KK) < 5.5 TeV`; code uses HARD ratio `5.5/m_KK`, so 6.0 TeV gives 0.917 pass and 3.0 TeV gives 1.833 fail. `CR001.yaml:82-109`, `CR001.py:64-68,400-411`.

4. NIT — Anchor numbers match snapshots: CMS2026 5.5 TeV, CMS2019 4.55 TeV, PDG Live 3.8 TeV; units are TeV after GeV→TeV conversion. `cms_2026_arxiv2603_23454.txt:16`, `hepdata_2026_ins3134005_record_extract.json:8`, `pdg_live_2026_s071kkg.txt:2-4`, `quarkConstraints/collider_resonance.py:109-127`.

5. NIT — Model-dependence caveat is correctly stated: full `sigma*BR`, width, interference, acceptance need a collider recast; current mass edge is explicitly marked NEEDS-HUMAN-PHYSICS and not presented as model-independent. `CR001.py:11-23,408-411`, `quarkConstraints/collider_resonance.py:3-15`.

PHYSICS-OK