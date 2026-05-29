1. NIT: ΔF=2 checklist is inapplicable: T004 is a top decay BR, not an M12 observable; no real/imag M12 issue and no `mu_had=2 GeV` QCD-running verdict path should be used. Verdict path is `t_to_q_gamma_from_couplings` at `T004.py:450`, using photon dipole BR in `quarkConstraints/top_fcnc.py:355`.

2. NIT: Correct photon-dipole amplitude is used: `Γ(t->qγ)=0.5*alpha_em*m_t*(|lambda_L|^2+|lambda_R|^2)` at `quarkConstraints/top_fcnc.py:215`, so phases enter by magnitudes, not an inappropriate imaginary part.

3. SHOULD-FIX: Budget is conservative but chirality-mismatched for generic L/R proxy points: active limit is always ATLAS LH `0.85e-5` at `T004.py:397`/`T004.py:403`, while YAML has RH `1.2e-5` and CMS observed `0.95e-5` at `T004.yaml:108`/`T004.yaml:119`/`T004.yaml:130`; pure RH points are too tight by up to `1.2/0.85=1.41`.

4. NIT: Anchor numbers match snapshots/YAML: PDG combined `9.5e-6`, ATLAS LH `8.5e-6`, ATLAS RH `1.2e-5`, CMS observed `9.5e-6`, SM context `~1e-15`; see `T004.yaml:96`, `:108`, `:119`, `:130`, `:160` and snapshots `atlas_2205_02537_arxiv_abs.txt:15`, `cms_2312_08229_public_page.txt:14`.

5. NIT: Severity/units/diagnostics are physics-consistent: HARD pure-NP collider-limit comparison with dimensionless branching fractions, SM omitted as negligible, RS proxy flagged NEEDS-HUMAN-PHYSICS at `T004.py:27`, `T004.py:457`, `T004.py:501`, `T004.py:510`.

PHYSICS-NEEDS-FIXES