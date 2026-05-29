1. NIT: ΔF=2 amplitude/running checks are not applicable to CR008: this is a collider mass lower-bound recast, with no `M12`, CP phase, Wilsons, or `mu_had` path. Verdict uses `collider_resonance`, not non-running `deltaf2`: `CR008.py:52`, `CR008.py:404`, `CR008.py:409`; core mass test is `m_pred >= m_limit` at `quarkConstraints/collider_resonance.py:166`.

2. NIT: Budget is defensible and channel-matched: active HARD limit is ATLAS singlet `T`, `m_T > 1.36 TeV`, for `B(T->Wb:Ht:Zt)=1/2:1/4:1/4`; not a bare residual. YAML/source match: `CR008.yaml:99`, `CR008.yaml:121`, `atlas_2024_arxiv2401_17165.txt:25`, `pdg2025_tprime_vlq_limits.txt:12`. CMS `1.48 TeV` is correctly diagnostic because it is a broader all-branching envelope: `CR008.yaml:125`, `CR008.py:435`.

3. NIT: Units and ratio are consistent: YAML stores TeV, code compares TeV, and only exposes GeV conversion as diagnostics. Example: `M_KK=1600 GeV -> m_T=1.6 TeV`, ratio `1.36/1.6=0.85`, pass. See `CR008.py:107`, `CR008.py:448`, `CR008.py:450`.

4. NIT: σ×BR recast is correctly marked `NEEDS-HUMAN-PHYSICS` in verdict diagnostics and notes: `CR008.py:438`, `CR008.py:455`. Minor wording issue only: `BR(...)^2` at `CR008.py:17` / adapter `collider_resonance.py:152` should eventually say branching-fraction products/mixtures for pair decays.

PHYSICS-OK