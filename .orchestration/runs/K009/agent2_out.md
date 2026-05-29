1. NIT: No ΔF=2/M12 issue: K009 is ΔS=1 semileptonic, not mixing. Verdict uses `Im(lambda_t*y7V/A)` for direct CP, not `M12`; see `rare_kaon_dilepton_muon.py:147-153` and `K009.py:786-795`.

2. NIT: QCD-running checklist is N/A for this K008-style semileptonic Q7V/Q7A proxy. Code records LO semileptonic factor 1.0 and effect 0%; see `rare_kaon_dilepton_muon.py:325-329` and core note `rare_kaon_dilepton.py:691-698`.

3. NIT: Budget is defensible and YAML-based: HARD budget is PDG/KTeV total-rate upper limit `3.8e-10`; SM direct is `2.02e-12`, SM total constructive is `1.53e-11`, both far below limit. See `K009.py:620-625`, `K009.py:792-795`.

4. NIT: Anchor spot-checks match snapshots: PDG/KTeV limit `3.8e-10`, KTeV `2` events and `0.87 +/- 0.15` background, ISU coefficients `3.7, 1.6, 1.0, 5.2`, `|a_S|=1.2`, SM `(1.5 +/- 0.3)e-11`; see `K009.yaml:92-109`, `118-151`, `185-254`, `313-331`.

PHYSICS-OK