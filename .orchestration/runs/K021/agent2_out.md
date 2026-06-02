1. NIT: ΔF=2 M12/Re/Im checks are not applicable here; K021 is a ΔS=1 LFV semileptonic BR. The code correctly compares pure-NP BR to an upper limit, using |YV|^2+|YA|^2 for the charge-summed LFV proxy. `flavor_catalog_constraints/secondary/kaon/K021.py:438`, `flavor_catalog_constraints/physics_adapters/rare_kaon_lfv_semileptonic.py:596`.

2. NIT: No QCD `*_with_running` path is required for this vector-current semileptonic proxy; the color-singlet quark vector current has unity QCD running at this level. Code routes through rare-kaon dilepton matching plus K→π f+(q²) phase space, not deltaf2. `rare_kaon_lfv_semileptonic.py:299`, `:370`, `:553`.

3. NIT: Budget is correct: HARD pure-NP upper bound = PDG YAML 7.6e-11 at 90% CL, with KTeV provenance 7.56e-11 rounded by PDG. `flavor_catalog/processes/secondary/kaon/K021.yaml:97`, `:113`; `K021.py:314`, `:457`.

4. NIT: Inputs match neutral mode: m_KL=0.497611 GeV, m_pi0=0.1349768 GeV, tau_KL=5.116e-8 s, f+(0)=0.9698, charge_state_factor=2 for summed e±mu∓. `rare_kaon_lfv_semileptonic.py:729`, `:744`.

5. NIT: Severity/diagnostics are appropriate: HARD, SM LFV rate set to 0, and the e-mu coupling plus K_L CP/orientation treatment are explicitly flagged NEEDS-HUMAN-PHYSICS. `K021.py:18`, `:30`, `:466`; `rare_kaon_lfv_semileptonic.py:712`.

PHYSICS-OK