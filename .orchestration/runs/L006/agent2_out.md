1. NIT: Amplitude check OK. L006 is not a neutral-meson Δm/CP observable; the correct proxy is dimensionless `P_MMbar ∝ |G_C/G_F|^2`. Code uses `abs(g_value) ** 2` and `P/P_limit`, not an imaginary part: `muonium_conversion.py:126-130`, `:144`, `:177`; returned prediction is probability at `L006.py:330`.

2. NIT: QCD running check N/A. This is a charged-lepton four-fermion proxy, not `Delta F=2`; no hadronic `mu_had=2 GeV` running is physically required. L006 imports only `muonium_conversion`, not `deltaf2`: `L006.py:43-47`. Running effect: N/A.

3. NIT: Budget OK. HARD veto uses the pure-NP MACS/PDG probability limit, not a central residual: YAML `P < 8.3e-11` at 90% CL, 0.1 T (`L006.yaml:85-96`); budget returns that value (`L006.py:93-97`); comparison is `P/P_limit` (`muonium_conversion.py:144`).

4. NIT: Anchors OK. `G_C/G_F < 0.0030` (`L006.yaml:97-108`) and original MACS `P <= 8.3e-11` (`L006.yaml:109-122`) match snapshots: PDG `pdg...txt:22-30`, Willmann `willmann...txt:23-31`. Calibration is `8.3e-11 / 0.0030^2 = 9.222e-6`.

5. NIT: Wording precision. Some implementation prose says generic `Delta L=2` (`L006.py:5`, `:18`, `:56`, `:339`; `muonium_conversion.py:4`, `:41`). Precise physics is `Delta L_mu = -Delta L_e = 2` family-number violation with total lepton number conserved; YAML/Tex already state this correctly (`L006.yaml:188`, `:213`; `L006.tex:16-17`, `:80-81`).

PHYSICS-OK