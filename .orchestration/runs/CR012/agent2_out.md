1. NIT: Δm/CP/QCD-running checklist is N/A for CR012; this is a collider mass-lower-bound proxy, not an \(M_{12}\) or Wilson-running constraint. Correct observable is \(m_{\rm limit}/m_{\rm proxy}\). `CR012.py:488`, `collider_resonance.py:144`.

2. NIT: Active budget is physics-consistent with the prompt: HVT model-B `W' -> WZ`, observed `m(W') > 4.4 TeV` at `95% CL`; GeV→TeV conversion and ratio `4.4 / m_proxy_TeV` are consistent. `CR012.yaml:87`, `CR012.py:419`, `CR012.py:544`, `CR012.py:613`.

3. NIT: Stronger YAML entries are handled appropriately as diagnostics: `V' -> VV = 4.5 TeV` assumes mass-degenerate HVT triplet, and `VV+VH = 4.8 TeV` is not pure diboson-only. Keeping active `4.4 TeV` is defensible for the single W'-like proxy. `CR012.yaml:107`, `CR012.yaml:127`, `CR012.py:579`.

4. NIT: NEEDS-HUMAN-PHYSICS is correctly documented for missing σ×BR, branching surface, width, production mix, interference, acceptance, and limit-curve recast. `CR012.py:17`, `collider_resonance.py:518`, `CR012.py:599`.

PHYSICS-OK