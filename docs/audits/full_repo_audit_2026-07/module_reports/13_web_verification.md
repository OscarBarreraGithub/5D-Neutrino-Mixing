# Report 13 — Online verification of dated experimental anchors (2026-07-05)

Checks performed against live web sources for the anchors the offline audit could not corroborate.

## Verdicts

### CR001 "CMS2026 5.5 TeV g_KK edge" — CONFIRMED GENUINE (audit NOTE refuted)
CMS-B2G-25-009 = arXiv:2603.23454, "Search for new particles decaying into top quark-antiquark pairs in pp collisions at √s = 13 TeV" (2016–2018, 138 fb⁻¹; 0/1/2-lepton combination). A Randall–Sundrum Kaluza–Klein gluon is excluded for masses **0.5–5.5 TeV** at 95% CL — the most stringent tt̄-final-state limits to date. The CR001 active HARD edge of 5.5 TeV is therefore a real observed exclusion, not a projection. (Z′ widths 1/10/30%: excluded to 4.8/6.2/7.4 TeV.)
→ Revises Report 11's low-confidence NOTE: anchor provenance is valid. The separate benchmark-model-mismatch concern for CR005/CR006 (SSM limits on √L-suppressed states) is unaffected.

### K004 "NA62 Moriond-2026 preliminary 9.6e-11" — CONFIRMED GENUINE (value real; policy point stands)
NA62 presented at Moriond 2026 (arXiv:2604.12649, 2605.02415): 2023–2024 data BR = (7.2+2.3−2.1)×10⁻¹¹; **2016–2024 combined BR = (9.6+1.9−1.8)×10⁻¹¹**, compatible with the SM within 1σ at <20% precision. The K004 anchor value is genuine. The audit's remaining point is policy only: whether a HARD veto should anchor on a preliminary combination vs the published JHEP 02(2025)191 observation (13.0e-11); with the (separately flagged) 1σ budget, the choice moves the pass/fail boundary by ~1.2σ.

### EW001 S,T anchors — AUDIT CLAIM SOFTENED (values plausibly current; floor mismatch stands)
The compendium's C-3 asserted the positive centrals (S=+0.026±0.075, T=+0.047±0.066, ρ=0.90) "match no published PDG U=0 fit". Live check: the Nov-2025 electroweak precision review (Reina & Silvestrini, arXiv:2511.16534) quotes a current U=0 fit **S = 0.05±0.09, T = 0.08±0.07, ρ = +0.91**, and PDG-2021-era U=0 fits give S=0.05±0.08, T=0.09±0.07, ρ=0.92 — positive centrals of this size are standard in the current literature. The reviewer's quoted "PDG 2024: S=−0.05, T=0.00" could not be corroborated and may itself have been wrong. C-3 therefore reduces to its numerically solid core:
- the shipped anchors + coefficients solve to a **15.96 TeV** M_KK floor, vs the documented "~18–20 TeV" — this doc-vs-code mismatch is independent of anchor provenance and still needs resolution;
- the exact triplet's citation ("PDG 2025 Table 10.8") remains unpinned to a checkable source, and the "Gfitter year: 2026" mislabel stands.

### NuFIT "6.1 (Nov 2025)" — EXISTS; centrals plausible, not fully pinned
A NuFIT 6.1 release exists; its published 3σ ranges (e.g. Δm²₂₁ ∈ [7.236, 7.823]×10⁻⁵, sin²θ₁₂ ∈ [0.289, 0.330]) contain the repo values (7.537e-5, 0.3088). Post-JUNO 2025 global-fit estimates (δm² = 7.48±0.10 ×10⁻⁵, sin²θ₁₂ = 0.3085±0.0073) sit ~0.5σ from the repo centrals. Verdict: the repo's attribution is plausible (Report 01's "possibly fabricated" is downgraded); the five centrals should still be pinned against the actual NuFIT 6.1 table (nu-fit.org blocked our fetch) before publication use.

### τ→3e limit 2.5e-8 — PLAUSIBLE (source exists)
Belle II published "Search for the lepton-flavor-violating τ⁻ → e∓ℓ±ℓ⁻ decays" (JHEP 12(2025)169, arXiv:2507.18236; 428 fb⁻¹, 2019–2022), which covers τ→eee and supersedes Hayasaka 2010 (2.7e-8). The repo's 2.5e-8 is consistent with being that paper's limit; exact digit not extracted from the abstract — pin when citing.

## Sources
- https://arxiv.org/abs/2603.23454 (CMS-B2G-25-009)
- https://arxiv.org/abs/2604.12649 , https://arxiv.org/pdf/2605.02415 (NA62 2026)
- https://arxiv.org/abs/2511.16534 (Reina–Silvestrini EW precision review)
- https://arxiv.org/pdf/2204.03796 , https://arxiv.org/pdf/2604.03382 (U=0 oblique fits)
- https://arxiv.org/abs/2410.05380 (NuFit 6.0), JUNO-era update https://arxiv.org/pdf/2511.21650
- https://arxiv.org/html/2507.18236 , https://link.springer.com/article/10.1007/JHEP12(2025)169 (Belle II τ→eℓℓ)
