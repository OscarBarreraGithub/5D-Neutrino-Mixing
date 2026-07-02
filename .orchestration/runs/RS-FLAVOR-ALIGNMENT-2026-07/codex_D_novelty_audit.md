**Analysis/Spec**

| Claim | Known? | Closest reference | What we add |
|---|---:|---|---|
| Anarchic RS `epsilon_K` dominated by KK-gluon LR `C_4` | Yes | Csaki-Falkowski-Weiler [0804.1954](https://arxiv.org/abs/0804.1954), Blanke-Buras-Duling-Gori-Weiler [0809.1073](https://arxiv.org/abs/0809.1073), Agashe-Azatov-Zhu [0810.1016](https://arxiv.org/abs/0810.1016) | Mostly re-derivation in this repo’s convention. |
| Generic anarchic RS has a CP/EDM problem | Yes | Agashe-Perez-Soni [hep-ph/0408134](https://arxiv.org/abs/hep-ph/0408134); Keren-Zur et al. [1205.5803](https://arxiv.org/abs/1205.5803) expects neutron EDM near sensitivity in partial compositeness | No novelty unless tied to a new scan/observable discriminator. |
| S2 / RH-down flavor protection suppresses kaon CP | Yes | Santiago minimal flavor protection [0806.1230](https://arxiv.org/abs/0806.1230); Csaki-Falkowski-Weiler simple protection [0806.3757](https://arxiv.org/abs/0806.3757) | Our distinction: these suppress down-sector flavor violation/magnitude, not just CP phase. |
| FPR / 5D-MFV solves RS flavor by tying flavor sources | Yes | Fitzpatrick-Perez-Randall [0710.1869](https://arxiv.org/abs/0710.1869); Cacciapaglia et al. GIM [0709.1714](https://arxiv.org/abs/0709.1714) | Our distinction: FPR-like models reduce/specialize flavor violation; they do not keep anarchic down `|C_4|`. |
| Flavor alignment via down-Yukawa/bulk-mass alignment | Yes | Csaki-Perez-Surujon-Weiler [0907.0474](https://arxiv.org/abs/0907.0474) | Again, known as magnitude/flavor alignment, not CP-only alignment. |
| Spontaneous CP / strong-CP sequestering in warped space | Yes | Cheung-Fitzpatrick-Randall [0711.4421](https://arxiv.org/abs/0711.4421); Harnik-Perez-Schwartz-Shirman split fermions [hep-ph/0411132](https://arxiv.org/abs/hep-ph/0411132) | Not new to use warped geometry plus spontaneous CP. Their kaon footprint is not “large anarchic `Delta m_K` but small phase.” |
| Nelson-Barr in warped/composite models | Yes | Girmohanta-Lee-Nakai-Suzuki [2203.09002](https://arxiv.org/abs/2203.09002); recent doublet-NB Alves-Nishi-Vecchi [2604.02506](https://arxiv.org/abs/2604.02506) | Strong prior art. The 2022 paper even notes CP phases can be moved into the up sector, but uses flavor symmetries/diagonal Yukawas rather than anarchic down magnitudes. |
| CP-only / “Minimal CP Violation” in composite/RS | Yes, very close | Redi-Weiler [1106.6357](https://arxiv.org/abs/1106.6357) | This is the sharpest prior art: real strong sector, CP from elementary mixings, third-generation phase, EDMs zero at leading order, Im `C_4^K` suppressed while Re remains. Do not claim this broad idea as new. |
| `(Delta m_K, nEDM)` plane as mechanism discriminator | Partly known conceptually | Blanke/Buras-style “DNA” flavor correlations; Barbieri et al. `U(2)^3` [1203.4218](https://arxiv.org/abs/1203.4218); Redi-Weiler EDM/flavor discussion | A focused RS-kaon discriminator separating magnitude-alignment vs CP-only alignment looks defensibly new if you show it quantitatively. |
| Empirical statement: anarchic `epsilon_K` survival is magnitude-dominated with codim-1 phase tail | I did not find this published in this form | Related scans/fine-tuning in Blanke et al. [0809.1073](https://arxiv.org/abs/0809.1073) | The scan result and survival-mechanism taxonomy look new. The `P_pass ~ 1/R_K` explanation is an elementary analytic reinterpretation, not a deep new mechanism. |

**Honest Novelty Boundary**

Do **not** claim novelty for “solve RS `epsilon_K` by CP alignment” in broad form. Redi-Weiler already did a CP-only/minimal-CP version in composite Higgs/RS language, including suppressed Im `C_4^K` and EDM relief.

Do **not** claim novelty for “warped Nelson-Barr” or “single-sector/up-sector CP transmission” without qualification. Girmohanta et al. have a warped NB construction with flavor symmetries and movable CP phases.

The defensible new lane is narrower:

1. **NB-compatible RS phase-only kaon alignment with anarchic down magnitudes**: unlike S2/FPR/shining, retain large LR `|C_4|` and hence potentially large `Delta m_K`; unlike Redi-Weiler, embed the CP-only idea in a Nelson-Barr-like strong-CP-safe construction; unlike Girmohanta et al., make anarchic RS kaon phenomenology central.

2. **Phenomenological taxonomy/diagnostic**: `Delta m_K` high + neutron EDM low is the proposed smoking-gun region for CP-only alignment, while S2/FPR-like magnitude alignment gives low `Delta m_K`, and generic anarchy gives high EDM / high `epsilon_K`.

**Sharpest Contributions To Stake**

1. A quantitative RS scan showing that anarchic `epsilon_K` survivors are mostly small-`|C_4|`, with a small codimension-one phase-aligned tail, and explaining this by the `P_pass(R_K) ~ 1/R_K` weighting.

2. A model/ansatz and observable discriminator for **CP-only down 1-2 alignment**: down-sector `G_L^{12}G_R^{12}` real enough to pass `epsilon_K`, while its magnitude remains anarchic enough to saturate `Delta m_K`, with neutron EDM parametrically suppressed by the same CP structure.

**Risks That Could Sink Claims**

Redi-Weiler is the main sink. If your claim reads like “Minimal CP Violation in partial compositeness suppresses `epsilon_K` and EDMs,” it is already published.

Girmohanta et al. is the main sink for “warped NB with phases moved to the up sector.” Your differentiator must be anarchic down magnitudes plus kaon/EDM phenomenology.

Radiative stability is nontrivial: CKM/up-sector CP can regenerate down 1-2 phases through brane kinetic terms, KK loops, Yukawa running, or NB-sector thresholds. This needs an explicit spurion proof.

`Delta m_K` is a messy discriminator because SM long-distance uncertainties and sign/cancellation issues are large. Treat it as a plane/region, not a precision observable.

The empirical 57%/7% numbers are prior- and cut-dependent. Present them as results for the specified Bauer S1 scan, not universal RS facts.