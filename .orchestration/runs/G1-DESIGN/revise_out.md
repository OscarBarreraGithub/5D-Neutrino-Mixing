Updated [.orchestration/rs_ew_sector_design_CONSENSUS.md](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/.orchestration/rs_ew_sector_design_CONSENSUS.md:1). No production code changed.

Changelog:
- codex #1: neutral contacts now use full product-minus-SM, including LFV `g_q^SM * delta_g_l`.
- Opus a: added explicit `KNOWN_EXTRA_KEYS` requirement for all new required/optional extras.
- Opus b + codex #6: rare Wilsons now bypass `_wilson_prefactor`; contacts are `GeV^-2` and already carry `1/m_Z^2` or `1/M_V^2`.
- codex #2-#5: added vector-path Zbb sign test, corrected universal-profile tests, fixed lepton Yukawa shapes, and added typed contact/Wilson key schemas.
- Arbitrations: numerical overlap first; closed-form `a(c)` deferred; gauge-only Zbb partial; `minimal_rs` default diagnostic; custodial/BKT and dipole loops human-gated.

Contact formula: `C_AB^Z = (g_Z^2/m_Z^2)[(g_qA^SM delta_ij + delta_g_qA_ij)(g_lB^SM delta_ab + delta_g_lB_ab) - g_qA^SM g_lB^SM delta_ij delta_ab]`, plus heavy `sum_V g_Vq g_Vl/M_V^2`.

Scope counts: FULL 26 after phase-6 fermion-KK Zbb, PARTIAL 17, HUMAN 7; before phase 6, FULL is 24.

HUMAN-INPUT items: custodial/BKT model choice; dipole-loop matcher/normalization; EDM basis, CP inputs, and matrix-element promotion criteria.