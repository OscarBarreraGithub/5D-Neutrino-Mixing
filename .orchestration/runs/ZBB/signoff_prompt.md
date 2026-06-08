OPUS SIGN-OFF (final gate before surfacing to the physicist). Work SYNCHRONOUSLY; end with a verdict line. Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing.

Two codex passes converged on a SURPRISING conclusion about the Z to bb constraint. Independently verify it before it goes to the physicist. Read `.orchestration/runs/ZBB/proposal_codex_out.md` and `.orchestration/runs/ZBB/review1_codex_out.md`. The claim:
1. In MINIMAL non-custodial RS, the dominant left-handed Zb_L correction is the GAUGE-KK profile term (Casagrande arXiv:0807.4937, "ZbbRS"): delta g_bL ~ g_bL^SM * [-m_Z^2/(2 M_KK^2) * F(c_Q3)^2/(3+2c_Q3) * (L - ...)]. Its SIZE is set by F(c_Q3)^2 (large because the third-gen doublet is IR-localized to get the top mass), so it is top-DRIVEN through compositeness, NOT a literal m_t^2 fermion term.
2. The m_f^2 B(c) fermion-admixture for Zbb is down-type (m_b^2 B(c_d)) and is already coded; m_t^2 B(c_u) belongs to Ztt, not Zb_L.
3. The literal large m_t^2/M_KK^2 top-partner Zb_L term is custodial/representation-dependent (Agashe-Contino-Da Rold-Pomarol hep-ph/0605341), hence deferred, not universal.
4. CODE: the gauge Zb_L profile term is ALREADY in z_delta_g_L_d[2,2] via the exact KK-tower overlap a(c) (rs_ew_couplings.py ~505-551; rs_ew_spectrum.py ~789-917), and the m_b^2 admixture is separately coded (rs_ew_couplings.py ~704-744), on in the scan. Therefore adding the analytic term would DOUBLE-COUNT.
5. RECOMMENDATION: nothing to add; minimal non-custodial Z to bb is ALREADY computed; the only action is to RE-TAG it from `partial` (non-vetoing) to a live minimal-complete vetoing constraint, keeping custodial/top-partner/BKT variants deferred.

INDEPENDENTLY CHECK:
- Is the physics right (is the minimal Zb_L correction gauge-driven via F(c_Q3)^2, and is the literal m_t^2 piece genuinely custodial/representation-dependent rather than a universal minimal term)? Spot-check against the cited equations.
- Is the CODE claim right: read the cited lines yourself and confirm (a) z_delta_g_L_d[2,2] already contains the third-gen-doublet gauge profile shift via a(c) (exact tower), (b) the m_b^2 admixture is separate and active, (c) adding the proposed analytic term WOULD double-count.
- Is the recommendation sound: is re-tagging (partial -> vetoing for the minimal piece, custodial deferred) the correct, honest action, or is there a real gap still missing for a minimal non-custodial Z to bb?
- Flag anything both codex passes may have gotten wrong or glossed.

OUTPUT (<=14 lines): per-point AGREE/DISAGREE with your own evidence (cite file:line for the code); the decisive answer to "is minimal Z to bb already computed and just mis-tagged?"; whether the right action is re-tag (no new physics) vs add-a-term; any caveat the physicist must know. End with EXACTLY ONE line: `ZBB-SIGNOFF: APPROVE` or `ZBB-SIGNOFF: NEEDS-WORK`.
