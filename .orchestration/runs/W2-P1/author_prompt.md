# W2 PHASE 1 — derivation pins (codex, gpt-5.x xhigh). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. Implement PHASE 1 of the APPROVED design `.orchestration/rs_ew_sector_design_CONSENSUS.md` (§5 item 1, and §6 normalization pin). Phase 1 is DERIVATION NOTES only — NO production code yet. First a short plan, then author. A codex reviewer and an Opus reviewer will independently verify the PHYSICS (both must APPROVE).

SCOPE (exactly as the approved design specifies for Phase 1):
Add derivation notes for: normalized gauge profiles in the RS slice; the zero-mode profile `g_0(c,z)`; the gauge-KK overlap integral that defines `a(c)` (the fermion-overlap form factor entering the Z-coupling shift `δg_f/g_Z ∝ (m_Z²/M_KK²)[a(c_f) − a_ref]`); the vector (gauge-boson KK) diagonalization; and the Z/W coupling normalization (g_Z, cosθ_W, the parentheses convention `L_Z = g_Z Z_μ f̄ γ^μ(g_L P_L + g_R P_R) f` used by `quarkConstraints/zpole.py`). These notes PIN the numerical overlap `a(c)` and the normalization that Phase 2 will implement.

GROUND IN: `derivations/` (match its existing LaTeX/markdown style + conventions — read `derivations/conventions.tex` and the kk_modes/zero_modes notes the design references), `warpConfig/wavefuncs.py` (f_IR/f_UV), `solvers/bessel.py` (KK roots; gauge NN root x_1≈2.405), `CLAUDE.md` (c=M_5/k, ε, overlap formulas), and `quarkConstraints/zpole.py` (the coupling convention you must match). Cite standard RS references (Agashe et al.; Casagrande et al. 0807.4937; Buras-Duling-Gori 0905.2318) where used.

REQUIREMENTS:
- Write the derivation to an appropriately-named file under `derivations/` (match repo conventions; e.g. `derivations/rs_ew_gauge_kk_coupling.tex` or the directory pattern already in use). Do NOT modify existing derivations; add new.
- The notes must be precise enough that Phase 2 can implement `a(c)` numerically and reproduce: (i) universal-c ⇒ all δg=0 (the a_ref subtraction), and (ii) IR-localized b_R ⇒ δg_R^b < 0, ~1e-3 at M_KK~3 TeV. State these explicitly as the checks Phase 2 must pass.
- Be honest about what is NOT derived here (closed-form a(c) is deferred to a future derivation per the design; Phase 1 grounds the NUMERICAL overlap + normalization only).
- No production code, no test changes in Phase 1.

OUTPUT (<=14 lines): short plan, the file path written, the key results (the a(c) overlap definition, the normalization convention, the a_ref subtraction, the two numerical checks Phase 2 must satisfy), and what is explicitly deferred. End with: PHASE1-DONE.
