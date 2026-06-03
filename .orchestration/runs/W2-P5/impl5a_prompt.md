# W2 PHASE 5 — SUB-STEP 5a (codex, gpt-5.x xhigh). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. Implement EXACTLY sub-step 5a of the DUAL-APPROVED plan `.orchestration/runs/W2-P5/plan.md`, building on the 3a/4a builder (committed b936c99: `quarkConstraints/rs_ew_spectrum.py`, `rs_ew_couplings.py`, `flavor_catalog_constraints/rs_ew_builder.py`). First a SHORT plan, then implement code+tests. Codex + Opus dual-review (both must APPROVE). **5a does NOT rewire any constraint** — charged-current builder + W mass-matrix/W' + δG_F + charged contacts + 5a gate tests ONLY.

BUILD (per the approved plan):
- New `quarkConstraints/rs_charged_current.py` with frozen `RSChargedCurrentCouplings`; add `rs_charged_current` to `KNOWN_EXTRA_KEYS`; `build_rs_ew_extras` populates it after rs_ew_spectrum + lepton inputs (raise loudly if c_L/lepton inputs absent while charged-current requested).
- W mass matrix on the gauge-NN tower (EXACT `solve_kk("gauge","NN",exact=True)` roots, x_1≈2.4505, INCLUDE the zero mode): `M²_{mn} = m_n² δ_{mn} + g₂² v² χ_m(1)χ_n(1)/4`; diagonalize → `m_w_gev` (zero eigenvalue→physical W), `m_wprime_gev` (first KK). State where the charged diagonalization lives.
- `eta_W` (the W-mixing sign coefficient — NOT named s_W) PINNED from the W mass-eigenvector (charged-current analog of neutral s_z=-1.0); include a test asserting its sign.
- `delta_g_W_ud_L = eta_W (m_W²/M_KK²) U_L_u† diag(a(c_Q) − a_ref) U_L_d` relative to V_CKM (left-handed); `delta_g_W_ud_R = 0` WITH a status string; `delta_g_W_lnu_L` from c_L (charged-lepton/PMNS basis), `delta_g_W_lnu_R=0`/status. Use the SHARED EW-universal `a_ref = spectrum.a(DEFAULT_A_REF_C=0.65)` (same object, do NOT recompute).
- `charged_contact_LL` (W' LL contact tensors, GeV^-2, from the charged tower mode sum, L_W=g₂/√2 convention).
- `delta_G_F_over_G_F`: total muon-decay amplitude shift (light-W e/μ vertex + W' exchange) in the same C_SM convention; raw components in diagnostics; the value to be subtracted EXACTLY ONCE by consumers (analog of a_ref absorption).
- Per-channel amplitude shift `epsilon_ij^a = delta_g_W_ud_L[i,j]/V_ij + delta_g_W_lnu_L[a,a] + charged_contact_LL[i,j,a,a]/C_SM − delta_G_F_over_G_F`. Field names per consensus §2.

5a GATE TESTS (recompute independently, NO rewiring):
- Exact W' root (≈2.4505·Λ_IR-ish from the mixed matrix); zero mode → physical W mass sane.
- eta_W sign pinned from eigenvector (assert sign).
- universal-c / a(c)=a_ref ⇒ all delta_g_W=0, all epsilon=0 (after δG_F subtraction) ⇒ SM recovered.
- delta_G_F/G_F magnitude/sign sanity at M_KK~3 TeV ((m_W/M_KK)²~7e-4 scale); subtracted exactly once (no double-count); no second 1/M_KK².
- Hermiticity/shape/finite/determinism; extras immutable. Independent recompute of one delta_g_W_ud_L entry + one epsilon.

CONSTRAINTS: ADD the charged-current module + extend builder + new tests; do NOT rewire constraints or touch the 103 constraint files (beyond KNOWN_EXTRA_KEYS). `python -m pytest tests/ -q` stays green + new tests pass. Numeric outputs real finite floats (complex matrices in extras/diagnostics).

OUTPUT (<=16 lines): short plan; new module + API; ACTUAL asserted numbers (m_w, m_wprime, eta_W sign, a sample delta_g_W_ud_L + epsilon, δG_F/G_F, universal-c=0 check); pytest counts. End with: P5A-DONE.
