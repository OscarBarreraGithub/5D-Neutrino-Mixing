Updated `.orchestration/runs/W2-P3/plan.md` only; no production code changed.

- Lepton side: Phase 3 sets charged-lepton `delta_g_l=0` but uses SM `g_{L,R}^{l,SM}` / `g_{V,A}^{l,SM}` in light-Z contacts, so stored `rs_semileptonic_wilsons` includes `z_delta_g_q * g_l^SM`; heavy-vector `g_Vl` terms defer to Phase 4.
- `z_delta_g`: pinned as `s_Z * g_A^{f,SM} * (m_Z^2/M_KK^2) * [U^\dagger diag(a(c_f)-a_ref) U]`, not raw overlap.
- 3a gate: full-subtraction SM-limit and IR-`b_R` sign/magnitude pins now run before rewiring.
- Rare B: B017/B018/B019 preserved as LFU ratio observables; B015 keeps existing C7 proxy, replacing only C9/C10/C9p/C10p.
- Rare-B files enumerated: `rare_b_meson.py`, `rare_b_electronic.py`, `rare_b_tauonic.py`, `rare_b_kstar_dilepton.py`, `rare_b_baryon.py`.