# Quark-only scan: Z→bb made live + collider KK-mass searches added (2026-06)

This note records two changes to the quark-only (`--quark-only`) constraint set,
what was done, why, the decisions taken, and how each was verified. Both were
implemented and reviewed under the dual-signoff gate (codex *gpt-5.x xhigh* +
Opus, both APPROVE) and committed in `4c73dd7`. Custodial work is explicitly
**deferred** (see end).

## Context

The quark-only mini-scan evaluates the Bucket-1 constraints that need **zero
swept-lepton data**. Before this change the allowlist was 37 constraints, Z→bb
was held non-vetoing, and the direct collider KK-mass searches were excluded.

## Change 1 — Z→bb is now a live (vetoing) constraint

**What was wrong (the reframing).** We initially thought Z→bb was missing its
dominant `m_t²/M_KK²` term and needed new code. That was incorrect. In **minimal
non-custodial RS** the dominant correction to the left-handed b coupling is the
**gauge-KK profile term** (Casagrande et al, arXiv:0807.4937), whose size is set
by the third-generation doublet overlap `F(c_Q3)²` — large because that doublet
is IR-localized to give the top its mass. So it is top-*driven* through
compositeness, not a literal `m_t²` fermion factor, and it was **already
computed** in the scan: it lives in `z_delta_g_L/R_d[2,2]` via the exact KK-tower
overlap `a(c)`, alongside the small `m_b²` Casagrande admixture. Adding the
analytic `m_t²` term would have **double-counted**.

**What was actually done.** T010/T011 already computed real `R_b`/`A_b` pulls but
were held non-vetoing solely by an unconditional custodial `needs_human` flag
(mapped to `partial` → advisory). The fix re-tags only: the minimal piece is now
reported as `minimal_rs_tree_veto_ready` → tagged **rigorous** and vetoing
(failing points route to `excluded_by_rigorous`), while the custodial /
top-partner `Zb_L` / brane-kinetic-term refinements are recorded as `*_deferred`
diagnostics that **do not** suppress the veto. **No coupling or pull math changed**
— the R_b/A_b predictions are byte-identical (regression-pinned: T010
predicted 0.21328 / ratio 4.47042; T011 0.93832 / 0.08827). Files: `T010.py`,
`T011.py`. Provenance of the verification: `.orchestration/runs/ZBB*`
(codex propose → codex review → Opus sign-off, all agreeing) and
`.orchestration/runs/ZBB-RETAG`, `.orchestration/runs/COLLIDER` (implementation +
combined dual review).

**Decisions taken (orchestrator, recorded):**
- `A_FB^{0,b}` (the observable carrying the long-standing 2.8σ SM tension) stays
  **context-only**, not a veto scalar — the existing, defensible design.
- The Casagrande chiral cross-assignment (singlet→g_L, doublet→g_R) is kept as
  documented.

## Change 2 — Lepton-free collider KK-mass searches added

The direct LHC KK-mass searches (`collider_rs`, CR001–CR014) were excluded from
the quark-only allowlist. The **lepton-free** ones are now included, so the scan
reflects the direct-search M_KK floor. Allowlist: **37 → 46**.

- **IN (proxy / HARD, 9):** CR001 (g_KK→tt̄), CR002/CR003/CR004/CR008/CR010
  (vector-like-quark pair searches), CR007 & CR013 (KK graviton → WW/ZZ, γγ),
  CR012 (KK spin-1 → WW/WZ/ZZ). Each verified to read only quark/geometry extras
  (`kk_gluon_mass_gev`, `kk_ew_mass_gev`, `quark_mass_basis_couplings`); the
  quark-only point already supplies these; no silent SM/zero default.
- **OUT (deferred, lepton-involving):** CR005 (dilepton resonance), CR006
  (W′→ℓν), CR009 (ℓℓqq contact), CR011 (VBS, INFO/non-vetoing), CR014 (four-top
  two-lepton).

These are σ×BR / mass-exclusion **proxy** recasts that act as a near-flat
direct-search floor on M_KK. Smoke check: at M_KK = 1 TeV CR001 (proxy) excludes;
by 30 TeV it passes.

## Net effect on the quark-only scan

At M_KK = 1 TeV the excluding constraints now include the rigorous ΔF=2 set
(ε_K = K001, Δm_s = B003, φ_s = B004) **plus the now-live Z→bb (T010/T011,
rigorous)** and the collider direct searches (CR001 proxy, …). By ~30 TeV all
clear. Full suite: 1716 passed / 1 skipped.

## Still deferred (not in this change)

- **Custodial relaxation** of both S,T,U (currently the non-custodial estimate)
  and Z→bb (the protected top-partner `Zb_L`). This needs a committed custodial
  structure / fermion representation and is a separate decision.
- The full lepton sector (μ→eγ etc.) remains out of the quark-only scan.
