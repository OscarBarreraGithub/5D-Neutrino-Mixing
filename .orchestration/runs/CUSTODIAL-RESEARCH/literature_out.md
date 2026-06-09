# Custodial protection in RS (Z->bb and oblique S,T) — web literature survey (Opus, a6c25146b9cb38d49)

## The two problems custodial solves
Minimal RS (SM gauge in bulk, brane Higgs) generates large tree-level **ΔT** (gauge wavefunction shapes + KK mixing; Csaki-Erlich-Terning hep-ph/0203034), scaling ~ (v^2/M_KK^2)·L with L=ln(M_Pl/TeV)~37 — the volume log makes T dangerous. And a sizable shift to **Zb_L** (third-gen doublet IR-localized for the top). Together: minimal M_KK ~ O(5-10 TeV).

## Mechanism
- **Taming T:** bulk SU(2)_L×SU(2)_R×U(1)_{B-L} (Agashe-Delgado-May-Sundrum hep-ph/0308036); SU(2)_R = gauged custodial isospin; Higgs is (2,2) of SU(2)_L×SU(2)_R; leading volume-enhanced ΔT cancels in the custodial limit.
- **Protecting Zb_L:** SU(2)_R alone does NOT protect Zb_L. Need discrete **P_LR** (L<->R parity). ACDRP hep-ph/0605341 protection theorem (Eq.3): Z coupling protected iff **T_L=T_R AND T^3_L=T^3_R**. b_L gets it via T_L=T_R=1/2, T^3_L=T^3_R=-1/2 inside a **(2,2)_{2/3} bidoublet** (Eq.7-8). Full custodial group SU(2)_L×SU(2)_R×U(1)_X×P_LR (Carena-Ponton-Santiago-Wagner hep-ph/0701055).

## Representations (3rd gen)
- Q_L=(t_L,b_L) in **(2,2)_{2/3}** bidoublet — load-bearing for Zb_L protection; brings custodial partners incl. charge-5/3 state.
- t_R: **(1,1)_{2/3}** singlet (minimal), or **(1,3)+(3,1)** (control ΔT/top couplings; putting singlet-top in same multiplet as doublet disfavored -> large negative one-loop ΔT, Carena et al §II).
- b_R: several embeddings (ACDRP Table 1) — choice sets sign/size of the UNPROTECTED δg_bR (can fix A_FB^b anomaly with δg_Rb~+0.02).
- KEY TRADEOFF (ACDRP §3): protecting Zb_L is INCOMPATIBLE with simultaneously protecting Zt_L and W_tb (t_L,b_L share the SU(2)_L doublet) -> residual Zt_L, W_tb corrections.

## With vs without custodial (KK-scale bounds)
- Minimal RS, T-only (light Higgs): M_KK >~ 4 TeV; full fit (Casagrande/Neubert): >~5 TeV; conservative/T-driven: ~O(10 TeV). NON-CUSTODIAL IS A RANGE ~4-10 TeV (Higgs-mass + fit dependent), not one number.
- **Custodial (RSc): M_KK ~ 2-3 TeV** (Carena et al: m^1_gauge >~ 2.5 TeV). Custodial lowers the EW bound to LHC reach.

## Residual / unprotected effects (must keep in the fit)
1. **δg_bR is NOT protected** (P_LR protects only g_bL). [Subtlety to NOT conflate: in MINIMAL RS the dominant tree shift is g_bL with g_bR small; in CUSTODIAL RS g_bL is killed and the residual concern is g_bR + loop g_bL.]
2. **Top sector (Zt_L, W_tb) unprotected** (the §3 incompatibility).
3. **S parameter custodially UNPROTECTED** (custodial fixes T, not S) — sets part of the residual ~few-TeV bound.
4. **Calculable, correlated one-loop ΔT and δg_bL** from top/custodian sector (Carena et al hep-ph/0701055 Eq.28-30; finite because symmetries forbid 5D counterterms; correlated per their Fig.1).

## Schematic formulas (as reported, not invented)
- Custodial T: T ~ -[pi v^2/(4 cos^2θ_W M_KK^2)]·(1/L) (review 0910.4876) — 1/L suppression vs minimal.
- δg_bL ∝ F^2(c_Q3)·(v^2/M_KK^2)·[log], δg_bR ∝ F^2(c_d3)·(v^2/M_KK^2); under P_LR the LH shift "basically vanishes."
- Protection: δg=0 iff T_L=T_R, T^3_L=T^3_R (ACDRP Eq.3-7).
- One-loop T, δg_bL: Carena et al Eq.28-30.

## Sources
hep-ph/0203034 (Csaki-Erlich-Terning); hep-ph/0308036 (Agashe-Delgado-May-Sundrum); hep-ph/0412089 (Agashe-Contino-Pomarol MCHM); **hep-ph/0605341 (ACDRP custodial Zbb — THE key paper)**; hep-ph/0701055 (Carena-Ponton-Santiago-Wagner EW custodial); 0807.4937 (Casagrande et al RS flavor I, §6 minimal Zbb); 0903.2415 (RSc EW+flavor w/ custodial); 0910.4876 (RS precision review); 0912.1625 (Bauer-Casagrande-Haisch-Neubert RS II); 1310.1070 (Snowmass).

## Caveats
Non-custodial bound is a RANGE ~4-10 TeV (not one number). Minimal-RS tree g_L^b/g_R^b from 0807.4937 §6 not pulled verbatim (render truncation) — canonical ref for exact expressions. Don't conflate "g_bR small in minimal RS" with "g_bR unprotected in custodial RS."
