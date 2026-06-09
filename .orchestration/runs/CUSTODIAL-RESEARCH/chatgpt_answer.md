# Expert (ChatGPT) answer to the custodial problem statement — TO BE VALIDATED

Use as the custodial branch: Q_L^i in (2,2)_{2/3}, b_L chosen as the P_LR eigenstate.
KEY IMPLEMENTATION CLAIM (contradicts our earlier Opus panel): do NOT keep the full
minimal-RS gauge-wavefunction contribution to delta g_L^b. In the custodial model the
leading Q_3-overlap-enhanced GAUGE piece is protected TOO. Keep only subleading
custodial-breaking residuals + loop effects.

## 1. Representations
- Q_L=(2,2)_{2/3}, b_L assigned as P_LR eigenstate => ACDRP T^3_L=T^3_R=-1/2 for b_L.
  In code, enforce by setting the protected Zb_L charge factor omega_Z^{b_L}=0 rather
  than relying on a sign convention. ACDRP: this protects Zb_Lb_L; Zt_Lt_L and Wtb NOT
  simultaneously protected.
- t_R in (1,1)_{2/3} (clean default), or SO(5) 5_{2/3}; triplets (1,3)+(3,1) more model-specific.
- b_R: conservative => mostly elementary, delta g_R^b ~= 0. SO(5) => b_R from 10_{2/3}
  containing (1,3)_{2/3} (Carena et al). A_FB^b-targeted => model delta g_R^b explicitly
  (positive ~+0.02 helps the anomaly, NOT automatic).

## 2. Scope: third gen vs all down-type left
- Minimal EW fix: protect Q_L^3 only.
- Anarchic flavor scan (cleaner): all Q_L^i in (2,2)_{2/3}, i=1,2,3, then rotate to mass
  basis — avoids reintroducing unprotected down-type left Z couplings via flavor rotations.
  Carena et al use the same rep for the first two gens.

## 3. delta g_L^b (THE KEY POINT)
Current: delta g_L^b = dgL_gaugeKK + dgL_fermionKK.
Custodial P_LR branch: ZERO THE LEADING PROTECTED PIECES OF BOTH TERMS. ACDRP separate
gauge-KK and fermion-KK corrections; for the P_LR eigenstate b_L the LEADING GAUGE piece
vanishes (the protected charge difference), and the fermion-mixing correction is also
killed when embedding+Yukawa respect the custodial condition.
So NOT: dgLb = minimal_gauge + 0. Instead:
  dgLb = dgL_custodial_residual + dgL_explicit_PLR_breaking + dgL_loop.
First-pass residual: delta g_L,residual^b ~ kappa_b (1/L) delta g_L,gauge,minimal^b,
kappa_b=O(1); or zero in the ideal P_LR limit with a flag for subleading irreducible
breaking. Leading L-enhanced term removed; subleading (BCs, imperfect symmetry) remain.
=> zeroing only the m_b^2/Lambda_IR^2 admixture while keeping the full minimal
F(c_Q3)-driven gauge piece would OVERCONSTRAIN the custodial scan.

## 4. Oblique
Replace Delta T_min = (pi L/2c_W^2) v^2/M_KK^2 by
  Delta T_cust = -(pi/(4 c_W^2 L)) v^2/M_KK^2,  Delta U=0.
=> suppressed by -1/(2L^2) (~ -1/2450 at L~35). S essentially unchanged:
  Delta S ~= 2 pi (v^2/M_KK^2)(1 - 1/L). So S, not T, is the dominant tree oblique
  constraint now. Check M_KK convention (1/R_IR vs first KK mass m_1 ~= 2.4 M_KK).

## 5. b_R, A_b, A_FB^b
Default: delta g_R^b = 0 (or small rep-computed residual) unless explicitly targeting
A_FB^b. A_b = (gL^2-gR^2)/(gL^2+gR^2). A_FB^b branch: separate param/rep model for
delta g_R^b (positive ~1e-2..2e-2 classic target, not automatic).

## 6. One-loop top-partner
Add the FLAG now; can defer numerics for first tree-level validation but do NOT claim
2-3 TeV viability without it. Carena et al: once tree T and Zb_L protected, one-loop
top-partner corrections relevant + correlated; viable gauge-KK ~2-3 TeV, custodians
sometimes few-hundred GeV. At M_KK~2-3 TeV: Delta T_loop ~ O(0.1)+, delta g_L^b|loop
~ 1e-3 (relevant for Zbb); sign model-dependent (bidoublet partners push T negative,
singlet-like positive).

## Bottom line
P_LR protects the leading GAUGE-WAVEFUNCTION contribution AS WELL AS the leading
fermion-mixing contribution. Zeroing only the m_b^2 admixture while keeping the full
minimal F(c_Q3) gauge piece would overconstrain and partially defeat the model.
