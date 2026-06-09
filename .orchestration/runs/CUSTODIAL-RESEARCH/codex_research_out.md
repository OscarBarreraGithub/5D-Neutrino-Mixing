**Custodial RS Design Doc**

**Scope**
This is research/design only. I did not change code.

The current repo implements a minimal non-custodial RS EW layer:

- [rs_ew_couplings.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/rs_ew_couplings.py:456): builds light-Z mass-basis shifts
  \[
  \delta g_A^q=s_Z g_A^{SM}\frac{m_Z^2}{M_{KK}^2}\,U_A^\dagger \mathrm{diag}(a(c)-a_{ref})U_A .
  \]
- [rs_ew_couplings.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/rs_ew_couplings.py:553): optionally adds minimal-RS Casagrande ZMA bottom fermion-KK admixture to `z_delta_g_L/R_d[2,2]`.
- [T010.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/top_higgs_ew/T010.py:522) and [T011.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/top_higgs_ew/T011.py:699): consume only `rs_ew_couplings.z_delta_g_L/R_d[2,2]`.
- [oblique_stu.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/oblique_stu.py:143): uses the aggressive minimal-RS proxy
  \[
  \Delta S=c_S v^2/M_{KK}^2,\quad
  \Delta T=\frac{\pi L}{2c_W^2}v^2/M_{KK}^2,\quad U=0.
  \]

**Literature Survey**
The canonical custodial RS gauge group is
\[
SU(2)_L\times SU(2)_R\times U(1)_X\times P_{LR},
\quad Y=T_R^3+X,
\]
up to sign conventions for \(T_R^3\). Agashe, Delgado, May, Sundrum introduced the enlarged RS electroweak gauge symmetry as a cure for the large \(T\) parameter: the Higgs sector preserves custodial \(SU(2)_V\), suppressing the isospin-breaking contribution to \(\rho\)/\(T\) [Agashe et al. 2003](https://arxiv.org/abs/hep-ph/0308036).

Agashe, Contino, Da Rold, Pomarol then showed that the same custodial structure can protect \(Zb_L\bar b_L\) if the BSM sector has \(O(4)\simeq SU(2)_L\times SU(2)_R\times P_{LR}\) and the fermion is a \(P_{LR}\) eigenstate with
\[
T_L=T_R,\qquad T_L^3=T_R^3 .
\]
For \(b_L\), this requires embedding the SM doublet in a bidoublet \(Q_L\in(2,2)_{2/3}\), with \(b_L\) carrying \(T_L^3=T_R^3=-1/2\). Then the leading non-universal \(Zb_L\bar b_L\) correction vanishes at zero momentum. This protection cannot simultaneously protect \(Zt_L\bar t_L\) and \(Wtb\); sizable top/charged-current deviations remain possible [ACDP 2006](https://arxiv.org/abs/hep-ph/0605341).

Carena, Ponton, Santiago, Wagner emphasized that custodial \(SU(2)_V\times P_{LR}\) makes the tree-level \(T\) and \(Zb_L\bar b_L\) problems small, but calculable one-loop top-partner effects remain and can be correlated and sizable. Their scans typically allow gauge KK masses around a few TeV, with potentially light exotic/vectorlike fermions [Carena et al. 2006](https://arxiv.org/abs/hep-ph/0607106), [Carena et al. 2007](https://arxiv.org/abs/hep-ph/0701055).

Casagrande, Goertz, Haisch, Neubert, Pfoh are the closest match to the repo’s current minimal-RS formulas. Their 2008 paper gives exact/zma minimal-RS gauge and fermion mixing contributions to \(Z\) couplings and compares minimal vs custodial oblique behavior. In their notation the custodial extended EW sector leaves \(S\) essentially of the same order but changes \(T\) from volume-enhanced to \(1/L\)-suppressed:
\[
T_{\rm cust}\propto -\frac{\pi}{4c_W^2 L}\frac{v^2}{M_{KK}^2}
\]
in the same relative normalization as the minimal expression [Casagrande et al. 2008](https://arxiv.org/abs/0807.4937).

Albrecht, Blanke, Buras, Duling, Gemmler give a full custodially protected warped flavour structure with all left-handed down-type \(Z\) couplings protected when the quark representations are chosen consistently across generations. This matters if the repo later wants custodial protection for FCNC \(Z d_L^i d_L^j\), not only diagonal \(Zb_L\) [Albrecht et al. 2009](https://arxiv.org/abs/0903.2415). Casagrande et al. 2010 then rederived custodial-RS gauge and fermion couplings in the mass basis, identifying irreducible \(P_{LR}\)-breaking residuals: the protection is real, but not “everything is exactly zero” once boundary conditions, fermion mixing, and loops are included [Casagrande et al. 2010](https://arxiv.org/abs/1005.4315).

**Representation Choices**
Recommended default for this framework:

| Field | Standard custodial choice | Protects | Tradeoff |
|---|---:|---|---|
| \(Q_L^i\) | \((2,2)_{2/3}\), ideally all generations | \(Z d_L^i d_L^j\), including \(Zb_L\) | Adds exotic \(Q=5/3,2/3\) partners; top/W couplings not protected |
| \(t_R\) | \((1,1)_{2/3}\), or SO(5) \(5_{2/3}\) realization | \(Zt_R\) via \(P_C\) | Singlet option is economical; loop spectrum still model-dependent |
| alternative \(t_R\) | \((1,3)_{2/3}\oplus(3,1)_{2/3}\), SO(5) \(10_{2/3}\) | also \(Zt_R\) at tree level | More custodians; loop \(T\)/\(Zb_L\) more sensitive |
| \(b_R\) | \((1,3)_{2/3}\) in common RSc constructions | permits bottom Yukawa with \(Q_L=(2,2)_{2/3}\) | \(Zb_R\) is not generically protected; affects \(A_b\) |
| alternative \(b_R\) | \((1,1)_{-1/3}\) plus separate bottom operator/multiplet | can make \(Zb_R\) small/protected | More model-building; can softly break \(Zb_L\) protection |
| light generations | same \(Q_L=(2,2)_{2/3}\) | protects left-handed down FCNC \(Z\) couplings | More exotics and \(SU(2)_R\) effects; cleaner for anarchic CKM |

The key physicist decision is whether this repo is modeling only the dominant diagonal \(Zb_L\) relaxation or the full custodial flavour structure. If only diagonal \(Zbb\), third-generation \(Q_L\) protection is enough. If the same `rs_ew_couplings` object will feed rare \(B/K\) semileptonic Wilsons, all-generation down-left protection is the more coherent choice.

**What Changes In This Repo**
For `Zbb`:

Current minimal path:
\[
\delta g_L^b =
(\delta g_L^b)_{\rm gauge}
+\frac{m_b^2}{2\Lambda_{\rm IR}^2}B_d,
\quad
\delta g_R^b =
(\delta g_R^b)_{\rm gauge}
-\frac{m_b^2}{2\Lambda_{\rm IR}^2}B_Q .
\]

Custodial \(P_{LR}\) path:
\[
\delta g_L^b \rightarrow \delta g_{L,\rm residual}^b \simeq 0
\]
for the ideal tree-level protected bidoublet embedding. The repo should not silently reuse the minimal Casagrande \(m_b^2 B_d\) term for protected \(b_L\); the custodial fermion spectrum is different. A strict proxy should set the protected `z_delta_g_L_d[2,2]` contribution to zero or to an explicit residual parameter, with diagnostics saying that top-partner loop and irreducible \(P_{LR}\)-breaking effects are omitted.

For \(b_R\), do not claim protection by default. Either set `delta_g_R_b` by a chosen `bR_strategy` or leave it as a separate human-input residual. The standard \((1,3)_{2/3}\) choice permits the bottom Yukawa but does not protect \(Zb_R\); this is directly relevant to \(A_b\) and \(A_{FB}^{0,b}\).

For `EW001`:
\[
\Delta S \text{ should remain unchanged by default}
\]
because custodial symmetry does not protect \(S\). The repo’s current \(c_S=30\) PDG-style aggressive proxy can stay unless the physicist supplies a custodial-specific \(S\) coefficient.

Replace the minimal \(T\) coefficient
\[
c_T^{\rm min}=\frac{\pi L}{2c_W^2}
\]
with a custodial tree-level coefficient preserving the literature suppression ratio,
\[
c_T^{\rm cust}=-\frac{\pi}{4c_W^2 L}.
\]
For \(L=35\), this changes \(c_T\) from about \(71.5\) to about \(-0.029\), a factor \(-1/(2L^2)\). Keep \(U=0\) at this proxy level.

**Ranked Implementation Options**
1. **Recommended: representation-aware custodial proxy switch.**
   Add a custodial config to the RS-EW builder: `ew_model="custodial_rs_plr"`, `qL_rep="bidoublet_2_2_X_2over3"`, `tR_rep`, `bR_rep`, `protect_scope=("bL_only"|"all_down_L")`, `bR_strategy`, and residual knobs. T010/T011 need no new physics code if `rs_ew_couplings.z_delta_g_L/R_d[2,2]` is populated correctly. EW001 gets a second oblique proxy with the custodial \(T\) coefficient. Captures the dominant intended relaxation. Missing: exact \(SU(2)_R\) gauge tower, custodial fermion masses, one-loop top partners, residual \(P_{LR}\) breaking. Cost: low-medium. Accuracy: good for scan-level “custodial protected vs not” classification, not final precision EW fit.

2. **Next-best: leading-log custodial neutral-current builder.**
   Implement custodial charge factors \(\omega_Z^q\), \(Z/Z_H/Z'\) mixing approximations, and representation-specific left/right down couplings following the custodial-RS mass-basis papers. This can protect all \(d_L\) entries consistently and improves rare-decay consistency. Inputs needed: \(g_R/g_L\), \(g_X\), \(Q_X\), representation tables, extra custodial bulk masses such as triplet partners, and residual \(P_{LR}\)-breaking choices. Missing: exact tower diagonalization and loops. Cost: medium-high. Accuracy: much better tree-level model.

3. **Gold standard: full custodial RS-EW sector.**
   Build \(SU(2)_L\times SU(2)_R\times U(1)_X\) gauge spectra, neutral/charged mixing, custodial fermion multiplet mass matrices, exact or truncated KK diagonalization, and one-loop top-partner \(T/Zb_L\) corrections. This is physically best but overlarge for the current repo inputs. Cost: high. Accuracy: publication-grade if validated.

**Open Physics Decisions**
The physicist must decide:

1. Whether \(Q_L\) protection applies only to the third generation or to all generations.
2. The third-generation embedding: \(Q_L=(2,2)_{2/3}\) is the standard \(Zb_L\) answer; choose \(t_R\) singlet vs triplet/SO(5) realization.
3. The \(b_R\) representation and whether \(A_b/A_{FB}^{0,b}\) should be treated as SM-like, fitted with a residual, or modeled explicitly.
4. Whether to include loop-level top-partner corrections to \(T\) and \(Zb_L\) now or explicitly defer them.
5. Whether the repo’s aggressive \(S\) coefficient should remain unchanged in custodial mode or be replaced by a model-specific custodial \(S\) input.

**Executive Recommendation**
1. Implement first a custodial proxy, not a fake full model: `ew_model="custodial_rs_plr"` with explicit representation metadata.
2. Use \(Q_L^i=(2,2)_{2/3}\) for all generations by default, \(t_R=(1,1)_{2/3}\), and \(b_R=(1,3)_{2/3}\) unless the physicist chooses a \(b_R\)-protected variant.
3. Set the protected `z_delta_g_L_d[2,2]` tree-level shift to zero or an explicit residual; do not reuse the minimal Casagrande \(b_L\) admixture.
4. Replace EW001 \(c_T=\pi L/(2c_W^2)\) with \(c_T=-\pi/(4c_W^2L)\); keep \(S\) and \(U=0\) unchanged by default.
5. Required new inputs: representation choice, protection scope, \(b_R\) treatment, \(L/M_{KK}/s_W^2\), and explicit flags for omitted top-partner loops/residual \(P_{LR}\) breaking.

CUSTODIAL-RESEARCH-DONE