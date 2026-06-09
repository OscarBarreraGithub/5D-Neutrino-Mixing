# Custodial protection for our warped (RS) flavor scan — problem statement / questions

## Setup
We run a large parameter scan of a Randall–Sundrum warped extra dimension with
anarchic 5D Yukawas and bulk fermions localized by dimensionless masses `c`,
with a brane-localized Higgs and an IR/KK scale `M_KK`. We impose flavor and
electroweak-precision constraints on each point.

We currently use the **minimal (non-custodial)** RS gauge sector (SM gauge group
in the bulk). In that setup two electroweak observables dominate and force the
KK scale high:
- **Z → b b̄**: the shift to the left-handed b coupling, δg_L^b (the RS Zb_L problem); and
- the **oblique T parameter** (we track S, T, U).

In our scan these currently push the rigorous floor to ~10 TeV (in Λ_IR units).
We want to add the standard **custodial** structure
`SU(2)_L × SU(2)_R × U(1)_X × P_LR` so these relax toward ~2–3 TeV, and we need
to know exactly how to model it at the level of our scan's granularity.

## What we compute now (so you know the granularity)
- **Z → bb:** shifts δg_L^b, δg_R^b to the b-quark Z couplings, built from **two
  distinct pieces**: (i) a **gauge–KK wavefunction** piece (Z zero-mode profile /
  KK mixing, set by the third-generation doublet IR overlap `F(c_Q3)`), and
  (ii) a **fermion–KK mixing admixture** (a Casagrande ZMA `∝ m_b²/Λ_IR²` term).
- **Oblique:** a leading-order proxy, ΔS = c_S · v²/M_KK², ΔT = (πL / 2c_W²) ·
  v²/M_KK² with `L = ln(M_Pl/TeV) ≈ 35`, and U = 0.

Standard references we're working from: Agashe–Contino–Da Rold–Pomarol
(hep-ph/0605341, the P_LR protection of Zb_L) and Carena–Pontón–Santiago–Wagner
(hep-ph/0701055, custodial electroweak fit + the one-loop T / δg_bL).

## Questions we need answered

1. **Third-generation representations.** Is the standard custodial choice the
   right one here: Q_L = (t_L, b_L) in the **(2,2)_{2/3} bidoublet** of
   SU(2)_L × SU(2)_R (with P_LR) to protect Zb_L? And which embeddings for
   **t_R** (singlet (1,1)_{2/3} vs (1,3)⊕(3,1)) and **b_R** (which of the
   standard choices), given we care about A_b / A_FB^{0,b} and the *unprotected*
   δg_R^b?

2. **Protection scope.** Should custodial protection apply only to the third
   generation, or to all down-type left-handed quarks?

3. **The Zb_L decomposition — the key one for us.** P_LR protects the b_L
   neutral-current vertex from **fermion mixing**. In our two-piece δg_L^b, do we
   (a) zero **only** the fermion-mixing admixture and **keep** the gauge-
   wavefunction piece (reduced but nonzero — an irreducible P_LR-breaking
   residual), or (b) is the gauge piece also protected/negligible? Concretely:
   which piece do we set to zero, and what residual should we retain at tree level?

4. **The T parameter.** With custodial, what is the correct leading replacement
   for our ΔT coefficient? Is the `~1/L²` suppression
   (`c_T : πL/2c_W² → −π/(4 c_W² L)`) the right leading form and sign, and does
   **S stay essentially unchanged** (custodial protects T, not S)?

5. **b_R / A_b.** How should we treat the unprotected right-handed coupling and
   the A_b / A_FB^{0,b} observable — SM-like, a fitted residual, or modeled from
   the chosen b_R representation?

6. **Loop level.** The custodial setup has calculable, correlated **one-loop
   top-partner** contributions to T and δg_bL (Carena et al., Eqs. 28–30). For a
   scan at M_KK ~ 2–3 TeV, should we include these now or defer them with an
   explicit flag — and roughly how large are they relative to the tree-level
   pieces?

## What we'll do with the answers
Wire a representation-aware custodial mode into the scan that: zeros/residualizes
the protected Zb_L piece per Q3; uses the custodial ΔT coefficient per Q4; keeps
S; handles b_R / A_b per Q5; and flags (or includes) the one-loop residuals per
Q6 — with the chosen representations (Q1) and scope (Q2) as explicit metadata.
