# Custodial — reproducible plots & equation attribution (for the constraints we implement)

**Purpose.** Answer two questions before building any ours-vs-literature figure:
1. **Which local papers have custodial vs non-custodial PLOTS** we could reproduce,
   for the constraints actually in our implementation (NOT the not-yet-wired ones)?
2. **Who derived / computed each custodial equation we implement**, how it affects the
   process, and where we get it from?

Restricted to the ~10 constraints we run: **ε_K, ΔM_K, ΔM_Bd, ΔM_Bs, D⁰ mixing, S/T/U,
Z→bb**, plus the M12-derived deviation channels (S_ψφ, C_Bd, C_Bs, D⁰ phase). The
not-yet-implemented pieces (h→γγ, KK-Z′, exact custodian spectrum, leptons, b→sγ loops)
are deliberately excluded. Built from [`CUSTODIAL_PROVENANCE.md`](CUSTODIAL_PROVENANCE.md)
(Phase 0) + [`CUSTODIAL_LITERATURE_CATALOG.md`](CUSTODIAL_LITERATURE_CATALOG.md) (Phase 1).

---

## The one-line answer (which of our constraints custodial even touches)

**Custodial changes ONLY two of our constraints: S/T/U (T protected, S not) and Z→bb
(g_L^b protected). It is BLIND to every flavor constraint we run — ε_K, ΔM_K, ΔM_Bd,
ΔM_Bs, D⁰ mixing, and the M12 deviation channels — because those are KK-gluon
(color-octet, EW-singlet) operators.** So a "custodial vs minimal" plot is only
meaningful for S/T/U and Z→bb; for the flavor constraints the meaningful plot is the
*null result* — the wall stays put (Blanke, in the custodial model, vs Bauer minimal).

> **Real vs proxy (read `CUSTODIAL_PROVENANCE.md` "CRITICAL" section).** The ε_K
> wall being custodial-blind is REAL physics (KK-gluon is EW-singlet; Blanke's *full*
> custodial ΔF=2 confirms M_KK≳18–30 TeV). But our code showing custodial touching
> *literally nothing* in flavor is ~10% a proxy limitation: we omit the subleading
> custodial flavor effects (EW KK-boson ΔF=2, custodian loops). The plots below use
> only the **leading** custodial physics, which we implement correctly.

---

## TABLE A — Custodial-relevant plots in our local papers, per constraint

Legend: ✅ custodial changes it · ❌ custodial-blind · "ours" = the figure we already have.

| Our constraint | Custodial? | CUSTODIAL plot (paper, fig, page) | MINIMAL/non-cust. plot (paper, fig, page) | Our existing fig | Reproduce target? |
|---|---|---|---|---|---|
| **S, T, U** (EW001) | ✅ (T yes, S no) | CPSW-2006 Figs. 2,3,5,6 (T vs c_t, custodial top-partner sign), Fig. 7 (ΔS_f); CPSW-2007 Fig. 2 (k̃ bound vs c), Figs. 3,4 (S sensitivity) | CGHNP §6.2 (minimal T → 4 TeV, analytic); ADMS analytic S≈2πv²z_v² (no data plot) | `solo_ST_recentered.png` (minimal+custodial trajectories) | **YES** — overlay literature k̃/S anchors on our S-T |
| **Z→bb** (T010/T011) | ✅ (g_L^b P_LR) | **CPSW-2006 Fig. 1** (δg_bL/g vs c_q: "gauge/doublets/b_R′/**bidoublet**" curves — the protected-vs-not plot); **CPSW-2007 Fig. 1** (δg_bL/g vs ΔT, loop) | **CGHNP Fig. 8** (g_L^b–g_R^b plane, 3000 minimal-RS pts + 68/95/99% ellipse) | `solo_Zbb_gLgR.png` (our g_L–g_R plane) | **YES** — our panel ≈ CGHNP Fig.8; add CPSW protected/loop curves |
| **ε_K** (K001) | ❌ blind | **Blanke Fig. 9** (ε_K fine-tuning vs M_KK, **in the custodial model** → wall stays 18/30 TeV); Blanke Fig. 4 | **Bauer II Fig. 4** (anarchic ε_K cloud vs M_KK); **CFW Fig. 1 & Fig. 7-right** (ε_K KK-gluon scan); Archer Fig. 3.5/3.6 | `solo_epsK_cloud.png` (our Lane-A cloud) | **YES (the key one)** — Blanke(cust) vs Bauer/CFW(min): both ~tens TeV |
| **ΔM_K** | ❌ blind | Blanke Fig. 5 (tuning, custodial model) | Bauer II (text; weak, hadronic) | (in existence/typical fig) | minor |
| **ΔM_Bd, ΔM_Bs** | ❌ mostly (partial EW KK) | Blanke (Z_H/Z′ contributions compete w/ gluon) | Bauer II C_Bd=0.89±0.17, C_Bs=0.93±0.19 (SM-like) | (our M12 ReIm panel) | optional |
| **D⁰ mixing** (C001/C002) | ❌ blind (up-sector KK-gluon) | — (no custodial-specific D plot) | **Gedalia-Grossman-Nir-Perez** funnel (the (sin2σ, x12) funnel) | `solo_D0_funnel.png` | already ours; no custodial change |
| **S_ψφ, C_Bd, C_Bs, D⁰ phase** (M12 channels) | ❌ blind | Blanke/Bauer II: S_ψφ∈[−0.5,0.5] (NP allowed) | Bauer II Figs (M12 planes) | `solo_sm_overlap_scoreboard.png` | custodial irrelevant (flavor) |

**Net:** the genuinely reproducible "custodial vs non-custodial" comparisons in our
local archive, for our constraints, are **(1) ε_K** [Blanke custodial Fig.9 vs Bauer
Fig.4 / CFW Fig.1 minimal — both tens of TeV], **(2) Z→bb** [CPSW Fig.1 protected vs
CGHNP Fig.8 minimal], and **(3) S/T/U** [CPSW T-vs-c / k̃-bound vs CGHNP minimal-T].
Everything else is custodial-blind by construction.

---

## TABLE B — Who derived / computed each custodial equation we implement ("who computed what")

The chain is consistently: **original mechanism (who first derived) → published closed
form (who distilled the equation we copy) → our code**. Full eq#/page in
[`CUSTODIAL_PROVENANCE.md`](CUSTODIAL_PROVENANCE.md).

| Process (our constraint) | What custodial does | Who DERIVED the mechanism | Who computed the FORM we use | Our code |
|---|---|---|---|---|
| **Oblique T** (EW001) | SU(2)_R kills the log(L≈30) T enhancement; ΔT ∝ +L → −1/L (sign flip) | **Agashe-Delgado-May-Sundrum 2003** (hep-ph/0308036) — custodial-for-T origin | **CGHNP 2008** (0807.4937) Eq.147 (minimal) / **Eq.153** (custodial, distilling ADMS) | `oblique_stu.py` minimal/custodial_rs_plr_t_coefficient |
| **Oblique S** (EW001) | NOT protected — stays + → sets the residual floor | (no protection; ADMS notes S tamed by UV fermion loc.) | **PDG-2025** value c_S=30 (CGHNP S-part ~36.6 consistent) | `oblique_stu.py` ΔS = c_S v²/M² |
| **Z→bb g_L^b** (T010/T011) | P_LR (L↔R) zeroes δg_L^b at tree, b_L in (2,2)_{2/3} | **Agashe-Contino-Da Rold-Pomarol 2006** (hep-ph/0605341) — P_LR-for-Zbb origin | ACDP Eqs.3–12 (zeroing); minimal baseline = **CGHNP Eq.170** | `rs_ew_couplings.py` _apply_custodial_rs_plr_proxy |
| **Z→bb g_R^b** (T010/T011) | b_R∈(1,1)_{−1/3} → δg_R^b=0 (P_C) | **ACDP 2006** Table 1 (rep→sign map) | ACDP Table 1 row 2 (our `elementary_zero` = this rep) | `rs_ew_couplings.py` bR_strategy='elementary_zero' |
| **Z→bb loop** (T010/T011) | top-partners re-introduce −δg_L^b & ΔT at 1-loop | **Carena-Pontón-Santiago-Wagner 2006** (hep-ph/0607106) | CPSW **2007** (hep-ph/0701055) Eqs.13–15 + T_top Eq.35 | `rs_ew_couplings.py` _build_custodial_top_partner_loop_proxy |
| **ε_K, ΔM_K** (K001) | NOTHING (KK-gluon LR C₄, EW-singlet) | C₄ mechanism: **Csáki-Falkowski-Weiler 2008**, **Bauer II 2009**; custodial-model: **Blanke 2008** | hadronic ME: GGMS (hep-ph/9604387) normalization in our `deltaf2.py` | `deltaf2.py` (no ew_model branch) |
| **ΔM_Bd, ΔM_Bs** (B-mixing) | NOTHING (leading); partial EW KK subleading | **Bauer II 2009**, **Blanke 2008** | same `deltaf2.py` hadronic path | `deltaf2.py` |
| **D⁰ mixing** (C001/C002) | NOTHING (up-sector KK-gluon) | **Gedalia-Grossman-Nir-Perez 2009** (the funnel) | our D funnel post-processing | `deltaf2.py` D + funnel render |

**Plain-language "who computed what":**
- **Agashe-Delgado-May-Sundrum (2003)** first showed custodial **SU(2)_R protects T**.
- **Agashe-Contino-Da Rold-Pomarol (2006)** first showed the discrete **P_LR protects Z b_L b̄_L**, and tabulated which **b_R representation** gives which δg_R^b (Table 1).
- **Carena-Pontón-Santiago-Wagner (2006/2007)** computed the **top-partner one-loop** corrections (the −T bidoublet / +T singlet, and the loop δg_L^b) — the equations our loop proxy copies.
- **CGHNP (2008)** computed the **minimal-RS** S/T/Zbb and **distilled the custodial T into the compact Eq.153** that our code actually uses.
- **Blanke-Buras-Duling-Gori-Weiler (2008)** computed the **full ΔF=2 in the custodial model** and showed **ε_K stays ≳18–30 TeV** — the proof custodial doesn't touch flavor.
- **Csáki-Falkowski-Weiler (2008)** and **Bauer II (2009)** computed the **anarchic ε_K wall** (~21–33 / ~10 TeV) we reproduce as Lane A.
- **We** (this repo) wired these published forms into `oblique_stu.py` / `rs_ew_couplings.py` / `deltaf2.py`, with **four proxy choices that are ours, not theirs** (κ_b/L residual; ξ/ρ loop-mass inputs; the b_R rep pinning; singlet-only ΔT) — flagged in CUSTODIAL_PROVENANCE.md §3.

---

## So, what to reproduce (the shortlist)

Three custodial-vs-minimal comparisons are well-supported by local plots **and** map to
our constraints:
1. **ε_K** — Blanke Fig.9 (custodial) vs Bauer Fig.4 / CFW Fig.1 (minimal): the headline
   "custodial leaves the wall at tens of TeV." Overlay our Lane-A cloud.
2. **Z→bb** — CGHNP Fig.8 (minimal g_L–g_R, = our `solo_Zbb_gLgR.png`) with CPSW Fig.1
   protected/loop curves layered on.
3. **S/T/U** — our `solo_ST_recentered.png` (minimal+custodial trajectories) with the
   ADMS/CPSW S-limited ~2.5–3 TeV anchor and CGHNP minimal-T ~4 TeV anchor marked.
The other constraints (ΔM_K, ΔM_B, D⁰, deviation channels) are custodial-blind — their
"comparison" is a deliberate null (the curve doesn't move), best shown as one row on the
decoupling ladder rather than separate panels.
