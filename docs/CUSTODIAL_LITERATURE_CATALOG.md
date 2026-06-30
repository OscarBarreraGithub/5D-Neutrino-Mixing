# Custodial RS — literature catalog (what custodial touches, and what others did)

**Purpose.** Phase 1 of the custodial program: a clean, sourced catalog of *every
place custodial symmetry appears in the RS literature*, organised by observable —
what custodial does to it, by what mechanism, the quantitative KK-scale floor each
paper reports, and which figures are the best "ours vs literature" comparison
targets. Built from a paper-by-paper sweep of the local `references/papers/` archive
(2026-06-30). Companion to [`docs/CUSTODIAL_PROVENANCE.md`](CUSTODIAL_PROVENANCE.md)
(Phase 0: where OUR equations come from).

> **Units caveat (read first).** Floors are quoted on different objects. `M_KK` = the
> gauge KK scale (≈ Λ_IR). `M_g(1) ≈ 2.45 M_KK` = the physical first KK gluon/gauge
> mass. A "~10 TeV M_g(1) wall" = "~4 TeV M_KK." Each entry flags which it uses.

---

## 0. The one-sentence headline (unanimous across 7 papers)

**Custodial symmetry (SU(2)_R × P_LR) lowers the ELECTROWEAK floor to ~2–3 TeV (set
by the unprotected S parameter) by protecting T and the Z b_L b̄_L vertex — but it
does NOTHING for ε_K, which stays at ~18–33 TeV in the anarchic lane because the
binding kaon operator is the color-octet KK-gluon LR (C₄/Q_LR), a *pure-QCD* object
blind to the electroweak custodial group.** Only *flavor* alignment (RH-down U(3)
degeneracy / 5D-MFV / V_5KM) moves the ε_K wall.

This is exactly our note's thesis, and the literature underwrites it quantitatively.

---

## 1. Master table — by observable

| Observable | Does custodial help? | Mechanism | Literature numbers (paper, floor, page) |
|---|---|---|---|
| **T (tree)** | ✅ YES | SU(2)_R kills the log(kπr_c≈30) enhancement | ADMS: protected, floor moves to ~3 TeV [hep-ph/0308036 §4.1]; CGHNP: cure #4, T∝1/L "tiny" [0807.4937 Eq.153]; CFW: "no tree-level T constraint" [0804.1954 §4.1] |
| **T (loop, top-partners)** | ◐ sign-dependent | bidoublet → **−T**; singlet → **+T** (compensation) | CPSW: ΔT≈−0.08 (c=+0.42) to −0.25 (c=−0.42) bidoublet [0607106 p.23]; singlet +T [Eq.43]; +T favored by +S |
| **S** | ❌ **NO** (unprotected) | no symmetry protects S | ADMS: S≈2πv²z_v²≈0.2 → ~3 TeV floor [Eq.7.1]; CPSW: S sets k̃≳1 TeV → M_g(1)≳2.5 TeV [Eq.17]; CFW: S<0.2 → M_KK≳3 TeV [Eq.4.25]; Archer: custodial S-bound M_KK≳1.5–2 TeV [p.56] |
| **U** | — | ≈0 both | negligible everywhere |
| **Z b_L b̄_L** (g_L^b) | ✅ YES (tree) | P_LR: b_L∈(2,2)_{2/3}, T³_L=T³_R | ACDP: δg_Lb→0 [hep-ph/0605341 Eq.12]; non-protected → M_KK≳3.5 TeV (M_g≳8.75 TeV) [CPSW p.11]; protected → δg/g≈0.8×10⁻³ |
| **Z b_L b̄_L (loop)** | ◐ partially undone | top-partners re-introduce −δg_Lb, correlated with +T | CPSW Fig.1: |δg_bL/g| grows, exits 2σ band by ΔT≈0.2–0.5 [0701055 p.11] |
| **Z b_R b̄_R** (g_R^b) | ◐ embedding-choice | b_R∈(1,1)_{−1/3} → δg_Rb=0 (P_C); other reps nonzero either sign | ACDP **Table 1** (rep→sign map, §2 below); minimal RS leaves g_R^b≈SM anyway [CGHNP p.50] |
| **A_FB^b** | ◐ only via tunable +δg_Rb | needs δg_Rb≈+0.02 (a Z-vertex fix, not flavor) | ACDP: b_R↪4 of SO(5) gives +δg_Rb [Eq.24]; minimal RS makes the −2.1σ anomaly WORSE [CGHNP Fig.8] |
| **Z t_L / W t_L b_L** | ❌ NO | P_LR can't protect top while protecting bottom | CPSW: δg_Ztt/g~−0.2, δg_Wtb/g~−0.07 [0607106 Eq.62] — a custodial *prediction*, not a fix |
| **ε_K (LR / C₄ / Q_LR)** | ❌ **NO (the crux)** | color-octet KK-gluon, EW-singlet → custodial-blind | **Blanke (custodial model!): M_KK≳18 TeV (Δ<20) / ≳30 TeV (Δ<10)** [0809.1073 Eq.6.3–6.4]; CFW: M_g≳21 TeV (RS)/33 TeV (GHU) [0804.1954]; Bauer II: M_g(1)≳10 TeV (S1,10%) [0912.1625 Eq.124]; Archer: M_KK≳10/22/30 TeV [1201.1561 p.37] |
| **ΔM_K** | ❌ no (QCD) | same KK-gluon sector | weaker than ε_K everywhere; large hadronic uncertainty |
| **ΔM_Bd, ΔM_Bs** | ◐ partial | EW KK bosons (Z_H, Z′) "compete with KK gluons" | Blanke abstract; Bauer II C_Bd=0.89±0.17, C_Bs=0.93±0.19 (SM-like) [0912.1625] |
| **S_ψφ, A_SL^s** | n/a | NP allowed, can exceed SM | Blanke/Bauer II: S_ψφ∈[−0.5,0.5] (unconstraining) |
| **b→sγ, b→sZ⁰** | ❌ NO | left unprotected | CGHNP: "leaving t→c(u)Z⁰ and b→sZ⁰ unprotected" [p.56]; Bauer II: b→sγ limits the large-Yukawa (S4) ε_K cure [p.55] |

---

## 2. ACDP Table 1 — b_R embedding → δg_R^b (resolves our provenance gap #3)

Source: Agashe-Contino-Da Rold-Pomarol hep-ph/0605341, Table 1, p.9. Defines the sign
and size of the RH-bottom Z-coupling shift per embedding. **Our production proxy's
`g_R^b = 0` ('elementary_zero') corresponds to the b_R∈(1,1)_{−1/3} row** — the
P_C-protected, literature-standard minimal custodial choice. Pinning this rep (and
citing this table) closes provenance gap #3.

| b_R rep | δg_R^b gauge (∝ −Q_A) | δg_R^b fermionic | note |
|---|---|---|---|
| (1,3)_{2/3} | −1 | −½sin²θ^{(2,2)} + sin²θ^{(2,4)} | negative |
| **(1,1)_{−1/3}** | **0** | **0** | **← our proxy (P_C protected)** |
| (1,3)_{−1/3} | 0 | 0 | also protected |
| (1,2)_{1/6} | −1/2 | −½sin²θ^{(2,1)} + ½sin²θ^{(2,3)} | negative |
| (1,2)_{−5/6} | +1/2 | +½sin²θ^{(2,1)} − ¼sin²θ^{(2,3)} | positive |
| (1,3)_{−4/3} | +1 | +½sin²θ^{(2,2)} − ⅓sin²θ^{(2,3)} | positive (A_FB^b cure) |

---

## 3. Per-paper floor summary (custodial EW vs ε_K)

| Paper | Custodial EW floor | ε_K floor | Key statement |
|---|---|---|---|
| **ADMS 2003** (hep-ph/0308036) | M_KK ~3 TeV (from S + Zb_L) | not addressed | origin of custodial-for-T; ε_K is a separate axis |
| **CPSW 2006/07** (0607106/0701055) | M_g(1) ~2.5–3 TeV (S-set) | absent | custodial lets KK be light; tree Zb_L protected, loop partially undoes it |
| **CGHNP I 2008** (0807.4937) | M_KK >2.4 TeV (~6 TeV M_1) | ~10–20 TeV M_1 | "four cures" = four cures for **T**, not ε_K |
| **CFW 2008** (0804.1954) | ~3 TeV (S) | **21 TeV (RS) / 33 TeV (GHU)** | most explicit EW-vs-flavor split |
| **Blanke 2008** (0809.1073) | **2–3 TeV** (EW precision) | **18 TeV (Δ<20) / 30 TeV (Δ<10)** | THE custodial-model ΔF=2 paper — ε_K wall survives custodial |
| **Bauer II 2009** (0912.1625) | Zb_L removable | M_g(1)≳10 TeV (S1) → **S2 align ~3 TeV** | C₄ is "a pure QCD effect... insensitive to the EW embedding" |
| **Archer 2012** (1201.1561) | M_KK ~1.5–2 TeV | M_KK ~10/22/30 TeV (tuning) | "custodial does not protect S or the gauge-fermion couplings" |

---

## 4. Ours vs literature — the mapping (Phase 2 targets)

Our floors (lane-tagged) line up with the literature as follows:

| Our quantity | Our value | Literature counterpart |
|---|---|---|
| S,T,U existence floor (non-custodial) | **18–20 TeV** | the *un-custodialized* oblique floor; literature custodial drops this to ~2–3 TeV (what a full custodial build would buy) |
| Z→bb floor (post-B1) | **~5 TeV** | P_LR-protected Zb_L (~few TeV; non-protected would be ~8.75 TeV [CPSW]) |
| ε_K Lane A (anarchic) | **~9–24 TeV** | anarchic ε_K wall: Blanke 18–30, CFW 21–33, Bauer II ~10 (M_g(1)), Archer 10–30 — **and custodial does NOT move it (Blanke)** |
| ε_K Lane C (FPR/aligned) | **~2 TeV** | custodial + alignment floor (Bauer S2 ~3 TeV; Blanke fine-tuned 2.5 TeV at Δ~700) |

**Best comparison figures to reproduce (Phase 2):**
- **ε_K wall:** Bauer II **Fig. 4** (ε_K vs M_KK quantile cloud — directly our Lane-A scatter); Blanke **Fig. 9** (ε_K fine-tuning vs M_KK — the custodial-model 18/30 TeV result); CFW **Fig. 1 / Fig. 7** (anarchic ε_K scans).
- **EW custodial vs minimal:** CPSW **Fig. 1** (δg_bL/g vs ΔT correlation); ADMS/CFW S→3 TeV curve; our recentered S-T + C1 floor table.
- **Zbb:** CGHNP **Fig. 8** (g_L^b–g_R^b plane with 68/95/99% ellipse — our `solo_Zbb_gLgR.png` analog).

---

## 5. Implications

**For the note.** §7 ("the cures and what they cost") is now fully literature-backed:
custodial → EW ~2–3 TeV (S-limited), ε_K untouched (Blanke, in the custodial model
itself). The §7 custodial paragraph can cite Blanke 0809.1073 Eq.(6.3–6.4) directly
for "custodial does not relieve ε_K," and CFW 0804.1954 for the explicit EW(3 TeV)-vs-
flavor(21–33 TeV) split. This is stronger than asserting it from the operator
structure alone.

**For the four provenance gaps** (CUSTODIAL_PROVENANCE.md §3):
- gap #3 (g_R rep) — **closable now**: declare b_R∈(1,1)_{−1/3}, cite ACDP Table 1 row 2.
- gap #4 (singlet-only ΔT) — the literature confirms the bidoublet ΔT is **negative and
  often dominant** (CPSW); our singlet-only truncation is anti-conservative and should
  carry the override, or be documented as an upper-T proxy.
- gaps #1 (κ_b/L residual) and #2 (ξ/ρ loop inputs) — still need a derivation; the
  literature (CPSW App. A) is where the wavefunction-correlated masses live, so the
  full custodial build can source them from there.

**For the full custodial build (roadmap D).** The literature gives the target: a
custodial implementation should reproduce (i) the S-limited ~2–3 TeV EW floor, (ii) the
−T bidoublet / +T singlet loop competition (CPSW Eqs.41–43), (iii) the unchanged ε_K
wall (Blanke), and (iv) the b_R-embedding-dependent A_FB^b. Reproducing Blanke Fig. 9
(custodial ε_K vs M_KK) would be the headline validation.

---

## References (all local in `references/papers/`)
- ADMS — Agashe, Delgado, May, Sundrum, hep-ph/0308036 (custodial-for-T origin).
- ACDP — Agashe, Contino, Da Rold, Pomarol, hep-ph/0605341 (P_LR-for-Zbb, Table 1).
- CPSW — Carena, Pontón, Santiago, Wagner, hep-ph/0607106 + hep-ph/0701055 (light KK, oblique).
- CGHNP — Casagrande, Goertz, Haisch, Neubert, Pfoh, 0807.4937 (minimal RS, four T-cures).
- CFW — Csáki, Falkowski, Weiler, 0804.1954 (EW-vs-flavor split, anarchic ε_K).
- Blanke — Blanke, Buras, Duling, Gori, Weiler, 0809.1073 (**custodial-model ΔF=2**).
- Bauer II — Bauer, Casagrande, Haisch, Neubert, 0912.1625 (anarchic ΔF=2; our Lane A).
- Archer — 1201.1561 (comprehensive warped constraints).
