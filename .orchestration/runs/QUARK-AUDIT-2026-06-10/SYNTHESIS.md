# Quark-sector audit — synthesis (2026-06-10)

Six read-only audit slices completed (see STATUS.md table; full evidence in
slice-*.md). No repo code was modified. Findings below are agent-derived;
per the program's dual-signoff gate, cross-review (e.g. Codex) is recommended
before any fix lands.

## Bottom line

The headline numbers — minimal RS floor 25–30 TeV, custodial strict 2–3 TeV,
custodial inclusive 7 TeV — are **not defensible as computed**. Three
independent implementation errors were found, pushing in BOTH directions:
the minimal floor is dominated by a mistranslated Zbb̄ term (over-excludes),
while ε_K is ~×1.7 too weak (under-excludes) and additionally
convention-dependent. The architecture (model construction, rotation
plumbing, scan harness, QCD running, anchors) is sound; the damage is
localized to specific formula transcriptions and one missing phase
convention, plus a test suite that pins code-vs-itself and so could never
catch them.

## Blockers

B1. **Zbb̄ fermion-KK admixture mistranslated from CGHNP 0807.4937**
    (`quarkConstraints/rs_ew_couplings.py:1835-1847, 957-971`; slice 3).
    Wrong per-generation 1/F² factor (uses light-generation c's instead of
    c_bR), wrong c-sign in the bracket, missing F²=2f² factor. At a
    scan-representative point: δg_L^b = −2.5e-3 (repo) vs +1.3e-5 (correct)
    — wrong sign, ~190×; δg_R^b ~2000× off. Slice 5 confirmed the minimal
    25–30 TeV floor is set 100% by T010 (694,123 vetoes; only rigorous veto
    above ~3 TeV) ⇒ **the minimal headline floor is an artifact of this
    term**. The gauge-KK Zbb piece (independently re-derived via Green's
    function) is exact; corrected minimal Zbb is gauge-dominated and ~25×
    weaker.

B2. **ε_K is convention-dependent: SVD column phases never fixed to a CKM
    convention** (`quarkConstraints/fit.py:245-258`; slice 2). Physically
    equivalent rephasing changed ε_K ratio-to-bound 17.34 → 105.1 (×6) with
    masses, |CKM|, Jarlskog bit-identical. Im(M12) is not rephasing-invariant
    as built; point-level pass/fail and any MFV phase-protection structure
    in r are undefined. Fix: rephase SVD factors to PDG convention (or use a
    rephasing-invariant ε_K vs arg[(V_ts*V_td)²]).

B3. **ΔF=2 matrix-element layer: O4/O5 chiral coefficients basis-swapped +
    state normalization 1/(2m_M) applied zero times**
    (`quarkConstraints/deltaf2.py:719-724, 908-913` + 3 vendored copies in
    `modern/phenomenology.py`; slice 1). Code uses (R/6+1/4)B4, (R/2+1/12)B5
    (BMU transcription) where the CFW-matched SUSY basis requires
    (R/4+1/24)B4, (R/12+1/8)B5; and all ⟨O⟩ are 2× the M12-ready values —
    contradicting the documented audit claim at STATE_OF_PROJECT.md:29.
    Net: ε_K underestimated ×1.72 ⇒ **ε_K M_KK floors ~31% too LOW
    (anti-conservative)**; B_d/B_s nearly cancel (×1.02–1.08), D0 ≈0.9.
    Slice 4 quantified at the gate level: kaon LR M12 understated 0.58×,
    VLL overstated 2.0×. The two bugs partially cancel in some channels —
    **fix together, not one at a time**.

## Majors

M1. **T010 gate**: vetoes at max|pull| ≤ 1.0σ while SM R_b pull is already
    −0.996 ⇒ 0.004σ NP headroom. With corrected couplings the floor is
    108 TeV (1σ, artifact) vs 6.8 TeV (2σ). Gate choice must be revisited
    (slice 3).
M2. **EW001 ΔT convention mix**: CGHNP geometric-Λ_IR coefficient used with
    physical M_KK ⇒ T ×6 underestimated (anti-conservative). Proxy lane, so
    strict floors unmoved; consistent treatment puts minimal EW001 at
    ~15–17 TeV (slice 3). Only scale-convention leak found anywhere
    (slice 5 swept all families).
M3. **"Custodial strict 2–3 TeV" violates the program's own caveat**:
    CORRECTED_PRESCRIPTION.md says do not claim it without top-partner loop
    numerics (δg_L^b|loop ~ 1e-3); loops deferred yet results tagged
    rigorous (slice 3). Also: custodial INCLUSIVE 7 TeV is a
    model-independent proxy wall — EW001/CR001/CR012/CR013/B012 veto counts
    are bit-identical between minimal and custodial runs (slice 5).
M4. **Floor prose vs grid**: grid {...,15,20,30,50} supports only
    "(20,30] TeV for r≥0.1; (15,20] for r=0.05"; "25–30" is the
    10–12 TeV ×2.45 literature conversion, not a scan measurement; the
    r=0.05 minimal floor ≤20 TeV is nowhere stated (slice 5).
M5. **tag_result substring bug**: T001/T002 evaluated+passing but tagged
    partial ("proxies" misses "proxy") → counted hard_not_evaluated on 100%
    of 878k rows, coverage_complete=False everywhere; a FAILING partial-HARD
    would silently never veto (slice 5).
M6. **Bag/Wilson scale pairing**: kaon B4/B5 quoted at 3 GeV used at the
    2 GeV Wilson point (~19% on O4); B-meson bags at m_b (~7%); mixed-scale
    quark masses in chiral factors (D0 ~17% anti-conservative) (slices 1/4/6).
M7. **Test suite structurally blind**: constraint tests pin adapter ≡ core
    and regression snapshots pin the buggy values (77/77 pass); no
    literature-anchored absolute pins anywhere in ΔF=2/Zbb (slices 1, 4).
M8. **qcd/ package**: d3 3-loop mass-decoupling constants wrong (≤2.2e-4
    effect, latent); `pdg_quark_masses_at_scale(mu)` crashes for mu < 4.18
    GeV with the PDG top input (violates documented contract) (slice 6).

## Verified sound (explicitly checked, non-exhaustive)

- MFV spurion construction, eigen-ordering, c-map, bulk rotations, SVD
  conventions (no Vh/V conflation), CKM = U_Lu†U_Ld, v=174 convention,
  f_IR formula and c→½ limit (slice 2).
- ΔF=2 tree matching C1/C4/C5 vs CFW incl. signs; LO RG: BMU→scalar ADM map
  re-derived, γ_VLL=4, LR factors closed-form exact; running applied once;
  κ_ε=0.94 and ε_K master formula (slices 1, 2).
- Gauge-KK Zbb shift re-derived independently (matches tower sum to 1e-4);
  T014 width factors; custodial PR1 mechanics per prescription;
  pr1_minimal_offdiag conservative and apples-to-apples (slice 3).
- All experimental anchors current; B_K(2 GeV)=0.5503 correctly paired (no
  RGI mixup); ps⁻¹↔GeV exact; Δm→M12 halving once everywhere; B002/B004
  phase logic (slice 4).
- Harness: no other scale leaks; drops are perturbativity skips at low r
  with correct denominators; minimal/custodial pairing verified at draw
  level; floors threshold-insensitive; reports match artifacts (slice 5).
- alpha_s + γ_m running agree with CRunDec to 1e-9 from 1 GeV to 50 TeV;
  g_s* uses n_f=6 MS-bar correctly (slice 6).

## Recommended fix order (after dual-signoff review of each)

1. B2 (CKM rephasing in fit.py) — prerequisite for any ε_K statement.
2. B3 (O4/O5 coefficients + 1/(2m_M), all 4 copies, together) + M6 scales.
3. B1 (Zbb admixture retranslation; corrected term is tiny — consider
   dropping with a documented bound).
4. M1/M2 (T010 gate policy; EW001 coefficient convention).
5. M5 (tag_result substring) + M7 (literature-anchored absolute tests that
   would have caught B1/B3).
6. Re-run the paired 1M scans, rebuild comparison + Scan Explorer, restate
   floors with grid-honest brackets (M4) and the custodial loop caveat (M3).

## Expected post-fix picture (rough, to be confirmed by re-scan)

Minimal: T010 collapses as a driver (B1); corrected ε_K strengthens ~×1.3
and likely becomes the (convention-stated) dominant floor; EW001 (proxy)
~15–17 TeV if kept consistently. Custodial: strict floor driven by corrected
ΔF=2, with the 2–3 TeV claim suspended pending top-partner loops. The
minimal-vs-custodial qualitative contrast survives; every headline number
changes.

## Agent ids (SendMessage to continue any slice)

slice 1 ae4caca93ff50d30e · slice 2 aa4a9ba5107a5fc73 · slice 3
a171609616c2ebce9 · slice 4 a2c74aeeae398ff88 · slice 5 a3fe349291685f7f6 ·
slice 6 af727ea41669273c3
