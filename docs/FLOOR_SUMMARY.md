# Minimal-RS Floor Summary (canonical current state)

**Status:** authoritative. Post-audit (June 2026). This is the one-page,
up-to-date summary of the minimal (non-custodial) Randall–Sundrum quark-sector
KK-scale floor. Where any older doc disagrees, this file and
[`reports/collaborator_2026-06/CONTENT.md`](../reports/collaborator_2026-06/CONTENT.md)
win.

All scales below are the **physical first gauge-KK mass** `M_KK = x₁·Λ_IR`
with `x₁ ≈ 2.4487`. Do not confuse with the geometric `Λ_IR ≡ 1/z_v`.

## The floor

| Constraint | Floor (physical M_KK) | Type | Tunable? | Driver |
|---|---|---|---|---|
| `epsilon_K` (K001) | **~30 TeV** | typical (median anarchic) | **yes** — align Im M12 → 0 | flavor / CP |
| oblique S,T,U (EW001) | **~18–20 TeV** | **existence** (best-tuned) | **no** — no Yukawa freedom | RS T-problem |
| Z→bb (T010/T011) | **~5 TeV** | — | no (gauge-dominated; c_Q3 pinned by m_t) | EW |
| Δm_s (B003), Δm_d, Δm_K, D⁰ | few TeV → <~1 TeV | tunable | yes | flavor |
| collider KK searches (CR*) | ~4 TeV | subleading | — | direct search |

Two distinct notions of "floor":

- **Typical (median anarchic) floor ≈ 30 TeV**, set by `epsilon_K`. This is
  where a *median* anarchic point survives. It is **tunable**: aligning the
  down-sector so Im M12 → 0 drives the NP contribution to zero, so the
  *existence* floor for any single flavor constraint is ≲1 TeV.
- **Existence (fine-tuned, irreducible) floor ≈ 18–20 TeV**, set by oblique
  **S,T,U** (the RS T-parameter problem). This carries **no Yukawa freedom**, so
  min ≈ median (spread only ~1.5×): even maximally fine-tuned Yukawas cannot put
  minimal RS below ~18–20 TeV. This is precisely the motivation for custodial RS.

**Z→bb is NOT the driver.** The old claim of a "25–30 TeV Z→bb-dominated floor"
(and "108 TeV at 1σ") was a **B1** bug (see below); corrected, Z→bb collapses to
~5 TeV and is gauge-dominated.

## Custodial caveat

Custodial RS (`SU(2)_L × SU(2)_R × U(1)_X × P_LR`) **fixes the oblique T
problem** and protects Zb_L, but does **not** relax `epsilon_K` (a flavor/CP
constraint). The custodial scan's strict floor is ~2–3 TeV from rigorous ΔF=2
and ~7 TeV inclusive once proxy EW/collider failures are counted — but
`epsilon_K` tunability, not Z→bb, is what governs the residual flavor floor.

## Audit fixes (June 2026, all dual-reviewed)

Six implementation errors found by an adversarial audit, now fixed:

- **B1** (blocker): Z→bb fermion-KK admixture mistranslated from Casagrande
  0807.4937 — wrong sign and ~190× too large. This faked the 25–30 TeV Z→bb
  floor. After the fix (+ residual-B1, the `1−2c → 1+2c` flavour-sum
  denominator sign), Z→bb is gauge-dominated at ~5 TeV.
- **B2**: `epsilon_K` SVD→PDG rephasing. Independently re-verified **CORRECT** —
  pinning V_ts as the 5th rephasing anchor is a harmless gauge choice; `epsilon_K`
  and all Δm observables are exactly rephasing-invariant (bit-identical numeric
  check). See the gauge note in `quarkConstraints/fit.py:_rephase_to_pdg_convention`.
- **B3**: ΔF=2 O4/O5 chiral coefficients + meson-state normalization.
- **M1**: Z→bb (T010) gate.
- **M2**: EW001 ΔT convention — oblique T was ~6× too weak. Correct minimal form
  is `ΔT_minimal = x₁² · (πL)/(2 cos²θ_W) · v²/M_KK²` (physical-M_KK convention).
- **M5**: a tagging bug (plural "proxies" missed by substring match).

## Validated literature reproductions

Run the fixed code in **anarchic-forward** mode (draw complex O(1) Yukawas,
forward-compute, keep mass+CKM-reproducing points):

- **Bauer 0912.1625**: `epsilon_K` 95%-quantile floor ≈ 10 TeV (paper-era
  inputs); ~106× Im/Re(M12) "ε_K problem" (≈ Bauer's ~100×). Script:
  `scripts/anarchic_bauer_s1.py` (+ `scripts/run_bauer_s1_deband.sbatch`,
  `scripts/run_bauer_scenarios.sbatch`).
- **Gedalia, Grossman, Nir, Perez 0906.1879**: D⁰ funnel containment rising to
  ~88% at 10 TeV. Script: `scripts/anarchic_complex_m12.py`.
- **Blanke 0809.1073**: Re/Im(M12) for K and B_s. Scripts:
  `scripts/anarchic_complex_m12.py`, `scripts/anarchic_reproduction_extract.py`.

The production scan **fits** Yukawas to masses+CKM (a thin locus), whereas the
literature **scatters** anarchic Yukawas (a fat ~6-decade cloud); the two
methods give very different spreads, which matters when comparing floors.

## Pointers

- Collaborator report (figures + per-figure physics):
  [`reports/collaborator_2026-06/CONTENT.md`](../reports/collaborator_2026-06/CONTENT.md)
- Full project state: [`STATE_OF_PROJECT.md`](STATE_OF_PROJECT.md)
- Known open issues: [`KNOWN_ISSUES.md`](KNOWN_ISSUES.md)
- Oblique convention: `quarkConstraints/oblique_stu.py`,
  `flavor_catalog_constraints/primary/top_higgs_ew/EW001.py`
- Legacy ΔF=2-only writeups (SUPERSEDED banners):
  `docs/quark_scan_methodology_note.tex`, `docs/quark_scan_consolidation_report.tex`,
  `docs/quark_scan_constraint_update_2026-06.md`
