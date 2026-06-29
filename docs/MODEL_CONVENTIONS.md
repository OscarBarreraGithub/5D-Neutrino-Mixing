# Model Conventions — quark-sector flavor structure (canonical)

**Status:** authoritative reference for *which flavor model we are running* and
*how it is parametrized*. Written to stop the recurring confusion between the
production (aligned, MFV) model and the anarchic literature-reproduction lane.
Date: 2026-06-29.

---

## 0. TL;DR — there are TWO lanes, do not conflate them

| Lane | Flavor structure | Purpose | Spread | epsilon_K floor |
|---|---|---|---|---|
| **Production / our model** | **RS-MFV (Fitzpatrick-Perez-Randall 0710.1869), aligned** | the floors, catalog, website, the model we advocate | thin locus | LOW (FPR ~2 TeV / production run ~7 TeV) — *to be re-verified, see sec 5* |
| **Anarchic reproduction** | **anarchy** (Bauer/Casagrande) | reproduce the published Bauer/Casagrande/Blanke plots, side-by-side | fat ~6-decade cloud | HIGH (~10 TeV paper-era, ~30 TeV current) |

**The anarchic ~10-30 TeV epsilon_K wall is the BASELINE/strawman the literature
criticizes — it is NOT our model.** Our model is the *aligned* FPR-MFV one, whose
whole point is to suppress that wall. When we quote "our floor," it must be the
FPR-aligned number, not the anarchic one.

The Bauer scenarios **S1/S2/S3/S4 live ONLY in the anarchic reproduction lane**
(`scripts/anarchic_bauer_s1.py`). They are Bauer's benchmarks, not our model.
"Are we S1 or S2?" -> for the production model, NEITHER; we are FPR-MFV.

---

## 1. The C_Q spurion combination — two equivalent parametrizations

The left-handed bulk-mass matrix is built from the up/down Yukawa **spurions**
`C_u = Y_u^dagger Y_u`, `C_d = Y_d^dagger Y_d`. Two ways to write the combination,
both in the repo:

**(a) Single-knob `r`** (simplified production lane, `quarkConstraints/model.py:237`):
```
C_Q = r * (Y_u Y_u^dagger) + (Y_d Y_d^dagger)
```
One coefficient `r` (up-weight); the down term is normalized to 1.

**(b) Two-coefficient `(a, r)` == `(gamma, alpha)`** (full FPR lane,
`quarkConstraints/paper_0710_1869/model.py:104`, Eq. (3) of 0710.1869):
```
diag(C_Q) = a * diag( r * V5KM^dagger C_u V5KM + C_d )
          = (a*r) * (V5KM^dagger C_u V5KM)  +  a * (C_d)
```
- up-term coefficient = `a*r`  (this is **gamma**)
- down-term coefficient = `a`  (this is **alpha**)
- `a` = overall scale of the LH bulk masses (default 0.8), `r` = up/down ratio
  (default 0.3). Stored separately on the FPR point.

So `(gamma, alpha) = (a*r, a)`, and the single-knob form is the special case with
the overall scale `a` absorbed elsewhere. **Keep `r` — it is the clean ratio — but
the canonical model is the explicit two-coefficient `(a, r)`/`(gamma, alpha)` form,
because it separates the overall scale from the up/down weighting.**

---

## 2. V5KM — the actual ALIGNMENT knob

`V5KM` in Eq. (3) is a 5D CKM-like unitary rotation `(theta12, theta23, theta13,
delta)` that rotates `C_u` before it is added to `C_d`. **This is what controls
alignment / flavor violation:**
- `V5KM = I` (or = the SM CKM in the appropriate basis): the up and down sectors are
  **aligned** -> `C_Q` is (near-)simultaneously diagonalizable with the down sector
  -> NP flavor violation is suppressed -> **this is the cure** (low KK floor).
- `V5KM` misaligned: off-diagonal structure in `C_Q` -> NP flavor violation turns on
  -> the epsilon_K problem returns.

This is the FPR "GIM-from-extra-dimensions" / next-to-minimal-flavor-violation
mechanism: it is a *principled* alignment, parametrized by how far `V5KM` sits from
the aligned point.

---

## 3. Bulk-mass map — USE THE EIGENVALUES DIRECTLY (decision 2026-06-29)

**Decision:** the bulk masses are an **affine function of the spurion eigenvalues**,
per sector (FPR mapping):
```
c_i^(x) = alpha_x * lambda_i(C_x) + beta_x        (x = Q, u, d)
```
where `lambda_i(C_x)` are the (ordered) eigenvalues of the spurion combination and
`(alpha_x, beta_x)` are fit to reproduce the measured quark masses + CKM.

**We do NOT use the `BulkMassMap` bounded-squash surrogate**
(`quarkConstraints/model.py:48`, `c = c_uv - (c_uv-c_ir)*lambda/(lambda+scale)`).
That nonlinear bounded map was an exploratory device; it obscures the direct
eigenvalue -> bulk-mass relation that the FPR mechanism is built on. The canonical
path takes the eigenvalues straight into the affine map above.

> **CODE TODO (not done in this doc):** the active production scan
> (`quarkConstraints/model.py` -> `fit.py`) still routes through `BulkMassMap`.
> Switching the production path to the direct affine-eigenvalue map `c_i = alpha_x
> lambda_i + beta_x` (the FPR `paper_0710_1869` lane already does this) is a separate
> implementation task; until then, production floors carry the BulkMassMap surrogate
> caveat.

---

## 4. Bauer's S2 vs our alignment (settled from 0912.1625 Sec. 5)

**What Bauer's S2 actually is.** Per 0912.1625 Sec. 5.2: S2 "aligned" keeps the same
Y_max=3 and warp L≈37 as S1, and only restricts to "common bulk masses c_{d_i} in the
sector of right-handed down-type quarks," achieved "by imposing a U(3) flavor symmetry"
on the RH down fields (Table 1: c_d1=c_d2=c_d3=0.60). It is genuinely an **alignment
scenario** — a U(3) DEGENERACY of the right-handed down bulk masses.

**It suppresses epsilon_K (directly).** Bauer Eq. (123): common c_d suppresses the
chirally/RG-enhanced LR operator C_4 by `(C_4)_aligned/(C_4)_hier ~ (Y_d^2 v^2/2M_KK^2)
(m_d/m_s) ~ 8e-3`, an O(100) cut. The epsilon_K-consistent fraction rises 19% -> 38%
and the 10%-consistency floor drops 3.6 -> 1.9 TeV. **Bauer's S2 cloud stays FAT (same
anarchic Y_max=3 draws) but gains far more epsilon_K-consistent (orange) points** — the
"improvement" is more consistency, not a thinner cloud.

**Is our anarchic "S2 = common c_d" faithful?** The DEFINITION is faithful and NOT
mislabeled: `scripts/anarchic_bauer_s1.py` SCENARIOS["S2"] keeps Y_max=3/c_max=2/L=geom
identical to S1 and only toggles `common_cd=True` (sets c_d to one common value) — exactly
Bauer's S2. The underlying S2 ensemble (`anarchic_bauer_S2.parquet`, 195k draws) is also
genuinely fat (~5.6-decade per-tile spread). **The problem is only with the rendered
PANEL** (`solo_epsK_cloud_S2.png`): its Z->bb overlay comes from the OLD sparse
`anarchic_bauer_s2_zbb.parquet` (~46 blue points; never densified like S1), so it *looks*
collapsed/empty next to Bauer's dense fat S2. This is a sampling/rendering artifact, not a
scenario or physics error. Fix = re-render S2 with the same hybrid treatment used for S1
(fat 195k grey + a densified S2 Z->bb overlay). [Correction to an earlier claim: S2 is
NOT a wrong/mislabeled scenario.]

**S2 alignment vs our FPR alignment — different sectors, and the key insight.**
- **Bauer S2** degenerates the **right-handed down** bulk masses (U(3) on RH singlets).
  This directly diagonalizes RH s-d mixing and chops the **LR / C_4** kaon operator — the
  single dominant, chirally-enhanced driver of the epsilon_K problem.
- **Our FPR (V5KM)** alignment is a **left-handed doublet** statement: V5KM ~ SM CKM
  rotates the up/down Yukawas into near-coincidence, suppressing LH s-d mixing and the
  **LL / C_1** operator. It only touches C_4 indirectly.

> **KEY INSIGHT:** epsilon_K is dominated by the LR (C_4) operator, which lives in the
> product (LH s-d mixing) x (RH s-d mixing). Bauer's S2 (RH-down U(3) degeneracy) is the
> **more direct, more targeted epsilon_K cure**; our FPR V5KM (LH) alignment suppresses
> the wrong (LL) sector for epsilon_K. If the goal is to kill the binding flavor
> constraint (epsilon_K), an RH-down-degeneracy alignment is arguably what we want MORE
> than (or in addition to) the FPR LH rotation. This is a concrete design lesson for the
> production model, not just a reproduction detail.

---

## 5. Floors by lane (to be stated separately, never blended)

- **Anarchic reproduction (Bauer S1):** epsilon_K typical (median) floor ~30 TeV
  current / 95%-quantile ~10 TeV paper-era; S,T,U existence ~18-20 TeV (irreducible,
  flavor-independent); Z->bb ~5 TeV. These are the numbers in the collaborator
  report / notes.pdf reproductions. They describe ANARCHY.
- **FPR-aligned production:** the aligned model suppresses the flavor wall. The old
  FPR lane quotes ~2 TeV; the QUARK-FIX production 100k recorded ~7 TeV epsilon_K and
  ~20 TeV S,T,U. **This number must be re-verified directly from the production scan
  output and stated next to the anarchic one** (open item).

**Whenever a floor is quoted anywhere (FLOOR_SUMMARY, the reassessment plan, the
report), it MUST be tagged with its lane.** The RS-vs-SM reassessment note is about
*anarchic* RS; the production model is *aligned* FPR-MFV — different models, different
floors.

---

## References
- Fitzpatrick, Perez, Randall, 0710.1869 — "A GIM mechanism from extra dimensions"
  (the production model's flavor structure; `quarkConstraints/paper_0710_1869/`).
- Bauer, Casagrande, Haisch, Neubert, 0912.1625 — anarchic scan (the reproduction
  lane; `scripts/anarchic_bauer_s1.py`).
- `docs/FLOOR_SUMMARY.md`, `docs/RS_vs_SM_flavor_reassessment_plan.md`.
