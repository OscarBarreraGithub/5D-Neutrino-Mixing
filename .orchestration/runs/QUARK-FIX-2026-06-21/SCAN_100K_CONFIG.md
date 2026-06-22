# 100k-point Minimal Quark-Only Scan — Run Config (READ-ONLY recon)

Date: 2026-06-22. Repo root: `/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing`.
This is a recon report only — nothing was launched or modified.

---

## 0. TL;DR — exact command(s)

The scan harness evaluates **all** quark-only catalog constraints on every draw and
records per-constraint veto fractions. You therefore do **not** restrict the constraint
set at scan time; you run the full quark-only scan and **filter to your 10 constraints
post-hoc** when reading floors. (The one exception, μ→eγ, is a lepton-sector constraint
that quark-only mode drops entirely — see §4.)

### Recommended: single 100k tile-batch on one compute node (no array needed)

100k points = 10 M_KK tiles × 10,000 draws at a single `r`. This is exactly what the
existing non-array wrapper `scripts/run_wq_quarkonly.sbatch` already does. Pick the
nominal `r = 0.25` (scan default `DEFAULT_QUARK_FIT_R = 0.25`):

```bash
cd /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing
# Option A — direct (login/compute node, ~8-15 min on 48 cores):
conda activate ising_bootstrap
python scripts/run_full_catalog_scan.py \
    --quark-only \
    --ew-model minimal_rs \
    --output-dir scan_outputs/wq_quarkonly_100k_minimal_$(date -u +%Y%m%dT%H%M%S) \
    --n-draws 10000 \
    --m-kk-tev 1,2,3,5,7,10,15,20,30,50 \
    --quark-fit-r 0.25 \
    --n-workers 48 \
    --smoke-report scan_outputs/.../run_report.json \
    --smoke-report-md scan_outputs/.../run_report.md
```

```bash
# Option B — via the existing sbatch wrapper (sets all of the above; N_DRAWS=10000 default):
sbatch scripts/run_wq_quarkonly.sbatch
#   wrapper hardcodes --quark-only, minimal EW (default), M_KK=1,2,3,5,7,10,15,20,30,50,
#   N_DRAWS=10000/tile -> 100k total, 48 cores, 1h walltime cap.
```

`--ew-model minimal_rs` is the default, so it can be omitted; it is shown explicitly to be
unambiguous that this is the **MINIMAL (non-custodial)** model.

### Read out per-constraint M_KK floors (post-hoc, no re-scan)

```bash
python scripts/analyze_wq_quarkonly.py --input scan_outputs/wq_quarkonly_100k_minimal_<ts>
# -> plots/ + analysis_report.md with per-constraint veto fraction vs M_KK for every ID.
```
Then read the floors of just your 10 IDs (table in §3) off the per-constraint curves /
the `constraint_tallies` in the tile summaries. The joint floor for the subset is the
lowest M_KK where the cell-wise **max** veto fraction over those IDs drops below threshold
(Scan Explorer uses 0.5).

---

## 1. Driver / entry point and how the prior minimal 1M run was launched

- **Driver:** `scripts/run_full_catalog_scan.py` (`main()` → `ScanConfig` → tiled scan).
- **Per-point physics:** `quarkConstraints.fit.fit_quark_sector` fits anarchic O(1)
  Yukawas to reproduce quark masses + CKM each draw; the draw is **rejected at the fit
  gate** if `require_fit_success` fails or `fit_score > max_fit_score`
  (`run_full_catalog_scan.py:835-839`), and `ckm_residuals` are recorded
  (`:800 max_ckm_residual`). **This fit gate IS constraints #1 (quark masses) and
  #2 (CKM)** — they are not catalog IDs.
- **Prior MINIMAL 1M run = SLURM array job `20128400`**, launched via
  `scripts/run_wq_quarkonly_1m_array.sbatch`:
  - grid: `r ∈ {0.05,0.10,0.25,0.50,1.00}` × `M_KK ∈ {1,2,3,5,7,10,15,20,30,50} TeV`
    × 20,000 draws/(r,M_KK) = **1,000,000 draws**.
  - 50 array tasks (`--array=0-49%50`); each task = one (r, draw-shard) running all
    10 M_KK tiles × 2,000 draws. Grid/seeds come from
    `scripts/wq_quarkonly_1m_plan.py` (`task-env` / `summary` / `seed-proof`).
  - output root: `scan_outputs/wq_quarkonly_1M_20128400/` (confirmed present, with
    `scan_plan.json`, `r0p05/ … r1/`, per-shard `tile-*.jsonl` + `tile-*.summary.json`).
- **Custodial counterpart 1M run = `20675555`**
  (`scan_outputs/wq_quarkonly_1M_custodial_20675555/`), launched via
  `run_wq_quarkonly_1m_custodial_array.sbatch` with `--ew-model custodial_rs_plr`.

---

## 2. Draw count = 100k and MINIMAL model selection

- **Draw count:** `--n-draws N` is draws **per M_KK tile**. Total = `N × (#M_KK tiles) × (#r)`.
  With the standard 10 M_KK tiles at a single `r`: `--n-draws 10000` → **100,000 points**.
- **MINIMAL (non-custodial) EW model:** `--ew-model minimal_rs`.
  - Defined in `run_full_catalog_scan.py:168-170`:
    `MINIMAL_RS_EW_MODEL = "minimal_rs"`, `CUSTODIAL_RS_PLR_EW_MODEL = "custodial_rs_plr"`,
    `SUPPORTED_EW_MODELS = (minimal_rs, custodial_rs_plr)`.
  - argparse: `--ew-model` choices = those two, **default = `minimal_rs`**
    (`:1627`). This flag is the W9 custodial toggle: minimal = no tree-level custodial
    protection; custodial = `custodial_rs_plr`. For the user's "MINIMAL, NO custodial"
    request, use `minimal_rs` (or just omit the flag).

---

## 3. Constraint-ID mapping for the requested 10 items

The quark-only allowlist (`QUARK_ONLY_ALLOWLIST_IDS`, `run_full_catalog_scan.py:53-72`)
evaluates ~47 constraints. The 10 the user wants map as follows (titles read from
`flavor_catalog_constraints/primary/...`):

| # | Requested constraint | Catalog ID(s) | Title / observable | In quark-only allowlist? |
|---|---|---|---|---|
| 1 | Reproduce SM quark masses | (fit gate, not an ID) | `fit_quark_sector` mass fit, gated by `--max-fit-score` | n/a — always enforced |
| 2 | Reproduce CKM matrix | (fit gate, not an ID) | `fit_quark_sector` CKM fit; `max_ckm_residual` recorded | n/a — always enforced |
| 3 | ε_K | **K001** | indirect CP violation in K0–K̄0 mixing | YES |
| 4 | D0–D̄0 mixing | **C001** (mass/mixing) [+ **C002** for CPV] | C001 = neutral charm mixing; C002 = CPV in charm mixing | YES |
| 5 | Δm_d | **B001** | neutral B_d mixing mass difference | YES |
| 6 | Δm_s | **B003** | B_s mixing mass splitting | YES |
| 7 | μ→eγ | **L001** | charged-LFV radiative μ→eγ | **NO — lepton sector, dropped in quark-only** (see §4) |
| 8 | Z→bb | **T010** [+ **T011**] | T010 = Z→bb̄ R_b, A_b; T011 = A_b, A_FB^{0,b} | YES |
| 9 | S, T, U (oblique) | **EW001** | Peskin–Takeuchi oblique S,T,U | YES |
| 10 | Collider M_KK cutoff (simple) | **a Λ_IR/M_KK floor on the grid** (see §5); native KK-resonance IDs = CR001 (KK-gluon ttbar), CR007 (KK-graviton), CR012 (diboson), CR013 (diphoton) | approximate via grid floor, **not** a single toggle | partial |

Notes:
- D0 mixing: the mixing amplitude itself is **C001**; **C002** is the CP-violation piece.
  Include both if "D0–D̄0 mixing" is meant to cover CPV; C001 alone for pure Δm/mixing.
- Δm_d / Δm_s: B001 / B003 are the mass-splitting observables. (B002/B004 are the
  associated CP phases φ_d, φ_s — exclude unless wanted.)
- Z→bb: T010 carries R_b/A_b; T011 carries the asymmetries. Use both for full Z→bb̄.
- The other ~40 allowlist IDs (B011/B012 radiative, B032-34 charmless, K002/K003/K013,
  C003, T001-T008 top-FCNC, T012, EW002/EW003, CR0xx VLQ/resonance limits, E004-E009
  EDMs) are **also evaluated and will independently veto** in a default scan. Because
  floors are per-ID (§6), you simply ignore those columns when reading your 10.

---

## 4. Gaps / caveats

1. **μ→eγ (item 7) is NOT available in quark-only mode.** `--quark-only` drops the entire
   lepton sector (`QUARK_ONLY_LEPTON_SECTOR_LABEL = "dropped (not rigorous)"`); no `L0xx`
   IDs (L001 = μ→eγ) are in `QUARK_ONLY_ALLOWLIST_IDS`. The quark scan **cannot** produce a
   μ→eγ floor because μ→eγ needs swept lepton bulk-mass / Yukawa data this harness does not
   generate. Options:
   - Drop μ→eγ from the 10 for this quark-sector run (recommended; it is lepton-sector and
     orthogonal to the quark M_KK reach), OR
   - Run the full (non-quark-only) catalog scan with a lepton sweep to get L001 — a
     different, much heavier harness, outside the quark-only 100k recipe here.
2. **"Simple collider M_KK cutoff" (item 10) is not a single native toggle.** The catalog's
   collider IDs (CR001/CR007/CR012/CR013, plus the VLQ pair-production limits CR002-CR010)
   are full recast/limit constraints, not a flat M_KK floor. To get "just a simple mass
   cutoff," apply a hard floor on the scan grid: keep only tiles with
   `M_KK ≥ M_KK^min` (e.g. drop the 1-3 TeV tiles, keep ≥ ~5 TeV per current di-jet/ttbar
   KK-gluon reach). Equivalently restrict `--m-kk-tev` at launch, or filter rows by
   `params.M_KK` post-hoc. No code change needed.
3. **No precomputed scalar "M_KK floor per constraint" field exists.** Every layer stores
   per-(r, M_KK, constraint) **veto fractions / vetoed-evaluated counts**; the floor is the
   derived lowest M_KK where the fraction crosses below threshold (Scan Explorer uses 0.5).

---

## 5. Grid, Λ_IR↔M_KK conversion, ξ_KK

- **Standard grid** (both 1M runs, and `run_wq_quarkonly.sbatch`):
  - `M_KK = {1, 2, 3, 5, 7, 10, 15, 20, 30, 50} TeV` (10 tiles) via `--m-kk-tev`.
  - `r = {0.05, 0.10, 0.25, 0.50, 1.00}` (the MFV up-spurion weight in
    `C_Q = r·Yu Yu† + Yd Yd†`); default single value `--quark-fit-r 0.25`.
  - A 100k run uses ONE `r` × 10 M_KK × 10k draws. To mirror the 1M's r-sweep at 100k
    total you would instead do 5 r × 10 M_KK × 2k draws (use the array plan with reduced
    draws), but the simplest 100k is single-r.
- **ξ_KK = `DEFAULT_XI_KK = 2.4487`** (`run_full_catalog_scan.py:145`).
  Conversion in `_build_tiles`: **`Lambda_IR = M_KK / xi_kk`** (grid is specified as the
  physical first-KK gauge mass M_KK; Λ_IR is derived). Equivalently `M_KK = 2.4487·Λ_IR`.
  The website rounds the convention note to `m1 ≈ 2.45·Λ_IR`.
- **Floor readout locations:**
  - per-tile `constraint_tallies` dict in `tile-*.summary.json` (keys = constraint IDs,
    each `{points, evaluated, active, failed, vetoed, severity, tag}`).
  - `scripts/analyze_wq_quarkonly.py` → `_plot_constraint_veto_by_mkk` (per-ID veto
    fraction vs M_KK, grouped by exact `(r, M_KK)`).
  - `scripts/build_wq_quarkonly_comparison.py` → `constraint_veto_by_r_mkk.csv`
    (cols: `ew_model, r, mkk_tev, constraint_id, tag, severity, …, veto_fraction`).
  - `flavor_catalog/website/scripts/build_scan_explorer.py` → `scan_explorer.json`
    `veto[model][r][cid]` arrays over the M_KK grid; `FLOOR_THRESHOLD = 0.5`;
    `_build_bare_floor` = lowest M_KK with any evaluated draws.

---

## 6. Per-constraint isolation (the key enabling fact)

**A single full quark-only scan suffices; no constraint-restricted re-run is needed.**
The shared veto predicate everywhere is:
`active AND evaluated AND severity=="HARD" AND passes is False AND tag in {rigorous, proxy}`.
Veto counts are keyed by constraint ID at every layer (tile `constraint_tallies`, the
analyze script's `(r,M_KK,pid)` stats, the comparison CSV, and the Scan Explorer veto
tree). To get the binding floor for the requested subset {K001, C001(+C002), B001, B003,
T010(+T011), EW001}, take the **per-(r,M_KK) cell-wise max veto fraction over those IDs**
and find the lowest M_KK where it drops below threshold. Everything else (the other ~40
allowlist IDs) is just ignored downstream. (μ→eγ/L001 is the lone exception — not produced
at all in quark-only mode; see §4.)

---

## 7. Compute cost and wall-clock (100k points)

Hard data from prior minimal array `20128400` (per task = 20k draws, 48 cores):
- `sacct 20128400_0`: Elapsed `00:07:51`, TotalCPU `01:18:08` (= 78.1 CPU-min) for 20k pts.
- ⇒ **~0.234 core-seconds/point** (consistent with the `run_wq_quarkonly.sbatch` note of
  ~0.28 s/point post-cache).

For **100,000 points**:
- CPU cost ≈ 100,000 × 0.234 s ≈ **23,400 core-seconds ≈ 6.5 core-hours**
  (≈ 7.7 core-hours at the conservative 0.28 s/pt).
- **Single-core wall:** ~6.5 h (not recommended).
- **48-core wall:** 23,400 / 48 ≈ **8.1 min** (~10–15 min with startup/cache warmup).
- A SLURM **array is unnecessary** for 100k: one 48-core node finishes in ~10 min.

### Does it need sbatch? — No.
A 100k run fits comfortably on a single interactive compute node (`salloc -p
serial_requeue -c 48 --mem 32G -t 0:30:00`) or via the one-shot
`sbatch scripts/run_wq_quarkonly.sbatch` (1h cap, 48 cores, 32G — already sized for exactly
this 10×10k=100k job). Avoid running heavy multi-worker jobs on the login node; use
`salloc` or `sbatch`.

---

## 8. Bottom line

1. Run `sbatch scripts/run_wq_quarkonly.sbatch` (or the explicit `python
   scripts/run_full_catalog_scan.py --quark-only --ew-model minimal_rs --n-draws 10000
   --m-kk-tev 1,2,3,5,7,10,15,20,30,50 --quark-fit-r 0.25 --n-workers 48 …`) — MINIMAL
   model, 100k points, ~10 min on 48 cores.
2. Quark masses (#1) + CKM (#2) are the per-draw fit gate (always on).
3. Items #3,5,6,8,9 = K001, B001, B003, T010(+T011), EW001 — already in the quark-only
   allowlist and reported per-ID. Item #4 (D0) = C001 (+C002 for CPV).
4. Item #7 (μ→eγ = L001) is **not available** in quark-only mode — drop it or run the
   lepton-sweep harness separately.
5. Item #10 (collider) — approximate with a hard `M_KK ≥ M_KK^min` grid floor (restrict
   `--m-kk-tev` or filter on `params.M_KK`); the native CR0xx resonance recasts are NOT a
   simple cutoff.
6. Isolate your 10 by reading per-ID veto fractions from
   `analyze_wq_quarkonly.py` / the comparison CSV / `scan_explorer.json` — no re-scan.

---

# ADDENDUM (2026-06-22): include μ→eγ + collider cutoff — decisive answer to the two follow-ups

This addendum supersedes §4 item 1 and §4 item 2 of the original recon where they
say μ→eγ is unavailable and collider is "not a toggle". After reading the harness
source (`scripts/run_full_catalog_scan.py`) and the collider review
(`review_local/collider_kk_review.tex`), the decisive findings are below.

## A. INCLUDING μ→eγ (L001): use FULL-CATALOG mode (drop `--quark-only`)

**Key fact:** `--quark-only` is the ONLY thing that drops the lepton sector.
Running WITHOUT `--quark-only` (full-catalog mode) already does everything needed:

1. It draws lepton parameters every point via `_draw_lepton_inputs`
   (`run_full_catalog_scan.py:744-761`): `c_L`, `c_E[3]`, `c_N` ~ U(c_min,c_max);
   `log10 M_N` ~ U(log10 m_n_min, log10 m_n_max); ordering by
   `normal_ordering_probability`; Majorana phases α,β ~ U(0,2π); lightest ν mass
   ~ U(lightest_nu_min, lightest_nu_max).
2. It computes the full lepton Yukawa set incl. PMNS-induced off-diagonal
   structure via `compute_all_yukawas(...)` (`:564-576`) and rejects
   non-perturbative draws (`_require_perturbative_leptons`).
3. It builds the point with `include_charged_current=True`,
   `include_fermion_kk_mixing=True`, `include_higgs_yukawas=True` (`:591-594`) —
   i.e. it generates exactly the `lepton_lmfv_parameters` /
   `lepton_mass_basis_couplings` / `rs_higgs_yukawas` extras that quark-only mode
   FORBIDS (`QUARK_ONLY_FORBIDDEN_EXTRAS`, `:75-83`).
4. It then calls `registry.evaluate_all(point)` (`:609`) over the **entire**
   registry — so **L001 (μ→eγ) is evaluated**, and so is the whole lepton sector
   L0xx, plus the held-out lepton-involving CR/EW IDs.

There is NO separate "lepton+quark combined flag." Full mode = the combined run.
**Dropping `--quark-only` is exactly the combined run that evaluates L001 while
still drawing quark params (the quark fit gate is identical in both modes,
`:540-551`) and evaluating K001, C001/C002, B001, B003, T010/T011, EW001.**

### Are the lepton draw defaults sensible / non-degenerate?  YES.

The defaults (`run_full_catalog_scan.py:154-159`) ARE the W7 μ→eγ LMFV draw config:
- `c_min=0.3, c_max=0.9`  → bulk-mass localization spread (`--c-min/--c-max`)
- `M_N ∈ [1e10, 1e15] GeV` log-uniform (`--m-n-min-gev/--m-n-max-gev`)
- `lightest_nu ∈ [0, 1e-2] eV` (`--lightest-nu-min-ev/--lightest-nu-max-ev`)
- random PMNS Majorana phases α,β ∈ [0,2π); ordering 50/50.

μ→eγ is **physically meaningful, not degenerate**, here because `c_E[3]` are drawn
*independently per generation* (non-universal RH charged-lepton localization) and
the neutrino Yukawa carries the PMNS rotation — both are required to generate a
nonzero (e,μ) dipole. (The degenerate "sanity" path
`_force_degenerate_neutrino_yukawas`, `:854`, is used ONLY in the universal-c
self-test, never in the production draw loop.) No extra flags are needed: the
bare full-mode command produces a meaningful L001.

### FULL combined command line (single 100k run, MINIMAL model)

```bash
cd /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing
conda activate ising_bootstrap
TS=$(date -u +%Y%m%dT%H%M%S)
OUT=scan_outputs/full_catalog_100k_minimal_${TS}
python scripts/run_full_catalog_scan.py \
    --ew-model minimal_rs \
    --output-dir ${OUT} \
    --n-draws 10000 \
    --m-kk-tev 1,2,3,5,7,10,15,20,30,50 \
    --quark-fit-r 0.25 \
    --n-workers 48 \
    --smoke-report ${OUT}/run_report.json \
    --smoke-report-md ${OUT}/run_report.md
```
(NO `--quark-only`. `--ew-model minimal_rs` is default; shown for clarity. All
lepton-draw flags omitted → W7 defaults above. 10 tiles × 10k = 100,000 points,
single r=0.25.) Run under `salloc -p serial_requeue -c 48 --mem 48G -t 1:00:00`
or wrap in a one-off sbatch; do NOT run 48-worker on the login node.

### Can we still isolate the user's 10 post-hoc in full mode?  YES, with one caveat.

- Every per-draw JSONL row carries a complete per-ID `constraints` dict
  (`passes/active/evaluated/severity/tag/ratio` per pid, `_classify_results`
  `:877-906`), so L001 + all 10 IDs are recoverable per point and per (r,M_KK)
  tile. Per-ID M_KK floors are obtained by aggregating those rows.
- **Caveat:** the convenient pre-aggregated `constraint_tallies` block in the tile
  `*.summary.json` is built **only in quark-only mode**
  (`constraint_tallies = … if cfg.quark_only else None`, `:310-311`). In full mode
  that pre-rolled tally is absent, and `scripts/analyze_wq_quarkonly.py` /
  `build_wq_quarkonly_comparison.py` are quark-only-shaped readers. So full-mode
  per-ID floors require a small custom aggregation pass over the raw `tile-*.jsonl`
  rows (group by (r, M_KK), count `severity==HARD & evaluated & !passes &
  tag∈{rigorous,proxy}` per pid). This is a ~30-line script, NOT a re-scan.

### Cost / wall-clock for full mode (re-estimate)

The per-tile RS-EW spline cache build and the per-draw quark fit (the cost
drivers) are identical to quark-only. Full mode adds: one lepton Yukawa solve
(cheap linear algebra) + ~56 extra constraint evals (arithmetic on precomputed
couplings) ⇒ estimate **~1.3–1.8× the quark-only per-point time**, i.e.
**~0.3–0.45 core-s/point** (vs 0.234 quark-only).
- 100k points ⇒ ~30k–45k core-s ≈ **8.5–12.5 core-hours**.
- **48-core wall: ~11–16 min** (+cache warmup ⇒ budget ~20 min). Still one node,
  no array. (No prior production full-mode timing exists at scale, so treat this
  as an estimate; the 1h salloc/sbatch cap is ample.)

## B. RECOMMENDATION: ONE combined full-catalog run (Approach 1), not two runs

A single full-catalog 100k run covers **all 10** items cleanly:
- #1 quark masses, #2 CKM = the per-draw quark fit gate (always on, both modes).
- #3 K001, #4 C001(+C002), #5 B001, #6 B003, #8 T010(+T011), #9 EW001 = quark/EW
  catalog IDs, evaluated identically to quark-only.
- #7 L001 (μ→eγ) = evaluated because full mode draws leptons (this is the whole
  point of NOT passing `--quark-only`).
- #10 collider = impose as a flat M_KK cutoff (see §C), no toggle needed.

The two-run alternative (quark-only 100k for the 7 quark/EW items + a separate
lepton μ→eγ floor harness) is **NOT cleaner**: it doubles bookkeeping and there is
no standalone "L001 M_KK floor" harness — μ→eγ only exists inside the full point
builder. **Recommend the single combined full-catalog run.** The only added work
vs. the quark-only recon is the ~30-line post-hoc tally aggregator for per-ID
floors (the quark-only analyze scripts won't read full-mode summaries).

## C. COLLIDER cutoff: the concrete number + convention

**The repo's scan grid variable `M_KK` (= `--m-kk-tev`) IS the physical first
gauge KK resonance mass m1.** Proof: `_build_tiles` sets
`lambda_ir = M_KK / xi_kk` (`:1608`), i.e. `M_KK = ξ_KK·Λ_IR` with
`ξ_KK = 2.4487` (`DEFAULT_XI_KK`, `:145`). The collider adapter compares the
experimental edge directly to `M_KK/1000` in TeV
(`kk_mass_tev_from_m_kk_gev`, `quarkConstraints/collider_resonance.py:109-112`:
`return m_kk_gev/1000.0`). So **there is no separate "ξ·M_KK"**: in this repo's
naming the grid `M_KK` already equals m1 = ξ·Λ_IR. The two conventions the
follow-up asks for therefore are:

| Convention | Symbol | Relation | Collider cutoff value |
|---|---|---|---|
| Physical first gauge KK resonance | m1 (= repo grid `M_KK`) | m1 = 2.4487·Λ_IR | **m1 ≥ 5.5 TeV** |
| IR/warp scale | Λ_IR | Λ_IR = m1/2.4487 | **Λ_IR ≥ 2.25 TeV** (5.5/2.4487) |

**Number to use: 5.5 TeV on the physical resonance mass.** Source: the KK-gluon→ttbar
all-topology search **CMS-B2G-25-009 (arXiv:2603.23454)**, which sets the active
HARD edge in `CR001.yaml`; the collider review states it explicitly:
`review_local/collider_kk_review.tex` lines 352–359 ("records the upper edge
$5.5$~TeV as the active budget … a model with $\MKK=1$~TeV is below the $5.5$~TeV
edge and is excluded by CR001; $\MKK\gtrsim 30$~TeV clears it"). 5.5 TeV is the
**highest single direct edge** in the catalog (review §floor, lines 466–469: KK
gluon ~5.5 TeV is the top of the near-flat direct floor), so using it as THE
cutoff is the conservative single-number choice.

Convention caveat for completeness (review lines 305, 311–314, 826–833): the KK
**graviton** edges (CR007/CR013) carry the spin-2 root x1^grav ≈ 3.83 vs the
gauge root x1^gauge ≈ 2.45, so a graviton resonance at a given Λ_IR is ~1.56×
heavier than the gauge KK. The 5.5 TeV cutoff is on the **gauge** resonance m1,
which is the right one for a single flat M_KK cutoff (KK gluon is the binding
direct channel). No graviton conversion is needed if you cut on the gauge m1.

### How to impose it (simple cutoff, no code change)

Two equivalent options:
1. **At launch — restrict the grid** to tiles at/above the cutoff:
   `--m-kk-tev 7,10,15,20,30,50` (drops 1,2,3,5 TeV; 7 TeV is the first grid tile
   ≥ 5.5 TeV). The exact-5.5 tile can be added: `--m-kk-tev 5.5,7,10,15,20,30,50`.
2. **Post-hoc — filter rows** by `params.M_KK ≥ 5500` GeV when reading floors.

**Reporting rule:** report the per-constraint physics floor for each of the user's
flavour/EW items, then the headline reach as
`floor = max(physics_floor, collider_cutoff)` with `collider_cutoff = 5.5 TeV`
(on m1) ⇔ `Λ_IR = 2.25 TeV`. In the anarchic minimal model the rigorous flavour
floor (~tens of TeV) sits far above 5.5 TeV, so the collider cutoff is rarely
the binding constraint — but it is the right floor for any point where flavour is
relaxed.

## D. BOTTOM LINE (the two follow-ups, answered)

1. **μ→eγ:** run the FULL-catalog scan (omit `--quark-only`). It draws leptons
   with the W7 defaults (meaningful, non-degenerate L001) and evaluates all 10
   items in one run. Per-ID floors via a ~30-line aggregation over the JSONL rows
   (the quark-only analyze scripts won't read full-mode summaries). ~11–16 min on
   48 cores (~8.5–12.5 core-h). This single run is the recommended approach.
2. **Collider cutoff:** impose **m1 ≥ 5.5 TeV** (physical first gauge KK
   resonance; repo grid `M_KK`), equivalently **Λ_IR ≥ 2.25 TeV**. Source:
   CR001 / CMS-B2G-25-009, recorded in the collider review
   (`review_local/collider_kk_review.tex` L352–359, L466–469). Apply by
   restricting `--m-kk-tev` to ≥5.5 TeV tiles or filtering `params.M_KK≥5500`,
   and report `floor = max(physics_floor, 5.5 TeV)`.
