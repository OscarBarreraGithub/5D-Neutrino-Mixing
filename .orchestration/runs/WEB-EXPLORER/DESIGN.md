# Scan Explorer tab — DESIGN SPEC (locked decisions)

A new interactive tab in the warpedflavor Astro site (`flavor_catalog/website/`) that visualizes
the quark-only 1M scan results: minimum M_KK to satisfy constraints (vs r), with per-constraint
toggles + custodial on/off, and the fitted Yukawas per r. LOCALHOST ONLY — do NOT deploy/push.

Data sources (already on disk, committed):
- Non-custodial: `scan_outputs/wq_quarkonly_1M_20128400/`
- Custodial:     `scan_outputs/wq_quarkonly_1M_custodial_20675555/`
- Their paired comparison: `scan_outputs/wq_quarkonly_1M_custodial_20675555/comparison/`
- Each raw row (r*/shard-*/tile-*.jsonl) has: `quark_fit_r`, `params.M_KK` (GeV), `seed`,
  `survives_all_HARD_strict/inclusive`, `skipped`, `constraints.<ID>{passes,ratio,tag,severity,active,evaluated}`,
  `fit_diagnostics.{fitted_up_yukawa_singular_values[3], fitted_down_yukawa_singular_values[3],
  bulk_c_Q[3],bulk_c_u[3],bulk_c_d[3], max_abs_quark_yukawa}`, `params.quark_yukawa_seed`.
- Grids: r ∈ {0.05,0.1,0.25,0.5,1.0}; M_KK ∈ {1,2,3,5,7,10,15,20,30,50} TeV (physical, = 2.45·Λ_IR).

## Part A — data artifact: `src/content/scan_explorer.json`
A Python extraction script (`flavor_catalog/website/scripts/build_scan_explorer.py`) reads BOTH scan
roots and emits compact JSON. Schema:

```
{
  "meta": {
    "r_grid": [0.05,0.1,0.25,0.5,1.0],
    "mkk_grid_tev": [1,2,3,5,7,10,15,20,30,50],
    "minimal_root": "...", "custodial_root": "...",
    "floor_threshold": 0.5,            // veto_fraction <= this counts as "satisfied" (median point)
    "mkk_convention": "physical first KK mass m1 = 2.45 * Lambda_IR",
    "n_draws_per_cell": 20000, "generated_from_git": "<sha>"
  },
  "constraints": [   // ONLY the quark constraints that actually bite somewhere (veto_fraction>threshold at some r,mkk)
    {"id":"T010","label":"Z → b b̄","tag":"rigorous","group":"electroweak","default_on":true},
    {"id":"EW001","label":"Oblique S,T,U","tag":"proxy","group":"electroweak","default_on":true},
    {"id":"K001","label":"ε_K","tag":"rigorous","group":"meson_mixing","default_on":true},
    {"id":"CR001","label":"KK gluon → t t̄","tag":"proxy","group":"collider","default_on":true},
    ...
  ],
  // veto_fraction over EVALUATED (non-skipped) rows, per ew_model/r/constraint/mkk:
  "veto": {
    "minimal":   { "<r>": { "<cid>": [vf@1,vf@2,...,vf@50] } },   // array aligned to mkk_grid_tev
    "custodial": { "<r>": { "<cid>": [...] } }
  },
  // bare baseline = masses+CKM+fit only (no flavor constraint): the min grid M_KK with a non-skipped
  // fit (i.e. fit succeeds). Per ew_model/r:
  "bare_floor_tev": { "minimal": {"<r>": <mkk>}, "custodial": {"<r>": <mkk>} },
  // fitted Yukawa singular-value summary among EVALUATED rows, per ew_model/r/mkk:
  // each of up/down is 3 entries [p25,p50,p75] for the 3 singular values (ascending).
  "yukawa": {
    "minimal":   { "<r>": { "<mkk>": { "up":[[p25,p50,p75]x3], "down":[[..]x3], "n":N } } },
    "custodial": { "<r>": { "<mkk>": {...} } }
  },
  // 3 representative SURVIVING points, full 3x3 |Y| reconstructed from seed (same draw path the scan
  // uses) and VERIFIED: svd(|Y|) must match the stored fitted_*_singular_values to <1e-6.
  // NOTE: Yu_abs/Yd_abs are ENTRY MAGNITUDES |Y_ij| of the complex anarchic spurion (standard
  // display). up_singular/down_singular are the singular values of the COMPLEX matrix Y (NOT of
  // |Y| — abs drops phase, so SVD(|Y|) != these). Faithfulness check the builder MUST pass:
  // SVD of the reconstructed COMPLEX Y_u/Y_d == stored fitted_*_singular_values to <1e-6.
  // The front-end must label the matrix as |Y_ij| and the singular values as "of the complex Y".
  "rep_matrices": [
    {"label":"low r (down-localized)","ew_model":"custodial","r":0.05,"mkk_tev":3,"seed":...,
     "Yu_abs":[[..3x3..]], "Yd_abs":[[..]], "up_singular":[..], "down_singular":[..]},
    {"label":"r = 0.25","...":...},
    {"label":"r = 1.0 (up-dominated)","...":...}
  ]
}
```
Notes:
- Floor for a constraint = first mkk in grid with veto_fraction <= meta.floor_threshold (front-end
  derives it from the `veto` arrays so the threshold is transparent; ship arrays, not just floors).
- A constraint that never reaches the threshold even at 50 TeV → floor = null (front-end shows ">50").
- A constraint with veto_fraction 0 everywhere → omit from `constraints` list (never bites).
- Only include constraints in `constraints[]` that bite for at least one (ew_model,r,mkk).
- Keep the JSON well under ~1 MB (summaries only, NO per-draw rows).

## Part B — front-end: `src/pages/explore.astro` + navbar link "Explore" in BaseLayout.astro
Static Astro page, vanilla client JS, hand-rolled inline SVG (NO new npm deps), match existing
global.css tokens + dark/light theming. KaTeX (already loaded by layout) for matrices.

### Controls (shared state for BOTH plots)
- A checklist of the `constraints[]` (label + a small rigorous/proxy chip), each a toggle.
  Group visually by `group` (electroweak / meson_mixing / collider / ...).
- "Toggle all" button (all on / all off).
- "Custodial" switch: off = minimal dataset, on = custodial dataset.
- State drives both plots; changing anything re-renders both immediately (no reload).

### Plot 1 — "Minimum M_KK vs r"
- X-axis: r (the 5 grid values; spaced on a log-ish/categorical axis, labeled 0.05…1.0).
- Y-axis: minimum M_KK [TeV], LOG scale, 1–50 (+ a ">50" cap row).
- Bold ENVELOPE line = per r, max over ACTIVE constraints of floor_c(r), floored at bare_floor.
  This is the headline "min M_KK to satisfy the selected set".
- Faint thin colored line per ACTIVE constraint = its own floor_c(r), SHOWN BY DEFAULT (the
  "many plots superimposed" view), with a legend. A "Show all constraint traces" checkbox (default
  CHECKED) can collapse to envelope-only for clarity. Hovering/focusing a legend row highlights that
  one constraint's trace (full opacity) and dims the rest. (This default was chosen after the CSS
  root-cause fix made the traces render as clean thin lines rather than a muddy mass.)
- Custodial vs minimal must visibly differ (e.g. minimal envelope ~20–30 TeV when Z→bb on,
  custodial ~2–3 TeV when Z→bb on). A short caption states the headline.

### Plot 2 — "Fitted Yukawas per r" (inherits the selection)
- A grid of boxes, ONE PER r (5 boxes; responsive — up to 3–4 across, wraps cleanly; NOT superimposed).
- Each box header: "r = <value>", and "M_KK = <envelope M_KK for this r>" (from Plot 1's current
  envelope + custodial state). As toggles/custodial change, the envelope M_KK per r changes and each
  box reads the precomputed `yukawa[ew_model][r][envelope_mkk]` summary.
- Box body: grouped mini bar chart of the 3 up + 3 down fitted Yukawa singular values (p50 bar with
  a p25–p75 whisker), log-y, clearly labeled y1<y2<y3 for up (warm) and down (cool). n shown.
- If envelope M_KK has no surviving sample (shouldn't happen at envelope), show "—".
- BELOW the grid: the 3 `rep_matrices`, each a nicely formatted 3×3 |Y| matrix (KaTeX bmatrix or a
  styled HTML grid), labeled with its (r, M_KK, ew_model) and its singular values. Show |Y_u| (and
  optionally |Y_d|) — pick a clean, readable layout. These are static illustrations of the anarchic
  structure (they do NOT need to react to toggles).

## Quality bar
- Clean, matches the site's typography/spacing/dark-light theme; no layout overflow; no console errors.
- Dependency-free (SVG by hand); deterministic; honest labels (tag chips, the floor-threshold note,
  the M_KK convention note).
- Dual sign-off: codex + Opus review BOTH the data extraction (faithful to the scans) AND the
  front-end (renders correctly, interactivity correct, clean). Render-verified on localhost.
- Do NOT push or deploy.
