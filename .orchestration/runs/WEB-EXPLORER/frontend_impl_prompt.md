# Scan Explorer — FRONT-END build (codex)

Read `.orchestration/runs/WEB-EXPLORER/DESIGN.md` (Part B is authoritative) and the produced
`flavor_catalog/website/src/content/scan_explorer.json` (the real data — match its exact shape).
Build the front-end. Dependency-free (hand-rolled inline SVG, vanilla client JS); match the existing
site (global.css tokens, dark/light theme); KaTeX is already loaded by the layout.

## Deliverables (only these files)
1. `flavor_catalog/website/src/pages/explore.astro` — the new page (title "Scan Explorer"), using
   `BaseLayout`, importing `../content/scan_explorer.json`, embedding it for the client script (e.g.
   `<script type="application/json" id="scan-data">` or `define:vars`), and rendering:
   - **Controls** (shared state): a grouped checklist of `constraints[]` (label + small rigorous/proxy
     chip, grouped by `group`), a "Toggle all" button, and a "Custodial" on/off switch. All changes
     re-render BOTH plots immediately (no reload).
   - **Plot 1 "Minimum M_KK vs r"**: x=r (5 grid values, log-ish/categorical), y=min M_KK [TeV] LOG
     1–50 (+ a ">50" cap). Bold ENVELOPE line = per r, max over ACTIVE constraints of floor_c(r)
     (floor_c = first mkk in meta.mkk_grid_tev with veto[ew][r][cid] <= meta.floor_threshold; null =>
     ">50"), floored at bare_floor_tev. Faint thin colored line per ACTIVE constraint = its floor_c(r),
     with a legend. Custodial vs minimal must visibly differ. Short caption with the headline.
   - **Plot 2 "Fitted Yukawas per r"**: a responsive grid, ONE box per r (NOT superimposed). Each box
     header shows "r=<v>" and "M_KK=<envelope M_KK for this r>" (from Plot 1's current envelope +
     custodial). Body: grouped mini bar chart of the 3 up + 3 down fitted singular values
     (yukawa[ew][r][envelope_mkk]: p50 bar + p25–p75 whisker), log-y, labeled, n shown. Reacts to
     toggles/custodial. If no sample, show "—".
   - **Below Plot 2**: the 3 `rep_matrices` as nicely formatted 3×3 |Y| matrices (KaTeX bmatrix or a
     styled HTML grid), labeled with (r, M_KK, ew_model) + singular values. Static (no toggle react).
   - A short notes line: floor-threshold (median point) + M_KK convention from meta.
2. Add a navbar link "Explore" (href `/explore/`) in `src/layouts/BaseLayout.astro` next to "Browse".
3. Any page-specific CSS: prefer a `<style>` block in explore.astro using existing CSS variables;
   only touch global.css if truly necessary.

## Quality + verification
- Must build cleanly: run `npm run build` in `flavor_catalog/website/` (astro build) — 0 errors.
  (A dev server may already be running on 4321; build is the reliable check.)
- No console errors; no layout overflow; works in BOTH light and dark theme; numbers/labels honest.
- The envelope math, the per-constraint floor derivation, and the Plot-2 "read yukawa at envelope
  M_KK" logic must be correct and match the DESIGN.
- Do NOT commit. Do NOT push/deploy. Do NOT modify scan_explorer.json or the data script.

Write `.orchestration/runs/WEB-EXPLORER/frontend_impl_summary.md` (files changed, build result,
any caveats) ending `FE-READY` (npm run build clean) or `FE-BLOCKED`.
