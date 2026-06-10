# Scan Explorer front-end implementation summary

## Files changed

- `flavor_catalog/website/src/pages/explore.astro`
  - Added the static Scan Explorer page.
  - Imports `src/content/scan_explorer.json` without modifying it.
  - Renders grouped constraint controls, all-on/all-off toggle, and minimal/custodial switch.
  - Computes per-constraint floors from `veto[ew_model][r][constraint_id]` using `meta.floor_threshold`.
  - Computes the envelope as the max active floor per `r`, floored by `bare_floor_tev`.
  - Draws the minimum `M_KK` envelope and active constraint lines with hand-rolled inline SVG.
  - Renders one Yukawa mini chart per `r`, reading `yukawa[ew_model][r][envelope_mkk]`.
  - Renders representative `|Y_u|` and `|Y_d|` matrices using KaTeX.
  - Uses page-scoped CSS and existing global theme tokens for light/dark mode.
- `flavor_catalog/website/src/layouts/BaseLayout.astro`
  - Added the navbar link `Explore` with `href="/explore/"` next to `Browse`.
- `.orchestration/runs/WEB-EXPLORER/frontend_impl_summary.md`
  - This summary.

## Verification

- Ran `npm run build` in `flavor_catalog/website`.
  - Astro build completed successfully.
  - 114 pages built.
  - Pagefind completed successfully and indexed 102 catalog entry pages.
- Ran a headless Chromium smoke test against `npm run preview` for:
  - `/explore/?theme=light`
  - `/explore/?theme=dark`
- The browser-rendered DOM contained populated envelope SVG paths, active legend rows, KaTeX-rendered labels/matrices, and Yukawa mini-chart SVGs. Chromium printed DBus warnings from the headless Linux environment, but the page rendered and the client script populated the expected DOM.

## Caveats

- No JSON or data-generation script changes were made.
- No commit, push, or deploy was performed.

FE-READY
