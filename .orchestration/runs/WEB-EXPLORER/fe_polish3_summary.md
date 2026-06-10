POLISH ITERATION 3 SUMMARY

Implemented the Astro scoping root-cause fix in `flavor_catalog/website/src/pages/explore.astro`.

Moved the JS-created SVG presentation rules into a `<style is:global>` block with container-id prefixes:

- `#minimum-plot`: grid lines, cap grid, axis lines/labels/titles, point labels, constraint paths/markers, envelope path, envelope points.
- `#yukawa-grid`: mini Yukawa SVG root, mini grid/axis elements, group labels, axis title, median bars, whiskers, whisker caps, up/down bar colors.

Kept static page layout/control/card CSS in the scoped style block. Also moved the JS-created legend child selectors under `#constraint-legend` because this implementation creates legend rows, swatches, ID chips, and labels with DOM APIs; leaving those scoped made swatches invisible and rows unaligned.

Verification:

- `npm run build` passes with 0 errors.
- Chromium CDP computed-style check on `http://127.0.0.1:4323/explore/?theme=light` and `?theme=dark` after client JS render:
  - Envelope path fill is `none`.
  - Envelope stroke equals `--color-accent` and stroke width is `3.4px`.
  - Light theme accent/envelope stroke: `rgb(108, 74, 216)`.
  - Dark theme accent/envelope stroke: `rgb(161, 140, 242)`.
  - With "Show all constraint traces" off, initial constraint path count is `0`.
  - Legend hover creates exactly one highlighted constraint path with `fill:none`, colored stroke, and `2.6px` stroke width.
  - "Show all constraint traces" creates 17 low-opacity colored paths with `fill:none` and `1.1px` stroke width.
  - Plot 2 has 30 median bars: up bars are orange (`rgb(210, 132, 34)`), down bars are blue (`rgb(47, 126, 188)`), with colored whiskers and `fill:none` on whisker elements.
  - SVG black-fill scan returned no default black fills in either theme.
  - Legend has 17 visible swatches and rows compute to a 3-column grid `[swatch][ID chip][description]`.

No JSON/data files were touched. No commit or push performed.

POLISH3-READY
