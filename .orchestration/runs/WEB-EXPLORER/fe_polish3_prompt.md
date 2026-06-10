# Scan Explorer front-end — POLISH ITERATION 3: ROOT-CAUSE FIX (codex)

## The bug (confirmed via built CSS + post-JS DOM)
The `<style>` block in `flavor_catalog/website/src/pages/explore.astro` is Astro-SCOPED. Astro
rewrites its selectors to require a `data-astro-cid-jsy7jxlt` attribute (e.g. built CSS contains
`.envelope-path[data-astro-cid-jsy7jxlt]`). But ALL the SVG elements in both plots are created at
runtime by the client script via `document.createElementNS` (`svgEl(...)`), so they NEVER receive
that `data-astro-cid` attribute. Result: every JS-created SVG element ignores its intended CSS and
falls back to the browser SVG default `fill:black; stroke:none`. That is why:
- the envelope renders as a BLACK FILLED wedge (its CSS `fill:none; stroke:var(--color-accent)`
  never applies — no inline stroke either),
- constraint lines showed colored strokes (set inline) but were ALSO black-filled (CSS fill:none
  never applied), producing the muddy mass,
- Plot 2 bars/whiskers and any class-based colors are similarly not receiving their CSS.

## The fix (do this cleanly)
Make the styling of all JS-created SVG elements NOT depend on Astro's scoped CSS. Preferred approach:
- Move EVERY CSS rule that targets a dynamically-created SVG element (envelope-path, constraint-path,
  constraint-marker, envelope-point, point-label, grid-line, axis-line, axis-label, axis-title, and
  all Plot-2 yukawa bar/whisker/axis/label classes, and any matrix SVG if applicable) into a
  SEPARATE `<style is:global>` block. To prevent leakage, PREFIX every moved selector with the plot
  container id (`#minimum-plot ...`, `#yukawa-grid ...`, etc.) so they only match inside those
  containers. Keep all the NON-SVG rules (controls, checklist, legend rows/chips/swatch, panels,
  captions, layout — these are real template elements) in the existing SCOPED `<style>` block.
- Equivalent acceptable alternative: set the presentation attributes (fill, stroke, stroke-width,
  opacity, font-size, text-anchor, etc.) INLINE on every svgEl(...) call. But the is:global +
  id-prefixed approach is cleaner — pick that unless you have reason not to.

## After the fix, VERIFY (this is the whole point)
- Envelope renders as a clean PURPLE (`--color-accent`) line ~3.4px, NO fill, with markers + labels.
- With "Show all constraint traces" OFF: only the envelope shows (no constraint lines).
- Hovering a legend row reveals exactly that one constraint's line in ITS color with `fill:none`
  (a clean line, not a black blob), markers where it bites; leaving removes it.
- With "Show all traces" ON: thin low-opacity colored lines (fill:none), no black fills.
- Plot 2 bars render in their intended colors (up vs down distinguishable), whiskers visible, NO
  stray black fills from the default.
- Legend swatches visible; rows left-aligned [swatch][ID chip][description].
- Both light AND dark theme correct. `npm run build` passes (0 errors).

Do NOT touch the JSON/data script. Do NOT commit/push. Write
`.orchestration/runs/WEB-EXPLORER/fe_polish3_summary.md` ending `POLISH3-READY` or `POLISH3-BLOCKED`,
listing what moved to is:global and confirming the envelope is now a clean purple line.
