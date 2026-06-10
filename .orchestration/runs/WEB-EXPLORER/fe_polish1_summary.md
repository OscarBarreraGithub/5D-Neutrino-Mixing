Scan Explorer front-end polish iteration 1 completed.

Changes made:
- Redesigned Plot 1 constraint rendering in `flavor_catalog/website/src/pages/explore.astro`: the envelope now uses the site accent color at a thinner hero weight, the y ticks are the clean log set with a separated `>50` cap label, and per-constraint traces only draw at r points where that constraint raises the bare floor or caps the scan.
- Added low-opacity per-constraint binding markers and legend hover/focus highlighting. Legend rows are keyboard-focusable buttons that highlight the matching trace while dimming the others.
- Reworked legend row layout so the swatch, mono ID chip, and KaTeX description are visually separated and scroll within the plot height.
- Added the representative-matrix note clarifying that displayed entries are `|Y_ij|` magnitudes while `s_u` and `s_d` are singular values of the complex matrices.
- Made Plot 2 p25-p75 whiskers render above the median bars so the quartile ticks remain visible.

Verification:
- `npm run build` in `flavor_catalog/website` completed successfully with 0 errors.

POLISH1-READY
