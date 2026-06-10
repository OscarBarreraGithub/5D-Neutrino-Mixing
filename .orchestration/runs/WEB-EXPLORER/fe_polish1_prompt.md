# Scan Explorer front-end — POLISH ITERATION 1 (codex)

The page (`flavor_catalog/website/src/pages/explore.astro`) builds and renders; data is faithful
(dual-reviewed). Fix the following readability issues found by render inspection (screenshots).
Do NOT change `scan_explorer.json` or the data script. Keep dependency-free, theme-correct. Re-run
`npm run build` (0 errors) when done.

## 1. Plot 1 "Minimum M_KK vs r" is a muddy near-black mass with many constraints active (MAIN FIX)
With 17 active constraints, the per-constraint lines (opacity 0.42, 1.8px) stack into a near-black
blob and the 4.4px pure-black envelope adds to it. Redesign for clarity:
- **Envelope**: keep it the clear hero line but use the site ACCENT color (a CSS var, e.g.
  `--color-accent`/link color), width ~3.2, with the existing M_KK point labels + markers. It must
  read clearly against the constraint lines.
- **Per-constraint lines**: thinner (~1.1px) and lower opacity (~0.28), distinct hues. To kill the
  muddy mass, only draw each constraint's line where it actually bites — i.e. draw the segment only
  where floor_c(r) > bare_floor (skip/flatten the baseline-at-1-TeV portions so they don't add ink),
  OR draw per-r markers at each constraint's floor instead of full baseline-spanning lines. Pick
  whichever is cleaner; the goal is the envelope + the few binding constraints stand out, not a blob.
- **Legend-hover interactivity**: hovering or focusing a legend row highlights THAT constraint's line
  (opacity 1, width ~2.4) and dims the others (opacity ~0.12). Mouse-leave restores. This makes 17
  constraints explorable without clutter. Keyboard-focusable for accessibility.
- Y-axis ticks: ensure the ">50" and "50" labels don't collide; clean log ticks (1,2,3,5,10,20,50).
- Keep the headline caption.

## 2. Legend rows concatenate the constraint ID and description (e.g. "EW001S, T, U")
Render each legend row as: [color swatch] · [mono ID chip, e.g. `EW001`] · [KaTeX description], with
clear spacing between the ID and the description. Compact; vertically scroll if it overflows the plot
height. The swatch color must match that constraint's line color.

## 3. Representative matrices — labeling honesty (from data review)
Add a short note under the 3 matrices: the matrices show entry magnitudes |Y_ij| of the complex
anarchic spurion, and s_u / s_d are the singular values of the COMPLEX matrix Y (so SVD of the
displayed |Y| will not equal them). Keep the clean KaTeX layout. (No data change.)

## 4. Plot 2 whiskers
Confirm the p25–p75 whiskers actually render on each singular-value bar; if not visible, make them
visible (thin tick at p25 and p75 around the p50 bar). Keep boxes readable.

## After
Run `npm run build` in `flavor_catalog/website` (0 errors). Write
`.orchestration/runs/WEB-EXPLORER/fe_polish1_summary.md` (what changed) ending `POLISH1-READY` or
`POLISH1-BLOCKED`. Do NOT commit/push.
