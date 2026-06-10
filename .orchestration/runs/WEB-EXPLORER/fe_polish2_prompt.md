# Scan Explorer front-end — POLISH ITERATION 2 (codex)

Render inspection of the post-JS DOM confirms: Plot 1 draws all 17 active constraint lines at low
opacity by default; in the 2–30 TeV band ~8 of them fan/overlap into a muddy near-black mass that
buries the purple envelope. The fix is to change the DEFAULT presentation. File:
`flavor_catalog/website/src/pages/explore.astro`. Dependency-free, theme-correct. `npm run build`
must pass (0 errors). Do NOT touch the JSON/data script. Do NOT commit/push.

## 1. Plot 1 default = ENVELOPE ONLY (the decisive change)
- By DEFAULT, draw ONLY: axes + log grid + the bold ENVELOPE line (var(--color-accent), ~3.4px) with
  its markers and M_KK value labels. Do NOT draw the per-constraint floor lines by default. The
  envelope already responds to the constraint checkboxes (toggling a constraint in/out raises/lowers
  it) — that IS the core "click constraints on/off" interaction, and it must read crystal-clear.
- The envelope must be unmistakably visible now (purple line, white plot area). Verify it is not
  hidden by anything.

## 2. Reveal ONE constraint on legend hover/focus (clean exploration)
- Hovering OR keyboard-focusing a legend row draws THAT single constraint's floor line on top of the
  envelope: its own color, full opacity, ~2.6px, small markers at each r where it bites. Mouse-leave
  / blur removes it. (Keep it to one line at a time — no stacking.)
- Optionally show a tiny inline readout of that constraint's floor values per r while hovered.

## 3. Opt-in "Show all traces" (preserve the superimposed view, default OFF)
- Add a small checkbox/toggle near Plot 1 labeled "Show all constraint traces", default OFF. When ON,
  render all active constraint lines (thin ~1.1px, opacity ~0.22, their colors) UNDER the envelope —
  the power-user superimposed view. When OFF (default), only the envelope (+ hover reveal) shows.

## 4. Legend row clarity
Each legend row must clearly read LEFT-ALIGNED as: [visible color swatch] [mono ID chip] [KaTeX
description], with obvious separation. Make the swatch a ~10x10px filled rounded square (or a ~14px
thick line sample) in the constraint color — the current 0.22rem hairline is invisible. Ensure the
ID chip has a visible subtle background/border and there's a clear gap before the description.
Left-align the row content (it currently looks centered/concatenated).

## 5. Keep intact
Plot 2 (per-r Yukawa boxes + whiskers), the 3 representative matrices + the magnitude/singular-value
note, dark+light theming. Re-verify all render.

## After
Run `npm run build` (0 errors). Write `.orchestration/runs/WEB-EXPLORER/fe_polish2_summary.md`
ending `POLISH2-READY` or `POLISH2-BLOCKED`.
