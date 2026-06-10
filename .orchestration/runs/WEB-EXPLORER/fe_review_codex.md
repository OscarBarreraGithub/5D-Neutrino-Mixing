**Findings**

- Medium: Plot 1 does not show the faint per-active-constraint traces by default. Part B requires a faint thin colored line per active constraint with a legend; the page initializes `showAllConstraintTraces = false` and only draws those traces when the extra checkbox is enabled. Runtime check confirmed `constraint-path-all` count is `0` on initial load with all 17 constraints active. See [explore.astro](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog/website/src/pages/explore.astro:203), [explore.astro](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog/website/src/pages/explore.astro:429), compared to [DESIGN.md](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/.orchestration/runs/WEB-EXPLORER/DESIGN.md:91).

**Verified**

`cd flavor_catalog/website && npm run build` completed with 0 errors.

Envelope logic matches the JSON/spec: all 17 on gives minimal `[20, 30, 30, 30, 30]` for `r=[0.05,0.1,0.25,0.5,1.0]`; custodial gives `[7, 7, 7, 7, 7]`; all-off falls to bare floor `[1,1,1,1,1]`; toggling a single constraint recomputed immediately.

Plot 2 reads `yukawa[model][r][envelope_mkk]`: live mini-chart geometry matched JSON p25/p50/p75 with max error `0` for checked states, and `n` readouts matched.

SVG root-cause fix looks good: dynamic SVG selectors are in `<style is:global>` and id-prefixed or inline; runtime check found no visible dynamic SVG nodes with default black fill. Legend rows are keyboard-focusable buttons. Representative matrices show `|Y_u|`/`|Y_d|` with the complex-singular-value note. No page console/runtime errors in light desktop or dark mobile checks.

VERDICT: NEEDS-FIXES - show all active constraint traces by default or otherwise satisfy Part B’s per-active-constraint trace requirement.