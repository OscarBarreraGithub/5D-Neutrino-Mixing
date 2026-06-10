# Front-end polish iteration 2 summary

Updated `flavor_catalog/website/src/pages/explore.astro` for Plot 1:

- Default Plot 1 rendering is now envelope-only: axes, log grid, purple envelope line, envelope markers, and M_KK labels.
- Added an opt-in `Show all constraint traces` checkbox, default off. When enabled, all active constraint traces render under the envelope at low opacity.
- Legend hover and keyboard focus now reveal one active constraint trace at a time above the envelope, with full-opacity line and markers.
- Legend rows now use a visible color swatch, a bordered mono ID chip, and left-aligned KaTeX description text with clearer spacing.

Kept the embedded JSON/data script unchanged. Plot 2, representative matrices, and theming code were left intact.

Verification:

- `npm run build` completed successfully with 0 errors.

POLISH2-READY
