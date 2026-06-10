Audit completed from raw JSONL, no files changed.

Raw stream totals: 1,000,000 rows per model; evaluated rows were `minimal=878,707`, `custodial=878,705`; every cell had 20,000 raw draws.

Checks that passed:
- Veto fractions matched the JSON arrays after 6-decimal rounding for all included constraints; spot checks included 5 cells across both models, e.g. `custodial r=0.1 M_KK=1 B003 = 11653/18249 = 0.638555537`, JSON `0.638556`.
- Constraint inclusion matched exactly: 46 constraints seen, 17 had `veto_fraction > 0.5` somewhere, JSON listed exactly those 17. No missing biting constraints, no extra non-biting constraints.
- Yukawa percentiles matched for `custodial r=0.25 M_KK=3`, `n=19989`, up/down p25/p50/p75 matching JSON to 6 significant digits.
- `bare_floor_tev` matched raw non-skipped minima: all r values in both models floor to `1 TeV`.
- Representative raw rows were found and all have `skipped=False`, `survives_all_HARD_strict=True`.

Issues:
- Shipped `rep_matrices` in [scan_explorer.json](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog/website/src/content/scan_explorer.json:5958) fail the literal `SVD(|Y_u|), SVD(|Y_d|) == stored singular values` check. Example: low-r `Yu_abs` SVD is `[0.091121, 0.290732, 4.945959]`, stored is `[0.436056, 1.02424, 4.82886]`, max diff `0.733508`. Raw complex reconstruction from seed does match raw stored singular values exactly, so the problem is the shipped rounded absolute matrices are not faithful under the design’s stated `SVD(|Y|)` check.
- The builder streams raw JSONL rows at [build_scan_explorer.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog/website/scripts/build_scan_explorer.py:274), but exact Yukawa percentiles store all singular values in arrays via `SingularCell` at [build_scan_explorer.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog/website/scripts/build_scan_explorer.py:67). That is compact and avoids loading raw rows, but it is not strictly bounded-memory with respect to row count.

VERDICT: NEEDS-FIXES - rep_matrices are not SVD-faithful as shipped, and percentile aggregation is compact streaming but not strictly bounded-memory.