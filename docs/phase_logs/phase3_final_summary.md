# Phase 3 Final Summary

All 10 ready-to-submit checklist items PASS for `quarkscan-paper-rc1`. The release gates completed on `paper/quark-scan-2026q2`: targeted quark-fit tests returned `14 passed`, the full suite returned `543 passed, 1 skipped, 0 xfails`, the methodology PDF rebuilt to 19 pages, the compact assumptions PDF rebuilt to 7 pages, artifact checksums spot-checked `OK`, publication figures are pruned to 9 PDFs, and the final acceptance log is committed at `1faba236ffa76df2a92784f6be199a7805b68cdf`.

## Headline numbers

| Quantity | Value | Convention |
|---|---|---|
| M_KK^min p50 (perturbative g_s) | 16.54 TeV | g_s ≈ 1.05 |
| M_KK^min p50 (g_s*=3 headline) | 47.26 +69.4 −25.0 TeV | g_s*=3, BGS 2020 + LO + factor-3 |
| M_KK^min p95 (g_s*=3) | 127.13 TeV | same |
| CFW reconciliation | factor 2.2 stronger at matched g_s*=3 | -- |
| Zero-pass UL (moreUV/moreIR) | p ≤ 2.3×10⁻⁶ | 95% Wilson |
| Zero-pass UL (Run C factor-1.5) | p ≤ 9.2×10⁻⁷ | 95% Wilson |
| ε_K binding fraction | 99.5% | post-audit |
| Test suite | 543 passed, 1 skipped, 0 xfails | full pytest |

## Branch and tag state

- RC1 tag commit count: `154` commits at `quarkscan-paper-rc1`.
- Paper branch final commit count after this summary commit: `155`.
- Tag confirmation: `git ls-remote --tags origin quarkscan-paper-rc1` returned tag object `f4d4fa837d9121c34dc1ca00d194fd89d160d0eb`.
- Peeled tag confirmation: `git ls-remote --tags origin 'quarkscan-paper-rc1^{}'` returned `1faba236ffa76df2a92784f6be199a7805b68cdf`.
- Branch push target: `origin/paper/quark-scan-2026q2`.

Phase 3 complete. `quarkscan-paper-rc1` is the release candidate.
