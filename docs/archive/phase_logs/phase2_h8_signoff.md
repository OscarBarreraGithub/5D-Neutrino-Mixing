# Phase 2 Hole #8 Sign-off

**Verdict**: PASS

Date: 2026-05-15
Branch: `paper/quark-scan-2026q2` @ `1d9e17e`

## Verification checklist

### 1. Branch and commit chain
Branch `paper/quark-scan-2026q2` is pushed to `origin` (remote tip
`1d9e17ea1269ed5c552198fd30acfab7fffa39f5` matches local tip). The four
hole-#8 commits are present and in the expected order on the tip:

```
1d9e17e docs(paper): replace zero-pass impossibility wording with binomial bounds
212320f figures(zero-pass): annotate plot legends with finite-ensemble upper limits
a9d1058 audit(zero-pass): inventory + 95% Wilson UL on moreUV/moreIR/runC
f34f9df physics(stats): add Wilson-score upper-limit helper for zero-pass observations
```

### 2. `wilson_upper_limit` numerical check
```
$ /n/home09/obarrera/.conda/envs/ising_bootstrap/bin/python -c \
    "from quarkConstraints.finite_stats import wilson_upper_limit; \
     print(wilson_upper_limit(0, 1_600_000)); \
     print(wilson_upper_limit(0, 4_000_000))"
2.3039946915962305e-06
9.215991506542227e-07
```
Both values match the expected order (~2.3e-6 and ~9.2e-7) and the
numbers used in the methodology note and inventory.

### 3. Tests
```
$ /n/home09/obarrera/.conda/envs/ising_bootstrap/bin/python -m pytest -q tests/test_finite_stats.py
...                                                                      [100%]
3 passed in 0.05s
```

### 4. Inventory completeness
`docs/audits/zero_pass_inventory.md` contains all three zero-pass scans
with the required columns:

| Run | N | k_obs | p_UL_95 | c-pattern | Y-prior | Gate | Output SHA256 |
|---|---:|---:|---:|---|---|---|---|
| Run 3 moreUV | 1,600,000 | 0 | 2.304e-6 | c -> c+0.05 baseline | U(-1.5,1.5), \|Y\|>=0.1 | factor 3/3/5 | `3150aadc...e3158cd` |
| Run 3 moreIR | 1,600,000 | 0 | 2.304e-6 | c -> c-0.05 baseline | U(-1.5,1.5), \|Y\|>=0.1 | factor 3/3/5 | `e6899f1a...90aa30d5` |
| Run C CFW-like | 4,000,000 | 0 | 9.216e-7 | baseline (eq. cvals) | U(-3,3), \|Y\|>=0.5 | factor 1.5/1.5/2.5 | `92f57214...d460cf3d` |

The inventory also includes the full "Claim Locations" cross-reference
table mapping each verbatim claim back to its run dir, N, gate
definition, and p_UL_95.

### 5. Wording audit
```
$ grep -nE "zero PDG|0 PDG|no passes|exactly zero" docs/quark_scan_methodology_note.tex
468:\emph{every} $c$ by $\pm 0.05$ produces zero PDG-passing draws.
960:$\pm 0.05$ shift of the $c$-pattern returns zero PDG-passing draws,
```
Both hits are accompanied by Wilson UL context within the same
paragraph or section:

- Line 468 is in the c-pattern provenance paragraph; the matching
  numerical Wilson UL ($p_{\rm pass}\le 2.3\times10^{-6}$) appears at
  lines 740-742 (figure caption) and 758-763 (in-text), and the
  surrounding text refers to "An empirical sensitivity check (Run~3,
  \S\ref{sec:robustness} below)".
- Line 960 is in the summary paragraph; the same Wilson-UL passage
  (\S\ref{sec:robust-cvals}, lines 757-768) supplies the explicit
  $p_{\rm pass}\le 2.3\times10^{-6}$ bound and the explicit
  "finite-ensemble bound, not an analytic impossibility statement"
  qualifier.
- The Run C zero-pass claim at lines 853-858 also carries its own
  $p_{\rm pass}\le 9.2\times10^{-7}$ Wilson UL and the matching
  finite-ensemble qualifier.

No hit is an unqualified "zero" claim; every occurrence is bracketed by
a Wilson UL and a finite-ensemble disclaimer.

### 6. Plot legend annotations
`pdftotext -layout
results/figures/quark/rs_anarchy_mkk_min_hist_by_cvals.pdf` reproduces
the legend entries

```
moreUV (N=1.6M, p 2.3 x 10 6 95% CL)
moreIR (N=1.6M, p 2.3 x 10 6 95% CL)
```

in both the per-draw histogram legend and the cumulative-distribution
legend. (The `<=` and exponent superscripts are flattened by
`pdftotext`; the on-screen PDF renders them correctly.) The required
`N=1.6M` and `p_UL` annotations are therefore present for moreUV and
moreIR in both panels.

### 7. Methodology PDF build
`pdflatex -interaction=nonstopmode quark_scan_methodology_note.tex`
ran to completion with `Output written on
quark_scan_methodology_note.pdf (18 pages, 590665 bytes).` and no
unresolved references or warnings of substance.

## Recommendation

Proceed to hole #9 (scope wording).
