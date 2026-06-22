# Heads-up: take a look at the oblique S,T,U fit anchors (possibly stale)

**Found 2026-06-17 (reference review + write-up-vs-code audit). Not yet fixed — flagging only.**

The active oblique-parameter fit anchors are **labelled "PDG 2025" but hold pre-2022
central values**, so the S,T,U cut is centred on the wrong ellipse (S even has the
wrong sign). Worth a look before relying on the oblique reach.

| | repo (stale) | current PDG 2024 RPP, U fixed (Eq. 10.99) |
|---|---|---|
| S | 0.026 +/- 0.075 | **-0.05 +/- 0.07** |
| T | 0.047 +/- 0.066 | **0.00 +/- 0.06** |
| rho(S,T) | 0.90 | **0.93** |
| chi^2 thresh (2 dof, 95%) | 5.991 | 5.991 (fine) |

Where it lives (would need updating in lockstep):
- Code source of truth: `flavor_catalog/processes/top_higgs_ew/EW001.yaml`, active
  anchors ~lines 86-118 (it already also carries HEPfit/de Blas 2022 + Gfitter
  alternatives that could be promoted instead of hand-editing).
- Write-ups quoting the same triple: `review_local/constraint_formulas.tex` (local);
  check `docs/quark_scan_methodology_note.tex` for any S/T statement too.

Bounded impact: oblique is a proxy (does NOT set the rigorous floor; Z->bbar at
~25-30 TeV dominates). The prediction-side formulas and constants
(dT = piL/2c_W^2 v^2/M_KK^2, dS = 30 v^2/M_KK^2, L=35, v=246.21965, sin^2thetaW=0.23122)
all check out; only the data-side fit anchors are stale.
