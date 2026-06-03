# RETROACTIVE CODE+PHYSICS REVIEW (codex, gpt-5.x xhigh). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. Do NOT modify code — produce findings + verdict. This reviews a CORE physics edit that was made OUTSIDE the per-constraint review gate: commit c6c949c fixed an m_μ⁵ normalization bug in `quarkConstraints/mu_e_conversion.py` (the μ→e conversion rate). It was "Opus re-verified" only — never codex-reviewed. It is load-bearing for L003/L004/L005.

REVIEW:
- quarkConstraints/mu_e_conversion.py — the μ→e conversion rate / overlap integrals; specifically the m_μ⁵ (muon-mass-to-the-fifth) normalization. Confirm the dimensional analysis: the conversion RATE ∝ G_F²·(relevant coupling)²·m_μ⁵ (× nuclear overlap factors), and the branching ratio normalizes by the capture rate. Verify the code's powers of m_μ are dimensionally and physically correct (the bug was treating KKO overlaps as dimensionless, off by ~m_μ⁵ ≈ 7.59e4 in some unit). 
- The consumers: flavor_catalog_constraints/primary/charged_lepton/L003.py, L004.py, L005.py (and tests) — confirm their pinned numbers are consistent with the fixed core and physically sane (BR(μ→e conv in nuclei) limits, e.g. SINDRUM-II Au ~7e-13, future Mu2e/COMET).
- Cross-check against a literature normalization (Kitano-Koike-Okada overlap integrals) at the level of dimensional consistency and order-of-magnitude.

CHECK (with evidence, run things):
1. DIMENSIONAL CORRECTNESS of the m_μ⁵ factor and all unit conversions (GeV vs natural units, ħ for rate→width).
2. Recompute the conversion rate / BR for a sample coupling INDEPENDENTLY (from the core formula, not by re-calling the same wrapper) and confirm it matches the code and the L003/L004/L005 pinned values.
3. Are the three consumers consistent with each other and the core (same normalization, correct target-nucleus overlap factors)?
4. Numeric fields real floats; deterministic; no NaN/Inf at sample points.
5. Run `python -m pytest tests/constraints/primary/charged_lepton/ -q` and `python -m pytest tests/constraints/ -q`; report counts.

OUTPUT (stdout, <=20 lines): numbered findings (BLOCKER/SHOULD-FIX/NIT) with file:line, the ACTUAL recomputed rate/BR vs code + the literature order-of-magnitude check, pytest counts. End with: MUE-OK or MUE-NEEDS-FIXES.
