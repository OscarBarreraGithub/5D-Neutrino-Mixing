# WQ NOTEBOOK — INDEPENDENT REVIEW (codex reviewer, gpt-5.x xhigh)

Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing (env: `source ~/.bashrc && conda activate ising_bootstrap`). An Opus agent built `notebooks/wq_quarkonly_explore.ipynb` (jupytext source `notebooks/wq_quarkonly_explore.py`; figs in `notebooks/figs/`) — clearer exploratory plots from the completed 1,000,000-point QUARK-ONLY scan in `scan_outputs/wq_quarkonly_1M_full/` (caches: `wq_quarkonly_cache.parquet` full 1M, `wq_quarkonly_fitted_subsample.parquet` 2100 fitted points). Do NOT trust the builder's report — RE-DERIVE from data + read the notebook code.

## VERIFY
1. **Executes clean:** re-run `jupyter nbconvert --to notebook --execute --inplace notebooks/wq_quarkonly_explore.ipynb` (or run the .py); confirm exit 0 and ZERO tracebacks. Report cell/fig counts.
2. **M_KK ≥ 4 TeV applied EVERYWHERE:** read the code — confirm every plot/aggregate filters to M_KK ≥ 4 TeV (i.e. the 5,7,10,15,20,30,50 TeV tiles; the 1/2/3 TeV tiles dropped) and axes floor at 4 TeV. Flag any plot that leaks sub-4-TeV data.
3. **No scale-mixing:** confirm every distribution/curve groups by EXACT (quark_fit_r, M_KK) or is faceted by one — never silently pooled across r or M_KK.
4. **Tags read per-row, not hardcoded:** confirm rigorous/proxy/partial come from each row's constraint `tag` field. Confirm `partial` INFO/SOFT advisories (e.g. T010 Z→bb with ratio≈52, B034, EW003, placeholder ratio=1 nonleptonic/EDM) are EXCLUDED from the carving/survival accounting and isolated into an "advisory-only, does NOT carve M_KK" panel — they must NOT be presented as real M_KK bounds.
5. **Fitted-Yukawa recompute is correct:** the notebook recomputes the full 3×3 fitted Yukawa by re-running `quarkConstraints.fit.fit_quark_sector` from the stored seed+r+M_KK. INDEPENDENTLY spot-check ≥3 points: re-run the fit and confirm the recomputed fitted singular values + bulk_c_Q match the values stored in the JSONL rows (to tight tolerance). Confirm Plot A labels clearly whether it shows fitted-physical vs input-anarchic Yukawas.
6. **Headline numbers match the raw data:** independently recompute from the parquet/JSONL and confirm the notebook's claims: (a) above 4 TeV, EW001 (oblique) is the strongest carver and B012 next; (b) binary veto: EW001 ~100% at 5 TeV, ~0 / nothing vetoes at ≥7 TeV; (c) strict survival = 1.0 at ≥5 TeV, inclusive survival goes 0→1 between 5 and 7 TeV; (d) c_Q descends gen1→gen3 crossing 0.5 (RS hierarchy); (e) the fitted j=3 (top) Yukawa column decreases with r. Flag any claim not supported by the data.
7. **Quark-only labeling:** confirm a cell lists the 37 allowlisted constraints by physics name (rigorous/proxy/partial) and clearly states lepton sector dropped / not fully rigorous.

## OUTPUT (≤14 lines)
nbconvert exit + fig count; the 4-TeV-floor + no-scale-mixing + per-row-tag confirmations; your independent fitted-recompute spot-check (values vs stored); your independent recomputation of ≥3 headline numbers vs the notebook; the advisory-not-a-bound isolation confirm; any discrepancy/fix. END with EXACTLY ONE line:
`WQ-NB-REVIEW: APPROVE`  or  `WQ-NB-REVIEW: NEEDS-FIXES`
