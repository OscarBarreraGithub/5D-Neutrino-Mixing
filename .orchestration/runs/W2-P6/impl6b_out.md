Implemented 6b code and tests, but I cannot honestly mark the dual gate complete: Codex reviewer approved; Opus review could not be run because this session has no `claude`/Opus CLI, Anthropic package/API key, Claude OAuth token, or installable Opus connector.

New API: `quarkConstraints/rs_higgs_yukawas.py` with frozen `RSHiggsYukawaCouplings`, builder wiring via `include_higgs_yukawas=True`, and `rs_higgs_yukawas` in `KNOWN_EXTRA_KEYS`.

T018/T019/T020 now require `rs_higgs_yukawas`, use the Higgs-LFV adapter/core formula, keep SM=0, and degrade absent input as non-vetoing `evaluated=False`.

Gate checks: diagonal-v1 offdiag `Y_h=0.0`; T018/T019/T020 predict `0`; LFV-live toy has `Y_mu_tau=-2.697e-05+4.163e-07j`, `BR=9.395e-07`, and 6 TeV/3 TeV Yukawa scaling `0.25`.

Tests: collection `1683 -> 1693` (+10); focused `48 passed`; full suite `1692 passed, 1 skipped`.

Codex review: APPROVE. Opus review: BLOCKED by unavailable tooling.