No remaining blockers found. I verified the revised plan against the real repo and the required fixes are now represented:

- Grid/seed reuse is byte-exact: plan uses unchanged `scripts/wq_quarkonly_1m_plan.py` constants and formulas, matching baseline `scan_plan.json`.
- Row pairing now normalizes raw `seed -> draw_seed`, `params.M_KK / 1000 -> mkk_tev`, and top-level `quark_fit_r`.
- Minimal hash protection is explicit: existing PR1 hash `45e21a07585f7489` stays pinned, minimal `ew_model` must be popped from `_config_payload`, and canonical WQ quark-only hash `c6939cc65d71f86a` is added.
- Both `_evaluate_draw` paths are called out for tests: quark-only and full-catalog.
- Worker round-trip `_config_from_payload(_config_payload(cfg))` is explicitly tested for both EW modes.
- Spectrum/model co-edit covers all real harness build sites, including universal-c sanity, and the mismatch `ValueError` is tested.
- Universal-c sanity policy is decided: it runs under `cfg.ew_model`.
- `paired_vetoes.parquet` enum/nullability now includes `not_evaluated` coherently.
- W7/W8 isolation now matches the real overlap and sequences W9 after both.

I also verified the current repo hashes directly: default/full hash is `45e21a07585f7489`, WQ `r=0.05` shard-00 hash is `c6939cc65d71f86a`, and WQ `r=0.25` shard-00 hash is `5de32abb4199300e`.

VERDICT: APPROVE