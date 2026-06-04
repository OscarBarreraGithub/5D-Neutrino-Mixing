Blocking findings: none.
Seed-disjointness: PASS; `max_offset = 1000003*9+1999 = 9002026 < 20000000`, and exhaustive 1,000,000-seed recomputation found 0 collisions; harness uses `np.random.default_rng(draw_seed)` per draw.
r/Yukawa pairing: keep unpaired for survival statistics; switch to paired only if the physicist wants fixed-input `c_Q` evolution across r.
Serialization: PASS; `quark_fit_r` is serialized in row/params/provenance, fitted u/d singular values are included, and full-mode config hash + row JSON were byte-identical to `HEAD`.
Analysis: PASS; smoke-1k read 4 JSONL/1000 rows, emitted all 11 PNGs and a non-empty quark-only report; old rows degrade via run-summary r and seed singular fallback.
SBATCH: PASS; `04:00:00`, 48 CPUs, 64G, `serial_requeue`/`randall_lab`, `0-49%50` are adequate; 48 CPUs is resource-heavy because only 10 tiles run per task.
Cosmetic: old smoke report is mislabeled W6b and shows skipped sanity badly; new analysis report is correctly labeled and does not present skipped sanity as failure.
Pytest: `1713 passed, 1 skipped in 776.08s`.
WQ-1M-REVIEW: APPROVE