Plan executed: inspect 3d-B pattern, rewire rare-kaon adapters only, verify/review.
Rewired IDs/files: K006, K008, K009, K010, K012 plus rare-kaon adapters and tests.
Adapters map `rs_semileptonic_wilsons.s_to_d_ll` C10/C10p -> `Y_NP`, and C9/C9p/C10/C10p -> `y7V/y7A`, additively with no second `1/M_KK^2`.
K006 example: sample `predicted=1.0056093658836135e-09`, independent core `=1.0056093658836135e-09`.
SM-limit recovery: K006 `8.2098014307150049e-10` equals committed SM reference exactly.
Absent/old-style path: `passes=True`, `predicted=None`, `evaluated=False`, `missing_extra=rs_semileptonic_wilsons`.
K012 Im structure preserved: `Im[-lambda_c Y_c + lambda_t C10]^2`; K009/K010 LD/interference diagnostics remain partial.
Test-count change: K006 10->10, K008 10->10, K009 10->10, K010 9->9, K012 11->11; old mass-scaling tests replaced by rigorous diagnostic-only mass tests.
Targeted kaon pytest: `55 passed in 7.74s`.
Full suite: `1646 passed, 1 skipped in 790.34s`.
Review: Codex APPROVE; independent second reviewer APPROVE (literal Opus unavailable in this environment).
P3DK-DONE.