1. No BLOCKER/SHOULD-FIX findings for B005 isolation, contract, anchors, determinism, or numerical behavior.
2. NIT `tests/constraints/primary/beauty/test_B005.py:131`: hardcoded budget duplicates the YAML-derived check just above; either remove it or comment that it is an intentional regression pin.

Isolation: `B005.py` imports physics only via `physics_adapters.rare_b_meson`; no tracked deltaf2 adapter/core edits; B005 additions are new files only.

Numerics: anchor exp `3.34e-09 +/- 2.7e-10`, SM anchor `3.64e-09 +/- 1.2e-10`, budget `5.95465734053883e-10`.  
SM constraint/manual/core: `3.6483515801441554e-09`, diff `0.0`.  
Safe point `left=1e-2`: pred `3.4867726848419893e-09`, NP shift `-1.615788953021661e-10`, ratio `0.27134877132584934`, pass `True`.  
Excluded point `left=-5e-1`: pred `1.639488920831644e-08`, NP shift `1.2746537628172286e-08`, ratio `21.405996851229173`, pass `False`.

Anchor probes: missing candidate raised `AnchorError`; unit mismatch raised `AnchorError`.  
Missing couplings: non-crashing pass result, `predicted=None`, `ratio=None`.  
Determinism: repeated evaluation equal; coupling arrays not mutated.  
Registry smoke: `B005 beauty HARD BR(B_s -> mu+ mu-)`.  
Pytest: `tests/constraints/primary/beauty/test_B005.py -q` = `10 passed`; `python -m pytest tests/constraints/ -q` = `157 passed`.

CODE-OK