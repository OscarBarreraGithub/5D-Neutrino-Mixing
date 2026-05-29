1. NIT `flavor_catalog_constraints/primary/beauty/B019.py:383`: temporary monkey-patch of `anchor_scaffold.load_pdg_block` is global import-time state. Fix: add a list-entry-aware anchor helper or injectable loader path later; current serial behavior passes.

Isolation: passed. B019 imports no other constraint, reaches physics only through `physics_adapters.rare_b_kstar_dilepton`, and B019-scoped code/test/core/adapter files are new; no tracked B019 diff modifies existing function bodies.

Contract/anchor: passed. Missing couplings returns non-crashing HARD result with `predicted=None`, `ratio=None`; numeric result fields are real floats; complex Wilsons stay in diagnostics. Missing anchor and mismatched `load_anchor` probes raise `AnchorError`.

Numerical cross-check: anchor exp=`1.028`, SM=`1.0`, budget=`0.112`, q2=`[1.1,6.0]`. SM manual BR=`3.676601751950680e-08`, constraint numerator=`3.676601736716339e-08`, R=`1.0`, ratio=`0.25`, pass=`True`.

NP independent check: manual BR=`3.567149411842634e-08`, constraint numerator=`3.567149397061819e-08`, R=`0.970230025580014`, ratio=`0.515803343035585`, pass=`True`.

Safe/excluded: left=`1e-2` gives R=`0.977059271664665`, ratio=`0.454827931565487`, pass=`True`; left=`5e-2` gives R=`0.890140625156602`, ratio=`1.23088727538748`, pass=`False`.

Determinism: same input produced equal `ConstraintResult`; coupling arrays unchanged. Registry smoke: `B019 beauty HARD R_K* central-q2`.

Pytest: `python -m pytest tests/constraints/ -q` -> `463 passed in 14.13s`.

CODE-OK