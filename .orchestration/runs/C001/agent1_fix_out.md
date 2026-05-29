Changed only `C001.py` and `test_C001.py`.

Override documented: C001 uses YAML `Delta_m_D=6.562e-15`, so budget is `3.281e-15`; this intentionally overrides the core’s stale `6.25e-15 / 2 = 3.125e-15`.

Corrected cross-check:
safe `pred=1.14472480657973238e-16`, ratio `3.48895094964868160e-02`, budget `3.28099999999999987e-15`.
excluded `pred=4.57889922631892903e-14`, ratio `1.39558037985947241e+01`, budget `3.28099999999999987e-15`.

Pytest: `test_C001.py -q` -> `8 passed`.
Pytest: `tests/constraints/ -q` -> `51 passed, 1 failed` on unrelated dirty-worktree `B003`.
Control: `tests/constraints/ -q --ignore=tests/constraints/primary/beauty/test_B003.py` -> `44 passed`.

`git diff --stat` is empty because these C001 files are untracked in this worktree; `--no-index` stat shows `C001.py` 231 insertions, `test_C001.py` 241 insertions.