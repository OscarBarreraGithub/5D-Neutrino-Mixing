1. BLOCKER: none found.
2. SHOULD-FIX: none found.
3. NIT: none found.
4. Isolation: OK. K021 imports no other constraint; adapter AST check vs `HEAD` found changed_existing_functions=`[]`, removed_functions=`[]`, added neutral helpers only.
5. Contract: OK. Numeric result fields are real floats; complex values stay in diagnostics; missing inputs return non-crashing unevaluated result.
6. Anchor: OK. Uses YAML/load_anchor; missing value-id and mismatched limit-type probes both raised `AnchorError`.
7. Numerical cross-check: predicted=`1.21931398124349708e-11`, independent=`1.21931398124349708e-11`, rel_diff=`0.0`, ratio=`0.16043605016361803`.
8. Safe/excluded: safe predicted=`9.29077450503687135e-12`, ratio=`0.12224703296101147`, passes=True; excluded predicted=`2.32269362625921776e-10`, ratio=`3.05617582402528676`, passes=False.
9. Tests: `test_K021.py`: 11 passed in 13.33s; full `python -m pytest tests/constraints/ -q`: 1033 passed in 33.82s. Registry smoke: `K021 SECONDARY kaon HARD`.

CODE-OK