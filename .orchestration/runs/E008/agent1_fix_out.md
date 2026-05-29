Changed E008 to virtualize dotted nested anchors through scaffold `load_anchor`, with a loud mismatch guard; honest INFO non-vetoing stub behavior remains unchanged.

Tests: `test_E008.py -q` -> 9 passed.  
Tests: `tests/constraints/ -q` -> 812 passed in 16.29s.