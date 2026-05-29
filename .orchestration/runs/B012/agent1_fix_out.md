Changed `B012.py` to raise `AnchorError` when neutral/charged anchor `block_key` is outside the expected candidate tuple; added the monkeypatch regression in `test_B012.py`.

`python -m pytest tests/constraints/primary/beauty/test_B012.py -q`: 12 passed  
`python -m pytest tests/constraints/ -q`: 194 passed