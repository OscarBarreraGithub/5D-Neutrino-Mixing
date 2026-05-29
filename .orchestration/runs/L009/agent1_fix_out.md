Implemented L009 flavor pin in `lfv_three_body_tau.py`: matrix/contact inputs are normalized to a pinned `tau -> mu` proxy before the shared core call, and `*_emu*` spurion aliases now raise instead of being reinterpreted.

Pure dipole now works without `m_kk`, using a zero-contact 3000 GeV reference placeholder.

Corrected matrix regression: tau-mu overlap `(0.003+0.004j)` gives `BR = 1.6896991450212318e-11`; the mu-e default `(0.09+0.01j)` is not used.

Tests: L009 `14 passed`; full `python -m pytest tests/constraints/ -q` `570 passed`.