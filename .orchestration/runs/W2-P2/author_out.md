Plan: read design/derivation, add standalone kernel, add independent tests, run targeted + full pytest.

Implemented [quarkConstraints/rs_ew_spectrum.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/rs_ew_spectrum.py:1) and tests in [tests/test_rs_ew_spectrum.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/test_rs_ew_spectrum.py:1). Added lazy export for `RSEWSpectrum`.

Public API: `RSEWSpectrum.build(...)`, `spectrum.kk_ew_mass_gev` / `kk_ew_mass`, `spectrum.a(c, a_ref=...)`, module `a(...)`, `kk_ew_mass(...)`, `g0`, `w0`.

Asserted numbers at `Lambda_IR=3000 GeV`, `epsilon=1e-15`: exact root `x1=2.450509663813736`, `kk_ew_mass_gev=7351.528991441208`.

Overlap samples: `a(0.2)=21.91654523590418` converges at `N=32`, rel `9.388952477967849e-05`; `a(0.65)=-1.4760974560250766`.

Endpoint identity max residual: `2.220446049250313e-15`. Universal subtraction check: `a(c)-a_ref=0.0`. Sign check: `a(0.2) > a(0.65)`.

Tests: `tests/test_rs_ew_spectrum.py` -> `6 passed`; full `python -m pytest tests/ -q` -> `1644 passed, 1 skipped`.

PHASE2-DONE.