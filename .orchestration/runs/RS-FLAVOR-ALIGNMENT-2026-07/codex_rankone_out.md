Implemented the rank-one/U(2) lane in [scripts/rankone_u2_lane.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/scripts/rankone_u2_lane.py), with the report written to [rankone_u2_lane_report.md](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/.orchestration/runs/RS-FLAVOR-ALIGNMENT-2026-07/rankone_u2_lane_report.md). I only touched the new target script/report; unrelated untracked files were left alone.

Validation run:
`python -m py_compile scripts/rankone_u2_lane.py`
and the requested conda/env run with `--n-draws 1000 --mkk-tev 2,3,5`.

Key results:
- Verdict: `PASS`
- U(2) PDG pass: `1000/1000` at each M_KK
- eps_K pass fraction: `51.2%` at 2 TeV, `93.5%` at 3 TeV, `99.4%` at 5 TeV
- Flat baseline eps_K pass: `0.0%`, `1.7%`, `6.1%`
- Median `|C4|` U2/flat: `0.002`, `0.00186`, `0.00229`
- Median `|sin Phi12|`: `0.722`, `0.692`, `0.716`, so this is not phase tuning
- Max `c_Q1-c_Q2` and `c_d1-c_d2`: `0.000e+00`, so the suppression is structural and the c-fit did not reintroduce light-family splitting
- Perturbativity reject fraction: `0.0%` in this default run

The generated numbers support the Codex-E2 prediction in this lane: the LR operator is suppressed at the median by more than two orders of magnitude relative to flat anarchy, and the eps_K floor reaches the 2-3 TeV range without forcing a small phase.