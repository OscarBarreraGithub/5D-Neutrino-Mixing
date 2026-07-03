# Rank-one/U(2) RS epsilon_K lane

Verdict: **PASS**

Configuration:
- draws per lane per M_KK: `1000`
- M_KK grid [TeV]: `2.0,3.0,5.0`
- xi_KK: `1.0`
- fixed U(2) doublet profiles: `c_Q1=c_Q2=0.5`, `c_Q3=0.25`
- right-down third-family leakage: `1.0` times CKM-sized angles
- profile-dressed left spurion boosts: `V23 x 4.0`, `V13 x 4.0`, `CP x 20.0`
- perturbativity: max |Y_ij| <= 3.0 with `reject` mode

Summary table, using PDG-passing draws for medians and pass fractions:

| M_KK [TeV] | U2 PDG pass | U2 eps_K pass | flat eps_K pass | med \|C4\| U2/flat | med \|GL12\| U2/flat | med \|GR12\| U2/flat | med \|sin Phi12\| U2 | c splits |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | 1000/1000 (100.0%) | 51.2% | 0.0% | 0.002 | 0.532 | 0.00416 | 0.722 | Q=0.000e+00, d=0.000e+00 |
| 3 | 1000/1000 (100.0%) | 93.5% | 1.7% | 0.00186 | 0.656 | 0.00352 | 0.692 | Q=0.000e+00, d=0.000e+00 |
| 5 | 1000/1000 (100.0%) | 99.4% | 6.1% | 0.00229 | 0.729 | 0.00383 | 0.716 | Q=0.000e+00, d=0.000e+00 |

Detailed medians:

| M_KK [TeV] | med \|C4\| U2 | med \|C4\| flat | med \|GL12\| U2 | med \|GL12\| flat | med \|GR12\| U2 | med \|GR12\| flat |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 6.692e-19 | 3.339e-16 | 8.255e-05 | 1.551e-04 | 3.254e-08 | 7.830e-06 |
| 3 | 3.049e-19 | 1.636e-16 | 8.192e-05 | 1.249e-04 | 3.376e-08 | 9.587e-06 |
| 5 | 1.074e-19 | 4.686e-17 | 8.116e-05 | 1.114e-04 | 3.306e-08 | 8.623e-06 |

Diagnostics:
- M_KK=2.0 TeV: rank-one perturbative reject fraction 0.0%; U2 all-DeltaF2 pass fraction 51.2% vs flat 0.0%; U2 fit/eval rejects 0/0, flat fit/eval rejects 0/0.
- M_KK=3.0 TeV: rank-one perturbative reject fraction 0.0%; U2 all-DeltaF2 pass fraction 93.5% vs flat 1.7%; U2 fit/eval rejects 0/0, flat fit/eval rejects 0/0.
- M_KK=5.0 TeV: rank-one perturbative reject fraction 0.0%; U2 all-DeltaF2 pass fraction 99.4% vs flat 6.1%; U2 fit/eval rejects 0/0, flat fit/eval rejects 0/0.

Interpretation:
- Median |C4| suppression >= 1e-2 at all measured M_KK: **True**.
- eps_K pass fraction exceeds 50% at 2-3 TeV in the U(2) lane: **True**.
- Phase tuning diagnostic, median |sin Phi12| among low-M_KK eps_K-pass draws < 0.1: **False**.
- Exact U(2) c-degeneracy retained in the fit: **True**.

The structural check is the last line: if `c_Q1-c_Q2` or `c_d1-c_d2` is nonzero, the lane has failed by construction because the c-fit reintroduced the anarchic source of the 1-2 coupling.
