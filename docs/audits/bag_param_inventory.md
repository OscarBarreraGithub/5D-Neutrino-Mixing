# Hadronic Bag-Parameter and Delta m/epsilon_K Input Inventory

Audit date: 2026-05-15  
Code target: `quarkConstraints/deltaf2.py` on branch `audit/bag-inputs`

The `Coded value` column records the pre-audit value in `deltaf2.py`.  The
relative drift is

```text
Delta = abs(coded - canonical) / abs(canonical)
```

Classification follows the Phase 2 rule: Green for `Delta < 5%` or exact
conventions, Yellow for `5% <= Delta < 10%`, and Red for `Delta >= 10%`.

| System | Symbol | Physical quantity | Coded value | Scale / scheme | Apparent source in code | Canonical source (year, eprint/page) | Canonical value | Delta | Class |
|---|---|---|---:|---|---|---|---:|---:|---|
| K | `F_K` | kaon decay constant | `0.1557` | GeV, physical | `PDG` comment | FLAG 2024, summary table, arXiv:2411.04268 | `0.1557` | `0` | Green |
| K | `M_K` | neutral kaon mass | `0.49761` | GeV, physical | `PDG` comment | PDG 2024, meson summary tables, Phys. Rev. D 110, 030001 | `0.497611` | `2.00960188e-06` | Green |
| K | `DELTA_M_K` | `K_L - K_S` mass difference | `3.484e-15` | GeV, physical | `PDG` comment | PDG 2024, `K_L-K_S = 3.484(6)e-12 MeV` | `3.484e-15` | `0` | Green |
| K | `M_S_2GEV` | strange-quark MS-bar mass | `0.0934` | GeV, MS-bar at 2 GeV | `FLAG` comment | FLAG 2024, summary table, arXiv:2411.04268 | `0.09346` | `0.000641985876` | Green |
| K | `M_D_2GEV` | down-quark MS-bar mass | `0.00467` | GeV, MS-bar at 2 GeV | `FLAG` comment | FLAG 2024, summary table, arXiv:2411.04268 | `0.00470` | `0.00638297872` | Green |
| K | `B_1_K` | SM kaon bag parameter `B_K` | `0.717` | Coded as 2 GeV matrix-element input; value matches RGI `hat B_K` | `FLAG average` comment | FLAG 2024, `B_K^MSbar(2 GeV)=0.5503(66)`, arXiv:2411.04268 | `0.5503` | `0.302925677` | Red |
| K | `B_4_K` | kaon LR bag parameter `B_4` | `0.78` | MS-bar, quoted lattice values at 3 GeV | `ETM 2013` comment | FLAG 2024 Nf=2+1 BSM average, `B_4=0.903(14)`, arXiv:2411.04268 | `0.903` | `0.136212625` | Red |
| K | `B_5_K` | kaon LR bag parameter `B_5` | `0.57` | MS-bar, quoted lattice values at 3 GeV | `ETM 2013` comment | FLAG 2024 Nf=2+1 BSM average, `B_5=0.691(14)`, arXiv:2411.04268 | `0.691` | `0.175108538` | Red |
| K | `KAPPA_EPSILON` | long-distance and phase correction to epsilon_K | `0.94` | dimensionless phenomenological factor | `Buras et al.` comment | Buras, Guadagnoli, Isidori 2010, arXiv:1002.3612 | `0.94` | `0` | Green |
| K | `EPSILON_K_EXP` | measured `|epsilon_K|` | `2.228e-3` | dimensionless physical observable | `PDG` comment | PDG 2024, neutral kaon summary | `0.002228` | `0` | Green |
| K | `EPSILON_K_SM` | SM central prediction for `|epsilon_K|` | `1.81e-3` | dimensionless theory central value | `CKMfitter central value` comment; numerically Brod-Gorbahn 2011 | Brod, Gorbahn, Stamou 2020, Eq. (20), arXiv:1911.06822 | `0.002161` | `0.162424803` | Red |
| B_d | `F_BD` | `B_d` decay constant | `0.1900` | GeV, physical | `FLAG 2024` comment | FLAG 2024 summary table, arXiv:2411.04268 | `0.1900` | `0` | Green |
| B_d | `M_BD` | `B_d` mass | `5.27972` | GeV, physical | `PDG` comment | PDG 2024, meson summary tables | `5.27972` | `0` | Green |
| B_d | `M_B_QUARK` | bottom-quark MS-bar mass | `4.18` | GeV, MS-bar at `m_b` | hardcoded comment | PDG 2024 quark summary, `m_b(m_b)=4.183 GeV` | `4.183` | `0.000717188621` | Green |
| B_d | `M_D_QUARK_BD` | down-quark MS-bar mass | `0.00467` | GeV, MS-bar at 2 GeV | hardcoded comment | FLAG 2024 summary table, arXiv:2411.04268 | `0.00470` | `0.00638297872` | Green |
| B_d | `B_1_BD` | `B_d` VLL bag parameter | `0.87` | MS-bar/NDR at `mu=m_b` | `FLAG 2024, renormalized` comment | HPQCD 2019 Table VI, arXiv:1907.01025 | `0.806` | `0.0794044665` | Yellow |
| B_d | `B_4_BD` | `B_d` LR bag parameter `B_4` | `1.02` | MS-bar/NDR at `mu=m_b` | `FLAG 2024` comment | HPQCD 2019 Table VI, arXiv:1907.01025 | `1.077` | `0.0529247911` | Yellow |
| B_d | `B_5_BD` | `B_d` LR bag parameter `B_5` | `0.96` | MS-bar/NDR at `mu=m_b` | `FLAG 2024` comment | HPQCD 2019 Table VI, arXiv:1907.01025 | `0.973` | `0.01336074` | Green |
| B_d | `DELTA_M_BD_EXP` | measured `Delta m_Bd` | `3.334e-13` | GeV, physical | `PDG` comment | PDG 2024, `0.5069e12 hbar s^-1` | `3.3364764095260997e-13` | `0.000742222999` | Green |
| B_d | `DELTA_M_BD_SM` | SM central prediction for `Delta m_Bd` | `3.6e-13` | GeV, theory central value | `CKMfitter` comment | HPQCD 2019 Eq. (17), `Delta M_d,SM = 0.555 ps^-1`, arXiv:1907.01025 | `3.653076360795e-13` | `0.0145292229` | Green |
| B_s | `F_BS` | `B_s` decay constant | `0.2303` | GeV, physical | `FLAG 2024` comment | FLAG 2024 summary table, arXiv:2411.04268 | `0.2303` | `0` | Green |
| B_s | `M_BS` | `B_s` mass | `5.36692` | GeV, physical | `PDG` comment | PDG 2024, meson summary tables | `5.36693` | `1.86326261e-06` | Green |
| B_s | `M_S_QUARK_BS` | strange-quark MS-bar mass | `0.0934` | GeV, MS-bar at 2 GeV | hardcoded comment | FLAG 2024 summary table, arXiv:2411.04268 | `0.09346` | `0.000641985876` | Green |
| B_s | `B_1_BS` | `B_s` VLL bag parameter | `0.87` | MS-bar/NDR at `mu=m_b` | `FLAG 2024` comment | HPQCD 2019 Table VI, arXiv:1907.01025 | `0.813` | `0.0701107011` | Yellow |
| B_s | `B_4_BS` | `B_s` LR bag parameter `B_4` | `1.02` | MS-bar/NDR at `mu=m_b` | `FLAG 2024` comment | HPQCD 2019 Table VI, arXiv:1907.01025 | `1.033` | `0.0125847047` | Green |
| B_s | `B_5_BS` | `B_s` LR bag parameter `B_5` | `0.96` | MS-bar/NDR at `mu=m_b` | `FLAG 2024` comment | HPQCD 2019 Table VI, arXiv:1907.01025 | `0.941` | `0.0201912859` | Green |
| B_s | `DELTA_M_BS_EXP` | measured `Delta m_Bs` | `1.1688e-11` | GeV, physical | `PDG` comment | PDG 2024, `17.765e12 hbar s^-1` | `1.16931354143285e-11` | `0.000439181977` | Green |
| B_s | `DELTA_M_BS_SM` | SM central prediction for `Delta m_Bs` | `1.17e-11` | GeV, theory central value | `CKMfitter` comment | HPQCD 2019 Eq. (17), `Delta M_s,SM = 17.59 ps^-1`, arXiv:1907.01025 | `1.1577948321871e-11` | `0.0105417363` | Green |
| D | `F_D` | `D` decay constant | `0.2120` | GeV, physical | `FLAG 2024` comment | FLAG 2024 summary table, arXiv:2411.04268 | `0.2120` | `0` | Green |
| D | `M_D0` | `D0` mass | `1.86484` | GeV, physical | `PDG` comment | PDG 2024, meson summary tables | `1.86484` | `0` | Green |
| D | `M_C_QUARK` | charm-quark MS-bar mass | `1.27` | GeV, MS-bar at `m_c` | hardcoded comment | PDG 2024 quark summary, `m_c(m_c)=1.2730 GeV` | `1.2730` | `0.00235663786` | Green |
| D | `M_U_QUARK` | up-quark MS-bar mass | `0.00216` | GeV, MS-bar at 2 GeV | hardcoded comment | PDG 2024 quark summary, `m_u(2 GeV)=2.16 MeV` | `0.00216` | `0` | Green |
| D | `B_1_D` | `D0` VLL bag parameter | `0.75` | MS-bar at 3 GeV in comparable lattice source | `less precise` comment | ETM 2015 Table 2, arXiv:1505.06639 | `0.757` | `0.00924702774` | Green |
| D | `B_4_D` | `D0` LR bag parameter `B_4` | `1.0` | MS-bar at 3 GeV in comparable lattice source | `estimated` comment | ETM 2015 Table 2, arXiv:1505.06639 | `0.91` | `0.0989010989` | Yellow |
| D | `B_5_D` | `D0` LR bag parameter `B_5` | `1.0` | MS-bar at 3 GeV in comparable lattice source | `estimated` comment | ETM 2015 Table 2, arXiv:1505.06639 | `0.97` | `0.0309278351` | Green |
| D | `DELTA_M_D_EXP` | measured `Delta m_D0` | `6.25e-15` | GeV, physical | `HFLAV` comment | PDG 2024/HFLAV average, `0.997e10 hbar s^-1` | `6.562373210293e-15` | `0.0476006469` | Green |

## Convention-fixed and out-of-scope constants

Post-B3 correction: the sentence that previously recorded the matrix-element
prefactors as inspected and unchanged was stale. The June-2026 B3 audit and the
2026-07 full-repo audit corrected the M12-ready GGMS forms used by
`_kaon_matrix_elements()` and `_meson_matrix_elements()`: `O1 = (1/3) f^2 m B1`,
`O4 = (r_chi/4 + 1/24) f^2 m B4`, and
`O5 = (r_chi/12 + 1/8) f^2 m B5`. The pre-B3 record used `2/3` for `O1` and
swapped and doubled the LR coefficients. This inventory remains an input
provenance document; the matrix-element correction is recorded here only to avoid
certifying the stale pre-B3 values.

The `DeltaF2Input` bounds and the legacy operator-weight bounds were also
inspected.  The hadronic path used by the scan computes budgets from the
physical constants above, while `bound`, `ll_weight`, `rr_weight`,
`lr1_weight`, `lr2_weight`, and `LEGACY_OPERATOR_WEIGHT_BOUNDS` encode the
repo-owned v1 gate/fallback convention rather than external hadronic inputs.

## Red-row amplitude propagation at M_KK = 3 TeV

For the red kaon bag rows, I used a one-operator estimate with
`Im C_i = 1 / (3 TeV)^2` and the current `epsilon_K` prefactor.  Amplitudes
below are dimensionless `|epsilon_K^NP|` estimates for a unit imaginary
coefficient in the named operator.

| Input | Operator estimate | Pre-audit amplitude | Canonical amplitude | Shift vs pre-audit |
|---|---|---:|---:|---:|
| `B_1_K` | `O1_VLL` | `122232.42308457832` | `93813.81091135768` | `0.23249651324965123` |
| `B_4_K` | `O4_LR` | `905734.5362109917` | `1048561.905382725` | `0.15769230769230766` |
| `B_5_K` | `O5_LR` | `1888476.57173552` | `2289363.703630253` | `0.21228070175438601` |

For `EPSILON_K_SM`, the allowed central-value NP budget changes from
`0.00041799999999999997` to `6.699999999999979e-05`, so all
`epsilon_K` ratios increase by a factor `6.238805970149273`.

Result: `INVALIDATION_GATE_TRIPPED` for scans whose final claims depend on the
pre-audit `epsilon_K` ratios, including RUNA.
