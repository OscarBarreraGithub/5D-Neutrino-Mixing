1. SHOULD-FIX: charged-B long-distance term is incorrectly scaled by NP. `quarkConstraints/rare_b_nunu.py:408-409` applies `BR = 5.58e-6 * R_K`, but HPQCD gives a separate `BR_LD = 6.09(53)e-7` and total `5.58(37)e-6`; correct proxy is `BR_LD + (BR_SM_total - BR_LD) * R_K`. For `R_K=4.002`, code gives `2.233e-5`; corrected is `2.050e-5`. See `flavor_catalog/references/B022/hpqcd2023_btoknunu_sm_arxiv.txt:34,67,70` and `B022.yaml:143,155`.

2. NIT: `X_t` default is close but not source-identical. Code computes `X_t=1.48123` from `m_t=163.5`, `eta_x=0.994` at `quarkConstraints/rare_b_nunu.py:75-78,226-230`; HPQCD snapshot quotes `X_t=1.469(17)` at `hpqcd2023_btoknunu_sm_arxiv.txt:17`. Difference is within `1σ`, but cite this convention or align it.

3. NIT: ΔF=2 amplitude/running checks are N/A for B022. This is not an `M12` observable; the rare-decay core correctly uses the full complex modulus squared `R_K=|C_L+C_R|^2/|C_L^SM|^2` at `quarkConstraints/rare_b_nunu.py:273-288`. No CP-imaginary vs magnitude confusion; no 2 GeV ΔF=2 running path is physically required for this semileptonic operator.

4. NIT: Budget anchors match YAML/snapshots and are uncertainty-aware: Belle II `2.30e-5`, exp `+7.071e-6/-6.403e-6`, HPQCD SM `5.58(0.37)e-6`, combined `+7.081e-6/-6.414e-6`; SM pull `2.716σ`, consistent with the YAML “2.7 sigma” note. Lines: `B022.yaml:111-120,137-142`, `flavor_catalog_constraints/primary/beauty/B022.py:233-255,280-330`.

PHYSICS-NEEDS-FIXES