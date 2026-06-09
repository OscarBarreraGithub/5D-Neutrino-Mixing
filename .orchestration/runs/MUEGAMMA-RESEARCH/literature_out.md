# mu->e gamma in anarchic RS — web literature survey (Opus, ae0c3df0421c8f03b)

## Bottom line
- **Dominant mechanism (all agree):** one-loop dipole from **KK-fermion + Higgs/Goldstone (or KK-gauge) exchange**, chirality flip via a KK-fermion mass insertion. KK-gluon/KK-Z loops are ~aligned with the SM Yukawa and subdominant; the KK-fermion+Higgs loop has a different flavor direction and dominates.
- **Flavor structure (all agree):** coefficient set by **anarchic Yukawa-cubed spurion ~Y Y† Y (= Y*^3)**, dressed by zero-mode overlaps f_L,f_E, warped by v/M_KK^2. In LMFV (Perez-Randall) the cubic is symmetry-removed and only **(Y_N Y_N†)_12** (two powers of the SMALL neutrino Yukawa) survives.
- **NDA vs explicit (the crux):** Perez-Randall C=0.02 is an NDA estimate of the brane-Higgs IR COUNTERTERM in their LMFV model. Warped Penguins (explicit 5D loop) c0 "matches the order of magnitude of the full 5D calculation" => **NDA agrees with explicit to O(1)**. Perez-Randall themselves: the NDA estimate is BIGGER than the explicit (log-enhanced) loop => **C=0.02 is CONSERVATIVE (safe)**. Explicit adds two things NDA can't: (i) alignment of the single-Yukawa term with the mass matrix (suppresses naive), (ii) a Y*-independent MISALIGNMENT spurion `b in (-0.03,0.03)` from bulk-mass-vs-Yukawa, which can enhance OR (tuned sign) cancel the anarchic `a` term. O(1)-O(few) modulations.

## Papers
- **Agashe-Blechman-Petriello hep-ph/0606021** (NB: prompt's 0606293 is wrong): one-loop dipole, two KK-fermion lines + internal Higgs (Eq.38-41), coeff ~ product of two anarchic Yukawa matrices summed over KK. BRANE Higgs: leading 1/M_KK cancels, residual masked by a divergent cutoff term => UNCALCULABLE (only param as a m_mu^2/Lambda5^2 U_12). BULK Higgs: finite (Eq.49-51). Bounds (anarchy Y=2, bulk Higgs, BR<1.2e-11): **mu->e gamma M_KK>15.8 TeV** (8.0 at Y=1); tau->mu gamma 1.6, tau->e gamma 0.9. Tree mu-e conversion 4.7-5.9 TeV. Tension: dipole ∝ Y*^2 (grows w/ Yukawa) vs tree conversion ∝ 1/Y* (shrinks) => irreducible.
- **Perez-Randall 0805.4652 (the NDA ref under test):** brane IR dipole O ~ (e m_mu/2·16pi^2 k^2)(Y_N Y_N†)_12 ... ; BR ~ 4e-8 (Y_N Y_N†)_12^2 (3TeV/M_KK)^4 => (Y_N Y_N†)_12 <~ 0.02 (3TeV/M_KK)^2 (Eq.31). Bulk-Higgs weaker (<~0.3). NDA conservative vs their explicit. LMFV + small Y_N collapses the bound to ~3 TeV.
- **Agashe-Csaki MFV-lepton 0804.2503:** MFV alignment (bulk masses diag with Y_e) kills leading anarchic dipole => M_KK ~3 TeV allowed.
- **Csaki-Grossman-Tanedo-Tsai "Warped Penguins" 1004.2037 (definitive explicit 5D loop):** spurion = aligned single-Yukawa `b` term (rotated away) + misaligned **Y*^3 `a` term**. C_L ~ (a Y*^3 + b Y* R0^2)(e/16pi^2)(v/sqrt2) 2 m_mu f_L^2 f_E^-1. Bound: |a Y*^2 + b|(3TeV/M_KK)^2 <= 0.015; b=0 => Y* <= 0.12|a|^-1/2 at 3 TeV. Anarchic consistency: **M_KK >~ 6 TeV (minimal), >~4.75 TeV (custodial)**; the `b` term alone forces >~4 TeV. Explicit matches NDA at order of magnitude; nontrivial chiral/gauge cancellations beyond naive power-counting.
- **Beneke-Dey-Rohrwild 1508.01705:** first complete dipole, 3 contributions w/ different Yukawa deps; regularizes finite-width brane Higgs => **M_KK >~ 2 TeV (minimal), ~4 TeV (custodial)** — much weaker than ABP's 15.8 TeV (cutoff-saturated).
- **MEG II 2025 arXiv:2504.15711: BR(mu->e gamma) < 1.5e-13** (repo default; ~80x tighter than the 1.2e-11 ABP/PR used; M_KK ∝ BR^-1/4).
- Csaki-Falkowski-Weiler 0804.1954 (quark-sector anarchy KK-gluon >~21-33 TeV) — context, not the lepton dipole.

## Bound spread (anarchic mu->e gamma, BR<1.2e-11): ~2-16 TeV depending on Higgs localization, custodial, coherent-sum. Central explicit ~ mid-single-digit TeV (~6 TeV minimal). With LMFV/MFV collapses to ~3 TeV.

## Disagreement to flag
Brane-Higgs calculability: ABP says UNCALCULABLE (UV-sensitive, 15.8 TeV cutoff-saturated); Beneke et al regularize finite-width => finite ~2 TeV. Not contradictory (singular limit) but very different numbers; the right number depends on brane-Higgs treatment. NDA (C=0.02) sits between, claims only order-of-magnitude/conservative.

## Net guidance
Treating Perez-Randall C=0.02 NDA as a conservative (safe) order-of-magnitude estimate agreeing with explicit 5D loops to O(1) is faithful. The only places flat NDA is genuinely incomplete: alignment suppression + the `b` misalignment spurion (O(1)-O(few)). In an LMFV lepton setup the dipole is suppressed by (Y_N Y_N†)_12 not the cubic Y*^3.
