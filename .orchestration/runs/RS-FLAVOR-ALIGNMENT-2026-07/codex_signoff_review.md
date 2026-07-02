**VERDICT: SOUND-WITH-FIXES**

1. **Claim 1:** Mostly sound. `C4_LR = -(left * right)/M_KK^2` matches [deltaf2.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/deltaf2.py:342). The `0.332` follows from `kappa_epsilon/(2 sqrt(2))` with `kappa=0.94`. Fix: repo `Delta m_K` uses `2|M12|`, not a signed `R_K cos theta_K` observable; if using `cos theta`, define theta as `arg M12`, not `arg(GL GR)`.

2. **F1:** Sound. I recomputed the S1 parquet: 474 PDG-passing survivors at `M_KK<=3 TeV`; 271/474 = 57.2% magnitude-small, 33/474 = 7.0% failer-typical `Delta m_K` tail, and `corr(log eps_K, log dm_K)=0.821`. The arcsin phase-volume explanation is correct for uniform phases.

3. **F3:** Core physics sound. Real `Y_d` makes down rotations/couplings real, so `sin theta_K=0`; `|C4|` remains anarchic; `I_d` remains nonzero from `Y_uY_u^\dagger`; fully real Yukawas make `I_d=0`. Required fix: write `I_d` as the mass-basis, mass-normalized invariant actually used in [instrument_epsK_phase.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/scripts/instrument_epsK_phase.py:151), not the bare flavor-basis `(11)` entry. Also call the draw “NB-like real-down/complex-up,” not genuine Nelson-Barr.

4. **EDM attribution:** Broadly correct. Tree KK-gluon exchange is not the leading quark EDM source; the RS CP problem is the loop dipole with Higgs/KK fermions / Yukawa misalignment. But the note should say `I_d` is an EDM proxy/invariant, not a calibrated neutron EDM prediction.

5. **Novelty boundary:** Mostly honest. Do not claim Redi-Weiler minimal CP or warped/NB CP sequestering as new. Add/cite missing adjacent prior art: Bauer-Malm-Neubert `1110.0471` as another RS flavor-protection cure. It does not sink the stated narrow claims.

6. **Worked target:** Logic is directionally right, but F3 alone proves only “CP-only real-down is insufficient.” Magnitude-only insufficiency needs the known RS EDM/CP problem or a separate scan/proof. Genuine NB must include a radiative/stability proof that suppresses the same `I_d` invariant.

**Required fixes:** tighten `Delta m_K` convention, correct the `I_d` definition, downgrade “NB eta=0” wording to “NB-like toy,” and avoid implying a computed absolute `d_n`.

Sources checked: APS [hep-ph/0408134](https://arxiv.org/abs/hep-ph/0408134), Redi-Weiler [1106.6357](https://arxiv.org/abs/1106.6357), Girmohanta et al. [2203.09002](https://arxiv.org/abs/2203.09002), Cheung-Fitzpatrick-Randall [0711.4421](https://arxiv.org/abs/0711.4421), FPR [0710.1869](https://arxiv.org/abs/0710.1869), Santiago [0806.1230](https://arxiv.org/abs/0806.1230).