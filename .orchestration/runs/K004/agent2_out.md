1. NIT (SM convention) quarkConstraints/rare_kaon_snd.py:8-15,200-214: SM structure is correct. Repo CKM gives λ=0.2254437, λc=-0.219313-0.000142i, λt=-3.249e-4+1.420e-4i, κ+=5.255e-11, Pc=0.404, Xt=1.481, BR_SM=8.498e-11. Exact Buras/BGS formula includes (1+ΔEM)=0.997; omission is a 0.3% convention issue, not a serious SM failure.

2. NIT (doc citation) quarkConstraints/rare_kaon_snd.py:16-18: wording overstates BGS-2021 compatibility. BGS 2021 quotes Xt=1.462 and BR(K+)=7.73(61)e-11; the 8.60±0.42e-11 anchor is Buras-Venturini 2022.

3. NIT/OK (NP proxy) quarkConstraints/rare_kaon_snd.py:20-34,53-57,267-296: proxy has defensible leading scaling and normalization: ΔX ∝ Δν(ΔsdL+ΔsdR)/(gSM² M²), with gSM²=1.781e-7 GeV^-2 and same-sign L/R for K→π vector matrix element. NEEDS-HUMAN-PHYSICS is honest; it clearly omits full EW KK/Z/Z′, chirality-specific couplings, and neutrino-sector matching.

4. OK (budget) flavor_catalog_constraints/primary/kaon/K004.py:230-299,302-419: NA62 9.6(+1.9/-1.8)e-11 plus SM 8.60±0.42e-11 gives σ+=1.946e-11, σ-=1.848e-11; dimensionless BR units and HARD severity are defensible.

Sources checked: arXiv:1503.02693, arXiv:2105.02868, arXiv:2203.10099, arXiv:2604.12649.

PHYSICS-OK