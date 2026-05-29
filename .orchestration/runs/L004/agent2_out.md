1. NIT: Check 1 N/A/OK. L004 is not a Δm or CP observable; code uses the coherent rate `CR = omega_conv/Gamma_capture` with modulus-squared amplitudes, not Re/Im `M_12`. Correct for μ→e conversion. `quarkConstraints/mu_e_conversion.py:6-14`, `:328-330`, `:367-374`.

2. NIT: Check 2 N/A. No ΔF=2/QCD `*_with_running` path applies here; inputs are already low-energy KKO nucleon proxies. No non-running deltaf2 verdict path is used. `flavor_catalog_constraints/primary/charged_lepton/L004.py:311-319`, `flavor_catalog_constraints/physics_adapters/mu_e_conversion.py:202-207`.

3. NIT: Budget OK. Uses observed SINDRUM II Au 90% C.L. limit `7.0e-13` as HARD pure-NP upper bound, with SM set to 0; not a projection or residual subtraction. `L004.yaml:72-87`, `L004.py:64-68`, `:263`, `:335-342`, `:364-372`.

4. NIT: Anchor/nuclear numbers OK. Au inputs match KKO: `D=0.189`, `Vp=0.0974`, `Vn=0.146`, `Sp=0.0614`, `Sn=0.0918`, `Gamma_capture=13.07e6 s^-1 = 8.60283e-18 GeV`. `quarkConstraints/mu_e_conversion.py:177-190`.

5. NIT: Severity/diagnostics OK. RS scalar/vector mapping is explicitly flagged `NEEDS-HUMAN-PHYSICS`; unknown dipole-contact phase uses destructive lower envelope for verdict and reports full interval. `L004.py:15-22`, `:373-380`; `mu_e_conversion.py:363-366`, `:463-467`.

PHYSICS-OK