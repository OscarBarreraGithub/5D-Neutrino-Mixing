# Implement L009 — BR(τ→3μ). family=charged_lepton. REUSE lfv_three_body module.
REUSE quarkConstraints/lfv_three_body.py (built for L002) with τ→μ flavors (initial=tau, final=mu). Dipole(+contact+interference, constructive envelope) as L002. Pin τ→3μ flavors. SM=0 → pure-NP vs Belle II/limit (from L009.yaml), HARD. RS proxy NEEDS-HUMAN-PHYSICS. Append-only (don't modify L002's μ→3e funcs).
