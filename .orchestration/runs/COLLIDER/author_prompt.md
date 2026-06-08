# Add lepton-free COLLIDER KK-mass searches to the quark-only scan (codex author, gpt-5.x xhigh)
Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. AUTHOR task: SHORT plan, then implement. A codex reviewer AND an Opus reviewer independently review; BOTH must APPROVE before commit.

## GOAL
The quark-only Bucket-1 allowlist (37 constraints, in `scripts/run_full_catalog_scan.py`, the `--quark-only` mode committed 070cdb6) currently EXCLUDES the entire `collider_rs` family (CR001..CR014), the direct LHC KK-mass searches. Add the LEPTON-FREE ones so the scan reflects the direct-search M_KK floor. These are pure-quark/boson resonance recasts and need no swept-lepton data.

## WHAT TO DO
1. For EVERY CR001..CR014, OPEN its module under `flavor_catalog_constraints/primary/collider_rs/` and determine what extras it reads. Build an IN/OUT table: a CR is IN only if it reads NO swept-lepton quantity (no z_delta_g_*_e, z_delta_g_L_nu, lepton couplings, lepton final-state object) AND the quark-only point already provides all its inputs (e.g. `kk_gluon_mass_gev`, `kk_ew_mass_gev`, `quark_mass_basis_couplings`). Likely IN: CR001 (g_KK->tt), CR002/003/004/008/010 (vector-like quark pair searches), CR007/013 (KK graviton->WW/ZZ/gamma gamma), CR012 and CR014 if purely bosonic/quark. Likely OUT (lepton final states/inputs): CR005 (KK->ll), CR006 (W_KK->l nu), CR009 (llqq contact), CR011 (jj W_L W_L). VERIFY each yourself; do not assume my list.
2. Add the verified lepton-free CRs to the quark-only allowlist. Confirm the quark-only point build already supplies their inputs (kk_gluon_mass_gev = tile.mkk_gev is set; kk_ew_mass derived from spectrum; quark couplings built). If any IN constraint needs an extra the quark-only point does not yet provide, either wire that quark-side extra in (if trivial and lepton-free) or move it to OUT with reason. NO constraint may silently degrade to an SM/zero default.
3. Tag them honestly: these are sigma*BR / mass-exclusion PROXY recasts -> tag `proxy`, HARD, vetoing in the inclusive accounting (they carve M_KK as a near-flat direct-search floor). Update the allowlist count + the candidate/deferred bookkeeping.

## GATE TESTS
- The added CRs evaluate active=True in `--quark-only` mode on a real point, read only quark/geometry inputs, and none degrade to a silent default (assert each reads its real mass/coupling input).
- A low-M_KK point (~1-2 TeV) is excluded by the direct-search CRs; a high-M_KK point (~30-50 TeV) passes them -> confirms the direct-search floor behaves monotonically.
- Allowlist count updates correctly; OUT (lepton) CRs remain excluded and tagged deferred_lepton_followup.
- `python -m pytest tests/ -q` stays green (+ a unit test that the quark-only allowlist now contains the lepton-free CRs and that they read no lepton extras).

## OUTPUT (<=16 lines)
Short plan; the IN/OUT table (CR id -> extras read -> lepton-free? -> IN/OUT -> why); new allowlist count; confirmation the quark-only point supplies their inputs (no silent defaults); the low-vs-high M_KK monotonicity check; pytest counts. End with: COLLIDER-AUTHOR-DONE.
