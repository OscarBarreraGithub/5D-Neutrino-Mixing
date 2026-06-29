# Known Issues

Tracked open issues that are **not** being fixed in the current documentation
pass. None of these affect any floor in [`FLOOR_SUMMARY.md`](FLOOR_SUMMARY.md)
or any validated literature reproduction.

---

## KI-1 — B002 / B004 CP-phase use a real SM box amplitude (convention-dependent)

**Severity:** latent bug; does NOT affect any floor or reproduction.
**Status:** open, not fixed. Found in the B2 rephasing audit (June 2026).
**Files:**
- `flavor_catalog_constraints/primary/beauty/B002.py` (`S_psiK_S`)
- `flavor_catalog_constraints/primary/beauty/B004.py` (`S_psiphi`)

### What is wrong

B002 and B004 form the new-physics phase as

```
phi_q_np = arg(1 + M12_NP / M12_SM),   with M12_SM = +Delta_m_SM / 2   (REAL)
```

and pass/fail on `S = sin(2β + phi_d_np)` (B002) / `S_psiphi` (B004). Using a
**real** `M12_SM = Δm_SM/2` is only correct in the convention where the SM box
amplitude `(V_tq* V_tb)²` is real-positive. It is not: `V_td` is irreducibly
complex in any standard CKM parameterization, so `(V_td* V_tb)²` is not real.
Consequently `arg(M12_NP)` measured against a real `M12_SM` is **not**
rephasing-invariant, and the predicted `S_psiKS` / `S_psiphi` shift (by up to
~2 rad in adversarial draws) under a change of the SVD→PDG rephasing convention.

The genuinely physical, rephasing-invariant quantity is

```
phi_q_np = arg(1 + M12_NP / M12_SM_box),
M12_SM_box ~ eta_B (V_tq* V_tb)^2 S0(x_t)   (COMPLEX SM box)
```

### Why it is NOT a B2 regression and NOT urgent

- This pre-existed the B2 (SVD→PDG rephasing) change. The B2 5th-anchor choice
  (V_ts vs Bauer's V_cs) is a harmless gauge fixing and is **not** the root
  cause: neither anchor makes `(V_td* V_tb)²` real. Switching to the Bauer
  V_cs anchor would NOT fix B002/B004 and would break nothing B2 got right.
- `epsilon_K` and every `|M12|` / Δm observable are exactly rephasing-invariant
  (bit-identical numeric check), so K001/B003/C001 and all floors are unaffected.
- **S_psiKS / S_psiphi (B002/B004) were NOT in the production quark-only scan**,
  so no scanned exclusion or published floor depends on them.

### Fix sketch (when acted on)

Replace the real `M12_SM = Δm_SM/2` in B002/B004 with the complex SM box
amplitude `M12_SM_box ∝ η_B (V_tq* V_tb)² S0(x_t)` and form
`phi_q_np = arg(1 + M12_NP / M12_SM_box)`. Then the predicted CP phases become
rephasing-invariant. Add a test that pins the prediction across two random
rephasings of the fit. See the B2 verdict notes for the analytic argument and
the 200-draw numeric demonstration.
