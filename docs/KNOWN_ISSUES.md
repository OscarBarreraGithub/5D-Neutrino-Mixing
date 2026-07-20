# Known Issues

Tracked open issues that are **not** being fixed in the current documentation
pass. None of these affect any floor in [`FLOOR_SUMMARY.md`](FLOOR_SUMMARY.md)
or any validated literature reproduction.

---

## KI-1 — B002 / B004 CP-phase use a real SM box amplitude (convention-dependent)

**Severity:** was a live bug in scanned constraints (see the scan-exposure
correction below).
**Status:** **FIXED in code 2026-07-20** (commit 29daa27): B002/B004 now use
the complex SM box amplitude via
`neutral_b_mixing_sm_amplitude()` with phase `arg((V_tq* V_tb)^2)`, and both
constraints carry random-rephasing invariance tests. **Existing scan
classifications that included B002/B004 predate the fix and remain
contaminated until the affected classifications are rerun** (audit P0-5).
Found in the B2 rephasing audit (June 2026).

> **Scan-exposure correction (2026-07-15 audit).** An earlier version of this
> entry claimed S_psiKS / S_psiphi (B002/B004) "were NOT in the production
> quark-only scan". That statement was **false**: both constraints are in
> `QUARK_ONLY_CONSTRAINT_IDS`, their required extras are wired, and the 1M-row
> analysis reports show nonzero B004 veto fractions (100% in several 1 TeV
> cells, 57.1% in one 2 TeV cell). Strict exclusion labels from those runs
> therefore contain a rephasing-dependent constraint.
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
- ~~S_psiKS / S_psiphi (B002/B004) were NOT in the production quark-only
  scan~~ **RETRACTED, see the scan-exposure correction banner above: they WERE
  scanned and did veto points.** Floors driven by K001/B003/C001 are unaffected
  (those are exactly rephasing-invariant), but strict survival labels from the
  June 1M runs include the rephasing-dependent B002/B004 verdicts.

### Fix applied (2026-07-20)

As sketched in the B2 verdict notes: the real `M12_SM = Δm_SM/2` was replaced
by the complex SM box amplitude with phase `arg((V_tq* V_tb)²)`
(`quarkConstraints/ckm_extraction.py::neutral_b_mixing_sm_amplitude`), the
magnitude still fixed by `Δm_SM/2`. `phi_q_np = arg(1 + M12_NP / M12_SM_box)`
is now rephasing-invariant; both constraints pin the prediction across 16
random quark-field rephasings in their test suites. Remaining work: rerun the
affected scan classifications (tracked as audit P0-5 / Gate 5).
