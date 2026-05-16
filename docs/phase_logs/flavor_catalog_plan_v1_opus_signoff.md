# Flavor Catalog Plan v1 — Opus Final Approval

**Reviewer**: Opus (independent final approval)
**Date**: 2026-05-16
**Plan**: `docs/phase_logs/flavor_catalog_plan_v1.md` (commit `744e70c`)

## Verdict
APPROVE-FOR-EXECUTION

## Summary
Plan v1 cleanly closes both Codex BLOCKERs (non-PKA role specs, hard iteration caps) and the status-semantics WARNING, and does so without weakening any of the v0 PASS items. The discovery-mode layout (`flavor_catalog/` on a new `flavor-catalog/2026q2` branch), the 128-entry initial process list, the writer/checker separation rule, the PKA/WA/CA/DA/Opus deliverables-and-success-criteria blocks, and the cap-hit status transition rules are all concrete enough for orchestration to begin. Residual gaps (notably oblique parameters S/T/U and CKM unitarity proper) are exactly the kind of thing the DA-3 scope is designed to catch and can be handled in-band; none rises to a BLOCKER.

## Item-by-item (1-6)

1. **Faithfulness to PI intent** — PASS.
   - Discovery-mode subdir separate from repo: explicit at `docs/phase_logs/flavor_catalog_plan_v1.md:22-79` and `:83`, with confirmed branch policy at `:90-93`.
   - PKA design: full spawn prompt at `:385-414`, deliverables at `:416-423`.
   - Per-process .tex template: required fields at `:110-144` (YAML) and `:163-195` (TeX).
   - Writer/Checker separation: enforced at `:453` (WA), `:484` (CA), and orchestration at `:547`.
   - Depth-not-breadth: PKA "one process by default … up to three only when … the same literature and PDG table covers all of them" at `:379`; WA "batches of 3-5 processes from the same family" at `:431`.
   - Code-coverage annotation: rubric and aliases at `:205-213` and rubric body at `:590-595`.
   - Implementation-difficulty rubric: explicit at `:590-595` with initial expectations at `:597-604`.
   - Multiple discovery rounds: DA at `:512` requires at least two rounds, capped at four at `:621`.

2. **Physics completeness of the initial process list (Section C)** — PASS-WITH-FOLLOWUP.
   - Kaon staples all present (K001-K022). The PI's ambiguous `K1 -> pi gamma` is correctly flagged for clarification at `:691` and mapped to K013/K014.
   - Charm: mixing parameters x_D, y_D folded into C001 at `:246`; |q/p|, phi_D at C002 (`:247`); Delta A_CP at C003 (`:248`); D -> mu mu, D -> gamma gamma, D -> rho gamma, D -> l nu all present (C004-C011). The `x_D, y_D` decomposition is hidden inside C001's "x, y" — acceptable but worth surfacing in the eventual PDG-snapshot file.
   - Beauty: full set including R_K, R_K*, R_D, R_D*, R_{J/psi}, B_s -> mu mu, b -> s gamma exclusive+inclusive, B -> K(*) nu nubar, Lambda_b -> Lambda ll, charmless non-leptonic, LFV B decays. Comprehensive.
   - Top/Higgs/EW: T001-T022 cover all FCNC top, Higgs LFV, Higgs quark FCNC, Z LFV, Z -> b/c pole observables.
   - LFV: L001-L020 cover all standard tau and mu modes including hadronic and mu-e conversion in three nuclei. L023-L025 add trident, PMNS, 0nbb.
   - EDM: e/mu/tau/n/p/d/Hg/Ra/Xe plus chromo-EDM and Weinberg. Comprehensive.
   - **Missing/under-specified** (NIT-level, picked up by DA-3 by design at `:492`):
     (a) Oblique electroweak parameters S/T/U — a direct RS-KK constraint historically; not in Section C.
     (b) CKM unitarity / first-row Cabibbo-angle anomaly — K018 covers V_us, but |V_ud| beta-decay + unitarity test is not its own entry.
     (c) g-2 of the muon — historically a flavor-adjacent dipole observable like mu->e gamma; arguable scope, worth listing as a DEFERRED-SCOPE candidate.
     (d) K_L-K_S width difference and K -> pi pi amplitude ratios are folded into K022; OK.
     (e) BR(B_s -> phi phi) CP phi_s^{sss} could be split out from B034 if needed.
   The plan's own "discovery agents should add missing baryon, hyperon, and nuclear processes" at `:373` and DA-3 scope at `:492` are the correct mechanism. Not a BLOCKER.

3. **Operational realism** — PASS.
   - Section E sequence at `:540-556` is orderable: branch creation → scaffold → PKA → WA → CA → iterate → DA round 1 (gated at 50% CA pass) → DA round 2 → Opus → master compile.
   - Concurrency budgets are concrete: 8-12 PKAs, 4 WAs, 4 CAs, 3 DAs (`:381`, `:429`, `:457`, `:488`).
   - Iteration caps are quantitative and bound to status transitions: PKA 2, WA/CA 3, DA 4, Opus 1 (`:619-622`); each cap forces a `WA-TRIAGE / OPUS-ARBITRATION / BLOCKED-PI / DEFERRED-SCOPE` transition.
   - Convergence rule at `:560-566` is testable.

4. **Risk hygiene** — PASS.
   - License: "Do not track publisher PDFs unless the license explicitly allows redistribution" at `:100`; arXiv PDFs and PDG text snapshots only.
   - Size caps: 200 MB total / 5 MB per file at `:99`.
   - Reproducibility: sha256 manifests at `:130-132` and `:421`.
   - Secrets: implicit (no API keys discussed); acceptable for catalog-only work.

5. **Internal consistency** — PASS.
   - Cross-references resolve: `quarkConstraints/deltaf2.py` line numbers in Section C aliases line up with the actual file; YAML schema fields used in Section B match those referenced in Sections D and G.
   - The `status_history` mechanism in B (`:110-160`) is consistently invoked in D (WA, CA, DA, Opus deliverables) and G (cap-hit transitions). No contradiction between the legal-transition graph and the convergence rule.
   - Process count statement "128 entries" at `:373` matches my count (K22 + C12 + B37 + T22 + L25 + E10 = 128).

6. **PI-facing open questions (Section I)** — PASS.
   - Six questions, each answerable in one sentence: `K1 -> pi gamma` disambiguation, scope confirmation, branch policy, PDF licensing posture, rc1.1/rc2 relationship, and `modern/` vs `deltaf2.py` integration target.
   - The set is the right set: it captures the only PI judgment calls the orchestrator cannot make safely.

## Specifically validated
- [x] Discovery-mode subdir separate from rest of repo (`:22-79`, `:83`, `:90-93`)
- [x] PKA / WA / CA / DA / Opus role specifications complete with deliverables + success criteria (`:385-538`)
- [x] Writer/Checker SEPARATION rule (`:427`, `:453`, `:456`, `:484`, `:547`)
- [x] Per-agent depth-not-breadth (`:379` PKA M=1 default; WA batch 3-5)
- [x] Iteration caps quantitative and tied to status transitions (`:619-622`, `:624-628`)
- [x] Initial process list >= 80 entries (128 entries verified)
- [x] Initial process list covers kaon, charm, B, B_s, top, LFV, EDM, EW, Higgs/top FCNC (Section C `:215-371`)
- [x] Code-coverage annotation rubric + initial guesses for at least the seed processes (`:205-213` aliases + seed mappings; `:597-604` initial difficulty)
- [x] License + size policies sensible (`:95-100`)

## Findings

### NIT 1 — Section C — Oblique parameters and CKM unitarity not pre-seeded
**Severity**: NIT. **Evidence**: `:215-371` lists no S/T/U entry and no standalone first-row CKM-unitarity / |V_ud| entry. **Fix**: instruct DA-3 (`:492`) to consider as proposed additions: `EW001 S/T/U oblique`, `EW002 first-row CKM unitarity (|V_ud|+|V_us|+|V_ub|)`, `EW003 |V_cb| vs |V_ub| inclusive-exclusive tension`. These can be admitted at DA round 1 with no impact on cap budgets.

### NIT 2 — Section C, charm — D mixing observables packed into one row
**Severity**: NIT. **Evidence**: `:246` C001 lists "`x`, `y`, `Delta m_D`" together. **Fix**: in the PKA execution for C001, require the YAML sidecar to emit a separate `pdg_or_equivalent` block per observable (x_D, y_D, Delta m_D, delta_y) rather than a single combined value. No plan edit needed; this is a PKA-prompt clarification at dispatch time.

### NIT 3 — Section D, Opus — Parallelism stated as "1-2 sign-offs at a time"
**Severity**: NIT. **Evidence**: `:516`. **Fix**: orchestrator should default to 1 for the first two Opus rounds (to lock global consistency), only allowing 2 when the catalog index passes a freeze-check. Not a plan defect; an execution-time policy.

### NIT 4 — Section E — DA round 1 gated at 50% CA pass
**Severity**: NIT. **Evidence**: `:550`. **Fix**: clarify in the DA-1 spawn message that "50% CA pass" is computed across the whole 128-entry initial list, not within a single family. Cosmetic clarification for the orchestrator.

### NIT 5 — Section H — Agent-hour total
**Severity**: NIT. **Evidence**: `:653-660` sums to ~466-768 agent-hours plus 30-50 orchestrator hours. This is correctly large but worth flagging to the PI in the kickoff message so the budget expectation is shared before execution.

No BLOCKERs. No WARNINGs.

## Recommendation
Ready for execution under the noted NITs. NITs 1, 2, 4 fold into DA-1 / PKA-dispatch policy; NITs 3, 5 are orchestrator-side. None requires another planner round.

## Suggested answers to Section I (your judgment, advisory)

1. `K1 -> pi gamma` — interpret as `K_L -> pi^0 gamma gamma` (K013). Belle-II/NA62-friendly, has a real branching-fraction limit in PDG, and matches the "(off-diagonal varepsilon')" framing the PI used next to it. K014 (`K_S -> pi^0 gamma gamma`) should still be carried as a secondary entry; flag for PI confirmation but proceed.
2. **Scope** — proceed with the full default scope (quark + charged-lepton + neutrino + EW + Higgs/top FCNC + EDM-adjacent). Pruning is cheaper after PKA than before.
3. **Branch** — execute on `flavor-catalog/2026q2` cut from `paper/quark-scan-2026q2`. Do not pollute the paper branch.
4. **PDG PDFs** — minimal text snapshots + URL + access date + sha256 only. Do not track publisher PDFs even if access is allowed. Cheap to upgrade later, expensive to retract.
5. **rc1.1/rc2 relationship** — companion artifact, not part of the rc paper. The quark-scan paper stays as-is; the catalog informs the next-paper roadmap.
6. **Integration target** — for any process that graduates from catalog to live constraint, target `quarkConstraints/modern/`. Treat `quarkConstraints/deltaf2.py` as legacy-stable; new constraints should not extend it.

===FLAVOR_CATALOG_PLAN_V1_OPUS_SIGNOFF_END===
