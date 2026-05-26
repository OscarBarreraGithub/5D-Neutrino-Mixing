# CLEANUP PLAN — Phase 2 (Planner Round 2)

Repo: `/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing`
Author: Claude Opus 4.7 (planner role, cleanup-phase round 2)
Status: REVISED — round-1 reviewer critique addressed (CR-1/CR-2/CR-3 + RC-1..7 + M-1..8)
Input: `.orchestration/ISSUES.md` (78 issues from Phase 1 review, all currently in `## Open`).
Round-1 inputs: `CLEANUP_PLAN.md` (r1), `CLEANUP_REVIEW_R1.md`.

The goal: design exactly how to close out every Phase-1 issue, scoped into atomic units
that an Opus agent (or in some cases the orchestrator's shell harness) can execute in
one run. Mirrors the structure of `MERGE_PLAN.md`. Per CLAUDE.md, Claude-only — no codex.

This plan does NOT execute. It is the blueprint the orchestrator will hand to
per-unit Opus agents (or run inline for orchestrator-deterministic units).

---

## A. Issue triage table

Severity totals across the 78 open issues: HIGH=1, MEDIUM=2, LOW=27, INFO=48 (49 minus
the closed `INFRA-1`). Tag distribution: physics=1, numerics=4, code=10, docs=51,
infra=4, provenance=4, tests=1, ux=1. (`physics/numerics/code/infra` counts include the
sub-tags the reviewers used; `tests`, `provenance`, `ux` are downgraded to the parent
tag in the bookkeeping below.)

Columns:
- **ID** — `R##-I#` from ISSUES.md
- **Sev** — H/M/L/I (HIGH/MEDIUM/LOW/INFO)
- **Tag** — physics/numerics/code/docs/infra
- **Effort** — T(<10min) / S(10-60min) / M(1-3h) / L(>3h)
- **Risk** — risk of touching prod code: lo/med/hi
- **Dep** — issues that must close first (or "—" if standalone)
- **Unit** — atomic cleanup unit (C##) that resolves it (assigned in §C)

| ID | Sev | Tag | Effort | Risk | Dep | Unit | One-liner |
|---|---|---|---|---|---|---|---|
| R01-I1 | L | code | S | lo | — | C12 | Dual source-of-truth M_CHARM/M_BOTTOM/M_TOP_MS vs PDG dict |
| R01-I2 | I | code | T | lo | — | C20 | `MASS_TOLERANCE_FLOOR=0.003` TODO — calibration deferred (see footnote) |
| R02-I1 | I | docs | T | lo | — | C18 | MERGE_PLAN R02 row says `pytest.ini`; landed in `pyproject.toml` |
| R02-I2 | I | code | T | lo | — | C13 | spurion-seed literals need provenance comment |
| **R03-I1** | **H** | **code** | **M** | **hi** | — | **C01** | **Duplicate pre-audit kaon constants in `modern/phenomenology.py`** |
| R03-I2 | M | numerics | S | lo | R03-I1 | C01 | Pin-test gap: `_KAON_*` not asserted vs `deltaf2` |
| R03-I3 | M | docs | L | med | R03-I1; hole-#6 done (R04 closed) | C02a-code (+ C02c deferred) | 5 hole-#5 follow-ups |
| R03-I4 | I | docs | T | lo | — | C18 | No top-level `REFERENCES.md` |
| R04-I1 | L | numerics | S | lo | — | C03 | Audit-script closed-form assumes upper-tri LR ADM — add assert+test |
| R04-I2 | L | code | T | lo | — | C03 | function-local import in `deltaf2.py:586` (lift unconditionally; see M-7) |
| R04-I3 | I | physics | T | lo | — | C03 | `paper_0710_1869` LR map sign consistency — verified, document |
| R04-I4 | I | docs | T | lo | — | C19 | Four R04 deferred follow-ups |
| R05-I1 | I | docs | T | lo | — | C18 | MERGE_PLAN "8 scans" vs 9 dirs (Run-3 moreUV/moreIR are zero-pass) |
| R05-I2 | L | infra | S | lo | — | C15 | `pytest_selection/` dir empty (24 missing files) |
| R05-I3 | I | docs | T | lo | — | C19 | 22.49x vs 4.47x cosmetic note in methodology |
| R06-I1 | L | code | T | lo | — | C20 | CFW constants live in audit driver; refactor only if 2nd driver added |
| R06-I2 | I | docs | T | lo | — | C18 | CFW PNG relocated to `exploratory/` — already harmless |
| R06-I3 | I | docs | T | lo | — | C19 | Methodology rounds 47.26 to ~47 TeV; cite audit doc |
| R07-I1 | L | physics | T | lo | — | C19 | Wilson z=1.92 vs standard z=1.96 — add footnote |
| R07-I2 | L | numerics | S | lo | — | C04 | `wilson_upper_limit` tests cover k=0 only — add k>0 parametrized test |
| R07-I3 | I | infra | T | lo | — | C20 | Manifest `code_sha` inferred from mtime — forward-looking SLURM fix |
| R07-I4 | I | docs | T | lo | — | C16 | `CLAUDE.md:12-15` paper-branch pointer (auto-resolves Phase 10) |
| R08-I1 | I | docs | T | lo | — | (no-op) | `bf6186c` no-op PDF rebuild — accept as commit-log noise |
| R08-I2 | I | numerics | S | lo | — | C05 | Spot-check `f_{u,2}=0.16` reproduces from documented c-pattern |
| R08-I3 | I | infra | M | lo | — | C20 | Historical figure dirs — defer (accept-as-risk) |
| R09-I1 | I | docs | T | lo | — | C18 | scaffold collapsed edm/+neutrino_universality/ into edm_neutrino/ |
| R09-I2 | I | docs | T | lo | — | C17 | 11 `.gitkeep` placeholders alongside real content |
| R10a-I1 | I | docs | T | lo | — | C06 | K001 `checker_agent_id: "CA"` stale wave label |
| R10a-I2 | I | docs | T | lo | R03-I1 | C06 | K001 `open_issues` flags same condition as R03-I1; cross-link |
| R10a-I3 | I | docs | T | lo | — | C06 | K003 PDG listing vs DataBlock canonical choice |
| R10b-I1 | I | docs | T | lo | — | C07 | B011 observable label typo "B_s gamma" → "B→X_s γ" |
| R10b-I2 | I | docs | T | lo | — | C07 | B009 status_history entry missing `state:` field |
| R10b-I3 | I | docs | S | lo | — | C07 | 5 beauty yamls have unresolved open_issues |
| R10c-I1 | I | docs | T | lo | — | C08 | T010 absorbs T011 — create stub T011 redirect OR doc mapping |
| R10c-I2 | I | docs | S | lo | — | C08 | T002/T010/E001 lack top-level `sha256sums.txt` |
| R10c-I3 | I | docs | T | lo | — | C08 | 4 R10c PKAs lack wave1-named worklog |
| R11-I1 | L | code | T | lo | — | C09 | B017 introduced inside WA-polish commit (workflow seam) |
| R11-I2 | I | docs | T | lo | — | C09 | DA-1 worklog: 4 PI escalations without resolution thread |
| R11-I3 | I | docs | T | lo | — | C09 | B025 Belle II CKM-2025 PDF-policy gap — monitor only |
| R12-I1 | I | docs | T | lo | — | C18 | R12 dispatch prompt mis-id'd L008 |
| R12-I2 | L | code | T | lo | — | C09 | PKA `7d3da08` bundles EW001+E004 (history-immutable) |
| R12-I3 | I | code | T | lo | — | C09 | `bcd8907` triple-bundles workstreams (history-immutable) |
| R12-I4 | I | docs | S | lo | — | C09 | DA-1 escalations partially resolved by W4 — back-annotate |
| R13-I1 | I | docs | T | lo | — | C18 | R13 dispatch prompt mis-id'd 3 process IDs |
| R13-I2 | L | docs | T | lo | — | C09 | `round_002_index.md` K009/K010 description column wrong |
| R13-I3 | I | docs | T | lo | — | C09 | DA-3 worklog: charm Wave-4 vs Wave-5a lineage |
| R14-I1 | I | docs | T | lo | — | C18 | MERGE_PLAN R14 row mis-attributes K012/K018 |
| R14-I2 | I | docs | T | lo | — | (no-op) | B001 x_d/chi_d open_issue — auto-resolved by R15 chain |
| R14-I3 | L | docs | T | lo | — | C16 | Codex quota pause cites gpt-5.5 vs config gpt-5.4 |
| R15-I1 | L | docs | T | lo | — | C18 | MERGE_PLAN R15 SHA typo `7500919` → `7500794` |
| R15-I2 | I | docs | S | lo | — | C09 | E009 fact-check verdict PARTIAL is JS-only artifact |
| R15-I3 | I | docs | T | lo | — | (no-op) | L023 trident family classification — accept |
| R15-I4 | I | docs | T | lo | — | C09 | DA-4 worklog: no round-2 closure addendum |
| R16-I1 | I | docs | T | lo | — | (no-op) | T012 absorbs T013 — accept |
| R16-I2 | I | docs | T | lo | — | (no-op) | R16 prompt anchor stale |
| R16-I3 | I | docs | T | lo | — | C09 | B012 promotion bookkeeping framing |
| R16-I4 | I | docs | T | lo | — | C18 | B001/B003 arbitration timestamp predates W7 PKA |
| R17-I1 | L | infra | S | lo | — | C14 | external_research PDFs lack MANIFEST + sha256sum |
| R17-I2 | L | docs | M | lo | — | C10 | factcheck_status.md drift v0.2 (Wave-7 5 rows missing at tag) |
| R17-I3 | I | docs | S | lo | — | (closed) | Catalog "v0.2" snapshot — **TAG ALREADY EXISTS** at `835cf48` (verified) |
| R18-I1 | L | docs | T | lo | — | C10 | v0.3 tag annotation "87V+1P" vs report "86V+2P" |
| R18-I2 | L | docs | T | lo | — | (no-op) | Commit `37beabb` message says K020 but diff is B013 — history-immutable |
| R18-I3 | L | docs | M | lo | — | C10 | factcheck_status.md not regenerated at v0.3 |
| R18-I4 | I | docs | T | lo | — | (no-op) | SECONDARY rationale in two places — desirable redundancy |
| R19-I1 | L | infra | T | lo | — | C18 | MERGE_PLAN R19 SHA typo `1cf8b57` → `1cd8b57` |
| R19-I2 | L | infra | T | lo | — | C18 | MERGE_PLAN R19 missing `82daa9b` (count 33 → 34) |
| R19-I3 | L | docs | T | lo | — | C10 | v0.4 tag annotation "100V+2P" vs report "101V+1P" |
| R19-I4 | L | docs | M | lo | — | C10 | factcheck_status.md not regenerated at v0.4 |
| R19-I5 | I | docs | T | lo | — | C18 | MERGE_PLAN `claim_level` should read `limit_type` |
| R20-I1 | L | docs | T | lo | — | C11 | SESSION_NOTES.md §1 tally drift |
| R20-I2 | I | docs | T | lo | — | C11 | AGENTIC_WORKFLOW.md cites gpt-5.5; global says gpt-5.4 |
| R20-I3 | I | docs | T | lo | — | C11 | SESSION_NOTES.md §8 lists 3 sign-off rounds; catalog has 5 |
| R21-I1 | L | docs | T | lo | — | C11 | WEBSITE_RUNBOOK.md tally drift |
| R21-I2 | I | infra | T | lo | — | C20 | cloudflare-pages.config.md is markdown not wrangler.toml — accept |
| R21-I3 | L | docs | T | lo | — | C11 | WEBSITE_RUNBOOK.md STOP-trigger wording |
| R22-I1 | L | code | M | lo | — | C04 | `notation.ts`/`prose.ts` no regression tests |
| R22-I2 | L | docs | T | lo | — | C11 | `/methodology/` page reachable but not in nav |
| R22-I3 | I | code | S | lo | — | C11 | Dead `<th>Tier</th>` / `data-tier` in EntryTable |

**Bookkeeping checkpoint** (RC-3 fix): 78 issue rows above. **7** rows are tagged
`(no-op)` (auto-resolved or accept-as-is, no code change): R08-I1, R14-I2, R15-I3,
R16-I1, R16-I2, R18-I2, R18-I4. One row (R17-I3) is now **(closed)** — verified by
`git tag -l 'flavor-catalog-v0.2'` returning a tag at `835cf48` (the actual
master-compile v0.2 commit, not the planner-r1-cited `42ac647`); this issue closes
without any cleanup action and is moved into `## Closed / Accepted-risk` by C21 with a
"PRE-EXISTING-TAG-VERIFIED" annotation. 70 rows route into 21 cleanup units (C01..C21).

---

## B. Priority tiers and serial-wall-time estimates

Tier ordering is by physics impact + dependency chain. Within a tier, units run in
the listed order. Tier total times are rough single-Opus-agent serial estimates;
**M-6 reclassification** moves four units (C15, C17, C18, C20) into the
**orchestrator-deterministic** lane where the orchestrator's shell harness executes
the action and an Opus agent only reviews the resulting commit for §E.2 compliance.

### Tier 1 — Physics-affecting fixes that need careful review (~5 issues, ~3-4 h Opus + SLURM wall-time deferred to C02c)

Issues: R03-I1, R03-I2, R03-I3 (CLI-flag portion only — see §D and §G1), R04-I1,
R04-I2, R04-I3.

This is the only tier where the user is asked to gate before merging. The R03-I1
fix changes a live numerical pipeline and the collaborator artifacts must be
re-exported (see §D). **Tier 1 ends with annotated checkpoint tag
`cleanup-tier1-complete`** (M-4), pushed to origin so the user can gate review.

- **C01** — R03-I1 (resolve duplicate kaon constants) + R03-I2 (pin-test extension).
  Split into two commits inside one Opus session: **C01a** (bugfix + tests),
  **C01b** (collaborator artifact regen).
- **C02a-code** — R03-I3 sub-tasks 1, 3, 4 (band paragraph using symbolic
  `1/sqrt(|Δε_K|)` scaling) plus the CLI flag (sub-task portion that does NOT
  require SLURM output) plus manifest placeholders. Single Opus atom.
- **C03** — R04-I1 (audit-script upper-tri assertion) + R04-I2 (function-local
  import — lift unconditionally per M-7) + R04-I3 (cross-check
  `paper_0710_1869` LR map; document-only)

Tier-1 serial estimate (Opus interactive): **3-4 hours**. SLURM RUNA reruns are
deferred to **C02c** (separate flow, see Tier-6 below).

### Tier 2 — Test/CI hardening (~3 issues, ~1-2 h)

- **C04** — R07-I2 (`wilson_upper_limit` k>0 parametrized test, scipy>=1.10 floor
  or closed-form Wilson UL fallback — see RC-1) + R22-I1 (`notation.ts`/`prose.ts`
  vitest snapshots for 5 entries, snapshot dir
  `flavor_catalog/website/src/lib/__tests__/__snapshots__/`).
- **C05** — R08-I2 (audit note that `f_{u,2}=0.16` reproduces from documented
  c-pattern).

Tier-2 serial estimate: **1-2 hours**.

### Tier 3 — Provenance + audit-doc consistency (~6 issues, ~2-3 h)

- **C10** — factcheck_status.md consolidation drift (R17-I2, R18-I1, R18-I3,
  R19-I3, R19-I4). Adopt lean (B) per G3: per-master-compile annotation, optional
  aggregate script.
- **C14** — R17-I1 external_research provenance manifest + sha256sums.

Tier-3 serial estimate: **2-3 hours**.

### Tier 4 — YAML/tex catalog metadata sweep (~14 issues, ~2-3 h)

- **C06** — Kaon yaml sweep (R10a-I1, R10a-I2, R10a-I3) — runs after C01 (dep)
- **C07** — Beauty yaml sweep (R10b-I1, R10b-I2, R10b-I3)
- **C08** — Top/Higgs/EW + EDM/charged_lepton yaml sweep (R10c-I1, R10c-I2, R10c-I3)
- **C09** — Worklog/signoff/factcheck text drift sweep (R11-I1, R11-I2, R11-I3,
  R12-I2, R12-I3, R12-I4, R13-I2, R13-I3, R15-I2, R15-I4, R16-I3)

Tier-4 serial estimate: **2-3 hours**.

### Tier 5 — Documentation drift + bookkeeping (~10 issues mixed-lane, ~0.5-1.5 h Opus + ~30 min orchestrator-shell)

- **C11** — SESSION_NOTES + HANDOFF + WEBSITE_RUNBOOK alignment (R20-I1, R20-I2,
  R20-I3, R21-I1, R21-I3, R22-I2, R22-I3) — **Opus**
- **C12** — Quark mass constants SSOT (R01-I1) — **Opus**
- **C13** — Spurion-seed provenance comment (R02-I2) — **Opus**
- **C15** — Backfill pytest_selection/*.txt (R05-I2; 24 files) —
  **orchestrator-deterministic** (shell loop; Opus reviews commit only)
- **C16** — CLAUDE.md + codex-version touchups (R07-I4, R14-I3). **R17-I3 dropped**
  per CR-1 (tag already exists at `835cf48`). — **Opus**
- **C17** — `.gitkeep` cleanup (R09-I2) — **orchestrator-deterministic**
  (`git ls-files` filter + `git rm`; Opus reviews commit only)
- **C18** — MERGE_PLAN.md retroactive corrections (R02-I1, R03-I4, R05-I1, R06-I2,
  R09-I1, R12-I1, R13-I1, R14-I1, R15-I1, R16-I4, R19-I1, R19-I2, R19-I5) —
  **orchestrator-deterministic** (shell `grep -n` for error patterns per RC-5; Opus
  reviews commit only)
- **C19** — Methodology-note polish (R04-I4, R05-I3, R06-I3, R07-I1) — **Opus**
- **C20** — Accept-as-risk / forward-only INFOs (R01-I2, R06-I1, R07-I3, R08-I3,
  R21-I2) — **orchestrator-deterministic** (5 stanza inserts into ISSUES.md; Opus
  reviews commit only)
- **C21** — Final ISSUES.md sweep + tag (`v2026q2-post-cleanup`) — **Opus**

Tier-5 serial estimate: **~0.5-1.5 hours of Opus** + **~30 min of orchestrator-shell**
(C15+C17+C18+C20 collectively reclassified per M-6).

### Tier 6 — DEFERRED-BUT-TRACKED (separate flow, not part of this cleanup wall-time)

- **C02b** — Optional sensitivity-band figure (R03-I3 task 5). Deferred to paper
  finalization phase (planner-r1 disposition retained).
- **C02c** — SLURM RUNA reruns at three ε_K-budget edges + ingest + paragraph
  refinement (R03-I3 task 2 substantive scans). Status: `deferred`.
  Tracked in `CLEANUP_QUEUE.md` so it does not get lost; executed when SLURM
  capacity + user attention allow. Not blocking tag `cleanup-tier1-complete` and
  not blocking the final `v2026q2-post-cleanup` tag.

**Grand total** (Tiers 1–5 only, Opus interactive): **~7-11 hours** of Opus agent
time. Orchestrator-deterministic units add ~30 min wall but no Opus budget.

With aggressive parallel dispatch (see §E), Tiers 2-5 could compress to ~3-4 wall-hours.

**Wall-time savings from M-6 reclassification**: ~1-2 hours of Opus budget moved to
shell scripts (formerly counted under Tier 5 Opus serial estimate).

---

## C. Atomic cleanup units

21 nominal units. The split of C02 into C02a-code/C02b/C02c yields **22 entries** in
the queue, but only **20 Opus-executed atoms** in cleanup proper (C02b + C02c are
deferred; four others are orchestrator-deterministic). Each unit carries the issue
IDs it resolves and the three-check criteria.

**Common preflight (M-2)** — every unit's first action is:
```
git status --porcelain
```
which MUST return empty before the agent begins edits. If non-empty, the unit
returns `BLOCK` and waits for orchestrator triage.

### C01 — Kaon-constants single-source-of-truth refactor (R03-I1 + R03-I2)

**Tier**: 1
**Execution lane**: Opus (single session, two commits).
**Files**: `quarkConstraints/modern/phenomenology.py:421-431`,
`tests/test_quark_deltaf2.py:111-117`, `tests/test_modern_phenomenology.py` (exists,
249 lines — see CR-2), `tests/test_modern_scan.py` (exists, 474 lines — see CR-2),
`artifacts/collaborator_5tev_points.csv`,
`artifacts/collaborator_5tev_points.provenance.json`,
`artifacts/collaborator_direct_affine_5_10tev_points.csv`,
`artifacts/collaborator_direct_affine_5_10tev_points.provenance.json`,
`docs/audits/epsilon_k_sm_decision.md`, `docs/phase_logs/phase2_h5_impl.md`.

**Pre-flight (M-1)**: run
`git grep -n "0.717\|2.228e-3\|1.81e-3\|F_K\s*=\s*0.155" quarkConstraints/paper_0710_1869/`
to confirm `paper_0710_1869/` does **not** carry duplicate kaon hadronic literals
(planner has verified at round-2: only `default_kaon` artifact IDs in
`paper_0710_1869/artifacts.py`, no `B_1/B_4/B_5/F_K/M_K` literal duplicates). If
the grep returns hits, expand C01 scope to include those files; otherwise proceed.

**Action — Commit C01a (bugfix + tests, RC-2)**:
Replace the 11 module-private `_KAON_*` literals at lines 421-431 with:
```python
from ..deltaf2 import (
    F_K as _KAON_F_K,
    M_K as _KAON_M_K,
    DELTA_M_K as _KAON_DELTA_M_K,
    M_S_2GEV as _KAON_M_S_2GEV,
    M_D_2GEV as _KAON_M_D_2GEV,
    B_1_K as _KAON_B_1,
    B_4_K as _KAON_B_4,
    B_5_K as _KAON_B_5,
    KAPPA_EPSILON as _KAON_KAPPA_EPSILON,
    EPSILON_K_EXP as _KAON_EPSILON_K_EXP,
    EPSILON_K_SM as _KAON_EPSILON_K_SM,
)
```
Place at top of file with other relative imports (not function-body) — same hygiene
as R04-I2. Extend
`test_audited_deltaf2_hadronic_constants_match_selected_sources` to assert byte-
identity (no `np.isclose`) for all 11 imported names.

**Action — Commit C01b (collaborator artifact regen, RC-2)**:
Re-run `python scripts/export_collaborator_5tev_points.py` and
`python scripts/export_collaborator_direct_affine_points.py`. Commit the resulting
CSV + provenance-JSON diffs. Append to each provenance JSON `changelog`:
"post-R03-I1 constant-correction; ε_K NP/budget ratios are ~6.24× tighter than
pre-fix CSV (4.18e-4 → 6.7e-5 budget)".

**Three checks**:
1. *Fix resolves issue*: `git grep "_KAON_B_1 = 0.717"` returns no results anywhere
   in `quarkConstraints/`.
2. *No regression* (**CR-2 explicit list**):
   ```
   pytest tests/test_quark_deltaf2.py tests/test_modern_phenomenology.py tests/test_modern_scan.py
   ```
   All three files MUST be in the selection. `tests/test_modern_phenomenology.py`
   and `tests/test_modern_scan.py` are both structural — they exercise
   `build_modern_point_phenomenology_artifact` and the monkeypatch stub respectively,
   but do **not** pin specific ε_K numerical ratios. Therefore **neither is expected
   to require numerical update**; both must still pass unchanged. If either fails
   numerically, escalate to BLOCK (the cleanup has uncovered an unexpected
   consumer). Also: `test_default_benchmark_point_has_stable_deltaf2_outputs` runs
   the LEGACY path (not modern) and must remain at `epsilon_k_ratio ≈ 1.9286`.
3. *Physics correctness*: re-derive `epsilon_K^NP / budget` shift; expected factor
   `(2.228e-3 − 1.81e-3) / (2.228e-3 − 2.161e-3) ≈ 6.24×` (closer to budget).
   Confirm `_evaluate_delta_mk_from_bridge` regression directions in §D1 footnote.

### C02a-code — Hole-#5 follow-ups (in-scope portion only)

**Tier**: 1
**Execution lane**: Opus (single session).
**Files**: `scripts/run_rs_anarchy.py` (add `--epsilon-k-budget {central,low,high}`
CLI flag; default = `central` so existing scans reproduce exactly),
`tests/test_run_rs_anarchy.py` or equivalent (extend with one parametrized test
across the three budget options),
`docs/quark_scan_methodology_note.tex` (insert band-quote paragraph using the
**symbolic** `1/sqrt(|Δε_K|)` scaling already documented in
`phase2_h5_signoff.md:100-101` — this avoids dependence on SLURM scan outputs),
`docs/artifact_manifest.md` (record three new scan-output **placeholder** stanzas
referencing C02c's planned output dirs).
**Scope (CR-3 + G1 hybrid)**: ONLY the CLI flag, tests, methodology paragraph,
and manifest placeholders. **No SLURM scans in this atom.** Substantive scan
output, ingest, and paragraph refinement are deferred to **C02c**.

**Three checks**:
1. *Fix resolves issue*: `phase2_h5_signoff.md:104-136` tasks 1, 3, 4 are marked
   `RESOLVED-BY: C02a-code`; task 2 is annotated `DEFERRED-TO: C02c`; task 5
   (figure) is annotated `DEFERRED-TO: C02b` (paper finalization).
2. *No regression*: `pytest tests/test_run_rs_anarchy.py` (or whatever test file
   exercises the CLI) passes with the new flag at default value; no existing scan
   output paths change.
3. *Physics correctness*: confirm the band paragraph cites the documented
   `M_KK^min(1e-5)/M_KK^min(6.7e-5) ≈ 2.59` and `M_KK^min(3e-4)/M_KK^min(6.7e-5) ≈ 0.47`
   ratios with the explicit `1/sqrt(6.7e-5 / Δε_K)` formula.

### C02b — Sensitivity-band figure (deferred-but-tracked)

**Tier**: 6 (deferred — paper finalization).
**Status**: `deferred`. Not executed during this cleanup. Recorded in
`CLEANUP_QUEUE.md` row so it is not lost.

### C02c — SLURM RUNA reruns + ingest + paragraph refinement (deferred-but-tracked)

**Tier**: 6 (deferred — SLURM queue / user attention).
**Status**: `deferred`.
**Files (when run)**: 3 new dated dirs under `scan_outputs/` (one per ε_K-budget
edge: central 6.7e-5, low ~1e-5, high ~3e-4); update
`docs/quark_scan_methodology_note.tex` band paragraph with **measured** M_KK^min
values replacing the symbolic estimates; update `docs/artifact_manifest.md` to
populate the three placeholder stanzas.
**Why deferred**: 3 SLURM submissions × 30–60 min queue+run each can leave an
Opus session idle for 2–6 wall-hours, burning budget. The hybrid disposition
(G1) leaves C02a-code completable in one Opus session and lets the user trigger
C02c on his timeline. C02c is NOT a blocker for `cleanup-tier1-complete` or
`v2026q2-post-cleanup`.
**Three checks (when run)**: same skeleton as C02a-code; check-3 hand-verifies
one band-edge against the symbolic formula.

### C03 — Wilson-RG audit hardening (R04-I1 + R04-I2 + R04-I3)

**Tier**: 1
**Execution lane**: Opus.
**Files**: `scripts/audit_wilson_rg.py:85-99`, `tests/test_wilson_rg_audit.py`
(new test); `quarkConstraints/deltaf2.py:586` (**unconditional** import lift,
per M-7 — verified non-circular with `qcd_running.py`);
`docs/audits/wilson_rg_inventory.md`.
**Action**: (1) Add `assert _GAMMA_LR[1, 0] == 0.0, "closed-form LR shortcut
requires upper-tri ADM"` at top of `scalar_lr_segment_matrix`, plus a test that
imports both and asserts the same. (2) **Lift `evolve_deltaf2_wilsons` import to
top of `deltaf2.py`**. Reviewer (M-7) confirmed `qcd_running.py` does not import
from `deltaf2.py`, so the lift is unconditionally safe. The previous "try to lift;
retain function-local if circular" wording is dropped. (3) Confirm
`paper_0710_1869/eft_deltaf2/rg.py` BMU LR map sign agrees with
`quarkConstraints/qcd_running.py`; record in `docs/audits/wilson_rg_inventory.md`.

**Three checks**:
1. *Fix resolves issue*: assertions present; test added; doc note added; import
   at top of file.
2. *No regression*: `pytest tests/test_wilson_rg_audit.py tests/test_qcd_running.py`
   passes; `python scripts/audit_wilson_rg.py` exit-0.
3. *Physics correctness*: re-execute audit, 22.49× cumulative factor still holds.

### C04 — Test hardening: Wilson UL + website notation (R07-I2 + R22-I1)

**Tier**: 2
**Execution lane**: Opus.
**Files**: `tests/test_finite_stats.py` (add parametrized k>0 cases),
`flavor_catalog/website/src/lib/notation.test.ts`,
`flavor_catalog/website/src/lib/prose.test.ts`,
`flavor_catalog/website/src/lib/__tests__/__snapshots__/` (vitest snapshot dir —
RC-1).
**Action — Python (RC-1)**: parametrize `k ∈ {1, 5, 10, 50}` and `n ∈ {100, 1000}`.
Compare against either (a) `scipy.stats.binomtest(k, n).proportion_ci(method='wilson')`
**requires scipy>=1.10** — pin in `pyproject.toml` if not already, OR (b) a
**closed-form Wilson UL fallback**:
```
z = 1.96
p_hat = k / n
denom = 1 + z**2 / n
center = (p_hat + z**2 / (2*n)) / denom
margin = z * sqrt(p_hat*(1-p_hat)/n + z**2/(4*n**2)) / denom
upper = center + margin
```
Prefer (a) if `scipy>=1.10` is already a dep; otherwise embed (b) in the test
file with a `# Wilson 95% CL, z=1.96, closed form` comment for the k>0 reference
and a separate `# Wilson 92.5% CL, z=1.92` for the codebase's choice (R07-I1).
**Action — TypeScript**: vitest snapshot tests for 5 entries (T010, CR002, K018,
B015, K020) at the snapshot dir cited above; commit snapshots.

**Three checks**: (1) tests added; (2) `pytest tests/test_finite_stats.py` and
`cd flavor_catalog/website && npx vitest run` both pass; (3) hand-verify one
k=10, n=100 case (closed-form gives upper ≈ 0.1740 at z=1.96).

### C05 — Spot-check f_{u,2}=0.16 reproduces from documented c-pattern (R08-I2)

**Tier**: 2
**Execution lane**: Opus.
**Files**: NEW `docs/audits/rc1p1_f_factor_check.md`.
**Action**: document `f_{u,2} = exp(c_{Q,2} k π R) exp(c_{u,2} k π R)` at the
canonical `ε = M_KK / M_Pl` from the methodology note's eq.(3) and the c-values
used in the WARNING-3 fix. Reproduce 0.16 to two sig figs.
**Three checks**: (1) note exists; (2) no code change → N/A; (3) physics
correctness IS the check (hand-calc reproduction).

### C06 — Wave-1 Kaon yaml polish (R10a-I1 + R10a-I2 + R10a-I3)

**Tier**: 4
**Execution lane**: Opus.
**Files**: `flavor_catalog/processes/kaon/K001.yaml`,
`flavor_catalog/processes/kaon/K003.yaml`.
**Dependency**: must run after C01 (R10a-I2 cross-links to R03-I1 closure).
**Three checks**: (1) yaml fields updated (K001 checker_agent_id →
`ca_w23_kaon_charm_edm`; K001 open_issues block replaced with
`tracked_in: cleanup-C01 ; resolved`; K003 PDG = canonical, DataBlock =
`supporting_value`); (2) all three files parse via `python -c "import yaml;
yaml.safe_load(open('...'))"`; (3) N/A.

### C07 — Wave-1 Beauty yaml polish (R10b-I1 + R10b-I2 + R10b-I3)

**Tier**: 4
**Execution lane**: Opus.
**Files**: `flavor_catalog/processes/beauty/B011.yaml`, `B009.yaml`, and the 5
beauty yamls cited in ISSUES.md R10b-I3.
**Three checks**: (1) typos / missing `state:` fields fixed and headlines chosen;
(2) all yamls parse; (3) N/A.

### C08 — Wave-1 Top/Higgs/EW + EDM + Charged-lepton yaml polish (R10c-I1/I2/I3)

**Tier**: 4
**Execution lane**: Opus.
**Files**: NEW `flavor_catalog/processes/top_higgs_ew/T011.{yaml,tex}` (stub
redirect to T010); NEW `flavor_catalog/references/{T002,T010,E001}/sha256sums.txt`;
worklog cross-references for T002/E001/C001/L001.
**Three checks**: (1) files added; (2) `(cd <dir> && sha256sum -c sha256sums.txt)`
passes in each of three reference dirs; (3) N/A.

### C09 — Worklog/signoff/factcheck text drift sweep (11 issues)

**Tier**: 4. Bundles: R11-I1, R11-I2, R11-I3, R12-I2, R12-I3, R12-I4, R13-I2,
R13-I3, R15-I2, R15-I4, R16-I3.
**Execution lane**: Opus.
**Files**: `flavor_catalog/worklogs/discovery/round_001_full_scope.md`,
`round_002_index.md`, `round_003_final_sweep.md`, `round_004_convergence.md`;
E009.yaml; B017/B025 worklog notes.
**Three checks**: (1) every cited file has a documented edit; (2) every yaml/md
still parses; (3) N/A.

### C10 — factcheck_status.md consolidation hardening (R17-I2 + R18-I1 + R18-I3 + R19-I3 + R19-I4)

**Tier**: 3
**Execution lane**: Opus.
**Disposition** (G3, lean B): per-master-compile annotation; optional aggregate
script. Both implemented.
**Files**: `flavor_catalog/audits/factcheck_status.md`,
`flavor_catalog/master_compile_v0{2,3,4,5+}_report.md` (add "consolidation status"
line per report naming addendum files for any counts not yet in consolidated
table); OPTIONAL `tools/aggregate_factchecks.py`.
**Three checks**: (1) HEAD `factcheck_status.md` row count matches wave totals
(102 entries through Wave-9 collider_rs); (2) per-family addenda counts sum to
the consolidated total; (3) N/A.
**Note**: v0.3 and v0.4 git tag annotations are immutable; C10 does NOT rewrite
history.

### C11 — SESSION_NOTES + HANDOFF + WEBSITE_RUNBOOK alignment (R20-I1/I2/I3 + R21-I1/I3 + R22-I2/I3)

**Tier**: 5
**Execution lane**: Opus.
**Files**: `flavor_catalog/SESSION_NOTES.md`, `flavor_catalog/HANDOFF_PROMPT.md`,
`flavor_catalog/AGENTIC_WORKFLOW.md`, `flavor_catalog/website/WEBSITE_RUNBOOK.md`,
`flavor_catalog/website/src/components/EntryTable.astro`,
`flavor_catalog/website/src/pages/browse.astro`,
`flavor_catalog/website/src/pages/methodology.astro`.
**Three checks**: (1) edits land; (2) `cd flavor_catalog/website && npm run build`
passes; (3) N/A.

### C12 — Quark mass constants single source-of-truth (R01-I1)

**Tier**: 5
**Execution lane**: Opus.
**Files**: `qcd/constants.py:17-19`, `qcd/mass_running.py:106-108`, relevant test.
**Three checks**: (1) `git grep "M_CHARM\s*="` returns only the PDG dict line;
(2) `pytest -k mass_running` passes; (3) document explicit drift bound (no
benchmark moves at < 1e-3 drift).

### C13 — Spurion-seed provenance comment (R02-I2)

**Tier**: 5
**Execution lane**: Opus.
**Files**: `quarkConstraints/benchmarks.py:189-223` (docstring only).
**Note** (RC-7): **no test re-run required** because no values change. Cleanup
agent skips Check-2 pytest invocation and explicitly notes "no functional change;
tests not re-run by design".
**Three checks**: (1) comment landed; (2) `python -c "import quarkConstraints.benchmarks"`
imports cleanly; (3) N/A (docs-only).

### C14 — external_research provenance manifest (R17-I1)

**Tier**: 3
**Execution lane**: Opus.
**Files**: NEW `flavor_catalog/external_research/MANIFEST.md`.
**Three checks**: (1) file exists; (2) `sha256sum -c` against recorded digests
OK; (3) N/A.

### C15 — Backfill `.orchestration/pytest_selection/*.txt` (R05-I2)

**Tier**: 5
**Execution lane**: **orchestrator-deterministic** (M-6).
**Action** (executed inline by orchestrator shell, NOT an Opus session):
```
for unit in R01 R02 R03 R04 R05 R06 R07 R08 R09 R10a R10b R10c \
            R11 R12 R13 R14 R15 R16 R17 R18 R19 R20 R21 R22 ; do
  git log --name-only <unit-commit-range> -- tests/ \
    | grep '^tests/' | sort -u \
    > .orchestration/pytest_selection/${unit}.txt
done
```
Opus reviews the resulting commit only for §E.2 compliance (24 non-empty files
written; spot-check by re-running `pytest $(cat .../R03.txt)`).
**Three checks**: (1) 24 files exist and are non-empty for code units; (2)
`pytest $(cat .orchestration/pytest_selection/R03.txt)` passes for the three
highest-risk units (R03, R04, R07); (3) N/A.

### C16 — CLAUDE.md + codex-model version touchups (R07-I4 + R14-I3) — **R17-I3 dropped**

**Tier**: 5
**Execution lane**: Opus.
**Files**: `CLAUDE.md:12-15` (paper-branch pointer → "main");
`docs/phase_logs/flavor_catalog_codex_quota_pause.md:46` (gpt-5.5 → gpt-5.4).

**CR-1 fix**: The v0.2 tag-creation task is **DROPPED**. Pre-flight verification
already performed by planner-round-2:
```
$ git tag -l 'flavor-catalog-v0*'
flavor-catalog-v0
flavor-catalog-v0.1
flavor-catalog-v0.2
flavor-catalog-v0.3
flavor-catalog-v0.4
$ git rev-list -n 1 flavor-catalog-v0.2
835cf48873e5addf89e58a8aec11cac8d9ee7837
```
The tag exists and points to the actual master-compile v0.2 commit (`835cf48`),
NOT the planner-r1-cited `42ac647` (which is unrelated). **R17-I3 is reclassified
as already-CLOSED**; C21 moves it to `## Closed / Accepted-risk` with rationale
"PRE-EXISTING-TAG-VERIFIED at 835cf48 on planner-r2 date".

**Three checks**: (1) `CLAUDE.md:12-15` and `flavor_catalog_codex_quota_pause.md:46`
edits land; (2) `git tag -l 'flavor-catalog-v0.*'` still shows v0.2 + v0.3 + v0.4
(unchanged); (3) N/A.

### C17 — `.gitkeep` cleanup (R09-I2)

**Tier**: 5
**Execution lane**: **orchestrator-deterministic** (M-6).
**Action** (executed inline by orchestrator shell): for each of the 11
`.gitkeep` files cited in ISSUES.md R09-I2, run
```
for f in <11 listed paths> ; do
  if [ -n "$(ls -A "$(dirname "$f")" 2>/dev/null | grep -v '^\.gitkeep$')" ]; then
    git rm "$f"
  fi
done
```
which preserves `.gitkeep` only in still-empty dirs (per G6 lean). Opus reviews
commit for: deletion count matches expected delta; remaining `.gitkeep`s are in
empty dirs.
**Three checks**: (1) `git ls-files | grep '\.gitkeep$'` count drops by expected
delta (likely 11 → 0 since the dirs cited have committed content); (2) directory
listings show real content; (3) N/A.

### C18 — MERGE_PLAN.md retroactive corrections (13 issues)

**Tier**: 5. Bundles: R02-I1, R03-I4, R05-I1, R06-I2, R09-I1, R12-I1, R13-I1,
R14-I1, R15-I1, R16-I4, R19-I1, R19-I2, R19-I5.
**Execution lane**: **orchestrator-deterministic** (M-6).
**Action** (RC-5): **never use absolute line numbers from ISSUES.md cite-text** —
MERGE_PLAN.md has been edited since those line refs were captured. Instead, for
each of the 13 issues, run a `grep -n` for the error pattern and then patch:
```
# Example for R15-I1 (SHA typo 7500919 → 7500794):
grep -n "7500919" .orchestration/MERGE_PLAN.md
# then sed -i 's/7500919/7500794/g' on the matching line(s)
```
The orchestrator scripts the 13 such `grep -n` + targeted-replace pairs. Opus
reviews the resulting diff for: (a) every cited error pattern is gone; (b) no
collateral edits; (c) `git diff` is restricted to `.orchestration/MERGE_PLAN.md`.
**Three checks**: (1) `grep -n` returns no hits for the 13 error patterns post-
fix; (2) pure docs → N/A; (3) N/A.

### C19 — Methodology-note polish (R04-I4, R05-I3, R06-I3, R07-I1)

**Tier**: 5
**Execution lane**: Opus.
**Files**: `docs/quark_scan_methodology_note.tex` (4 inserts), rebuild PDF.
**Three checks**: (1) `pdflatex` builds clean; (2) sha256 of PDF recorded in
`artifacts/checksums.sha256`; (3) N/A.

### C20 — Accept-as-risk / forward-only INFOs (R01-I2, R06-I1, R07-I3, R08-I3, R21-I2)

**Tier**: 5
**Execution lane**: **orchestrator-deterministic** (M-6).
**Action**: move 5 issues to `## Closed / Accepted-risk` with the rationales
below. Orchestrator script is a single `cat <<EOF >> ISSUES.md` block (Opus
reviews the resulting commit, doesn't write the edit).

**Per-issue rationale**:
- **R01-I2** (RC-6 footnote): `MASS_TOLERANCE_FLOOR=0.003` calibration is
  deferred. `scripts/calibrate_phase0.py` **exists** in the repo but has
  **not been re-run** since the post-audit kaon-constants update (i.e. since
  R03-I1 was first noted). Current floor 0.003 is conservative enough that no
  downstream benchmark fails; revisit next time Phase-0 is re-run.
- **R06-I1**: refactor CFW constants only when a 2nd-comparison driver is added
  (none planned).
- **R07-I3**: manifest `code_sha` is inferred from mtime; SLURM-side change;
  implemented next-time the dispatcher is touched.
- **R08-I3**: historical figure dirs (`quark_pre_audit_constants/`,
  `quark_baseline_800k/`) kept as paper-archive companion (G7 lean B).
- **R21-I2**: Cloudflare Pages dashboard-driven deployment; `wrangler.toml`
  not required.

**Three checks**: (1) ISSUES.md updated, 5 issues moved; (2) N/A; (3) N/A.

### C21 — Final ISSUES.md sweep + post-cleanup tag

**Tier**: 5
**Execution lane**: Opus.
**Files**: `.orchestration/ISSUES.md` (verify `## Open` empty or only Accepted-risk
+ the 7 no-ops + R17-I3 PRE-EXISTING-TAG-VERIFIED); add a "Cleanup phase summary"
section at top; new annotated git tag `v2026q2-post-cleanup` on the final cleanup
commit; update `progress.json` to add a `cleanup_phase` block; update `CLAUDE.md`
with a Phase-2 reference.
**Tag chronology (M-8)**: confirmed that `v2026q2-catalog-complete` exists today
(verified: `git tag -l 'v2026q2*'` returns only `v2026q2-catalog-complete`).
`v2026q2-post-cleanup` will sit chronologically AFTER `v2026q2-catalog-complete`
in the commit DAG (post-cleanup HEAD is descendant). User to confirm tag-signing
convention (annotated, unsigned — same as `v2026q2-catalog-complete`).
**Three checks**: (1) `## Open` count is the explicitly-Accepted set only; (2)
`git tag -l v2026q2-post-cleanup` returns the tag; (3) N/A.

---

## D. Detailed designs for the Tier-1 (high-stakes) fixes

### D1 — R03-I1: Duplicate pre-audit kaon constants (the only HIGH)

**What needs to change** (file:line, specific edits):

In `quarkConstraints/modern/phenomenology.py`, replace lines 421-431 with the
named import block in §C1 above. Place at top of file with other relative
imports. Verify no circular import:
`grep "from.*modern" quarkConstraints/deltaf2.py` returns empty.

**Tests that need updating**:

In `tests/test_quark_deltaf2.py`, extend
`test_audited_deltaf2_hadronic_constants_match_selected_sources`:
```python
from quarkConstraints.modern import phenomenology
assert phenomenology._KAON_B_1 is deltaf2.B_1_K     # identity, not equality
assert phenomenology._KAON_B_4 is deltaf2.B_4_K
assert phenomenology._KAON_B_5 is deltaf2.B_5_K
assert phenomenology._KAON_EPSILON_K_SM is deltaf2.EPSILON_K_SM
# ... plus the other 7 constants
```
`is` (identity) suffices because the import binds the same Python object.

**Tests already on disk that MUST stay green (CR-2)**:
- `tests/test_modern_phenomenology.py` (249 lines, exercises
  `build_modern_point_phenomenology_artifact`): structural assertions only — no
  ε_K ratio pins. **No update expected**, but MUST be in Check-2 selection.
- `tests/test_modern_scan.py` (474 lines, monkeypatch stub): structural. **No
  update expected**, but MUST be in Check-2 selection.
- `test_default_benchmark_point_has_stable_deltaf2_outputs` (LEGACY path):
  `epsilon_k_ratio ≈ 1.9286` MUST be unchanged.

**Numerical regression to expect** (with footnote per M-5):

`epsilon_K^NP / budget` (modern backend) shifts by factor:
```
budget_pre  = |EPSILON_K_EXP − EPSILON_K_SM_pre|  = |2.228e-3 − 1.81e-3|  = 4.18e-4
budget_post = |EPSILON_K_EXP − EPSILON_K_SM_post| = |2.228e-3 − 2.161e-3| = 6.7e-5
ratio_shift = budget_pre / budget_post            ≈ 6.24
```
i.e. the post-fix `ε_K^NP / budget` ratio is **~6.24× tighter** for the same NP
bridge.

**M-5 one-line derivation footnote** for the `_evaluate_delta_mk_from_bridge`
direction estimates:
```
ΔM_K ratio shift ≈ (B_4_post / B_4_pre) · (something) − 1 ;
  using B_1_pre=0.717 → B_1_post=0.5503; B_4_pre=0.78 → B_4_post=0.903;
  B_5_pre=0.57 → B_5_post=0.691.
For VLL-dominated bridges (Q1_VLL dominant), ΔM_K ∝ B_1, so
  shift ≈ 0.5503 / 0.717 − 1 ≈ −0.232  → "~−23%"
For LR-dominated bridges (Q4_LR + Q5_LR dominant, ~equal weight by κ_ε):
  shift ≈ ((B_4_post + B_5_post) / (B_4_pre + B_5_pre)) − 1
        ≈ (0.903 + 0.691) / (0.78 + 0.57) − 1
        ≈ 1.594 / 1.350 − 1 ≈ +0.181 → "~+15-18%"  (planner-r1 quoted +15%)
```
The footnote goes in the commit body for C01a so a reader of the diff sees the
expected magnitudes immediately.

**Collaborator artifacts** (RC-2 split):

`scripts/export_collaborator_5tev_points.py` and
`scripts/export_collaborator_direct_affine_points.py` are re-run in **commit
C01b** (separate from C01a bugfix). Provenance JSON gets a changelog entry as in
§C1.

**Documentation that needs updating**:
- `docs/audits/epsilon_k_sm_decision.md`: append "modern backend deduplicated;
  see C01 commit <sha>".
- `docs/phase_logs/phase2_h5_impl.md:71`: append
  "RESOLVED-BY: C01 ; previous out-of-scope declaration superseded".
- `flavor_catalog/processes/kaon/K001.yaml:153-155` open_issues block: handled
  in C06 (depends on C01).

**Risk + mitigation**:
- Risk: any consumer of `_evaluate_epsilon_k_from_bridge` outside the scan path
  changes. *Mitigation*: `git grep -l "_evaluate_epsilon_k_from_bridge"` + audit
  every caller before commit. Run extended pytest selection (CR-2).
- Risk: collaborator artifacts already shared by email. *Mitigation*: changelog
  in provenance JSON; ping Oscar for downstream notification.

### D2 — R03-I2: Pin-test extension

Same commit as C01a; the pin-test is the guard against R03-I1 re-introduction.

### D3 — R03-I3: Five hole-#5 follow-ups (hybrid scope per G1)

**In-scope (C02a-code)**:
1. Methodology-note band paragraph (using symbolic `1/sqrt(|Δε_K|)` scaling).
3. CLI flag (`--epsilon-k-budget`).
4. Methodology-note paragraph describing band construction.

**Deferred to C02c (G1 hybrid)**:
2. RUNA reruns at three ε_K-budget edges. (RC-4 reconciliation: there is NO
   "1 existing" — all three budget edges are NEW SLURM scans. The planner-r1
   "3 SLURM reruns" wording stands; the planner-r1 ambiguity is resolved here.)

**Deferred to C02b (paper finalization)**:
5. Optional sensitivity-band figure.

**Per-task numerical estimates** (used by C02a-code in the methodology paragraph):
```
M_KK^min(1e-5)   / M_KK^min(6.7e-5) ≈ sqrt(6.7e-5 / 1e-5)  ≈ 2.59
M_KK^min(3e-4)   / M_KK^min(6.7e-5) ≈ sqrt(6.7e-5 / 3e-4)  ≈ 0.47
```
So the band-text is approximately "M_KK^min ≈ X^{+1.6X}_{-0.5X} TeV" with X the
central value (~16.5 TeV per the p50 pert quote).

---

## E. Per-unit execution protocol

### E.1 Cleanup-agent contract

**Two lanes**:
1. **Opus-executed** units (default): C01, C02a-code, C03, C04, C05, C06, C07,
   C08, C09, C10, C11, C12, C13, C14, C16, C19, C21.
2. **Orchestrator-deterministic** units (M-6): C15, C17, C18, C20. Orchestrator
   shell executes the action; Opus reviews the resulting commit only.

**Inputs provided by orchestrator to Opus agents**:
1. Unit ID + title + §C section
2. Issue IDs the unit resolves
3. Current state of `.orchestration/ISSUES.md`
4. Pointer to this plan
5. Path for the report: `.orchestration/cleanup_reports/C##.md`
6. Reminder of write access

**Inputs provided to orchestrator for deterministic units**:
1. The exact shell script for the unit (one of §C15/C17/C18/C20 above)
2. Path for the report (same)

**Outputs the cleanup agent must produce**:
- One git commit per atomic unit (C01 is the only two-commit unit: C01a + C01b)
- One report file `.orchestration/cleanup_reports/C##.md`
- An update to `.orchestration/ISSUES.md` moving resolved issues to
  `## Closed / Accepted-risk`
- An update to `.orchestration/CLEANUP_QUEUE.md` flipping unit status

### E.2 The three checks

Same semantics as Phase-1 reviewer rounds, with cleanup-scoped wording.

#### E.2.1 Did the fix actually resolve the issue?

Re-read original issue from ISSUES.md. Verdict: PASS | PARTIAL | FAIL.

#### E.2.2 Did the fix break anything else?

For Opus code-touching units: run per-unit pytest selection (e.g. for C01 the
selection is the **explicit** triple `tests/test_quark_deltaf2.py
tests/test_modern_phenomenology.py tests/test_modern_scan.py` per CR-2). Document
wall time. PARTIAL if > 10 min (same MERGE_PLAN policy).

For docs-only units: `python -c "import yaml; yaml.safe_load(open('...'))"` on
every edited yaml + `pdflatex` on every edited `.tex` if feasible in < 5 min.

For C13 (RC-7): explicitly note "no test re-run required; pure-docstring change
with no value drift".

For orchestrator-deterministic units: Opus's Check-2 is reviewing the diff
against the script's intent (no collateral edits; expected files only).

#### E.2.3 Does the fix maintain physics correctness?

Skip for pure-docs. Mandatory for code change. C01's Check-3 cites the §D1
"Numerical regression to expect" calculation (including the M-5 footnote).

### E.3 Report template

Same as planner-r1's template, additionally noting the unit's execution lane
("Opus" or "orchestrator-deterministic; reviewed by Opus").

### E.4 Commit policy

One git commit per atomic unit C##, with **C01 as the explicit two-commit
exception** (C01a bugfix + C01b artifact regen per RC-2). C02a-code is a single
commit (CLI flag + tests + methodology paragraph + manifest placeholders).

C18 is a single commit even though it touches 13 issues (all in MERGE_PLAN.md).

**No `--amend`** of existing history. **No force-push**.

Commit-message format:
```
cleanup(C##): <one-line title>

Resolves: <R##-I# list>
Plan: .orchestration/CLEANUP_PLAN.md §C C##
Report: .orchestration/cleanup_reports/C##.md

<2-4 line body; for C01a, include the M-5 derivation footnote.>
```

### E.5 ISSUES.md update protocol

When a unit closes an issue (M-3 organization):
1. Locate issue under `## Open`.
2. Move the entire issue block under `## Closed / Accepted-risk`, **organized by
   cleanup unit C##** (not chronological). Within each unit subsection, issues
   are listed in their original R##-I# order. This gives reviewers a clean
   "what did C04 close?" view by scrolling to the C04 subsection.
3. Prepend the unit subsection header on first close into it:
   ```
   ### Closed by C04
   ```
4. Prepend each issue block with:
   ```
   - CLOSED <YYYY-MM-DD> by C## (commit <short-sha>):
   ```
   For ACCEPTED-RISK (C20 only):
   ```
   - ACCEPTED-RISK <YYYY-MM-DD> rationale: <one-line>
   ```
5. Verify `## Open` count decreases by expected delta.
6. The cleanup agent does NOT touch issues outside its own unit.

### E.6 Tier-1 checkpoint (M-4)

After C01, C02a-code, and C03 land, the orchestrator (NOT an Opus agent) creates
an annotated tag:
```
git tag -a cleanup-tier1-complete -m "Cleanup Tier 1 (physics-affecting) complete: C01, C02a-code, C03"
git push origin cleanup-tier1-complete
```
This is the only tag created mid-cleanup. **Tier-2..5 do NOT need user gating**.
After Tier-5 lands, C21 creates `v2026q2-post-cleanup`.

---

## F. Tracking and resumability

`.orchestration/` layout, post-cleanup-phase (unchanged from r1 except as noted):

```
.orchestration/
├── MERGE_PLAN.md
├── PLAN_REVIEW_R{1,2}.md
├── PLAN_CHANGES_R2.md
├── PRE_MERGE_*.md
├── REFERENCES.md
├── ROLLBACK.md
├── REVIEW_QUEUE.md
├── COMMIT_INDEX.csv
├── reviews/
├── pytest_selection/           # backfilled by C15 (24 files)
├── ISSUES.md                   # SHARED — Phase 1 + Phase 2; M-3 organized by C##
├── progress.json               # Phase 1 v2 + Phase 2 cleanup_phase block
├── CLEANUP_PLAN.md             # THIS FILE (r2)
├── CLEANUP_REVIEW_R{1,2}.md
├── CLEANUP_CHANGES_R2.md       # NEW r2 (this revision's changelog)
├── CLEANUP_QUEUE.md            # NEW — status board incl. C02b/C02c deferred
├── cleanup_progress.json       # NEW — Phase-2 progress mirror
└── cleanup_reports/
    ├── C01a.md / C01b.md       # C01 has two commits → two reports OR one report with both SHAs
    ├── C02a-code.md
    ├── C03.md
    ├── ...
    └── C21.md
```

### F.1 `CLEANUP_QUEUE.md` format

```
| Unit | Title | Tier | Lane | Status | Started | Done | Commit | Report |
|------|-------|------|------|--------|---------|------|--------|--------|
| C01 | Kaon-constants SSOT + collab artifacts | 1 | opus | pending | — | — | — | — |
| C02a-code | Hole-#5 CLI + symbolic band paragraph | 1 | opus | pending | — | — | — | — |
| C02b | Sensitivity-band figure | 6 | deferred | deferred | — | — | — | — |
| C02c | SLURM RUNA reruns + ingest | 6 | deferred | deferred | — | — | — | — |
| C03 | Wilson-RG audit hardening | 1 | opus | pending | — | — | — | — |
| C04 | Wilson UL + website notation tests | 2 | opus | pending | — | — | — | — |
| C05 | f_u,2 spot-check | 2 | opus | pending | — | — | — | — |
| C06 | Kaon yaml polish | 4 | opus | pending (dep:C01) | — | — | — | — |
| C07 | Beauty yaml polish | 4 | opus | pending | — | — | — | — |
| C08 | Top/EW/EDM/CL yaml polish | 4 | opus | pending | — | — | — | — |
| C09 | Worklog text drift sweep | 4 | opus | pending | — | — | — | — |
| C10 | factcheck_status consolidation | 3 | opus | pending | — | — | — | — |
| C11 | SESSION_NOTES + RUNBOOK | 5 | opus | pending | — | — | — | — |
| C12 | Quark mass SSOT | 5 | opus | pending | — | — | — | — |
| C13 | Spurion-seed comment | 5 | opus | pending | — | — | — | — |
| C14 | external_research MANIFEST | 3 | opus | pending | — | — | — | — |
| C15 | pytest_selection backfill | 5 | shell | pending | — | — | — | — |
| C16 | CLAUDE.md + codex version | 5 | opus | pending | — | — | — | — |
| C17 | .gitkeep cleanup | 5 | shell | pending | — | — | — | — |
| C18 | MERGE_PLAN retro corrections | 5 | shell | pending | — | — | — | — |
| C19 | Methodology-note polish | 5 | opus | pending | — | — | — | — |
| C20 | Accept-as-risk INFOs | 5 | shell | pending | — | — | — | — |
| C21 | Final sweep + tag | 5 | opus | pending (last) | — | — | — | — |
```

### F.2 `cleanup_progress.json` schema

```
{
  "schema_version": "v1",
  "phase": "cleanup",
  "queue": [
    {"unit": "C01", "status": "pending", "tier": 1, "lane": "opus",
     "resolves": ["R03-I1", "R03-I2"], "commits": 2},
    {"unit": "C02a-code", "status": "pending", "tier": 1, "lane": "opus",
     "resolves": ["R03-I3"], "partial_close": true},
    {"unit": "C02b", "status": "deferred", "tier": 6, "lane": "deferred",
     "resolves": ["R03-I3 task 5"]},
    {"unit": "C02c", "status": "deferred", "tier": 6, "lane": "deferred",
     "resolves": ["R03-I3 task 2"]},
    {"unit": "C03", "status": "pending", "tier": 1, "lane": "opus",
     "resolves": ["R04-I1", "R04-I2", "R04-I3"]},
    ...
    {"unit": "C15", "status": "pending", "tier": 5, "lane": "shell", ...},
    {"unit": "C17", "status": "pending", "tier": 5, "lane": "shell", ...},
    {"unit": "C18", "status": "pending", "tier": 5, "lane": "shell", ...},
    {"unit": "C20", "status": "pending", "tier": 5, "lane": "shell", ...},
    {"unit": "C21", "status": "pending", "tier": 5, "lane": "opus", ...}
  ],
  "tier1_checkpoint_tag": "cleanup-tier1-complete",
  "final_tag": "v2026q2-post-cleanup",
  "dropped_tasks": ["C16 v0.2 tag-creation (R17-I3): tag already exists at 835cf48"]
}
```

### F.3 ISSUES.md reuse (M-3 organization)

Per user requirement, ISSUES.md is the single source of truth. The `## Closed /
Accepted-risk` section is organized **by cleanup unit** (subsection per `### Closed
by C##`), with issues in original ordering within each subsection. A separate
`### Auto-resolved / no-op` subsection holds the 7 no-op issues + R17-I3
PRE-EXISTING-TAG-VERIFIED.

---

## G. Open-question dispositions (reviewer's calls, adopted)

The reviewer-r1 round adjudicated planner-r1's open questions. The dispositions
below are **the planner-r2 decisions** — open questions are now closed.

### G1. Scope of R03-I3 paper-touching tasks
**Adopted: hybrid (A').** CLI flag + methodology band paragraph (symbolic
`1/sqrt(|Δε_K|)` scaling, citable from `phase2_h5_signoff.md:100-101`) land in
**C02a-code** (single Opus session). Actual SLURM RUNA reruns deferred to
**C02c** (tracked as `deferred` in the queue; not blocking
`cleanup-tier1-complete` or `v2026q2-post-cleanup`).

### G2. Re-running + committing collaborator export artifacts
**Adopted: lean (A) split into two commits.** **C01a** = bugfix + tests; **C01b**
= collaborator artifact regen. Both inside the same Opus session, but separate
commits for clean revert/cherry-pick.

### G3. Audit-doc count-drift issues (C10 scope)
**Adopted: lean (B).** Per-master-compile annotation in each
`master_compile_v0{N}_report.md`, plus optional `tools/aggregate_factchecks.py`
that the next master-compile may use to regenerate the table automatically.

### G4. The (no-op) issues
**Adopted: lean (A) with corrected count.** **7 no-op issues** (RC-3 off-by-one
fix): R08-I1, R14-I2, R15-I3, R16-I1, R16-I2, R18-I2, R18-I4. All moved to
ISSUES.md `### Auto-resolved / no-op` subsection by C21 with one-line rationale
each.

### G5. Should C16 retroactively tag `flavor-catalog-v0.2`?
**Adopted: dropped (CR-1).** Tag already exists at `835cf48` (verified planner-r2).
R17-I3 reclassified as closed without any cleanup action; C21 moves it into
`### Auto-resolved / no-op` with rationale "PRE-EXISTING-TAG-VERIFIED at 835cf48".

### G6. `.gitkeep` cleanup scope
**Adopted: planner's lean (drop only redundant ones).** C17 keeps `.gitkeep` in
still-empty dirs.

### G7. Historical figure dir cleanup (R08-I3)
**Adopted: planner's lean (B).** Accept as paper-archive companion; closure
recorded in C20.

---

## H. Final repo state after cleanup

After C01..C21 (excluding C02b/C02c which are deferred) land, the expected
end-state is:

1. **`.orchestration/ISSUES.md`**:
   - `## Open`: empty.
   - `## Closed / Accepted-risk`: organized by `### Closed by C##` subsections;
     78 issues accounted for (70 cleanup-closed across 21 units; 5 accepted-risk
     under C20; 7 no-op + 1 PRE-EXISTING-TAG-VERIFIED under `### Auto-resolved /
     no-op`).
   - `## Infra follow-ups`: unchanged (`INFRA-1` already closed).

2. **Physics-affecting code changes** (C01, C02a-code, C03):
   - `quarkConstraints/modern/phenomenology.py` uses canonical post-audit kaon
     constants via a single import.
   - `tests/test_quark_deltaf2.py` pins modern duplicates against canonical.
   - `tests/test_modern_phenomenology.py` and `tests/test_modern_scan.py`
     unchanged (structural assertions only) — verified still green post-fix.
   - `scripts/run_rs_anarchy.py` has `--epsilon-k-budget {central,low,high}`
     CLI flag at default `central`.
   - `docs/quark_scan_methodology_note.tex` has a symbolic-band paragraph
     (citing `1/sqrt(|Δε_K|)`) and a Wilson z=1.92 footnote; PDF rebuilt.
   - `scripts/audit_wilson_rg.py` carries upper-tri ADM assertion;
     `tests/test_wilson_rg_audit.py` gates it.
   - Collaborator artifacts re-exported in C01b.

3. **Tests**: `tests/test_finite_stats.py` has k>0 parametrized cases (scipy>=1.10
   or closed-form Wilson UL fallback); website vitest snapshots at
   `flavor_catalog/website/src/lib/__tests__/__snapshots__/`.

4. **Documentation**: `SESSION_NOTES.md`, `HANDOFF_PROMPT.md`,
   `AGENTIC_WORKFLOW.md`, `WEBSITE_RUNBOOK.md`, `CLAUDE.md`, `MERGE_PLAN.md` (via
   C18 shell script using `grep -n` patterns rather than absolute line numbers),
   all aligned to v0.4+ state.

5. **Catalog**: yaml metadata cleanups across all six families;
   `references/{T002,T010,E001}/sha256sums.txt` exist;
   `flavor_catalog/external_research/MANIFEST.md` exists; `.gitkeep` cleanup
   complete (orchestrator-deterministic).

6. **Audit doc state**: `factcheck_status.md` consolidated table includes all
   102 entries; per-master-compile annotations added (G3 lean B); optional
   `tools/aggregate_factchecks.py` script in repo.

7. **Tags** (M-8):
   - `flavor-catalog-v0.2`, `v0.3`, `v0.4`: all already exist; untouched.
   - `cleanup-tier1-complete`: annotated, created post-C03, pushed to origin
     (M-4).
   - `v2026q2-post-cleanup`: annotated, created by C21, chronologically after
     `v2026q2-catalog-complete` (M-8).
   - DEFERRED `v2026q2-paper-finalized` for whenever C02b + C02c close.

8. **Tracking files**: `.orchestration/CLEANUP_QUEUE.md` shows 22 units (20 done,
   2 deferred); `cleanup_progress.json` mirrors; reports under
   `.orchestration/cleanup_reports/` (20 per-unit reports; C01 has two reports
   for C01a and C01b).

---

## Appendix: full unit → issue mapping (cross-check)

| Unit | Lane | Issues resolved | Count |
|------|------|-----------------|-------|
| C01  | opus | R03-I1, R03-I2 | 2 |
| C02a-code | opus | R03-I3 (tasks 1, 3, 4 + CLI) | 1 (partial) |
| C02b | deferred | R03-I3 (task 5) | 0 in cleanup |
| C02c | deferred | R03-I3 (task 2 SLURM) | 0 in cleanup |
| C03  | opus | R04-I1, R04-I2, R04-I3 | 3 |
| C04  | opus | R07-I2, R22-I1 | 2 |
| C05  | opus | R08-I2 | 1 |
| C06  | opus | R10a-I1, R10a-I2, R10a-I3 | 3 |
| C07  | opus | R10b-I1, R10b-I2, R10b-I3 | 3 |
| C08  | opus | R10c-I1, R10c-I2, R10c-I3 | 3 |
| C09  | opus | R11-I1, R11-I2, R11-I3, R12-I2, R12-I3, R12-I4, R13-I2, R13-I3, R15-I2, R15-I4, R16-I3 | 11 |
| C10  | opus | R17-I2, R18-I1, R18-I3, R19-I3, R19-I4 | 5 |
| C11  | opus | R20-I1, R20-I2, R20-I3, R21-I1, R21-I3, R22-I2, R22-I3 | 7 |
| C12  | opus | R01-I1 | 1 |
| C13  | opus | R02-I2 | 1 |
| C14  | opus | R17-I1 | 1 |
| C15  | shell | R05-I2 | 1 |
| C16  | opus | R07-I4, R14-I3 | 2 (was 3 in r1; R17-I3 dropped per CR-1) |
| C17  | shell | R09-I2 | 1 |
| C18  | shell | R02-I1, R03-I4, R05-I1, R06-I2, R09-I1, R12-I1, R13-I1, R14-I1, R15-I1, R16-I4, R19-I1, R19-I2, R19-I5 | 13 |
| C19  | opus | R04-I4, R05-I3, R06-I3, R07-I1 | 4 |
| C20  | shell | R01-I2, R06-I1, R07-I3, R08-I3, R21-I2 | 5 |
| C21  | opus | R17-I3 (PRE-EXISTING-TAG closure) + final sweep + tag | 1 (bookkeeping) |
| no-op | — | R08-I1, R14-I2, R15-I3, R16-I1, R16-I2, R18-I2, R18-I4 | 7 |

**Total**: 2+1+3+2+1+3+3+3+11+5+7+1+1+1+1+2+1+13+4+5+1+7 = **78** ✓

Cross-check matches the 78 open issues in ISSUES.md.
