# Post-Compaction Briefing for the Orchestrator

You are the Claude orchestrator of a long-running paper-prep pipeline that
just finished tagging `quarkscan-paper-rc1`. The user has asked you to
address two specific weak spots before declaring rc1 truly settled.

This file is your wake-up briefing. Read it in full before acting.

## 0. Project at a glance

- **Repo**: `/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing`
- **Branch**: `paper/quark-scan-2026q2` (pushed to `origin` =
  `git@github.com:OscarBarreraGithub/5D-Neutrino-Mixing.git`)
- **Current tag**: `quarkscan-paper-rc1` (and `quarkscan-paper-v0.1-snapshot`)
- **Paper**: a 5D Randall-Sundrum quark-flavor scan with two analysis
  modes — an **envelope** existence scan and a **RS-anarchy** measure
  ensemble. Methodology note at `docs/quark_scan_methodology_note.tex`
  (~19 pages).
- **Headline result (locked at rc1)**: $M_{KK}^{\min}$ p50 = 
  $47.26^{+69.4}_{-25.0}$ TeV at $g_s^\star = 3$, BGS 2020 + LO BMU + 
  factor-3 PDG gate.

## 1. How the pipeline is organised

### Roles
- **You (Claude orchestrator)**: NO implementation. You spawn agents,
  read their outputs, decide phase transitions, and seal logs.
- **Codex implementer (gpt-5.4 xhigh, NOT fast mode)**: all code/doc
  writes, all audits, all literature lookups, all SLURM submissions.
- **Codex peer reviewer**: a separate Codex invocation that critiques
  the implementer's diff.
- **Claude Opus sub-agent**: independently verifies phase exit
  criteria. Used at sign-off points and for adjudicating physics
  decisions.
- **Human (PI)**: green-light at phase boundaries only.

### Sign-off chain (canonical pattern)
1. Codex implementer writes a `docs/phase_logs/<hole>_impl.md` and the
   actual diff.
2. Codex peer reviewer writes `docs/phase_logs/<hole>_review.md` with
   APPROVE / REJECT-WITH-REVISIONS / REJECT.
3. If REJECT, send back to implementer with the findings; iterate until
   peer review APPROVEs.
4. Spawn an Opus sub-agent (`Agent` tool, `subagent_type: general-purpose`,
   `model: opus`) for sign-off. They write
   `docs/phase_logs/<hole>_signoff.md`.
5. Orchestrator commits all three docs to `paper/quark-scan-2026q2`
   and pushes.

### How to invoke Codex
- Wrapper: `/n/home09/obarrera/bin/codex_worker.sh`
- Form: `bash ~/bin/codex_worker.sh -C /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing "<prompt>"`
- **CRITICAL**: do NOT use the wrapper's `-o` flag. It overwrites the
  output file with Codex's summary message and clobbers any longer
  content. Instead, in the prompt instruct Codex to write outputs via
  shell (`cat > /tmp/foo.md <<'EOF' ... EOF`).
- Run in background via `run_in_background: true`. You'll be notified
  on completion. Don't poll.
- Codex's gpt-5.4 xhigh takes 5-30 min per task depending on scope.
  Plan accordingly.

### Communication discipline (impose on every Codex prompt)
- Every Codex task must print `<TASK_NAME>_DONE` or
  `<TASK_NAME>_FAILED: <reason>` to stdout at the end.
- Every Codex task must write its report to a specific file via shell.
- Every Codex task must specify branch hygiene (which branch to work
  on, how to merge back).

## 2. State at compaction time

### Phase status (all complete)
- Phase 1: snapshot bundle + branch + topic commits + v0.1 tag
- Phase 2 holes 4, 5, 6, 7, 8, 9, 2: all complete with sign-off chain
- Invalidation gate re-run: 9 SLURM jobs re-run, plots regenerated,
  methodology note updated
- Phase 3 hole 10 (figure prune): complete
- Phase 3 final acceptance: complete, rc1 tag pushed

### Sign-off chain commit SHAs (for reference)
- v0.1-snapshot: `6071dea`
- rc1: `1faba23` (peeled), `f4d4fa8` (annotated)

### Headline numbers (locked at rc1)
| Quantity | Value | Convention |
|---|---|---|
| $M_{KK}^{\min}$ p50 | 47.26 +69.4 -25.0 TeV | $g_s^\star=3$ + BGS 2020 + LO + factor-3 PDG |
| $M_{KK}^{\min}$ p95 | 127.13 TeV | same |
| CFW comparison | factor 2.2 stronger at matched conventions | -- |
| Zero-pass UL (moreUV/moreIR) | $p \le 2.3\times10^{-6}$ | Wilson 95% |
| Zero-pass UL (Run C) | $p \le 9.2\times10^{-7}$ | Wilson 95% |
| $\varepsilon_K$ binding fraction | 99.5% | -- |

### Decision log
`docs/paper_execution_decisions.md` (read this; it records every
adjudication made by Opus + orchestrator).

### Roadmap
`docs/paper_execution_roadmap.md` (the converged plan v2).

## 3. What you must do after compaction

The PI has acknowledged the rc1 state but flagged that two weak spots
need closure before rc1 is "truly settled":

### Task A — End-to-end Opus physics read of the methodology note

Spawn an Opus sub-agent (`Agent` tool, `subagent_type: general-purpose`,
`model: opus`) with this brief:

- Read `docs/quark_scan_methodology_note.tex` end-to-end as if you've
  never seen the repo before. ~19 pages, 9 figures.
- Read `docs/audits/*.md` and `docs/phase_logs/*signoff*.md` for the
  audit trail.
- Look for: physics errors not caught by intra-Codex review,
  unsupported claims, missing systematics, sloppy phrasing of
  load-bearing statements (the factor-2.2 CFW claim, the BGS choice,
  the BMU sign convention, etc.).
- Write `docs/phase_logs/phase3_final_opus_signoff.md` with PASS /
  PASS-WITH-FOLLOWUP / FAIL and findings.
- If FAIL or PASS-WITH-FOLLOWUP, list specific fixes. Then dispatch
  Codex to apply them, then re-run sign-off until clean.

### Task B — Re-execute the 4 paper-grade notebooks under post-audit constants

Spawn a Codex implementer (background, gpt-5.4 xhigh) with this brief:

- The 4 notebooks live under `notebooks/`:
  - `dense_scan_2sigma_vs_1sigma_comparison.ipynb`
  - `dense_scan_mkk_constraints_pdg2024.ipynb`
  - `pdg_quark_target_fix_verification.ipynb`
  - `rs_anarchy_analysis.ipynb`
- Re-execute each one under the post-audit `deltaf2.py` (current HEAD).
  Use `jupyter nbconvert --to notebook --execute --inplace` or
  equivalent. The conda env is `/n/home09/obarrera/.conda/envs/ising_bootstrap`.
- Verify the output cells now show the post-audit numbers (47 TeV at
  p50, etc.) rather than the pre-audit numbers (10 TeV).
- If any notebook errors during execution, fix the bug (likely a stale
  import or removed file path).
- Commit each notebook on the paper branch:
  `notebooks: re-execute <notebook name> under post-audit ΔF=2 constants`.
- After all 4: dispatch Codex peer review on the diff.
- After peer-review APPROVE: dispatch Opus sign-off.
- Seal logs at `docs/phase_logs/phase3_notebook_rerun_*.md`.

### After A and B complete

1. Check whether either task produced a substantive content change
   (numbers shifted, doc wording changed, etc.). If yes:
   - Update the methodology PDF (rebuild).
   - Update the headline summary if anything moved.
   - Refresh checksums in `artifacts/checksums.sha256`.
2. Tag a NEW release: `quarkscan-paper-rc1.1` (or `rc2` if substantive
   changes). Push to origin.
3. Report to the user: "rc1 settled. Task A and B both PASS. rc1.1
   tagged. The methodology note is now ready for the next stage
   (introduction/abstract drafting, external collaborator review,
   eventual submission)."

If A and B both PASS without changes, just confirm rc1 is sealed and
no rc1.1 is needed.

## 4. Conventions that survived compaction

- Headline `g_s^\star = 3` (CFW no-UV-boundary-term)
- Perturbative $g_s \approx 1.05$ reported alongside
- BGS 2020 $\varepsilon_K^{\rm SM} = 2.161 \times 10^{-3}$ central +
  asymmetric uncertainty band
- FLAG 2024 bag parameters (`B_K = 0.5503`, `B_4^K = 0.903`,
  `B_5^K = 0.691`)
- BMU LO Wilson running with `Q_1^{LR}_{\rm BMU} = -2\,O_5^{LR}` and
  scalar-LR ADM `[[-16, -6], [0, 2]]`
- Factor-3 PDG gate (mass + CKM); factor-5 on $J$
- Anarchic Yukawa prior $\mathcal{U}(-1.5, 1.5)$ on Re + Im,
  rejection floor $|Y| \ge 0.1$
- Fixed c-pattern: $c_Q = (0.63, 0.57, 0.20)$,
  $c_u = (0.66, 0.50, -0.50)$, $c_d = (0.66, 0.61, 0.55)$
- All M_KK^min via $\Mkk^{\rm tile}\sqrt{\max\text{ratio}}$ rescaling

## 5. Deferred items (NOT to address now; documented for rc1.1+)

1. Endpoint mismatch (B/D bag params at $m_b$/3 GeV vs Wilson endpoint
   at 2 GeV) — 10-30% systematic
2. NLO BMU Wilson running — 10-30% systematic
3. VLL/VRR NLO sanity check
4. Quantitative BMU/BBL convention table
5. Run C higher-statistics rerun (current $n=217$ at CFW gate)

## 6. Files / paths the orchestrator MUST know

- Plan: `docs/paper_execution_roadmap.md`
- Decisions: `docs/paper_execution_decisions.md`
- Artifact manifest: `docs/artifact_manifest.md`
- Audit decisions: `docs/audits/`
- Phase logs: `docs/phase_logs/`
- All canonical scan dirs: `scan_outputs/rs_anarchy_*_20260515T085*`
- All canonical figures: `results/figures/quark/*.{pdf,png}` (9 PDFs)
- Exploratory figures: `results/figures/quark/exploratory/`
- Pre-audit figure snapshot: `results/figures/quark_pre_audit_constants/`
- Codex wrapper: `/n/home09/obarrera/bin/codex_worker.sh`

## 7. Operating principles

- **Delegate everything**. The orchestrator does NOT write code, does
  NOT run pytest, does NOT WebFetch. The orchestrator reads agent
  outputs, makes go/no-go decisions, and commits logs.
- **Cross-review is non-negotiable**. Two BLOCKERs were caught at
  peer-review during Phase 2 (BMU sign, CFW g_s* convention). The
  pattern works. Use it.
- **Honest physics framing > impressive numbers**. The "factor-2.2
  stronger than CFW" is more credible than a "11% agreement" claim
  would have been.
- **Long-running jobs in background**. Codex and Opus calls always
  go via `run_in_background: true`; you'll get a notification on
  completion.

## 8. Final sanity check before you start

Run these in a single Bash call:

```bash
git status -sb && git log --oneline -5 && git ls-remote --tags origin | tail -5 && ls docs/phase_logs/phase3_*signoff* 2>&1
```

You should see:
- Branch `paper/quark-scan-2026q2`, in sync with origin (nothing ahead/behind)
- Last commit is the Phase 3 final summary
- Tags include `quarkscan-paper-rc1` AND `quarkscan-paper-v0.1-snapshot`
- A Phase 3 sign-off doc (`phase3_h10_signoff.md`) is present

If any of those don't match, STOP and investigate before delegating.
The pipeline state is the ground truth, not this briefing.

---

# Your first action after compaction

Read this briefing. Then dispatch Task A (Opus end-to-end physics read)
and Task B (Codex notebook re-execution) **in parallel** (single message,
two tool calls).

Report back to the user only after both have completed AND their
sign-off chain (Opus for A; Codex peer review + Opus for B) has
APPROVE / PASS verdicts.
