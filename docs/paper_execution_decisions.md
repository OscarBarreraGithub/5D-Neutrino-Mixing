# Paper Execution: Decisions Log

This document records orchestrator-level decisions made on the user's behalf
during execution of the roadmap at `/tmp/codex_plan_v2.md`. All agents
(Codex implementers, Codex reviewers, Claude Opus sign-off agents) MUST read
this file before acting.

Last updated: 2026-05-15.

## Scope

**Quark-sector paper only.** The lepton sector (`neutrinos/`,
`flavorConstraints/`, `yukawa/`) is OUT OF SCOPE for this submission. The
methodology doc, the RS-anarchy framework, and all 10 prioritized holes
from the plan are quark-sector work driven by the recent collaborator
critique. Lepton-sector work is deferred to a follow-up paper. Where the
methodology doc currently is mute on this, add a one-sentence scope note
("This paper addresses the quark sector; the lepton sector of this repo is
the subject of follow-up work.").

## External scan-output snapshot

**Location:** `/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/artifacts/quarkscan_paper_<release_tag>/`,
a tracked manifest (in `docs/artifact_manifest.md`) pointing to it. No
public DOI yet; defer Zenodo / Dataverse registration to the submission
moment. The intermediate snapshots live alongside `scan_outputs/` on the
cluster; they will be tar+gzipped, sha256summed, and recorded in the
manifest by file path + checksum + size + run command + code SHA.

## Branch and tag scheme

- Main paper branch: `paper/quark-scan-2026q2`, off verified `origin/main`.
- Per-hole branches off the paper branch: `fix/pdg-benchmarks`,
  `audit/bag-inputs`, `audit/wilson-rg`, `audit/cfw-comparison`,
  `scan/zero-pass-statistics`, `docs/paper-readiness`. Merge by
  fast-forward after a Codex peer review + Claude Opus sign-off.
- Tags: `quarkscan-paper-v0.1-snapshot` (Phase 1 capture),
  `quarkscan-paper-rc1` (audits passed), `quarkscan-paper-v1.0`
  (submission).
- Push every branch to `origin` once it lands a meaningful commit (the
  remote is the user's GitHub fork; pushing is safe).

## Ownership

- **Codex (gpt-5.4 xhigh)**: all implementation, audits, literature
  search, code changes, doc edits, scan re-runs.
- **Codex peer review** (a separate Codex invocation): reviews each
  Codex-implementer deliverable before sign-off.
- **Claude Opus subagent**: signs off after each phase. Validates that the
  phase exit criteria are met against the plan's checklist.
- **Claude main (orchestrator)**: makes phase transitions, records
  decisions, monitors background jobs. No implementation work.
- **Human (PI)**: green-light only at phase boundaries and at the final
  release tag.

## Audit philosophy

For the bag-parameter, Wilson-RG, and CFW audits: **prefer narrowing
claims over silently changing numbers.** If the audit finds the existing
code uses a non-standard convention, the first option is to document the
convention and state explicitly what was assumed; only if the convention
is clearly wrong (i.e. contradicts the literature consensus) do we change
the code and re-run scans.

## Re-run policy

If an audit changes a numerical input by more than 10% in any
$\Delta F=2$ Wilson coefficient at $M_{KK} = 3$ TeV, the affected
RS-anarchy scan(s) (RUNA at minimum) must be re-run before the final
manifest is signed. Below 10%, the existing scan outputs may be
retained with a footnote in the methodology doc.

## Stopping conditions for the orchestrator

The execution loop stops only when:
1. All 10 holes from the plan are marked completed in `TodoWrite`.
2. Phase 3 checklist (10 items) is fully checked.
3. The release tag `quarkscan-paper-rc1` exists on origin and the PDF
   builds cleanly from a fresh clone.

The orchestrator does NOT stop for:
- An audit returning "no change needed" (continue to next hole).
- A Codex implementer needing literature lookup (provide WebFetch
  permission and continue).
- A peer-review iteration (loop back to implementer with findings).

The orchestrator DOES stop and ask the human for:
- A previously undocumented physics decision (e.g., a new convention).
- An audit that finds a >50% bound shift (rare; would require user
  acknowledgement of new headline numbers).
- A blocker that 3 consecutive Codex iterations cannot resolve.

## Bookkeeping requirements for every agent task

Every Codex task must produce:
1. A diff (git or patch format) of code/doc changes.
2. A markdown task report at `docs/phase_logs/<task_id>.md` summarising
   what was done, what was verified, and what remains open.
3. A unit test or scripted verification that the change works (when
   applicable).

Every Codex peer-review task must produce:
1. APPROVE / REJECT-WITH-REVISIONS / REJECT verdict.
2. A numbered findings list at `docs/phase_logs/<task_id>_review.md`.

Every Claude Opus sign-off must produce:
1. PASS / FAIL with reasoning at `docs/phase_logs/<phase>_signoff.md`.
2. A list of any unresolved deviations from the plan.
