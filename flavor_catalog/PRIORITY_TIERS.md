# Flavor Catalog Priority Tiers

## 1. Policy

The existing 80 v0.2 catalog entries from Waves 1-7 are the **PRIMARY** tier.
They are the highest-leverage RS-flavor constraints, OPUS-APPROVED, and
fact-checked. They are the canonical set for any downstream code
implementation.

Any new process added in Wave-8 or later is **SECONDARY** by default. These
processes were deferred by DA-1 through DA-4 plus external review on leverage
grounds: long-distance dominated, redundant with stronger active entries,
future-projection-only, or relevant only under specific scope extensions such
as lepton-bulk RS.

The PRIMARY default for v0.2 entries is implicit. Absence of a `priority_tier`
field in a sidecar means **PRIMARY**. A retroactive PRIMARY tag pass is clean
and optional, but it is not required.

## 2. Repo-Structure Separation

This split is binding architecture.

PRIMARY entries live at `flavor_catalog/processes/<family>/<id>.tex` and
`flavor_catalog/processes/<family>/<id>.yaml`.

SECONDARY entries live at
`flavor_catalog/processes/secondary/<family>/<id>.tex` and
`flavor_catalog/processes/secondary/<family>/<id>.yaml`.

References (`references/<id>/`), worklogs
(`worklogs/{pka,writer,checker}/<id>.md`), and per-process signoff
(`signoff/by_process/<id>.md`) stay flat for both tiers. Process IDs are
globally unique, so nesting these support artifacts by tier adds no
information.

`catalog_master.tex` includes a dedicated top-level
`\section{Secondary Entries (Wave-8+, lower implementation priority)}` with
`\subsection*{<Family>}` and
`\input{processes/secondary/<family>/index}` per affected family.

Rationale: the PI requested a structural distinction visible from a bare `ls`
and easy to separate later. `git mv` of the `secondary/` subtree, whole or
per-entry, is the trivial promotion mechanism if a SECONDARY entry is later
judged worth promoting to PRIMARY. Tree-deletion is the trivial purge
mechanism.

## 3. YAML Sidecar Contract for SECONDARY

Every SECONDARY sidecar carries these fields near the top, after
`process_id` and `family`:

```yaml
priority_tier: SECONDARY
priority_rationale: <one line: why originally deferred + why being added now>
promoted_in_wave: <integer wave number>
```

Do not retroactively edit any PRIMARY sidecar in this round.

## 4. Per-Wave Rationale

Wave-8 uses the top-tier candidate list from `flavor_catalog/SESSION_NOTES.md`
Section 7 and the deferred-scope rationale in
`flavor_catalog/worklogs/discovery/round_004_addendum_deferred_scope.md`.

| ID | Process | Family | Why originally deferred | Why being promoted in Wave-8 |
|---|---|---|---|---|
| K019 | `K_L → e^± μ^∓` | kaon | DA-4: LFV-kaon trio deferred as lepton-bulk RS scope extension | Top-tier candidate; explicit lepton-bulk-RS handle |
| K020 | `K^+ → π^+ e^± μ^∓` | kaon | DA-4: same trio | Same; companion to K019 |
| K021 | `K_L → π^0 e^± μ^∓` | kaon | DA-4: same trio | Same; neutral-kaon LFV companion |
| B007 | `B_{s,d} → e^+ e^-` | beauty | DA-4: rare leptonic e tail | Top-tier; electron-specific scalar operator handle |
| B008 | `B_{s,d} → τ^+ τ^-` | beauty | DA-4: rare leptonic τ tail | Top-tier; tau-specific scalar operator handle |
| B013 | `B_s → φ γ` | beauty | DA-4: photon-helicity / time-dependent CP variant | Top-tier; helicity companion to active B011/B012 |
| B014 | `B → ρ γ`, `B → ω γ` | beauty | DA-4: b→dγ exclusive tail | Top-tier; down-dipole chirality probe |
| T014 | `Z → bs, bd, sd` | top_higgs_ew | DA-4: weak FCNC-Z tails | Top-tier; only dedicated FCNC-Z analysis row |

## 5. Implementation-Priority Rule

When the "first constraint into code" task arrives, currently scoped to
`quarkConstraints/modern/`, the implementation queue **MUST** process PRIMARY
entries first.

No SECONDARY entry may be pulled into code without explicit PI confirmation,
even if its physics looks superficially simpler. The SECONDARY tag means the
catalog authors and DA judged it lower-leverage for the discovery-mode
pipeline. Bypassing that judgment requires a PI override on the record.

## 6. Cross-References

This document is policy. The live wave state is in
`flavor_catalog/worklogs/orchestration/wave_008_runbook.md`.

Discovery and deferral context is in
`flavor_catalog/worklogs/discovery/round_004_addendum_deferred_scope.md`.

The PI handoff that introduced the tier system is
`flavor_catalog/HANDOFF_PROMPT.md`, under the "Priority tagging" section.

Orchestrator decisions on the PI's behalf are recorded in
`docs/phase_logs/flavor_catalog_orchestrator_decisions.md`.

## 7. How To Extend

A future Wave-N (`N ≥ 9`) follows the same pattern. Add a new row to the
per-wave rationale table with the new ID, why it was deferred, and why it is
being promoted. Increment `promoted_in_wave` to `N`. Use the same
`flavor_catalog/processes/secondary/<family>/` repo path. Update
`catalog_master.tex` only if a new family appears.
