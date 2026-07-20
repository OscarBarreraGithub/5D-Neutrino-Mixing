# Flavor Catalog Priority Tiers

## 1. Policy

The existing 80 v0.2 catalog entries from Waves 1-7 are the **PRIMARY** tier.
They are the highest-leverage RS-flavor constraints, OPUS-APPROVED, and
fact-checked. They are the canonical set for any downstream code
implementation.

A new process added in a later wave is **SECONDARY** if it came from the
DA-deferred pool (the 51 plan-v1 Section C rows DA-4 explicitly deferred).
These were deferred by DA-1 through DA-4 plus external review on leverage
grounds: long-distance dominated, redundant with stronger active entries,
future-projection-only, or relevant only under specific scope extensions such
as lepton-bulk RS.

A new process added in a later wave is **PRIMARY** (no `priority_tier`
field) if it opens a **new scope class** at PI directive — i.e., the entry
was not in the plan-v1 Section C deferred list, but represents a category
the catalog didn't previously cover. The canonical example is Wave-9's
`collider_rs` family (direct-collider RS resonance / heavy-partner searches),
which was deliberately OUT-OF-SCOPE for the low-energy catalog and was opened
as a sibling tier by explicit PI directive.

The PRIMARY default is implicit. Absence of a `priority_tier` field in a
sidecar means **PRIMARY**, regardless of which wave introduced the entry.
A retroactive PRIMARY tag pass is clean and optional, but it is not
required.

## 2. Repo-Structure Separation

This split is binding architecture.

**PRIMARY** entries live at `flavor_catalog/processes/<family>/<id>.tex`
and `flavor_catalog/processes/<family>/<id>.yaml`. Each family is a
top-level subdirectory:
- Original 6 low-energy families (Waves 1-7): `kaon/`, `charm/`,
  `beauty/`, `top_higgs_ew/`, `charged_lepton/`, `edm_neutrino/`.
- Wave-9 new-scope family: `collider_rs/` (KK gauge resonances,
  custodial top/bottom partners, VLQ, KK-graviton, EW-precision tail).
- Future new-scope families: each gets its own top-level directory.

**SECONDARY** entries live at
`flavor_catalog/processes/secondary/<family>/<id>.tex` and
`flavor_catalog/processes/secondary/<family>/<id>.yaml`, mirroring the
PRIMARY family layout under a single `secondary/` parent.

References (`references/<id>/`), worklogs
(`worklogs/{pka,writer,checker}/<id>.md`), and per-process signoff
(`signoff/by_process/<id>.md`) stay flat for both tiers. Process IDs are
globally unique, so nesting these support artifacts by tier adds no
information.

`catalog_master.tex` has:
- One top-level `\section{<Family Name>}` per PRIMARY family (e.g.,
  `\section{Kaon}`, `\section{Beauty}`, …, `\section{Collider RS
  Resonances}`).
- A single top-level `\section{Secondary Entries (Wave-8+, lower
  implementation priority)}` containing `\subsection*{<Family>}` +
  `\input{processes/secondary/<family>/index}` per affected SECONDARY
  family.

Rationale: the PI requested a structural distinction visible from a bare
`ls` and easy to separate later. `git mv` of the `secondary/` subtree,
whole or per-entry, is the trivial promotion mechanism if a SECONDARY
entry is later judged worth promoting to PRIMARY. Tree-deletion is the
trivial purge mechanism. For PRIMARY new-family scope-class extensions
(Pattern B in §7), each new family is structurally peer to the existing
six low-energy families — equally easy to add, equally easy to remove.

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
`docs/archive/phase_logs/flavor_catalog_orchestrator_decisions.md`.

## 7. How To Extend

Future waves fall into one of two patterns. Both use the same
PKA → WA → CA → fact-check → Opus pipeline, the same cycle caps, the
same arbitration precedents. The differences are tier label, sidecar
fields, and repo structure.

### Pattern A — SECONDARY re-promotion (Wave-8 was the canonical example)

DA-deferred entries from the plan-v1 Section C list (`round_004_addendum_deferred_scope.md`)
get promoted with explicit rationale.

- **Path**: `processes/secondary/<existing-family>/<id>.{tex,yaml}`
- **YAML required fields**:
  ```yaml
  priority_tier: SECONDARY
  priority_rationale: <one line: why originally deferred + why being added now>
  promoted_in_wave: <integer wave number>
  ```
- **Master TeX**: add a `\subsection*{<Family>}` line under the existing
  `\section{Secondary Entries (Wave-8+, …)}` (or, if the family is
  already represented, the entry's `\input` is auto-picked up by the
  existing per-family secondary `index.tex`).
- **Runbook**: `worklogs/orchestration/wave_<NNN>_runbook.md` (sibling
  of the Wave-8 runbook).
- **Per-wave rationale table** in §4 gets a new row per added entry.

### Pattern B — PRIMARY new family (Wave-9 was the canonical example)

PI opens a new scope class (e.g., direct-collider RS in Wave-9; could be
a fermion-mass-extension family, a leptogenesis observables family, a
neutrino oscillation-fit family, etc.).

- **Path**: `processes/<new-family>/<id>.{tex,yaml}`
- **YAML**: do NOT add `priority_tier` field (implicit PRIMARY).
- **Master TeX**: add a new top-level `\section{<New Family Name>}`
  placed appropriately (Wave-9 placed `\section{Collider RS Resonances}`
  between EDM/Neutrino and the SECONDARY section).
- **Runbook**: `worklogs/orchestration/wave_<NNN>_runbook.md` with
  decisions D-1..D-N documenting tier choice, family name,
  catalog_master.tex placement, batch grouping.
- **New audit file**: `audits/factcheck_<new-family>.md` (mirrors the
  Wave-1 audit format used by existing families).
- **Per-wave rationale table** in §4 gets a header noting the family is
  new and listing the entries with brief sub-class labels.

### Which to use

If the entry is in the plan-v1 deferred list → SECONDARY (Pattern A).
If the entry represents a new scope class not previously cataloged →
PRIMARY new family (Pattern B). When in doubt, ask the PI.
