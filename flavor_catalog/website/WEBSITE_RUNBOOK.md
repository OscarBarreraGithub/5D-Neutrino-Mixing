# WEBSITE_RUNBOOK — Flavor Catalog Static Site Build

**Mission**: collaborator-facing static site for 102-entry catalog at tag `flavor-catalog-v0.4`. Cloudflare Pages target. Cold-boot brief: [`../WEBSITE_BUILD_PROMPT.md`](../WEBSITE_BUILD_PROMPT.md).

**Branch**: `flavor-catalog-website/2026q2` off `flavor-catalog/2026q2` @ `6d99b17`. Never push to `flavor-catalog/2026q2`, `main`, or `paper/quark-scan-2026q2`.

**Orchestrator stance**: dispatch only. No inline code, no inline CSS, no per-entry content authored by orchestrator. Sample-audit only. Updates this ledger after every dispatch.

**Stack (provisional, scaffold subagent confirms or overrides)**: Astro 5.x + KaTeX + Pagefind, static output, Cloudflare-Pages-friendly. Local toolchain: node 25.6.1, python 3.14.3, npx, codex 0.130.0. No hugo, no pnpm, no codex wrappers.

---

## Required reading completed

- `flavor_catalog/README.md`, `CATALOG_METHODOLOGY.tex`, `PRIORITY_TIERS.md`
- `processes/kaon/K001.{tex,yaml}`, `processes/collider_rs/CR002.{tex,yaml}`, `processes/secondary/kaon/K019.{tex,yaml}`
- `references/K001/source_manifest.yaml` + representative snapshot
- `audits/factcheck_kaon.md`
- Skim of `AGENTIC_WORKFLOW.md`

Catalog facts: 102 entries (94 PRIMARY + 8 SECONDARY) across 8 families. 100 VERIFIED / 2 PARTIAL / 0 MISMATCH. Per-entry YAML schema documented in Explore agent report. M_KK^min,p50 = 47.26 TeV at g*=3 (rc1.1 quark scan) — load-bearing for honest collider-tier framing.

---

## Dispatch ledger

| # | Phase | Type | Model | Deliverable | Bg ID | Status | Output paths |
|---|---|---|---|---|---|---|---|
| 0 | Pre-flight | direct | — | branch cut + runbook | — | DONE | `flavor_catalog/website/WEBSITE_RUNBOOK.md` |
| 0a | Required reading | Explore | sonnet | structured catalog summary | — | DONE | (consumed inline) |
| 1 | Scaffold | subagent | opus | Astro 5.x + KaTeX scaffold, ingest pipeline, K001 template, working localhost | — | DONE | `flavor_catalog/website/{package.json,astro.config.mjs,scripts/ingest_catalog.py,src/**}` |

---

## STOP-and-ping triggers (escalate to PI)

- Codex citation-anchor failure > 20% on any family
- Catalog value disagrees with linked source after careful check
- Load-bearing framework choice (e.g., Pages vs Workers, static vs SSR)
- Out-of-scope creep (live API, comments, editing UI, etc.)

---

## Phases

1. **Scaffold** (opus subagent): stack + landing + 1 example entry template → working localhost.
2. **Codex extraction** (~6 parallel codex agents): 102 normalized entries + citation anchors in `_data/`.
3. **Data-driven rendering** (opus subagent): wire data, family pages, search, badges, dark mode, mobile.
4. **Citation audit** (direct, orchestrator): sample 5-10 entries, confirm anchors land correctly.
5. **Methodology page + Cloudflare config** (opus subagent).
6. **PI visual review at localhost**.
