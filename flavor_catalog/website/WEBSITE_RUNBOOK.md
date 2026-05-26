# WEBSITE_RUNBOOK, Flavor Catalog Static Site Build

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

Catalog facts: 102 entries (94 PRIMARY + 8 SECONDARY) across 8 families. 101 VERIFIED / 1 PARTIAL / 0 MISMATCH (post-v0.4-tag canonical; K020 cleared by `b5c2375` before the v0.4 tag — matches `master_compile_v04_report.md` §"Consolidation status" and the website's `catalog_index.json` `verdict_counts`). Per-entry YAML schema documented in Explore agent report. M_KK^min,p50 = 47.26 TeV at g*=3 (rc1.1 quark scan), load-bearing for honest collider-tier framing.

---

## Dispatch ledger

| # | Phase | Type | Model | Deliverable | Bg ID | Status | Output paths |
|---|---|---|---|---|---|---|---|
| 0 | Pre-flight | direct |, | branch cut + runbook |, | DONE | `flavor_catalog/website/WEBSITE_RUNBOOK.md` |
| 0a | Required reading | Explore | sonnet | structured catalog summary |, | DONE | (consumed inline) |
| 1 | Scaffold | subagent | opus | Astro 5.x + KaTeX scaffold, ingest pipeline, K001 template, working localhost |, | DONE (commit 460ece9) | `flavor_catalog/website/{package.json,astro.config.mjs,scripts/ingest_catalog.py,src/**}`; dev `http://localhost:4321/` |
| 2a | Codex anchor batch | codex via codex-delegate | gpt-5.5 xhigh | kaon (11) + secondary/kaon (3) = 14 entries citation anchors (K007/K011 didn't exist; K017/K018 backfilled) | wrapper a32bcb9d4d8d4ba64 | DONE (commit pending) | `_data/citation_anchors/K*.yaml` |
| 2b | Codex anchor batch | codex via codex-delegate | gpt-5.5 xhigh | beauty (22) + secondary/beauty (4) = 26 entries | wrapper a5e3a526fbe422358 | DONE | `_data/citation_anchors/B*.yaml` |
| 2c | Codex anchor batch | codex via codex-delegate | gpt-5.5 xhigh | top_higgs_ew (19) + secondary/top_higgs_ew (1) = 20 entries | wrapper adfffbb9315974f8e | DONE | `_data/citation_anchors/T*.yaml`,`EW*.yaml` |
| 2d | Codex anchor batch | codex via codex-delegate | gpt-5.5 xhigh | charged_lepton (11) + charm (8) = 19 entries | wrapper a70c8ceabba7c4244 | DONE | `_data/citation_anchors/L*.yaml`, `C*.yaml` |
| 2e | Codex anchor batch | codex via codex-delegate | gpt-5.5 xhigh | edm_neutrino (7) + collider_rs (14) = 21 entries | wrapper a940497068cec94e1 | DONE | `_data/citation_anchors/E*.yaml`, `CR*.yaml` |
| 2f | Codex backfill | codex via codex-delegate | gpt-5.5 xhigh | K017 + K018 (missed in 2a) | aec3284bcf3f2b973 | DONE | `_data/citation_anchors/K017.yaml`, `K018.yaml` |
| 3 | Data-driven rendering | subagent | opus | generalize K001 template → all 102; citation modal (RES/AMB/UNR); family pages; `/browse/` with Pagefind + 5 filters; dark-mode + 375px polish | a58db0da4350eb940 | DONE (commit pending) | new: `src/components/CitationModal.astro`,`EntryTable.astro`,`src/pages/browse.astro`,`scripts/screenshot_cdp.mjs`; modified: ingest, [id].astro, [family].astro, index.astro, global.css. Build: 112 pages, Pagefind index 102 entries/6907 words. |
| 4 | Citation audit | direct |, | sample 10 entries × 2 anchors each, verify line N of snapshot contains anchor value |, | DONE | 20/20 sample anchors PASS (4 "FAILs" in script were normalization gap, manually eyeballed all 4, snapshots contain cited values). Zero false positives in sample. |
| 5 | Methodology page + CF Pages config + README | subagent | opus | `/methodology/` rendering CATALOG_METHODOLOGY.tex + AGENTIC_WORKFLOW + wave timeline + signoff rounds + fact-check table + Phase-2 anchor stats; CF Pages config (build cmd, output dir, Node, `_redirects`, `_headers`); website README | a5ee872c8783e7233 | DONE (commit 5f087fd) | new: `src/pages/methodology.astro`, `cloudflare-pages.config.md`, `README.md`, `.node-version`, `public/_redirects`, `public/_headers`; modified: BaseLayout (nav), index.astro (hero CTA). Build: 113 pages. |
| 6a | LaTeX delimiter fix | codex via codex-delegate | gpt-5.5 xhigh | enable `$...$` (and keep `$$...$$`, `\(...\)`, `\[...\]`) as KaTeX auto-render delimiters site-wide; expose `window.fcatRenderMath()` helper for dynamic content | aabc0699b243dc0b9 | DONE (commit pending) | modified: `src/layouts/BaseLayout.astro`, `src/pages/entries/[id].astro`. |
| 6b | UI polish | subagent | sonnet | row-clickable tables (no Open column); drop "DIFF:" prefix on difficulty badge; sleek landing hero (move 94/8 stats out of hero, scope-note below family cards) | a610a0e5d3afc1de6 | DONE (commit pending) | modified: `Badges.astro`, `EntryTable.astro`, `browse.astro`, `index.astro`, `global.css`. |
| 6c | Methodology rewrite | subagent | opus | cut 660 → 224 lines; remove "Wave N"/"PI seed"/"PKA" jargon, signoff-round row enumeration, arbitration trivia; add 4-card pipeline (Identify→Draft→Review→Verify); keep fact-check table + honest framing + PARTIAL note | a56ddb9582e90fcfa | DONE (commit pending) | modified: `src/pages/methodology.astro` (-436 lines). |
| 10 | UI cleanup: tier/anchors | direct | - | Remove Primary/Secondary tier badges everywhere; remove ANCHORS pill from status row; fix PDG table observable column (was wrapped in `$...$` collapsing prose spaces); fix TeV unit duplication when `display` already encodes unit | - | DONE (commit 7053fb7) | modified: `[id].astro`, `EntryTable.astro`, `browse.astro`, `global.css`. |
| 11 | LaTeX rendering sweep | subagent | opus | Screenshot all page types; fix Greek letter boundary matching (nu_e, _Lambda, etc.); add capital Greek; add Kbar/sbar/dbar bar variants; multi-char subscript brace-wrap (_FB -> _{FB}); prose-word text-wrappers; tighten looksLikePhysics() | a7221cd51f8480fe7 | DONE (commit 4eba297) | modified: `src/lib/notation.ts`, `src/lib/prose.ts`. |
| 12 | Methodology rewrite v2 | subagent | opus | Trust-first hero; 3-card pipeline (Identify/Verify/Check) in plain English; T020+K020 error-catch cards with KaTeX; per-family fact-check table; honest scope-limit; no jargon | a448288d3102e40ee | DONE (commit 4eba297) | modified: `src/pages/methodology.astro`. |

## Operational notes

- **Dev server**: `npm run dev -- --host 127.0.0.1 --port 4321`. Stays alive between phase commits.
- **Pagefind**: `public/pagefind/` is gitignored (build artifact). Production path = `astro build && pagefind --site dist`. For local dev on a fresh clone, run `npm run build` once to populate `public/pagefind/`; otherwise search falls back to substring match.
- **Screenshot helper** at `scripts/screenshot_cdp.mjs` (Node 22+ built-in WebSocket) for true-375px mobile captures bypassing Chrome headless's 500px-default viewport.

**Phase 2 totals**: 102/102 entries · 644 anchors · 476 RESOLVED (73.9%) · 121 AMBIGUOUS (18.8%) · 47 UNRESOLVED (7.3%).
- Per-family UNRESOLVED %: B 1.0, C 23.7 (**>20% soft trigger**), CR 15.7, E 0.0, EW 9.1, K 2.1, L 12.1, T 11.5.
- 3 entries with >50% UNRESOLVED: C002 (3/5), CR004 (2/3), CR010 (4/7). UI marks UNRESOLVED clearly.
- Charm exceeds the brief's 20% soft-STOP trigger by 3.7 pp. Per PI directive "do not stop", continuing; surfaced for end-of-build decision on targeted re-resolve.

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
