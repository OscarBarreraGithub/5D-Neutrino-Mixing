# Cold-Boot Prompt: Build the Flavor Catalog Website

Paste the fenced block below into a fresh local Claude Code session. It's
self-contained — points the new Claude at the catalog source, the
verification trail, and the split between Claude (UI) and codex (per-entry
content extraction).

---

```
You are taking over a NEW task: orchestrate the build of a clean,
professional static website that surfaces the 5D-Neutrino-Mixing flavor
catalog (currently 102 OPUS-APPROVED entries at tag
`flavor-catalog-v0.4`). The catalog itself is built and verified — your
job is purely to make it collaborator-readable on the web.

**You are a pure orchestrator.** This is multi-day, multi-phase work
involving 102 per-entry content extractions, site scaffolding, layout,
CSS, math rendering, citation-anchor resolution, visual review, and
deployment config. If you try to do any of that yourself inline, your
context will fill up and you'll lose continuity. Instead, you dispatch
subagents and review short verdict lines.

The split is:
- **Codex GPT-5.5 xhigh** (via `~/bin/codex_worker.sh` or `codex exec`)
  for the high-volume repetitive per-entry work: content extraction,
  citation-anchor resolution, batch JSON/YAML generation.
- **Claude subagents** (via the Agent tool with `subagent_type:
  general-purpose` and `model: opus` for design-quality work, or
  `model: sonnet` for routine implementation) for UI work: framework
  scaffolding, per-page templates, CSS / design system, layout
  components, search UX, methodology-page authoring, Cloudflare-Pages
  config, visual review reports.
- **You (orchestrator)** dispatch, read short verdicts, sample-audit
  outputs, and decide next dispatch. You write no code, no CSS, and
  no per-entry content yourself. The only things you author directly
  are orchestration meta-docs: a `WEBSITE_RUNBOOK.md` (mirrors the
  catalog's wave runbooks) and short commit messages.

This will be hosted on Cloudflare Pages via the standard "Cloudflare
Pages connects to a GitHub repo and rebuilds on push" workflow. For
now you only need to confirm it runs locally (localhost dev server);
deployment config is the final-step deliverable, not the immediate
goal.

# Working environment

- Repo: `/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing` (or
  the local clone path)
- Cut a NEW branch off `flavor-catalog/2026q2` named
  `flavor-catalog-website/2026q2`. ALL website work goes there.
  Do not contaminate `flavor-catalog/2026q2` (the catalog content
  branch, frozen at v0.4) or `main` or `paper/quark-scan-2026q2`
  (frozen at rc1.1).
- Confirm starting commit is `2bda5f1` or later on
  `flavor-catalog/2026q2`. The catalog tag `flavor-catalog-v0.4` pins
  the canonical content state.

# Required reading (do in order, ~30 minutes)

1. `flavor_catalog/README.md` — layout overview and what's where.
2. `flavor_catalog/CATALOG_METHODOLOGY.tex` — the 1-page collaborator
   pitch. This is the substance the website needs to convey.
3. `flavor_catalog/PRIORITY_TIERS.md` — PRIMARY vs SECONDARY tiering.
   The website must surface this distinction visibly.
4. Three representative entries to understand schema and content:
   - `flavor_catalog/processes/kaon/K001.tex` and `K001.yaml`
     (epsilon_K — PRIMARY, in-code, well-cited)
   - `flavor_catalog/processes/collider_rs/CR002.tex` and `CR002.yaml`
     (T_{5/3} custodial — PRIMARY new-family from Wave-9)
   - `flavor_catalog/processes/secondary/kaon/K019.tex` and `K019.yaml`
     (LFV kaon, SECONDARY)
5. One reference manifest to see source-snapshot structure:
   `flavor_catalog/references/K001/source_manifest.yaml` and a
   snapshot file in that directory.
6. One fact-check report:
   `flavor_catalog/audits/factcheck_kaon.md` — to understand what
   per-claim verification looks like.

# Goal: what the website must do

The catalog is structured data + LaTeX prose. Your job is to render
it on the web with these properties:

## 1. Quality bar (you decide HOW; this is just the bar)

- Audience: physicists who will read this on a real screen and form
  an opinion about whether the underlying catalog is credible. The
  website is part of that credibility argument.
- Pick whatever stack, framework, design system, typography, color
  palette, density, layout, and theming approach you think serves
  the audience best. You're better positioned than this prompt to
  judge what looks credible to a 2026 physics collaborator. Make
  the call.
- Hard quality requirements (HOW is up to you):
  - Math rendering must work cleanly for LaTeX physics notation
    (entries are full of $\varepsilon_K$, $\mathcal{B}(\cdot)$,
    $M_{KK}$, integrals, etc.).
  - Both desktop and mobile must be usable. Collaborators will
    open this on phones.
  - Reading long-form prose must be comfortable (the per-entry
    Relevance + Post-2008 / Validity sections are several
    paragraphs).
  - The site must build to static files for Cloudflare Pages.
- Soft quality signals (do at least one well; doing all is great
  but not required):
  - Dark / light mode (if you ship it, do it well in both — partial
    dark mode is worse than no dark mode).
  - Accessibility / screen-reader basics (semantic HTML, alt text,
    contrast).
  - Loading performance (static = fast by default).
  - Search and / or faceted filtering across the 102 entries.

## 2. Catalog index — content requirements (design is up to you)

What the collaborator must be able to do, somehow:

- Land on the home page and, within ~10 seconds, understand:
  - how many entries there are (102) and what tag pins them (v0.4),
  - the fact-check status (100 VERIFIED / 2 PARTIAL / 0 MISMATCH),
  - which families exist and how many entries each has,
  - that there are two tiers (PRIMARY + SECONDARY) and roughly what
    each means,
  - how to drill down to a family or a specific entry.
- Navigate by family. Each family view shows its entries with the
  per-entry signals (status, tier, difficulty) at a glance.
- Find entries by search and/or filter (search box minimum; faceted
  filters are a nice-to-have).
- See SECONDARY entries — they must be visible but visually
  distinguishable from PRIMARY. They are NOT hidden; the point is
  the collaborator can see the full set and instantly know which
  tier each entry is in.

Honest framing the design must respect:
- PRIMARY > SECONDARY (the tier matters; collaborators may want to
  filter SECONDARY out before reading).
- Within PRIMARY, the low-energy flavor families (Waves 1–7) imply
  $M_{KK}^{\min,p50}\sim 47$ TeV at $g_*=3$. The `collider_rs`
  Wave-9 family is also PRIMARY but its LHC direct-search reach is
  $\sim$ few TeV. The site should reflect this honestly — don't
  oversell the collider tier as equally constraining.
- Implementation difficulty (LOW / MEDIUM / HIGH / BLOCKED) and
  existing in-code coverage (YES / PARTIAL / NO) are reasonable
  proxies for "where to start when implementing" — surface these
  somewhere so a collaborator scanning for a starting point can
  find them.

Family ordering on the landing page is your call; physics tradition
is kaon / charm / beauty / top–Higgs–EW / charged-lepton / EDM /
neutrino / collider, but feel free to deviate if you find a better
organization (e.g., grouping by sensitivity class).

## 3. Per-entry detail pages

Each of the 102 entries needs a page rendering the full content
cleanly:

- Header: process ID, standard notation (LaTeX rendered), family,
  tier badge.
- Status row: OPUS-APPROVED, fact-check verdict (VERIFIED / PARTIAL),
  cycle count, implementation difficulty (LOW/MEDIUM/HIGH/BLOCKED),
  code coverage status (YES/PARTIAL/NO with file:line if available).
- **The PDG/equivalent values block** as a structured table: each row
  is value + uncertainty/CL + year + units + source. Each value's
  source is a CLICKABLE LINK (see §4 below for the exact-location
  requirement).
- **Why this is relevant to RS / anarchic flavor** — render the .tex
  "Relevance" section.
- **Post-2008 / Post-2010 developments** — render the corresponding
  .tex section.
- **Constraint validity / model dependence** — render that section.
- **Code coverage** — render with hyperlinks to repo file:line
  evidence (links can be `https://github.com/.../blob/...#L<line>`
  format).
- **Implementation difficulty** — render with rationale.
- **Key references** — clickable list (see §4).
- **Provenance footer**: cycle count, who approved, fact-check
  agent, last update, sha256 of the local snapshot, link to the
  per-process worklog chain.

## 4. CITATION-EXACT REFERENCE LINKING (the critical requirement)

For every cited reference where a numerical claim is drawn from a
specific source, the website must let a collaborator click through
and see the EXACT spot the value comes from. Not just "this came from
arXiv:1234.5678" — but "this `< 1.46 TeV` mass-exclusion limit comes
from this specific line of the abstract / Table 1 / Figure 3".

Implementation:

- Each value entry in the YAML sidecar carries `source_url` (live
  page, usually arXiv abs or PDG live URL), `snapshot_path` (local
  text snapshot), `access_date`, and `sha256_of_local_snapshot`.
- Build a "Reference" view per cited source that shows:
  1. **The live URL** (clickable; opens in new tab).
  2. **The local snapshot** with the cited value highlighted in
     context. Render the snapshot as a code block / preformatted
     text; highlight the line(s) where the value appears with a
     subtle background color and an anchor link.
  3. **The sha256** displayed (for provenance integrity).
  4. **Access date** displayed (so collaborators know the snapshot
     vintage).

This means a codex agent must process each entry's snapshots and
find, for every claimed value, the EXACT character offsets (or line
numbers) where the value appears in the local text snapshot. Codex
must be EXTREMELY rigorous about this — false positives (linking to
the wrong line) destroy the whole credibility argument. Codex must:

- For each numerical claim in `pdg_or_equivalent.values`, search the
  cited `snapshot_path` for the value string (with sensible
  variants: `1.46 TeV` matches `1.46~\mathrm{TeV}`, etc.).
- If found unambiguously: record the line number(s) and a ±3-line
  context block.
- If found ambiguously (multiple matches): record ALL matches with
  context; flag for human review.
- If NOT found: do NOT silently link to a default line; flag the
  entry with `citation_resolution: UNRESOLVED` and exclude that
  specific value from the auto-linking. The website should display
  the live URL but mark the in-snapshot anchor as "unverified —
  manual confirmation needed."

The codex pipeline writes a JSON or YAML file per entry (e.g.,
`flavor_catalog/website/_data/citation_anchors/<id>.yaml`) recording
the anchor mappings. The website reads this and renders the
inline-source-view feature.

The Opus / Claude review step at the end checks a sample for false
positives.

## 5. Methodology page

A "How the catalog was built" page that:
- Embeds (or directly reproduces) `CATALOG_METHODOLOGY.tex` content.
- Links to the underlying playbook (`AGENTIC_WORKFLOW.md`).
- Briefly shows the 9-wave timeline and the 5 Opus sign-off rounds.
- Cites the 3 arbitration precedents (L001, B001+B003, B021+B023).
- Includes the fact-check summary table (per-family VERIFIED counts).

## 6. Static + Cloudflare-Pages-compatible (the only stack constraint)

- The site must build to static files. No backend, no runtime
  database.
- The build pipeline must read from the catalog (`processes/`,
  `references/`, the data files codex produces) and produce a static
  output directory Cloudflare Pages can serve.
- Picking a SSG framework, a custom build script, or anything else
  that meets the static-output bar is fine. Pick whatever ships best
  given the local toolchain you find on the machine.
- A README in the website source directory must document:
  - Local dev (whatever command starts the dev server).
  - Production build (whatever command produces the static output).
  - Cloudflare Pages config: build command, output directory, env
    vars if any.
- The build must be re-runnable: if a new entry lands on the catalog
  branch, a re-build of the website should pick it up without manual
  intervention.

# Tool split: orchestrator → codex (content) + Claude subagents (UI)

You (orchestrator Claude) write zero code and zero CSS. You dispatch.

## Codex GPT-5.5 xhigh (NOT fast) — high-volume repetitive work

Owns:
- For each of 102 entries: parse the YAML sidecar + `.tex`, extract
  structured content (prose sections, value blocks, references),
  output a normalized JSON/YAML file the website consumes.
- Citation-anchor resolution: for each value claim in
  `pdg_or_equivalent.values`, find exact line numbers + ±3-line
  context in the local snapshot at `references/<id>/`.
- Per-reference snapshot indexing: list all references + their
  snapshot paths + their live URLs + sha256s, per entry.

Dispatch pattern (matches the catalog's AGENTIC_WORKFLOW):
- Batch by family (kaon: ~16; beauty: ~26; collider_rs: ~14; etc.).
- One codex worker per batch, in parallel (max ~6 concurrent).
- Codex writes output into `flavor_catalog/website/_data/`.
- After a batch lands, you sample-audit (read 1-2 of the produced
  JSON files); if obvious errors, dispatch a corrective codex pass;
  do NOT fix inline.

Invocation:
- If `~/bin/codex_worker.sh` exists (cluster setup), use it.
- If not (local machine without wrapper), invoke `codex exec
  --dangerously-bypass-approvals-and-sandbox` directly.
- Always pipe `< /dev/null`.
- Run `~/bin/codex_check_usage.sh` if it exists, else eyeball the
  codex quota before each dispatch batch.

## Claude subagents (via Agent tool) — UI / design / authoring

Owns everything the website looks like, including framework choice,
scaffolding, components, CSS, math rendering, search UX, layout,
responsive behavior, dark mode, methodology page, Cloudflare config.

Use the Agent tool: `subagent_type: general-purpose`, with explicit
`model:` choice per task complexity:
- `opus` for high-design-judgment tasks: framework decision +
  scaffold, the per-entry detail template, the index/landing layout,
  the citation-anchor snapshot-view component, the methodology page,
  the final visual-polish pass.
- `sonnet` for routine implementation: adding a CSS rule, wiring
  one more page, adopting an existing pattern to a new component,
  fixing a specific layout issue called out in a screenshot audit.

Phase a Claude-subagent dispatch by deliverable, not by file. Each
subagent dispatch should produce a coherent, testable artifact
(e.g., "scaffold the site + landing page using catalog index
data," or "build the per-entry detail page template using one
example entry as a fixture").

Anti-patterns to avoid:
- Dispatching a subagent with a 50-line bullet list of CSS tweaks —
  too granular; the subagent loses the design vision. Group into a
  coherent "polish landing page" task instead.
- Dispatching a single subagent for "build the whole site" — too
  broad; that's what phases are for.
- Doing the work yourself "just this one small thing" — drift starts
  here; resist.

## What you (orchestrator) actually do

- Read the required-reading files (once, at start).
- Dispatch subagents per the phase plan below.
- After each subagent returns: read its short summary, sample-audit
  1-2 specific files it produced, verify on disk, commit + push.
- Maintain `flavor_catalog/website/WEBSITE_RUNBOOK.md` as the live
  dispatch ledger (mirror of the catalog's
  `wave_NNN_runbook.md` pattern): list each dispatch with
  agent type, model, deliverable, background ID, output paths,
  status. Survives compaction.
- Decide when to escalate to PI (see STOP-and-ping section).
- Visual verification: open localhost via screenshot tools (see
  Visual verification section). This IS orchestrator work because
  it's an integration judgment call, not authoring. Brief: take a
  screenshot, look at it, dispatch a subagent to fix what's wrong.

# Visual verification — Claude must actually LOOK at the rendered UI

CSS that compiles is not the same as CSS that looks good. Before
declaring any phase done, load the dev server at `http://localhost:<port>`
and inspect the rendered output, not just the source.

Concrete loop, on every meaningful UI change:

1. Start the dev server (whichever command the chosen stack uses).
   Note the port.
2. Use whatever browser-automation / screenshot tool your local
   Claude Code session has available. Common options, in order of
   preference:
   - **Playwright MCP** (if configured): take screenshots at multiple
     viewports (desktop 1440x900, tablet 768x1024, mobile 375x812).
   - **chrome-devtools-mcp** or any MCP-attached headless Chrome:
     same — screenshot at multiple viewports, inspect rendered DOM.
   - **Bash + headless Chrome**:
     `google-chrome --headless --disable-gpu --screenshot=out.png
      --window-size=1440,900 http://localhost:<port>/<path>`
     then `Read out.png` to view it.
   - **Last resort**: `curl http://localhost:<port>/<path>` to inspect
     raw HTML — useful for sanity but doesn't catch visual issues.
3. Take screenshots of at minimum: landing page, one family page,
   one PRIMARY entry detail, one SECONDARY entry detail, the
   methodology page. If you shipped dark mode, include both modes.
4. Actually look at the screenshots. Things to check:
   - Alignment: nothing visually drifting (headers not flush,
     columns misaligned, badges overlapping prose).
   - Typography: long-form content reads comfortably; hierarchy
     is clear.
   - Whitespace: not cramped, not vast.
   - Color contrast: usable for someone with imperfect vision.
   - Math rendering: LaTeX expressions actually render — no raw
     `\\varepsilon_K` or `$M_{KK}$` source showing through.
   - Citation-anchor view: when you click a reference, the snapshot
     view actually highlights the right line.
   - Mobile: nothing overflows; tap targets are reasonable;
     navigation collapses sensibly.
5. If something looks off, iterate. CSS / layout adjustments should
   loop through the visual check, not skip it.
6. Before final commit of a phase, do a full multi-page screenshot
   sweep at desktop + mobile + dark + light. If any look broken,
   fix before push.

Specific anti-patterns to catch by looking (regardless of stack):
- Cards or list items with mismatched heights breaking visual rhythm.
- Status badges wrapping awkwardly because content is too long.
- Math expressions overflowing horizontally on mobile.
- Code blocks (snapshot views) overflowing without horizontal scroll
  affordance.
- Insufficient contrast on muted SECONDARY-tier styling — muting
  shouldn't make text unreadable, just visually deprioritized.
- Theming gaps (if you ship dark mode): commonly forgotten elements
  include dropdowns, modals, code-block backgrounds, hover states.

You are the design quality bar. Codex extracts content; you make sure
the rendered result is professional. If you can't make it look good,
say so explicitly and ask for design input from the PI rather than
shipping something embarrassing.

# Deliverables

In the `flavor-catalog-website/2026q2` branch:

- `website/` (or `flavor_catalog/website/`) — full static site source.
- `website/README.md` — local dev + Cloudflare Pages deployment
  instructions.
- `flavor_catalog/website/_data/` — normalized JSON/YAML extracted by
  codex from the catalog.
- `flavor_catalog/website/_data/citation_anchors/<id>.yaml` — per-entry
  citation-anchor mappings.
- Optionally: a `Makefile` or `scripts/build-data.sh` that re-runs the
  codex data-extraction pipeline if the catalog source changes.

A working build at `https://localhost:<port>` showing:
- Landing page with all 7 family cards + SECONDARY card.
- 102 per-entry detail pages with full content.
- Citation-anchor "view in source" working for at least 90% of
  claimed values across the catalog (the other 10% can be UNRESOLVED
  with a manual-review flag).
- Light/dark mode toggle working.
- Methodology page rendering the 1-pager content.

# Quality bar

The website is collaborator-facing. Treat it as production-quality:
- No placeholder content. No "lorem ipsum." No broken links.
- All 102 entries render without errors.
- All citation anchors that ARE resolved actually point at the right
  line in the snapshot (random-sample audit by Opus before declaring
  done).
- Mobile + desktop both work.
- Lighthouse / accessibility audit reasonable (a11y matters for
  external collaborators using screen readers).

# Pre-flight when you start

```bash
cd /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing
git status -sb
git log --oneline -3
git ls-remote --tags origin | grep "flavor-catalog"
# Confirm flavor-catalog-v0.4 tag exists and you're on
# flavor-catalog/2026q2 at commit 2bda5f1 or later.

# Cut the website branch
git checkout -b flavor-catalog-website/2026q2

# Verify codex available (if on cluster):
~/bin/codex_check_usage.sh 2>/dev/null && echo "codex OK" || echo "no codex wrapper (use codex exec directly)"

# Optional: check available toolchains (you decide what to use)
node --version 2>/dev/null || echo "no node"
python3 --version 2>/dev/null || echo "no python3"
which hugo 2>/dev/null || echo "no hugo"
```

# Iteration (each phase = one or more delegated dispatches)

This is multi-day, multi-phase work. The orchestrator dispatches and
reviews; the orchestrator does not author. Each phase below names
WHO is dispatched and what they produce — your role is to set up
each dispatch, wait for the notification, sample-audit the output,
commit + push, then move on.

1. **Phase 1 — Scaffold (Claude subagent, opus, ~1-2 hr wall)**:
   one Claude subagent picks the stack (their call — any static-site
   approach is fine), scaffolds the site, hand-builds ONE example
   entry's detail page to define the per-entry template, and
   produces the landing layout. Returns a working dev-server
   localhost. Orchestrator action: dispatch subagent → review
   screenshot → commit + push.

2. **Phase 2 — Codex extraction (codex, ~3-6 hr wall, parallel)**:
   dispatch ~6 codex agents in parallel (one per family + the
   secondary tree). Each produces normalized JSON/YAML per entry
   plus citation-anchor data files. Total: 102 entries' content
   + ~600 citation anchors resolved.
   Orchestrator action: dispatch codex batches → wait for
   notifications → sample-audit 2 entries per family → commit.

3. **Phase 3 — Data-driven rendering (Claude subagent, opus,
   ~2-3 hr)**: one Claude subagent reads the codex data and wires
   it into the template from Phase 1. Implements family pages,
   search, filters, status badges, citation-anchor view, dark
   mode, mobile responsive layout. Returns a complete 102-page
   localhost site.
   Orchestrator action: dispatch → screenshot multi-page audit
   → if issues, dispatch a `sonnet` subagent to fix specific
   issues → commit.

4. **Phase 4 — Citation audit (orchestrator, direct, ~30 min)**:
   orchestrator does this directly because it's pure judgment
   work, no authoring. Pick 5-10 random entries; for each,
   click through to the cited reference; confirm the highlighted
   line in the local snapshot actually contains the value.
   If false positives found: dispatch codex to re-resolve anchors
   for that family.
   Orchestrator action: small bounded review.

5. **Phase 5 — Methodology page + deployment config (Claude
   subagent, opus, ~1 hr)**: one Claude subagent authors the
   methodology page (from `CATALOG_METHODOLOGY.tex` + the relevant
   sections of `AGENTIC_WORKFLOW.md`), writes the website README,
   and adds the Cloudflare Pages configuration (build command,
   output dir, env vars if any).
   Orchestrator action: dispatch → review → commit + push.

6. **Phase 6 — PI review**: PI visually reviews at localhost,
   approves the Cloudflare deployment.

After each phase, orchestrator commits + pushes to
`flavor-catalog-website/2026q2`. Update `WEBSITE_RUNBOOK.md` with
the dispatch row + outcome. Do NOT push to `flavor-catalog/2026q2`
or any frozen branch.

# When to STOP and ping the PI

- If codex citation-anchor resolution fails on > 20% of values across
  any family — that's a signal the snapshot quality varies and the
  PI may want to commission targeted fixes before website build.
- If you find a value in the catalog that doesn't match the linked
  source after careful checking — that's a fact-check escape; flag
  to PI before silently fixing.
- If a framework decision feels load-bearing (e.g., choosing
  Cloudflare Pages vs Cloudflare Workers, or static vs SSR), confirm
  with the PI before committing.
- If you want to add scope not in this prompt (live API, comments,
  editing UI, etc.), confirm with PI first. This is read-only.

# Read the required-reading docs first, then summarize

Before scaffolding anything: read the 6 required-reading files above
and summarize back to the PI in 5-10 sentences what you understand:
- What the catalog contains
- What the credibility chain looks like (the 5-pass pipeline)
- What "exact reference linking" means in practice
- Your recommended framework choice and why
- A rough Phase-1 plan

Then wait for the PI's go-ahead before scaffolding.
```

---

## What I added beyond the PI's brief

A few things the PI didn't explicitly specify but that materially affect
the deliverable:

1. **Pure-orchestrator role** for the new Claude. They dispatch codex
   (for repetitive per-entry work) and Claude subagents via the Agent
   tool (for UI / design / authoring) but write no code, CSS, or
   per-entry content themselves. Context-window preservation across
   the multi-day, multi-phase build.
2. **Dedicated branch `flavor-catalog-website/2026q2`** off the catalog
   branch — keeps website code from polluting the (frozen-at-v0.4)
   catalog content branch.
3. **Citation-anchor data format** — codex writes
   `_data/citation_anchors/<id>.yaml` files; UI subagents render. Clean
   split between "find exact line" (codex) and "show it nicely"
   (UI subagent).
4. **UNRESOLVED handling** — codex must NOT silently link to a wrong
   line if it can't find the value unambiguously. Mark as unresolved
   with manual-review flag. False-positive citations would destroy the
   whole credibility argument.
5. **Status badges + fact-check verdicts** — surface OPUS-APPROVED,
   fact-check verdict, cycle count, code coverage, implementation
   difficulty per entry. These signals already exist in the YAML;
   collaborators should see them at a glance.
6. **Tier visualization requirement** — SECONDARY entries visible but
   visually distinct. Honest framing per the PI's directive on the
   catalog itself.
7. **Methodology page** — surfaces `CATALOG_METHODOLOGY.tex` content
   as a web page so collaborators can read the verification pitch
   without downloading a PDF.
8. **Provenance display** — sha256 and access date shown per cited
   source. Visible provenance is part of the credibility argument.
9. **Code-coverage hyperlinks** — file:line evidence in entries should
   be GitHub-blob-anchor URLs, so a collaborator can verify the
   coverage claim by clicking through.
10. **Random-sample citation audit** — a final orchestrator pass on
    5–10 entries to confirm anchors land correctly. Sampling is cheap;
    false positives are expensive.
11. **6-phase iteration plan** with rough time estimates and an explicit
    "PI visual review at localhost before pushing" gate.
12. **STOP-and-ping conditions** — when to escalate vs proceed,
    mirroring the catalog's existing orchestration protocol.
13. **Pre-flight checklist** — confirms environment before scaffolding
    starts.
14. **Visual-verification loop** (per PI follow-up): orchestrator must
    actually open `http://localhost:<port>` via browser-automation /
    screenshot tools, view the rendered output at multiple viewports,
    and iterate based on what it sees. CSS that compiles ≠ CSS that
    looks good. Specific anti-patterns listed.

## What I deliberately did NOT specify

The PI explicitly asked me to leave design and stack choices to the
new Claude. The following are NOT specified in the prompt:

- **Framework / stack**. Astro, Next.js, Hugo, Eleventy, plain Vite,
  custom Python+Jinja, anything else — the new Claude picks. The only
  hard constraint is static output for Cloudflare Pages.
- **Math rendering library**. KaTeX, MathJax, or something else — pick
  what fits the chosen stack.
- **CSS approach**. Tailwind, vanilla CSS, CSS modules, a UI kit,
  anything — the new Claude decides.
- **Color palette, typography, density**. All design choices.
- **Dark mode**. Encouraged but optional; if it ships, must be done
  well in both modes.
- **Search vs facets**. Search box is the minimum; faceted filters
  are an option.
- **Navigation structure**. Family-based, sensitivity-based,
  hybrid — the new Claude picks.
- **Landing-page layout**. Cards, list, table, mixed — design call.

What IS specified is content + structure + quality bar:
- Which fields per entry must surface (process ID, notation, status
  badges, value blocks, RS-relevance prose, post-2008 dev,
  validity / model dependence, code coverage, key references,
  provenance footer).
- Citation-exact reference linking with UNRESOLVED handling.
- PRIMARY vs SECONDARY tier visible.
- Methodology page reproducing `CATALOG_METHODOLOGY.tex`.
- Static build deployable to Cloudflare Pages.
- Visual verification loop (orchestrator screenshots + iterates).
