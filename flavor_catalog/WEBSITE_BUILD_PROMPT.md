# Cold-Boot Prompt: Build the Flavor Catalog Website

Paste the fenced block below into a fresh local Claude Code session. It's
self-contained — points the new Claude at the catalog source, the
verification trail, and the split between Claude (UI) and codex (per-entry
content extraction).

---

```
You are taking over a NEW task: build a clean, professional static
website that surfaces the 5D-Neutrino-Mixing flavor catalog (currently
102 OPUS-APPROVED entries at tag `flavor-catalog-v0.4`). The catalog
itself is built and verified — your job is purely to make it
collaborator-readable on the web.

This will be hosted on Cloudflare Pages via the standard "Cloudflare
Pages connects to a GitHub repo and rebuilds on push" workflow. For
now you only need to run it locally (localhost dev server); deployment
config is the final-step deliverable, not the immediate goal.

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

## 1. Clean, minimalist, professional UI
- Collaborator audience (physicists). Looks credible without being
  flashy.
- Typography-driven. Generous whitespace. No marketing graphics.
- Dark / light mode toggle, system-preference default.
- Mobile responsive (collaborators read on phones).
- Math rendering via KaTeX or MathJax (KaTeX preferred — faster, no
  network dependency once bundled).

## 2. Catalog index with tiered organization

Surface the priority tiering visibly. Suggested layout (you decide
the best UX):

- **Landing page**: brief catalog overview (count, tag, fact-check
  status); 6 family cards for PRIMARY low-energy; 1 family card for
  PRIMARY new-family `collider_rs`; 1 muted card for SECONDARY
  re-promotions. Each card lists its entries with a one-line process
  description.
- **Family page**: lists all entries in that family sorted by
  process ID; each entry shows ID, standard notation, brief
  description, status badges (OPUS-APPROVED, fact-check verdict,
  cycle count if >1, implementation difficulty).
- **Filtering / search**: at minimum, a search box that matches on
  process ID, standard notation, plain-name, and key references.
  Optionally faceted filters: family, tier, implementation
  difficulty, code coverage status, fact-check verdict.
- **Tier visualization**: SECONDARY entries should be visually
  distinct (e.g., muted color or a "SECONDARY" badge) but NOT hidden.
  The point is collaborators can see the full set and immediately
  know which tier each entry is in.

You decide the precise visual hierarchy (highest-impact → lower-impact)
based on the data available. Reasonable defaults that should inform
your layout:
- PRIMARY > SECONDARY (tier).
- Within PRIMARY: low-energy flavor (Waves 1-7) is the
  highest-leverage tier (~47 TeV M_KK reach); `collider_rs` Wave-9 is
  PRIMARY but lower-reach (~few TeV). The website should reflect
  this honestly — don't oversell the collider tier.
- Within a family: implementation difficulty LOW > MEDIUM > HIGH > BLOCKED
  is a reasonable impact proxy.
- Within a family: existing in-code coverage (YES > PARTIAL > NO) is
  another impact proxy.
- Family ordering on the landing page can follow physics tradition:
  kaon, charm, beauty, top/Higgs/EW, charged_lepton, edm_neutrino,
  collider_rs.

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

## 6. Static site, regeneratable, deployable to Cloudflare Pages

- Pure static site (no backend, no database).
- Build command produces a `dist/` or `build/` or `out/` directory
  that Cloudflare Pages can serve.
- Source files in a top-level `website/` or `flavor_catalog/website/`
  subdirectory.
- README in that directory documenting:
  - Local dev: `npm run dev` (or equivalent) starts on localhost.
  - Production build: `npm run build`.
  - Cloudflare Pages config: build command, output directory, env
    vars if any.
- The build script must read the catalog's YAML sidecars + .tex
  prose + citation_anchors data and render the site. If a new
  entry lands on the catalog branch, the website rebuilds
  automatically on next deploy.

# Tool split: Claude UI / codex content

Use Codex GPT-5.5 xhigh (NOT fast) for repetitive per-entry work.
Use Claude (you) for UI/design/framework/integration. Specifically:

## Claude (you) owns
- Framework selection (recommend: Astro for static site with
  content-heavy pages; alternative: Next.js with `output: 'export'`;
  alternative: Hugo if you want zero-JS-framework simplicity).
- Site scaffolding, layout, CSS / Tailwind / vanilla styling.
- Light/dark mode, mobile responsive design, navigation.
- Search / filter UX.
- Methodology page authoring.
- Per-entry detail page TEMPLATE (one Astro/Next component, used 102
  times).
- Citation-anchor RENDERING (read the JSON codex produces, render the
  snapshot-view component).
- Code-review of the codex outputs before integration.
- Cloudflare Pages config.

## Codex (delegated) owns
- For each of 102 entries: parse the YAML sidecar + .tex, extract
  structured content (prose sections, value blocks, references), and
  output a normalized JSON/YAML file the website can consume.
- Citation-anchor resolution: for each value claim, find exact line
  numbers + ±3-line context in the local snapshot.
- Per-reference snapshot indexing: produce a list of all references
  + their snapshot paths + their live URLs + sha256s.

Dispatch pattern (matches the catalog's existing AGENTIC_WORKFLOW):
- Batch the entries by family (kaon: ~16 batched; beauty: ~26; etc.).
- One codex worker per batch, in parallel (max ~6 concurrent).
- Codex writes its output files into `flavor_catalog/website/_data/`.
- After codex batches land, Claude reads the data, fixes
  obvious errors, re-dispatches codex for outliers.
- Final Opus or Claude sample-audit on 5-10 random entries to confirm
  citation anchors point at the right lines.

If `~/bin/codex_worker.sh` exists (cluster setup), use it. If not
(local machine without the wrapper), invoke `codex exec
--dangerously-bypass-approvals-and-sandbox` directly. Always pipe
`< /dev/null`. Always run `~/bin/codex_check_usage.sh` (if it exists)
or otherwise eyeball the codex quota before each dispatch batch.

# Visual verification — Claude must actually LOOK at the rendered UI

CSS that compiles is not the same as CSS that looks good. Before
declaring any phase done, load the dev server at `http://localhost:<port>`
and inspect the rendered output, not just the source.

Concrete loop, on every meaningful UI change:

1. Run the dev server (`npm run dev` or equivalent). Note the port.
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
   methodology page. Both light mode and dark mode.
4. Actually look at the screenshots. Check for:
   - Alignment: nothing visually drifting (e.g., headers
     not flush, columns misaligned, badges overlapping prose).
   - Typography: line lengths comfortable (45-75ch for body text);
     hierarchy clear (H1 > H2 > H3 visually distinct);
     monospace where appropriate (file paths, sha256, IDs).
   - Whitespace: not cramped, not vast.
   - Color contrast: WCAG AA at minimum in both modes (light AND
     dark — dark mode often fails this).
   - Math rendering: KaTeX expressions actually render, no raw
     `\\varepsilon_K` showing through.
   - Citation-anchor view: when you click a reference, the snapshot
     view actually highlights the right line.
   - Mobile: nothing overflows; tap targets are ≥ 44px;
     navigation collapses sensibly.
5. If something looks off, iterate. CSS / layout adjustments should
   loop through the visual check, not skip it.
6. Before final commit of a phase, do a full multi-page screenshot
   sweep at desktop + mobile + dark + light. If any look broken,
   fix before push.

Specific anti-patterns to catch by looking:
- Cards with mismatched heights breaking the grid.
- Status badges wrapping awkwardly because content is too long.
- Math expressions overflowing horizontally on mobile.
- Code blocks (snapshot views) overflowing without horizontal scroll
  affordance.
- Insufficient contrast on muted SECONDARY-tier styling (the muting
  shouldn't make text unreadable, just visually deprioritized).
- Dark mode forgetting any element — common: dropdowns, modals,
  code-block backgrounds, hover states.

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

# Confirm node/npm if you plan to use Astro/Next:
node --version
npm --version
```

# Iteration

This is multi-session work, not one-shot. Reasonable phasing:

1. **Phase 1 (Claude, ~1-2 hours)**: framework decision, scaffold,
   one example entry hand-coded to design the per-entry template.
2. **Phase 2 (codex, ~3-6 hours)**: dispatch codex batches to extract
   structured content + citation anchors for all 102 entries. Run
   in parallel by family.
3. **Phase 3 (Claude, ~2-3 hours)**: data-driven rendering of the
   full site from codex output. Search, filters, dark mode,
   responsive design.
4. **Phase 4 (Opus / Claude, ~30 min)**: sample audit of 5-10
   citation anchors. Fix any false positives.
5. **Phase 5 (Claude, ~1 hour)**: methodology page, Cloudflare Pages
   config, README.
6. **Phase 6 (PI)**: visual review at localhost, then push and
   approve Cloudflare deployment.

After each phase, commit + push to `flavor-catalog-website/2026q2`.
Do NOT push to `flavor-catalog/2026q2` or any frozen branch.

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

A few things the PI didn't specify that I baked in because they materially
affect the deliverable:

1. **Dedicated branch `flavor-catalog-website/2026q2`** off the catalog
   branch — keeps website code from polluting the (frozen-at-v0.4) catalog
   content branch.
2. **Framework recommendation** (Astro primary, Next.js export or Hugo as
   fallbacks). Astro is the right pick for content-heavy static sites and
   has first-class Cloudflare Pages support.
3. **Citation-anchor data format** — codex writes
   `_data/citation_anchors/<id>.yaml` files Claude renders. Separates the
   "find exact line" work (codex) from the "render it nicely" work (Claude).
4. **UNRESOLVED handling** — codex must NOT silently link to a wrong line
   if it can't find the value unambiguously. Mark as unresolved + manual
   review instead. False-positive citations would destroy the whole
   credibility argument.
5. **Status badges + fact-check verdicts** — surface OPUS-APPROVED, fact-
   check verdict, cycle count, code coverage, implementation difficulty per
   entry. These signals already exist in the YAML; collaborators should see
   them at a glance.
6. **Tier visualization rule** — SECONDARY entries are muted but NOT hidden.
   Honest framing per the PI's earlier directive on the catalog itself.
7. **Methodology page** — surfaces `CATALOG_METHODOLOGY.tex` content as a
   web page so collaborators can read the verification pitch without
   downloading a PDF.
8. **Provenance display** — sha256 and access date shown per cited source.
   Visible provenance is part of the credibility argument.
9. **Math rendering** — KaTeX recommended (faster than MathJax, no network
   dep). Essential for the LaTeX physics notation throughout the catalog.
10. **Code-coverage hyperlinks** — file:line evidence in entries should be
    GitHub-blob-anchor URLs, so a collaborator can verify the coverage claim
    by clicking through.
11. **Random-sample citation audit** — a final Opus/Claude pass on 5-10
    entries to confirm anchors land correctly. Sampling is cheap, false
    positives are expensive.
12. **6-phase iteration plan** with rough time estimates and an explicit
    "PI visual review at localhost before pushing" gate.
13. **STOP-and-ping conditions** — when to escalate vs proceed, mirroring
    the catalog's existing orchestration protocol.
14. **Mobile + a11y requirements** — collaborators read on phones and
    accessibility matters for screen readers.
15. **Light/dark mode** — small UX win that signals professionalism.
16. **Pre-flight checklist** — confirms environment before scaffolding
    starts. Matches the catalog's existing handoff pattern.
17. **Visual-verification loop** (added per PI follow-up) — Claude must
    actually open `http://localhost:<port>` via browser-automation /
    screenshot tools, view the rendered output at multiple viewports
    in both light + dark mode, and iterate based on what it sees. CSS
    that compiles ≠ CSS that looks good. Anti-patterns to look for
    (mismatched card heights, badge wrapping, math overflow on mobile,
    dark-mode contrast failures, etc.) are listed explicitly.
