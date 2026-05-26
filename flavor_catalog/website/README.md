# Flavor Catalog -- Static Website

The Astro 5.x static site that surfaces the 102-entry flavor-physics
catalog (`../`) as a browsable, filterable web artifact for collaborators.
It renders one page per process, per-family indices, a Pagefind-powered
browse/filter view, and a "How this catalog was built" methodology page.

The catalog itself (YAML sidecars + `.tex` prose + reference manifests +
text snapshots) is the source of truth; this site is a read-only
projection. See the parent repo's
[`flavor_catalog/`](../) for the source data and
[`WEBSITE_BUILD_PROMPT.md`](../WEBSITE_BUILD_PROMPT.md) for the multi-phase
build history.

---

## Local development

### Prerequisites

- **Node 22+** (Astro 5 requires `^18.20.8 || ^20.3.0 || >=22.0.0`;
  Cloudflare Pages production target is 22 LTS, pinned in `.node-version`).
- **Python 3** (the ingest script has no third-party requirements; PyYAML
  is preferred if present, otherwise a minimal hand-rolled YAML parser is
  used).

### One-time setup

```bash
cd flavor_catalog/website
npm install
```

### Regenerate the content layer

The Astro pages read three JSON files generated from the catalog source:

- `src/content/entries/<process_id>.json` -- one per entry
- `src/content/catalog_index.json` -- summary used by landing + browse
- `src/content/families.json` -- per-family counts/labels

Re-run the ingest whenever a catalog entry changes (`.yaml` sidecar,
`.tex` body, or `source_manifest.yaml`):

```bash
python3 scripts/ingest_catalog.py
```

The ingest also stamps **citation anchors** from
`_data/citation_anchors/*.yaml` onto each entry JSON, so the entry
page's "view in source" modal works without an extra client fetch.

### Regenerate citation anchors (rare)

Citation anchors are normally resolved by a codex (gpt-5.4 xhigh) Phase-2
job that reads the YAML manifests, fetches snapshots, and writes one
YAML file per family under `_data/citation_anchors/`. If you need to
re-resolve a family locally (after a snapshot edit, for example), the
generic resolver lives at `scripts/resolve_citation_anchors.py`:

```bash
python3 scripts/resolve_citation_anchors.py --family beauty
```

For a fresh codex-driven re-run, the canonical incantation is in
[`WEBSITE_RUNBOOK.md`](./WEBSITE_RUNBOOK.md) (Phase 2 dispatch).

### Populate the Pagefind index (one-time, for dev-server search)

The dev server does not run Pagefind on the fly. Run a production build
once to populate `public/pagefind/`:

```bash
npm run build
```

Subsequent `npm run dev` invocations reuse the populated `public/pagefind/`
copy, so the browse-page search works locally.

### Start the dev server

```bash
npm run dev
# or, matching the orchestrator setup exactly:
npm run dev -- --host 127.0.0.1 --port 4321
```

Visit <http://localhost:4321/>.

---

## Production build

```bash
npm run build
```

This runs `astro build && pagefind --site dist`, producing:

- `dist/index.html` -- landing
- `dist/entries/<id>/index.html` -- 102 entry pages
- `dist/families/<key>/index.html` -- 8 family indices
- `dist/browse/index.html` -- filterable table + Pagefind search
- `dist/methodology/index.html` -- "How this catalog was built"
- `dist/pagefind/` -- static search index (fragments, bundle, UI)
- `dist/_astro/` -- hashed JS/CSS bundles
- `dist/_redirects` -- Cloudflare redirects (trailing-slash fix)
- `dist/_headers` -- Cloudflare cache headers

Verify locally:

```bash
npm run preview            # Astro preview server on :4322
# or
python3 -m http.server --directory dist 8080
```

---

## Cloudflare Pages deployment

Full configuration lives in [`cloudflare-pages.config.md`](./cloudflare-pages.config.md).
Key values:

| Field | Value |
|---|---|
| Repository | `OscarBarreraGithub/5D-Neutrino-Mixing` |
| Production branch | `flavor-catalog-website/2026q2` (move to `main` after merge) |
| Build command | `npm install && npm run build` |
| Build output directory | `dist` |
| Root directory | `flavor_catalog/website` |
| Node version | `22.x` (pinned in `.node-version`) |

`public/_redirects` and `public/_headers` are picked up automatically by
Cloudflare Pages from the build output.

---

## Architecture

### Data flow

```
flavor_catalog/processes/**/*.{yaml,tex}     (catalog source of truth)
              │
              ▼
scripts/ingest_catalog.py                    (Python, no third-party deps)
              │
              ├── reads _data/citation_anchors/<family>.yaml
              │   (stamps anchor metadata onto each entry)
              │
              ▼
src/content/entries/<id>.json                (one JSON per entry)
src/content/catalog_index.json               (rollup for landing/browse)
src/content/families.json                    (family chips)
              │
              ▼
src/pages/{index,browse,methodology}.astro   (Astro static SSG)
src/pages/entries/[id].astro                 (dynamic route per entry)
src/pages/families/[family].astro            (dynamic route per family)
              │
              ▼
dist/                                        (static output)
              │
              ▼
pagefind --site dist                         (full-text index)
              │
              ▼
dist/pagefind/                               (search bundle, content-addressed)
```

### Citation anchors

A per-value provenance trigger appears on every entry-detail page. The
flow:

1. **codex (Phase 2)** reads each family's `source_manifest.yaml` and the
   tracked text snapshots, then writes
   `_data/citation_anchors/<family>.yaml` with per-anchor RESOLVED /
   AMBIGUOUS / UNRESOLVED status plus context lines.
2. **`ingest_catalog.py`** loads these YAMLs and stamps anchor metadata
   onto each entry JSON.
3. **`src/pages/entries/[id].astro`** renders the trigger button next to
   each value; the modal in `src/components/CitationModal.astro`
   displays the snapshot lines with the hit row highlighted.

Phase 2 produced 644 anchors across 102 entries: 73.9% RESOLVED, 18.8%
AMBIGUOUS, 7.3% UNRESOLVED. See `/methodology/` for the per-family
UNRESOLVED breakdown.

### Per-phase dispatch ledger

The website was built in 6 phases dispatched by an orchestrator. See
[`WEBSITE_RUNBOOK.md`](./WEBSITE_RUNBOOK.md) for the per-phase decisions,
audit trail, and recovery instructions.

---

## Re-running citation anchors

If a catalog snapshot changes and an anchor goes stale:

1. **Local quick path** (single family, no codex):
   ```bash
   python3 scripts/resolve_citation_anchors.py --family <family>
   ```
   Writes `_data/citation_anchors/<family>.yaml`. Then re-run
   `python3 scripts/ingest_catalog.py` to re-stamp.

2. **Canonical codex path** (multi-family, fresh resolution): see
   [`WEBSITE_RUNBOOK.md`](./WEBSITE_RUNBOOK.md) for the gpt-5.4 xhigh
   dispatch template and the family-aware prompt with the orchestrator's
   audit-trail conventions.

After either path, run `npm run build` to confirm 112+ pages still build
cleanly.

---

## Tested on

- Node `25.6.1` (local dev)
- Node `22.x` (Cloudflare Pages deployment target)
- Python `3.14.3` (ingest)
- Astro `5.18.1`
- Pagefind `1.5.2`
- macOS `14`
