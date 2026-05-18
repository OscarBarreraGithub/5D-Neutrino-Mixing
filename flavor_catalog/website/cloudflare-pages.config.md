# Cloudflare Pages deployment configuration

This document captures the exact values to enter in the Cloudflare Pages
dashboard when connecting the GitHub repo for build-on-push deployment.
The repo is a Python physics codebase; the website is a subdirectory.

## Connect repository

- **Repository**: `OscarBarreraGithub/5D-Neutrino-Mixing`
- **Production branch (current)**: `flavor-catalog-website/2026q2`
  - Will change to `main` after merge.
- **Preview branches**: any other branch (default Cloudflare Pages behavior).

## Build configuration

| Field | Value |
|---|---|
| Framework preset | `None` (Astro preset works too, but explicit values below are authoritative) |
| Build command | `npm install && npm run build` |
| Build output directory | `dist` |
| Root directory (advanced) | `flavor_catalog/website` |
| Node version | `22.x` |

Notes:

- `npm run build` is defined in `package.json` as
  `astro build && pagefind --site dist`. It runs both the static-site
  generator and the Pagefind index builder in one shot.
- The Node version is pinned via `.node-version` (`22`) in the website root.
  Cloudflare Pages reads this automatically. You may also set the
  `NODE_VERSION` environment variable to `22` as a belt-and-braces fallback.
- Astro requires Node `^18.20.8 || ^20.3.0 || >=22.0.0`. Local development
  has been run on Node 25.6.1; 22 LTS is the recommended deployment target.

## Environment variables

| Name | Value | Scope |
|---|---|---|
| `NODE_VERSION` | `22` | Production + Preview (optional belt-and-braces) |

No secrets are required: the site is fully static and contains no runtime
API calls.

## Output structure

After `npm run build`, the `dist/` directory contains:

- `index.html`, landing page
- `entries/<process_id>/index.html`, one page per entry (102 total)
- `families/<family>/index.html`, per-family index pages
- `browse/index.html`, filterable table with Pagefind search
- `methodology/index.html`, "How this catalog was built"
- `pagefind/`, Pagefind static search index (fragments, bundle, UI)
- `_astro/`, hashed JS/CSS bundles
- `_redirects`, Cloudflare Pages static-redirects file (copied from
  `website/_redirects` via Astro's `public/` semantics; see below)
- `_headers`, Cloudflare Pages static-headers file

Total pages: **112+** (102 entries + 8 family indices + browse + methodology
+ landing).

## `_redirects` and `_headers`

These two files live in `flavor_catalog/website/public/` and Astro copies
them verbatim into `dist/` at build time (Astro's `public/` semantics).
Cloudflare Pages reads them at deploy time.

- `public/_redirects` ensures bare `/entries/<id>` and `/families/<f>`
  URLs redirect to their trailing-slash canonical forms (HTTP 301).
  Astro generates trailing-slash URLs by default; without this redirect,
  deep links missing the slash would 404.
- `public/_headers` sets long-cache headers on hashed `/_astro/*` and
  `/pagefind/*` bundles, no-cache on HTML, and sane security headers
  site-wide.

Verify after `npm run build`:

```bash
ls dist/_redirects dist/_headers
```

## Verifying deployment locally

```bash
cd flavor_catalog/website
npm install
npm run build
ls dist/                    # expect index.html, entries/, families/, methodology/, browse/, pagefind/, _astro/, _redirects, _headers
python3 -m http.server --directory dist 8080
open http://localhost:8080/methodology/
```

Pagefind search will not work under `python3 -m http.server` for paths
that need MIME-type sniffing; use `npx http-server dist -p 8080` for a
production-equivalent check, or `npm run preview` which uses Astro's
preview server.

## Custom domain (optional)

The Pages project exposes a `*.pages.dev` URL by default. To attach a
custom domain (e.g. `flavor-catalog.example.org`), use the dashboard's
Custom domains tab; Cloudflare will issue a Universal SSL cert
automatically.

## Branch protection reminder

- Production branch (`flavor-catalog-website/2026q2`) is currently the
  active branch. After merging to `main`, update the Cloudflare Pages
  project settings to point at `main`.
- Do **not** allow direct pushes to the production branch without the
  pre-merge build check passing.
