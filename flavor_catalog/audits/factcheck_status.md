# Flavor Catalog — Fact-Check Status (v0 → v0.1)

**Goal**: independently verify that each cited paper EXISTS *and* actually
CONTAINS the numerical value(s) the catalog entry claims to take from it.
Not just URL HTTP-200; the paper content must match.

**Method**: per-family codex (gpt-5.5 xhigh) fact-checkers with WebFetch.
Each agent reads its family's `.yaml` source_manifest, fetches the cited
arXiv / PDG / HFLAV / FLAG / experiment-collaboration pages, and verifies the
claimed numerical values appear at the cited source.

**Prior review credit**: the original PKA/WA/CA cycles + Opus round-1/round-2
sign-offs verified internal consistency (every TeX claim resolves to a manifest
entry; sha256 of snapshot matches; status_history transitions clean). They did
NOT independently re-fetch every source. This pass closes that gap.

**Status legend**:
- `PENDING` — not yet fact-checked.
- `VERIFIED` — all cited sources fetched; all claimed values match.
- `PARTIAL` — most sources verified; some unresolvable (e.g. paywalled or 404).
- `MISMATCH` — at least one claimed value does not match the cited source.
- `FAILED` — multiple claims unverifiable.

Each per-family fact-checker writes a detailed report at
`flavor_catalog/audits/factcheck_<family>.md`. Then an aggregator compiles
the rows below.

## Per-process status table

(Populated by aggregator after all 6 family fact-checkers land.)

| process_id | family | status | date | fact-checker | mismatches/notes |
|---|---|---|---|---|---|

## Summary

(Populated by aggregator.)
