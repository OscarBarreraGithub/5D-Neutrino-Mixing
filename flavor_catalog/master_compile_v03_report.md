# Flavor Catalog v0.3 Master Compile Report

Date: 2026-05-17
Branch: flavor-catalog/2026q2

## Scope

- Catalog version: flavor-catalog-v0.3
- Process coverage: PRIMARY = 80 (unchanged), SECONDARY = 8 (Wave-8 promotions), total = 88
- Family counts:
  - kaon: PRIMARY 13, SECONDARY 3, total 16
  - charm: PRIMARY 8, SECONDARY 0, total 8
  - beauty: PRIMARY 22, SECONDARY 4, total 26
  - top/Higgs/EW: PRIMARY 19, SECONDARY 1, total 20
  - charged lepton (LFV): PRIMARY 11, SECONDARY 0, total 11
  - EDM/neutrino: PRIMARY 7, SECONDARY 0, total 7
- Wave-8 SECONDARY additions included in the master indexes: K019, K020, K021, B007, B008, B013, B014, T014
- Per-process `.tex` and `.yaml` content was not modified during this master rebuild.
- `catalog_master.tex` was not modified; the existing SECONDARY section was used as wired.

## Index Rebuild

Each `flavor_catalog/processes/<family>/index.tex` and each `flavor_catalog/processes/secondary/<family>/index.tex` was regenerated from the sorted `<ID>.tex` files present in that family directory. Each entry uses an `\IfFileExists` guard around the corresponding `\input`.

All nine regenerated indexes matched the pre-existing checked-in content, so no index file content changed in this rebuild.

## Build Command

```bash
cd /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog
rm -f *.aux *.log *.out *.toc *.synctex.gz
pdflatex -interaction=nonstopmode catalog_master.tex
pdflatex -interaction=nonstopmode catalog_master.tex
pdfinfo catalog_master.pdf | head -20
rm -f *.aux *.log *.out *.toc *.synctex.gz
```

## Build Result

- Result: success
- Output PDF: `flavor_catalog/catalog_master.pdf`
- Page count: 148
- PDF size reported by `pdfinfo`: 796017 bytes
- LaTeX notes: build completed with overfull/underfull hbox layout warnings from long monospace paths and references, plus hyperref PDF-string warnings for math tokens in section/bookmark strings. No process inputs were skipped.

## Approval and Fact-Check Summary

- PRIMARY approval status unchanged: 80/80 OPUS-APPROVED.
- PRIMARY fact-check status unchanged: 79 VERIFIED, 1 PARTIAL (E009 INSPIRE JS-only; APS journal cross-check confirms content), 0 MISMATCH, 0 FAILED.
- SECONDARY approval status: 8/8 OPUS-APPROVED per `flavor_catalog/signoff/round_004_index.md`.
- SECONDARY fact-check status from Wave-8 addendum sections in `flavor_catalog/audits/factcheck_*.md`: 7 VERIFIED, 1 PARTIAL (K020 metadata-only author-line convention drift), 0 MISMATCH, 0 FAILED.
- Total catalog state after Wave-8: 88/88 OPUS-APPROVED, with 86 VERIFIED and 2 PARTIAL across PRIMARY plus SECONDARY.

## Priority Tier Cross-Reference

SECONDARY entries are lower implementation priority by policy. See `flavor_catalog/PRIORITY_TIERS.md`, especially the Wave-8 rationale table and implementation-priority rule.

## Consolidation status (added by C07 cleanup, 2026-05-25)

Reconciliation of the v0.3 totals (`88/88 OPUS-APPROVED`, `86 VERIFIED + 2 PARTIAL`) against the per-family fact-check audits at the v0.3 boundary (tag `flavor-catalog-v0.3`):

| Source | Entries | VERIFIED | PARTIAL | Notes |
|---|---:|---:|---:|---|
| Consolidated table in `audits/factcheck_status.md` (rows `:32-108`) | 75 | 74 | 1 | Unchanged from v0.1; PARTIAL row is E009. |
| Wave-7 PRIMARY addenda in `audits/factcheck_top_higgs_ew.md:294-352` + `factcheck_beauty.md` (B012 row) | 5 | 5 | 0 | T003, T004, T008, T012, B012 — all VERIFIED. |
| Wave-8 SECONDARY addenda in `audits/factcheck_kaon.md:254-313`, `factcheck_beauty.md:243-310`, `factcheck_top_higgs_ew.md:391-427` | 8 | 7 | 1 | K019/K021/B007/B008/B013/B014/T014 VERIFIED; K020 PARTIAL at the v0.3 tag (metadata-only NA62 author-line convention drift; cleared to VERIFIED by post-tag cycle-3 cleanup commit `b5c2375`). |
| **Total at v0.3 boundary** | **88** | **86** | **2** | Matches headline `86 VERIFIED + 2 PARTIAL` in §"Approval and Fact-Check Summary" above. |

The v0.3 consolidated `audits/factcheck_status.md` table was inherited unchanged from v0.1 (75 rows) and was not regenerated to fold in the 5 Wave-7 PRIMARY + 8 Wave-8 SECONDARY additions. The Wave-7 and Wave-8 verdicts are individually authoritative in the per-family addendum sections cited above; the headline `86 VERIFIED + 2 PARTIAL` figure is the sum of the consolidated table (74 VERIFIED + 1 PARTIAL) plus the Wave-7 addenda (5 VERIFIED) plus the Wave-8 addenda (7 VERIFIED + 1 PARTIAL).

Note on the v0.3 tag annotation: `git cat-file -p flavor-catalog-v0.3` reads `87 VERIFIED + 1 PARTIAL`, which projects the post-tag K020 cleanup state. The numbers in this report and in `worklogs/orchestration/wave_008_runbook.md` §8 reflect the at-tag state and remain canonical (see R18-I1 for context).

This block closes `R18-I3` (LOW, docs) — see `.orchestration/ISSUES.md`. Forward-looking convention (Wave-10+): regenerate `audits/factcheck_status.md` at each master-compile bump, optionally using `tools/aggregate_factchecks.py`.
