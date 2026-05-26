# Flavor Catalog v0.2 Master Compile Report

Date: 2026-05-16
Branch: flavor-catalog/2026q2

## Scope

- Catalog version: flavor-catalog-v0.2
- Process coverage: 80 OPUS-APPROVED processes across 6 families
- Family counts:
  - kaon: 13
  - charm: 8
  - beauty: 22
  - top/Higgs/EW: 19
  - charged lepton (LFV): 11
  - EDM/neutrino: 7
- Wave-7 additions included in the master indexes: T003, T004, T008, T012, B012
- Per-process `.tex` content was not modified during this master rebuild.

## Index Rebuild

Each `flavor_catalog/processes/<family>/index.tex` was regenerated from the sorted `<ID>.tex` files present in that family directory. Each entry uses an `\IfFileExists` guard around the corresponding `\input`.

## Build Command

```bash
cd /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog
rm -f *.aux *.log *.out *.toc *.synctex.gz
pdflatex -interaction=nonstopmode catalog_master.tex
pdflatex -interaction=nonstopmode catalog_master.tex
rm -f *.aux *.log *.out *.toc *.synctex.gz
pdfinfo catalog_master.pdf | head -20
```

## Build Result

- Result: success
- Output PDF: `flavor_catalog/catalog_master.pdf`
- Page count: 133
- PDF size reported by `pdfinfo`: 733862 bytes
- LaTeX notes: build completed with layout warnings from long monospace references and hyperref PDF-string warnings, with no process inputs skipped.

## Approval and Fact-Check Summary

- Opus round 1: 50/50 APPROVE
- Opus round 2: 21/21 APPROVE
- Opus round 3: 5/5 APPROVE plus subtleties APPROVE plus addendum APPROVE
- Opus arbitrations: L001, B021, B023, B001, B003
- Total approval status: 80/80 OPUS-APPROVED
- Fact-check status: 79 VERIFIED, 1 PARTIAL (E009 INSPIRE JS-only; APS journal cross-check confirms content), 0 MISMATCH, 0 FAILED

## Consolidation status (added by C07 cleanup, 2026-05-25)

Reconciliation of the headline `79 VERIFIED, 1 PARTIAL` count against the per-family fact-check audits at the v0.2 boundary (tag `flavor-catalog-v0.2`):

| Source | Entries | VERIFIED | PARTIAL | Notes |
|---|---:|---:|---:|---|
| Consolidated table in `audits/factcheck_status.md` (rows `:32-108`) | 75 | 74 | 1 | DA-4-converged set; PARTIAL row is E009. |
| Wave-7 PRIMARY addendum in `audits/factcheck_top_higgs_ew.md:294-352` | 4 | 4 | 0 | T003, T004, T008, T012 — all VERIFIED. |
| Wave-7 PRIMARY addendum in `audits/factcheck_beauty.md` (B012 row) | 1 | 1 | 0 | B012 VERIFIED. |
| **Total at v0.2 boundary** | **80** | **79** | **1** | Matches headline above. |

The v0.2 consolidated `audits/factcheck_status.md` table was inherited unchanged from the v0.1 tag (DA-4 converged at 75) and was not regenerated to fold in the 5 Wave-7 PRIMARY additions. The Wave-7 verdicts are individually authoritative in the per-family addendum sections cited above; the headline `79/80 VERIFIED` figure is the sum of the consolidated table (74 VERIFIED + 1 PARTIAL) and the two addenda (4 + 1 VERIFIED).

This block closes `R17-I2` (LOW, docs) — see `.orchestration/ISSUES.md`. Forward-looking convention (Wave-10+): regenerate `audits/factcheck_status.md` at each master-compile bump, optionally using `tools/aggregate_factchecks.py`.
