# Flavor Catalog v0.4 Master Compile Report

Date: 2026-05-17
Branch: flavor-catalog/2026q2

## Scope

- Catalog version: flavor-catalog-v0.4
- Process coverage: PRIMARY = 80 + 14 = 94, SECONDARY = 8, total = 102
- Family counts:
  - kaon: PRIMARY 13, SECONDARY 3, total 16
  - charm: PRIMARY 8, SECONDARY 0, total 8
  - beauty: PRIMARY 22, SECONDARY 4, total 26
  - top/Higgs/EW: PRIMARY 19, SECONDARY 1, total 20
  - charged lepton (LFV): PRIMARY 11, SECONDARY 0, total 11
  - EDM/neutrino: PRIMARY 7, SECONDARY 0, total 7
  - collider_rs: PRIMARY 14, SECONDARY 0, total 14
- Wave-8 SECONDARY additions included in the master indexes: K019, K020, K021, B007, B008, B013, B014, T014
- Wave-9 collider_rs PRIMARY additions included in the master indexes: CR001, CR002, CR003, CR004, CR005, CR006, CR007, CR008, CR009, CR010, CR011, CR012, CR013, CR014
- Per-process `.tex` and `.yaml` content was not modified during this master rebuild.
- `catalog_master.tex` was not modified; the existing Collider RS and SECONDARY sections were used as wired.

## Index Rebuild

Each `flavor_catalog/processes/<family>/index.tex`, each `flavor_catalog/processes/secondary/<family>/index.tex`, and `flavor_catalog/processes/collider_rs/index.tex` was regenerated from the sorted `<ID>.tex` files present in that family directory. Each entry uses an `\IfFileExists` guard around the corresponding `\input`.

The guarded input lists matched the checked-in sorted process inventory. The only index-file content change was normalizing `flavor_catalog/processes/collider_rs/index.tex` to the primary master-compile header.

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
- Page count: 179
- PDF size reported by `pdfinfo`: 913614 bytes
- LaTeX notes: build completed with overfull/underfull hbox layout warnings from long monospace paths and dense reference strings, plus hyperref PDF-string warnings for math tokens in section/bookmark strings. Captured two-pass output contained 998 overfull hbox warnings, 68 underfull hbox warnings, and 628 hyperref PDF-string warnings. No process inputs were skipped and no fatal errors occurred.

## Approval and Fact-Check Summary

- Existing PRIMARY approval status (Waves 1-7) unchanged: 80/80 OPUS-APPROVED.
- Existing PRIMARY fact-check status from `flavor_catalog/audits/factcheck_*.md`: 79 VERIFIED, 1 PARTIAL (E009 remains PARTIAL in `factcheck_edm_neutrino.md`; accepted per the v0.2 documentation), 0 MISMATCH, 0 FAILED.
- Wave-8 SECONDARY approval status: 8/8 OPUS-APPROVED per `flavor_catalog/signoff/round_004_index.md`.
- Wave-8 SECONDARY fact-check status: 8 VERIFIED, 0 PARTIAL, 0 MISMATCH, 0 FAILED. K020's earlier PARTIAL status was cleared by the cycle-3 cleanup.
- Wave-9 collider_rs approval status: 14/14 OPUS-APPROVED per `flavor_catalog/signoff/round_005_index.md`.
- Wave-9 collider_rs fact-check status from `flavor_catalog/audits/factcheck_collider_rs.md`: 14 VERIFIED, 0 PARTIAL, 0 MISMATCH, 0 FAILED.
- Total catalog state after Wave-9: 102/102 OPUS-APPROVED, with 101 VERIFIED and 1 PARTIAL across PRIMARY plus SECONDARY.

## Priority Tier Cross-Reference

The PRIMARY/SECONDARY split and implementation-priority rule are documented in `flavor_catalog/PRIORITY_TIERS.md`. Existing Waves 1-7 flavor entries and Wave-9 `collider_rs` entries are treated as PRIMARY for the v0.4 master catalog; Wave-8 additions remain SECONDARY and should not be pulled into code ahead of PRIMARY entries without explicit PI approval.

## Consolidation status (added by C07 cleanup, 2026-05-25)

Reconciliation of the v0.4 totals (`102/102 OPUS-APPROVED`, `101 VERIFIED + 1 PARTIAL`) against the per-family fact-check audits at the v0.4 boundary (tag `flavor-catalog-v0.4`):

| Source | Entries | VERIFIED | PARTIAL | Notes |
|---|---:|---:|---:|---|
| Consolidated table in `audits/factcheck_status.md` (Waves 1-7 PRIMARY base, rows `:32-108`) | 75 | 74 | 1 | DA-4-converged set; PARTIAL row is E009 (INSPIRE JS-only metadata caveat). |
| Wave-7 PRIMARY addenda in `audits/factcheck_top_higgs_ew.md:294-352` + `factcheck_beauty.md` (B012 row) | 5 | 5 | 0 | T003, T004, T008, T012, B012 — all VERIFIED. |
| Wave-8 SECONDARY addenda in `audits/factcheck_kaon.md:254-313`, `factcheck_beauty.md:243-310`, `factcheck_top_higgs_ew.md:391-427` | 8 | 8 | 0 | K019/K020/K021/B007/B008/B013/B014/T014 — all VERIFIED. K020 was PARTIAL at v0.3 and cleared by post-v0.3 cycle-3 cleanup commit `b5c2375` (before the v0.4 tag). |
| Wave-9 PRIMARY rows in consolidated table (`audits/factcheck_status.md`, CR001-CR014) | 14 | 14 | 0 | Added directly to the consolidated table by commit `0c5dacc` (folded in at v0.4 master-compile time). |
| **Total at v0.4 boundary** | **102** | **101** | **1** | Matches headline `101 VERIFIED + 1 PARTIAL` in §"Approval and Fact-Check Summary" above. |

At the v0.4 tag the consolidated `audits/factcheck_status.md` table has 89 rows (75 DA-4 base + 14 Wave-9 collider_rs); the 13 Wave-7 PRIMARY + Wave-8 SECONDARY verdicts remain in their per-family addendum sections rather than in the consolidated table. The Wave-7/Wave-8 verdicts are individually authoritative in the per-family addenda cited above; the headline `101 VERIFIED + 1 PARTIAL` figure is the sum of the consolidated table (88 VERIFIED + 1 PARTIAL) plus the Wave-7 addenda (5 VERIFIED) plus the Wave-8 addenda (8 VERIFIED).

Note on the v0.4 tag annotation: `git cat-file -p flavor-catalog-v0.4` reads `100 VERIFIED + 2 PARTIAL`, which projects the v0.3-tagged historical K020 PARTIAL state. The numbers in this report (101V + 1P) reflect the at-tag state of v0.4 (K020 cleared pre-v0.4) and are canonical (see R19-I3 for context).

This block closes `R19-I4` (LOW, docs) — see `.orchestration/ISSUES.md`. Forward-looking convention (Wave-10+): regenerate `audits/factcheck_status.md` at each master-compile bump, optionally using `tools/aggregate_factchecks.py`.
