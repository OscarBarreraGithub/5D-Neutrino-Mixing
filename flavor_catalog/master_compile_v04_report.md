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
