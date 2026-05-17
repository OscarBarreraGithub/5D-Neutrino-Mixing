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
