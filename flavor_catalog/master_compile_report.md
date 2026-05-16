# Flavor Catalog Master Compile Report

Last update timestamp: 2026-05-16

## Build Inputs

- Branch: `flavor-catalog/2026q2`
- Master file: `catalog_master.tex`
- Output PDF: `catalog_master.pdf`
- Process files indexed: 75
- Families covered: 6

Family input counts generated from `processes/<family>/*.tex`:

- `kaon`: 13
- `charm`: 8
- `beauty`: 21
- `top_higgs_ew`: 15
- `charged_lepton`: 11
- `edm_neutrino`: 7

## Sign-Off Summary

- 4 discovery rounds converged.
- 2 Opus sign-off rounds: round 1 had 50 APPROVE; round 2 had 21 APPROVE.
- 5 individual arbitrations completed: L001, B021, B023, B001, B003.
- Total release state: 75/75 OPUS-APPROVED.

## Build Commands

```bash
cd /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog
rm -f *.aux *.log *.out *.toc *.synctex.gz
pdflatex -interaction=nonstopmode catalog_master.tex
pdflatex -interaction=nonstopmode catalog_master.tex
rm -f *.aux *.log *.out *.toc *.synctex.gz
pdfinfo catalog_master.pdf | head -20
```

## Result

- Build status: succeeded
- Page count: 123
- PDF size: 699241 bytes
- PDF producer: pdfTeX-1.40.19
- Guarded inputs required: none

`pdfinfo catalog_master.pdf | head -20` reported:

```text
Title:
Subject:
Keywords:
Author:
Creator:        LaTeX with hyperref package
Producer:       pdfTeX-1.40.19
CreationDate:   Sat May 16 16:44:21 2026 EDT
ModDate:        Sat May 16 16:44:21 2026 EDT
Tagged:         no
UserProperties: no
Suspects:       no
Form:           none
JavaScript:     no
Pages:          123
Encrypted:      no
Page size:      612 x 792 pts (letter)
Page rot:       0
File size:      699241 bytes
Optimized:      no
PDF version:    1.5
```
