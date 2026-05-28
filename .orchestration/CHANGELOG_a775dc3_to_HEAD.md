# Structured changelog: a775dc3..HEAD

Source commands: `git log --oneline a775dc3..HEAD`; per-entry files from `git show --stat <sha>`.
Total commits read: 412.

Phase counts:
- Phase-0 (merge prep): 6 commits
- Phase-1 (24-unit review tracking): 383 commits
- Phase-2 (21-unit cleanup): 23 commits

Phase assignment used here: Phase-1 includes the 382 reviewed payload commits (`6ccf6d8^..5f31f2d`) plus the review-report commit `a05832e`; Phase-0 contains the merge-prep/consolidation commits around the reviewed payload; Phase-2 contains the cleanup commits `3ab1f8f^..71a736b`.

## Phase-0 (merge prep)

### 1. `c9d6921` — chore(merge-prep): commit collaborator exports + orchestration plan
- SHA: `c9d6921ee33e619555199709c791b3619c17b3cc`
- Message: chore(merge-prep): commit collaborator exports + orchestration plan
- Physical/numerical summary: Captured pre-merge decisions, merge plan, planner/reviewer notes, collaborator export scripts, and two small collaborator point CSVs/provenance files.
- Files touched (`git show --stat`):
```text
 .orchestration/MERGE_PLAN.md                                        | 672 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 .orchestration/PLAN_CHANGES_R2.md                                   |  82 ++++++++++++++
 .orchestration/PLAN_REVIEW_R1.md                                    | 178 +++++++++++++++++++++++++++++
 .orchestration/PLAN_REVIEW_R2.md                                    |  43 ++++++++
 .orchestration/PRE_MERGE_DECISIONS.md                               |  14 +++
 artifacts/collaborator_5tev_points.csv                              |   6 +
 artifacts/collaborator_5tev_points.provenance.json                  | 331 ++++++++++++++++++++++++++++++++++++++++++++++++++++++
 artifacts/collaborator_direct_affine_5_10tev_points.csv             |   6 +
 artifacts/collaborator_direct_affine_5_10tev_points.provenance.json | 149 +++++++++++++++++++++++++
 scripts/export_collaborator_5tev_points.py                          | 631 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 scripts/export_collaborator_direct_affine_points.py                 | 878 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 11 files changed, 2990 insertions(+)
```

### 2. `4d33bfc` — merge: trunk consolidation 2026-05-25 (post-review)
- SHA: `4d33bfc10928cf1735ccd987c1dc85c26f0e246d`
- Message: merge: trunk consolidation 2026-05-25 (post-review)
- Physical/numerical summary: Merged the reviewed trunk payload plus orchestration artifacts into main, preserving 24 review reports and collaborator exports.
- Files touched (`git show --stat`):
```text
 .orchestration/ISSUES.md                                            | 484 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 .orchestration/MERGE_PLAN.md                                        | 672 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 .orchestration/PLAN_CHANGES_R2.md                                   |  82 ++++++++++++++
 .orchestration/PLAN_REVIEW_R1.md                                    | 178 +++++++++++++++++++++++++++++
 .orchestration/PLAN_REVIEW_R2.md                                    |  43 ++++++++
 .orchestration/PRE_MERGE_DECISIONS.md                               |  14 +++
 .orchestration/PRE_MERGE_STATE.md                                   |  56 ++++++++++
 .orchestration/REVIEW_QUEUE.md                                      |  30 +++++
 .orchestration/progress.json                                        | 103 +++++++++++++++++
 .orchestration/reviews/R01.md                                       |  57 ++++++++++
 .orchestration/reviews/R02.md                                       |  66 +++++++++++
 .orchestration/reviews/R03.md                                       |  80 ++++++++++++++
 .orchestration/reviews/R04.md                                       | 134 ++++++++++++++++++++++
 .orchestration/reviews/R05.md                                       | 107 ++++++++++++++++++
 .orchestration/reviews/R06.md                                       | 100 +++++++++++++++++
 .orchestration/reviews/R07.md                                       |  64 +++++++++++
 .orchestration/reviews/R08.md                                       |  75 +++++++++++++
 .orchestration/reviews/R09.md                                       |  91 +++++++++++++++
 .orchestration/reviews/R10a.md                                      | 134 ++++++++++++++++++++++
 .orchestration/reviews/R10b.md                                      | 157 ++++++++++++++++++++++++++
 .orchestration/reviews/R10c.md                                      | 243 ++++++++++++++++++++++++++++++++++++++++
 .orchestration/reviews/R11.md                                       |  71 ++++++++++++
 .orchestration/reviews/R12.md                                       | 118 ++++++++++++++++++++
 .orchestration/reviews/R13.md                                       | 115 +++++++++++++++++++
 .orchestration/reviews/R14.md                                       | 106 ++++++++++++++++++
 .orchestration/reviews/R15.md                                       | 115 +++++++++++++++++++
 .orchestration/reviews/R16.md                                       | 125 +++++++++++++++++++++
 .orchestration/reviews/R17.md                                       | 125 +++++++++++++++++++++
 .orchestration/reviews/R18.md                                       |  98 ++++++++++++++++
 .orchestration/reviews/R19.md                                       | 115 +++++++++++++++++++
 .orchestration/reviews/R20.md                                       |  88 +++++++++++++++
 .orchestration/reviews/R21.md                                       |  84 ++++++++++++++
 .orchestration/reviews/R22.md                                       |  83 ++++++++++++++
 artifacts/collaborator_5tev_points.csv                              |   6 +
 artifacts/collaborator_5tev_points.provenance.json                  | 331 ++++++++++++++++++++++++++++++++++++++++++++++++++++++
 artifacts/collaborator_direct_affine_5_10tev_points.csv             |   6 +
 artifacts/collaborator_direct_affine_5_10tev_points.provenance.json | 149 +++++++++++++++++++++++++
 scripts/export_collaborator_5tev_points.py                          | 631 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 scripts/export_collaborator_direct_affine_points.py                 | 878 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 39 files changed, 6214 insertions(+)
```

### 3. `27f4e82` — docs(claude): point active-paper-branch reference at main after consolidation
- SHA: `27f4e826769b42012e3e838477f5523205fa5477`
- Message: docs(claude): point active-paper-branch reference at main after consolidation
- Physical/numerical summary: Repointed Claude/docs branch references from the retired paper branch to main after consolidation.
- Files touched (`git show --stat`):
```text
 CLAUDE.md | 8 +++++---
 1 file changed, 5 insertions(+), 3 deletions(-)
```

### 4. `cb58a36` — merge: fold website branch into main (consolidate to single trunk)
- SHA: `cb58a36cbad932e5db23168b3ac5bce92f4fac7e`
- Message: merge: fold website branch into main (consolidate to single trunk)
- Physical/numerical summary: Folded the website branch into main, adding the full Astro site, 102 citation-anchor datasets, priority data, and generated entry JSON.
- Files touched (`git show --stat`):
```text
 flavor_catalog/website/.gitignore                                 |    8 +
 flavor_catalog/website/.node-version                              |    1 +
 flavor_catalog/website/README.md                                  |  227 ++++++
 flavor_catalog/website/WEBSITE_RUNBOOK.md                         |   77 ++
 flavor_catalog/website/_data/citation_anchors/B001.yaml           |   77 ++
 flavor_catalog/website/_data/citation_anchors/B002.yaml           |  127 ++++
 flavor_catalog/website/_data/citation_anchors/B003.yaml           |  211 ++++++
 flavor_catalog/website/_data/citation_anchors/B004.yaml           |   80 +++
 flavor_catalog/website/_data/citation_anchors/B005.yaml           |  254 +++++++
 flavor_catalog/website/_data/citation_anchors/B006.yaml           |  412 +++++++++++
 flavor_catalog/website/_data/citation_anchors/B007.yaml           |  488 +++++++++++++
 flavor_catalog/website/_data/citation_anchors/B008.yaml           |   93 +++
 flavor_catalog/website/_data/citation_anchors/B009.yaml           |  157 ++++
 flavor_catalog/website/_data/citation_anchors/B011.yaml           |   61 ++
 flavor_catalog/website/_data/citation_anchors/B012.yaml           |  154 ++++
 flavor_catalog/website/_data/citation_anchors/B013.yaml           |  284 ++++++++
 flavor_catalog/website/_data/citation_anchors/B014.yaml           |  338 +++++++++
 flavor_catalog/website/_data/citation_anchors/B015.yaml           |  197 +++++
 flavor_catalog/website/_data/citation_anchors/B016.yaml           |   42 ++
 flavor_catalog/website/_data/citation_anchors/B017.yaml           |  155 ++++
 flavor_catalog/website/_data/citation_anchors/B018.yaml           |   98 +++
 flavor_catalog/website/_data/citation_anchors/B019.yaml           |  163 +++++
 flavor_catalog/website/_data/citation_anchors/B021.yaml           |  176 +++++
 flavor_catalog/website/_data/citation_anchors/B022.yaml           |  180 +++++
 flavor_catalog/website/_data/citation_anchors/B023.yaml           |  176 +++++
 flavor_catalog/website/_data/citation_anchors/B025.yaml           |  203 ++++++
 flavor_catalog/website/_data/citation_anchors/B026.yaml           |  310 ++++++++
 flavor_catalog/website/_data/citation_anchors/B032.yaml           |  445 ++++++++++++
 flavor_catalog/website/_data/citation_anchors/B033.yaml           |  154 ++++
 flavor_catalog/website/_data/citation_anchors/B034.yaml           |  271 +++++++
 flavor_catalog/website/_data/citation_anchors/C001.yaml           |  115 +++
 flavor_catalog/website/_data/citation_anchors/C002.yaml           |   70 ++
 flavor_catalog/website/_data/citation_anchors/C003.yaml           |  119 +++
 flavor_catalog/website/_data/citation_anchors/C004.yaml           |   79 ++
 flavor_catalog/website/_data/citation_anchors/C005.yaml           |   48 ++
 flavor_catalog/website/_data/citation_anchors/C006.yaml           |   36 +
 flavor_catalog/website/_data/citation_anchors/C007.yaml           |  102 +++
 flavor_catalog/website/_data/citation_anchors/C008.yaml           |   60 ++
 flavor_catalog/website/_data/citation_anchors/CR001.yaml          |  135 ++++
 flavor_catalog/website/_data/citation_anchors/CR002.yaml          |  175 +++++
 flavor_catalog/website/_data/citation_anchors/CR003.yaml          |   80 +++
 flavor_catalog/website/_data/citation_anchors/CR004.yaml          |   42 ++
 flavor_catalog/website/_data/citation_anchors/CR005.yaml          |   97 +++
 flavor_catalog/website/_data/citation_anchors/CR006.yaml          |   96 +++
 flavor_catalog/website/_data/citation_anchors/CR007.yaml          |   79 ++
 flavor_catalog/website/_data/citation_anchors/CR008.yaml          |   37 +
 flavor_catalog/website/_data/citation_anchors/CR009.yaml          |  156 ++++
 flavor_catalog/website/_data/citation_anchors/CR010.yaml          |  100 +++
 flavor_catalog/website/_data/citation_anchors/CR011.yaml          |   41 ++
 flavor_catalog/website/_data/citation_anchors/CR012.yaml          |   76 ++
 flavor_catalog/website/_data/citation_anchors/CR013.yaml          |  154 ++++
 flavor_catalog/website/_data/citation_anchors/CR014.yaml          |   71 ++
 flavor_catalog/website/_data/citation_anchors/E001.yaml           |   59 ++
 flavor_catalog/website/_data/citation_anchors/E002.yaml           |  194 +++++
 flavor_catalog/website/_data/citation_anchors/E004.yaml           |   61 ++
 flavor_catalog/website/_data/citation_anchors/E006.yaml           |   82 +++
 flavor_catalog/website/_data/citation_anchors/E007.yaml           |  187 +++++
 flavor_catalog/website/_data/citation_anchors/E008.yaml           |  142 ++++
 flavor_catalog/website/_data/citation_anchors/E009.yaml           |   97 +++
 flavor_catalog/website/_data/citation_anchors/EW001.yaml          |  390 ++++++++++
 flavor_catalog/website/_data/citation_anchors/EW002.yaml          |  262 +++++++
 flavor_catalog/website/_data/citation_anchors/EW003.yaml          |  232 ++++++
 flavor_catalog/website/_data/citation_anchors/K001.yaml           |  145 ++++
 flavor_catalog/website/_data/citation_anchors/K002.yaml           |   79 ++
 flavor_catalog/website/_data/citation_anchors/K003.yaml           |   23 +
 flavor_catalog/website/_data/citation_anchors/K004.yaml           |   91 +++
 flavor_catalog/website/_data/citation_anchors/K005.yaml           |   20 +
 flavor_catalog/website/_data/citation_anchors/K006.yaml           |   23 +
 flavor_catalog/website/_data/citation_anchors/K008.yaml           |  877 +++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/K009.yaml           |  445 ++++++++++++
 flavor_catalog/website/_data/citation_anchors/K010.yaml           |  348 +++++++++
 flavor_catalog/website/_data/citation_anchors/K012.yaml           |  197 +++++
 flavor_catalog/website/_data/citation_anchors/K013.yaml           |   23 +
 flavor_catalog/website/_data/citation_anchors/K017.yaml           |   80 +++
 flavor_catalog/website/_data/citation_anchors/K018.yaml           |  259 +++++++
 flavor_catalog/website/_data/citation_anchors/K019.yaml           |   55 ++
 flavor_catalog/website/_data/citation_anchors/K020.yaml           |   96 +++
 flavor_catalog/website/_data/citation_anchors/K021.yaml           |  149 ++++
 flavor_catalog/website/_data/citation_anchors/L001.yaml           |   48 ++
 flavor_catalog/website/_data/citation_anchors/L002.yaml           |   33 +
 flavor_catalog/website/_data/citation_anchors/L003.yaml           |   18 +
 flavor_catalog/website/_data/citation_anchors/L004.yaml           |   18 +
 flavor_catalog/website/_data/citation_anchors/L005.yaml           |   55 ++
 flavor_catalog/website/_data/citation_anchors/L006.yaml           |   66 ++
 flavor_catalog/website/_data/citation_anchors/L007.yaml           |   55 ++
 flavor_catalog/website/_data/citation_anchors/L008.yaml           |   55 ++
 flavor_catalog/website/_data/citation_anchors/L009.yaml           |   42 ++
 flavor_catalog/website/_data/citation_anchors/L010.yaml           |   51 ++
 flavor_catalog/website/_data/citation_anchors/L023.yaml           |  142 ++++
 flavor_catalog/website/_data/citation_anchors/T001.yaml           |  190 +++++
 flavor_catalog/website/_data/citation_anchors/T002.yaml           |  181 +++++
 flavor_catalog/website/_data/citation_anchors/T003.yaml           |  137 ++++
 flavor_catalog/website/_data/citation_anchors/T004.yaml           |   80 +++
 flavor_catalog/website/_data/citation_anchors/T005.yaml           |  210 ++++++
 flavor_catalog/website/_data/citation_anchors/T006.yaml           |  172 +++++
 flavor_catalog/website/_data/citation_anchors/T007.yaml           |  116 +++
 flavor_catalog/website/_data/citation_anchors/T008.yaml           |  154 ++++
 flavor_catalog/website/_data/citation_anchors/T010.yaml           |   90 +++
 flavor_catalog/website/_data/citation_anchors/T012.yaml           |   61 ++
 flavor_catalog/website/_data/citation_anchors/T014.yaml           |   61 ++
 flavor_catalog/website/_data/citation_anchors/T015.yaml           |  144 ++++
 flavor_catalog/website/_data/citation_anchors/T016.yaml           |  164 +++++
 flavor_catalog/website/_data/citation_anchors/T017.yaml           |  106 +++
 flavor_catalog/website/_data/citation_anchors/T018.yaml           |  221 ++++++
 flavor_catalog/website/_data/citation_anchors/T019.yaml           |  168 +++++
 flavor_catalog/website/_data/citation_anchors/T020.yaml           |  138 ++++
 flavor_catalog/website/_data/priority/B001.yaml                   |    8 +
 flavor_catalog/website/_data/priority/B002.yaml                   |    8 +
 flavor_catalog/website/_data/priority/B003.yaml                   |    8 +
 flavor_catalog/website/_data/priority/B004.yaml                   |    8 +
 flavor_catalog/website/_data/priority/B005.yaml                   |    8 +
 flavor_catalog/website/_data/priority/B006.yaml                   |    8 +
 flavor_catalog/website/_data/priority/B007.yaml                   |    8 +
 flavor_catalog/website/_data/priority/B008.yaml                   |    8 +
 flavor_catalog/website/_data/priority/B009.yaml                   |    8 +
 flavor_catalog/website/_data/priority/B011.yaml                   |    8 +
 flavor_catalog/website/_data/priority/B012.yaml                   |    8 +
 flavor_catalog/website/_data/priority/B013.yaml                   |    8 +
 flavor_catalog/website/_data/priority/B014.yaml                   |    8 +
 flavor_catalog/website/_data/priority/B015.yaml                   |    8 +
 flavor_catalog/website/_data/priority/B016.yaml                   |    8 +
 flavor_catalog/website/_data/priority/B017.yaml                   |    8 +
 flavor_catalog/website/_data/priority/B018.yaml                   |    8 +
 flavor_catalog/website/_data/priority/B019.yaml                   |    8 +
 flavor_catalog/website/_data/priority/B021.yaml                   |    8 +
 flavor_catalog/website/_data/priority/B022.yaml                   |    8 +
 flavor_catalog/website/_data/priority/B023.yaml                   |    8 +
 flavor_catalog/website/_data/priority/B025.yaml                   |    8 +
 flavor_catalog/website/_data/priority/B026.yaml                   |    8 +
 flavor_catalog/website/_data/priority/B032.yaml                   |    8 +
 flavor_catalog/website/_data/priority/B033.yaml                   |    8 +
 flavor_catalog/website/_data/priority/B034.yaml                   |    8 +
 flavor_catalog/website/_data/priority/C001.yaml                   |   12 +
 flavor_catalog/website/_data/priority/C002.yaml                   |   13 +
 flavor_catalog/website/_data/priority/C003.yaml                   |   12 +
 flavor_catalog/website/_data/priority/C004.yaml                   |   12 +
 flavor_catalog/website/_data/priority/C005.yaml                   |   12 +
 flavor_catalog/website/_data/priority/C006.yaml                   |   12 +
 flavor_catalog/website/_data/priority/C007.yaml                   |   11 +
 flavor_catalog/website/_data/priority/C008.yaml                   |   12 +
 flavor_catalog/website/_data/priority/CR001.yaml                  |   12 +
 flavor_catalog/website/_data/priority/CR002.yaml                  |   11 +
 flavor_catalog/website/_data/priority/CR003.yaml                  |   12 +
 flavor_catalog/website/_data/priority/CR004.yaml                  |   12 +
 flavor_catalog/website/_data/priority/CR005.yaml                  |   12 +
 flavor_catalog/website/_data/priority/CR006.yaml                  |   12 +
 flavor_catalog/website/_data/priority/CR007.yaml                  |   12 +
 flavor_catalog/website/_data/priority/CR008.yaml                  |   11 +
 flavor_catalog/website/_data/priority/CR009.yaml                  |   12 +
 flavor_catalog/website/_data/priority/CR010.yaml                  |   12 +
 flavor_catalog/website/_data/priority/CR011.yaml                  |   11 +
 flavor_catalog/website/_data/priority/CR012.yaml                  |   11 +
 flavor_catalog/website/_data/priority/CR013.yaml                  |   11 +
 flavor_catalog/website/_data/priority/CR014.yaml                  |   12 +
 flavor_catalog/website/_data/priority/E001.yaml                   |   12 +
 flavor_catalog/website/_data/priority/E002.yaml                   |   11 +
 flavor_catalog/website/_data/priority/E004.yaml                   |   12 +
 flavor_catalog/website/_data/priority/E006.yaml                   |   13 +
 flavor_catalog/website/_data/priority/E007.yaml                   |   12 +
 flavor_catalog/website/_data/priority/E008.yaml                   |   13 +
 flavor_catalog/website/_data/priority/E009.yaml                   |   12 +
 flavor_catalog/website/_data/priority/EW001.yaml                  |    8 +
 flavor_catalog/website/_data/priority/EW002.yaml                  |    8 +
 flavor_catalog/website/_data/priority/EW003.yaml                  |    8 +
 flavor_catalog/website/_data/priority/K001.yaml                   |   11 +
 flavor_catalog/website/_data/priority/K002.yaml                   |   12 +
 flavor_catalog/website/_data/priority/K003.yaml                   |   12 +
 flavor_catalog/website/_data/priority/K004.yaml                   |   12 +
 flavor_catalog/website/_data/priority/K005.yaml                   |   12 +
 flavor_catalog/website/_data/priority/K006.yaml                   |   12 +
 flavor_catalog/website/_data/priority/K008.yaml                   |   12 +
 flavor_catalog/website/_data/priority/K009.yaml                   |   12 +
 flavor_catalog/website/_data/priority/K010.yaml                   |   11 +
 flavor_catalog/website/_data/priority/K012.yaml                   |   12 +
 flavor_catalog/website/_data/priority/K013.yaml                   |   12 +
 flavor_catalog/website/_data/priority/K017.yaml                   |   12 +
 flavor_catalog/website/_data/priority/K018.yaml                   |   11 +
 flavor_catalog/website/_data/priority/K019.yaml                   |   12 +
 flavor_catalog/website/_data/priority/K020.yaml                   |   12 +
 flavor_catalog/website/_data/priority/K021.yaml                   |   11 +
 flavor_catalog/website/_data/priority/L001.yaml                   |   12 +
 flavor_catalog/website/_data/priority/L002.yaml                   |   12 +
 flavor_catalog/website/_data/priority/L003.yaml                   |   13 +
 flavor_catalog/website/_data/priority/L004.yaml                   |   12 +
 flavor_catalog/website/_data/priority/L005.yaml                   |   12 +
 flavor_catalog/website/_data/priority/L006.yaml                   |   12 +
 flavor_catalog/website/_data/priority/L007.yaml                   |   11 +
 flavor_catalog/website/_data/priority/L008.yaml                   |   11 +
 flavor_catalog/website/_data/priority/L009.yaml                   |   11 +
 flavor_catalog/website/_data/priority/L010.yaml                   |   11 +
 flavor_catalog/website/_data/priority/L023.yaml                   |   12 +
 flavor_catalog/website/_data/priority/T001.yaml                   |    8 +
 flavor_catalog/website/_data/priority/T002.yaml                   |    8 +
 flavor_catalog/website/_data/priority/T003.yaml                   |    8 +
 flavor_catalog/website/_data/priority/T004.yaml                   |    8 +
 flavor_catalog/website/_data/priority/T005.yaml                   |    8 +
 flavor_catalog/website/_data/priority/T006.yaml                   |    8 +
 flavor_catalog/website/_data/priority/T007.yaml                   |    8 +
 flavor_catalog/website/_data/priority/T008.yaml                   |    8 +
 flavor_catalog/website/_data/priority/T010.yaml                   |    8 +
 flavor_catalog/website/_data/priority/T012.yaml                   |    8 +
 flavor_catalog/website/_data/priority/T014.yaml                   |    8 +
 flavor_catalog/website/_data/priority/T015.yaml                   |    8 +
 flavor_catalog/website/_data/priority/T016.yaml                   |    8 +
 flavor_catalog/website/_data/priority/T017.yaml                   |    8 +
 flavor_catalog/website/_data/priority/T018.yaml                   |    8 +
 flavor_catalog/website/_data/priority/T019.yaml                   |    8 +
 flavor_catalog/website/_data/priority/T020.yaml                   |    8 +
 flavor_catalog/website/astro.config.mjs                           |   20 +
 flavor_catalog/website/cloudflare-pages.config.md                 |  111 +++
 flavor_catalog/website/package-lock.json                          | 5719 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/website/package.json                               |   22 +
 flavor_catalog/website/public/_headers                            |   19 +
 flavor_catalog/website/public/_redirects                          |   11 +
 flavor_catalog/website/scripts/ingest_catalog.py                  |  710 ++++++++++++++++++
 flavor_catalog/website/scripts/resolve_beauty_citation_anchors.py |  917 ++++++++++++++++++++++++
 flavor_catalog/website/scripts/resolve_citation_anchors.py        |  790 ++++++++++++++++++++
 flavor_catalog/website/scripts/screenshot_cdp.mjs                 |   99 +++
 flavor_catalog/website/src/components/Badges.astro                |   75 ++
 flavor_catalog/website/src/components/CitationModal.astro         |  211 ++++++
 flavor_catalog/website/src/components/EntryTable.astro            |  145 ++++
 flavor_catalog/website/src/content.config.ts                      |  141 ++++
 flavor_catalog/website/src/content/catalog_index.json             | 1867 ++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/website/src/content/entries/B001.json              |  314 ++++++++
 flavor_catalog/website/src/content/entries/B002.json              |  356 +++++++++
 flavor_catalog/website/src/content/entries/B003.json              |  527 ++++++++++++++
 flavor_catalog/website/src/content/entries/B004.json              |  270 +++++++
 flavor_catalog/website/src/content/entries/B005.json              |  364 ++++++++++
 flavor_catalog/website/src/content/entries/B006.json              |  475 ++++++++++++
 flavor_catalog/website/src/content/entries/B007.json              |  662 +++++++++++++++++
 flavor_catalog/website/src/content/entries/B008.json              |  339 +++++++++
 flavor_catalog/website/src/content/entries/B009.json              |  396 ++++++++++
 flavor_catalog/website/src/content/entries/B011.json              |  267 +++++++
 flavor_catalog/website/src/content/entries/B012.json              |  407 +++++++++++
 flavor_catalog/website/src/content/entries/B013.json              |  658 +++++++++++++++++
 flavor_catalog/website/src/content/entries/B014.json              |  798 +++++++++++++++++++++
 flavor_catalog/website/src/content/entries/B015.json              |  350 +++++++++
 flavor_catalog/website/src/content/entries/B016.json              |  196 +++++
 flavor_catalog/website/src/content/entries/B017.json              |  305 ++++++++
 flavor_catalog/website/src/content/entries/B018.json              |  246 +++++++
 flavor_catalog/website/src/content/entries/B019.json              |  286 ++++++++
 flavor_catalog/website/src/content/entries/B021.json              |  449 ++++++++++++
 flavor_catalog/website/src/content/entries/B022.json              |  306 ++++++++
 flavor_catalog/website/src/content/entries/B023.json              |  322 +++++++++
 flavor_catalog/website/src/content/entries/B025.json              |  384 ++++++++++
 flavor_catalog/website/src/content/entries/B026.json              |  460 ++++++++++++
 flavor_catalog/website/src/content/entries/B032.json              |  511 +++++++++++++
 flavor_catalog/website/src/content/entries/B033.json              |  272 +++++++
 flavor_catalog/website/src/content/entries/B034.json              |  362 ++++++++++
 flavor_catalog/website/src/content/entries/C001.json              |  436 +++++++++++
 flavor_catalog/website/src/content/entries/C002.json              |  332 +++++++++
 flavor_catalog/website/src/content/entries/C003.json              |  439 ++++++++++++
 flavor_catalog/website/src/content/entries/C004.json              |  322 +++++++++
 flavor_catalog/website/src/content/entries/C005.json              |  264 +++++++
 flavor_catalog/website/src/content/entries/C006.json              |  217 ++++++
 flavor_catalog/website/src/content/entries/C007.json              |  362 ++++++++++
 flavor_catalog/website/src/content/entries/C008.json              |  225 ++++++
 flavor_catalog/website/src/content/entries/CR001.json             |  450 ++++++++++++
 flavor_catalog/website/src/content/entries/CR002.json             |  658 +++++++++++++++++
 flavor_catalog/website/src/content/entries/CR003.json             |  472 ++++++++++++
 flavor_catalog/website/src/content/entries/CR004.json             |  316 ++++++++
 flavor_catalog/website/src/content/entries/CR005.json             |  494 +++++++++++++
 flavor_catalog/website/src/content/entries/CR006.json             |  509 +++++++++++++
 flavor_catalog/website/src/content/entries/CR007.json             |  380 ++++++++++
 flavor_catalog/website/src/content/entries/CR008.json             |  317 ++++++++
 flavor_catalog/website/src/content/entries/CR009.json             |  536 ++++++++++++++
 flavor_catalog/website/src/content/entries/CR010.json             |  451 ++++++++++++
 flavor_catalog/website/src/content/entries/CR011.json             |  393 ++++++++++
 flavor_catalog/website/src/content/entries/CR012.json             |  458 ++++++++++++
 flavor_catalog/website/src/content/entries/CR013.json             |  509 +++++++++++++
 flavor_catalog/website/src/content/entries/CR014.json             |  542 ++++++++++++++
 flavor_catalog/website/src/content/entries/E001.json              |  291 ++++++++
 flavor_catalog/website/src/content/entries/E002.json              |  491 +++++++++++++
 flavor_catalog/website/src/content/entries/E004.json              |  266 +++++++
 flavor_catalog/website/src/content/entries/E006.json              |  325 +++++++++
 flavor_catalog/website/src/content/entries/E007.json              |  357 +++++++++
 flavor_catalog/website/src/content/entries/E008.json              |  297 ++++++++
 flavor_catalog/website/src/content/entries/E009.json              |  442 ++++++++++++
 flavor_catalog/website/src/content/entries/EW001.json             |  527 ++++++++++++++
 flavor_catalog/website/src/content/entries/EW002.json             |  454 ++++++++++++
 flavor_catalog/website/src/content/entries/EW003.json             |  405 +++++++++++
 flavor_catalog/website/src/content/entries/K001.json              |  318 +++++++++
 flavor_catalog/website/src/content/entries/K002.json              |  321 +++++++++
 flavor_catalog/website/src/content/entries/K003.json              |  181 +++++
 flavor_catalog/website/src/content/entries/K004.json              |  307 ++++++++
 flavor_catalog/website/src/content/entries/K005.json              |  191 +++++
 flavor_catalog/website/src/content/entries/K006.json              |  163 +++++
 flavor_catalog/website/src/content/entries/K008.json              | 1097 ++++++++++++++++++++++++++++
 flavor_catalog/website/src/content/entries/K009.json              |  757 ++++++++++++++++++++
 flavor_catalog/website/src/content/entries/K010.json              |  468 ++++++++++++
 flavor_catalog/website/src/content/entries/K012.json              |  315 ++++++++
 flavor_catalog/website/src/content/entries/K013.json              |  166 +++++
 flavor_catalog/website/src/content/entries/K017.json              |  322 +++++++++
 flavor_catalog/website/src/content/entries/K018.json              |  457 ++++++++++++
 flavor_catalog/website/src/content/entries/K019.json              |  326 +++++++++
 flavor_catalog/website/src/content/entries/K020.json              |  383 ++++++++++
 flavor_catalog/website/src/content/entries/K021.json              |  441 ++++++++++++
 flavor_catalog/website/src/content/entries/L001.json              |  336 +++++++++
 flavor_catalog/website/src/content/entries/L002.json              |  240 +++++++
 flavor_catalog/website/src/content/entries/L003.json              |  187 +++++
 flavor_catalog/website/src/content/entries/L004.json              |  186 +++++
 flavor_catalog/website/src/content/entries/L005.json              |  328 +++++++++
 flavor_catalog/website/src/content/entries/L006.json              |  270 +++++++
 flavor_catalog/website/src/content/entries/L007.json              |  299 ++++++++
 flavor_catalog/website/src/content/entries/L008.json              |  294 ++++++++
 flavor_catalog/website/src/content/entries/L009.json              |  256 +++++++
 flavor_catalog/website/src/content/entries/L010.json              |  265 +++++++
 flavor_catalog/website/src/content/entries/L023.json              |  425 +++++++++++
 flavor_catalog/website/src/content/entries/T001.json              |  519 ++++++++++++++
 flavor_catalog/website/src/content/entries/T002.json              |  514 +++++++++++++
 flavor_catalog/website/src/content/entries/T003.json              |  468 ++++++++++++
 flavor_catalog/website/src/content/entries/T004.json              |  325 +++++++++
 flavor_catalog/website/src/content/entries/T005.json              |  579 +++++++++++++++
 flavor_catalog/website/src/content/entries/T006.json              |  516 +++++++++++++
 flavor_catalog/website/src/content/entries/T007.json              |  382 ++++++++++
 flavor_catalog/website/src/content/entries/T008.json              |  451 ++++++++++++
 flavor_catalog/website/src/content/entries/T010.json              |  262 +++++++
 flavor_catalog/website/src/content/entries/T012.json              |  231 ++++++
 flavor_catalog/website/src/content/entries/T014.json              |  290 ++++++++
 flavor_catalog/website/src/content/entries/T015.json              |  417 +++++++++++
 flavor_catalog/website/src/content/entries/T016.json              |  451 ++++++++++++
 flavor_catalog/website/src/content/entries/T017.json              |  343 +++++++++
 flavor_catalog/website/src/content/entries/T018.json              |  609 ++++++++++++++++
 flavor_catalog/website/src/content/entries/T019.json              |  484 +++++++++++++
 flavor_catalog/website/src/content/entries/T020.json              |  444 ++++++++++++
 flavor_catalog/website/src/content/families.json                  |   50 ++
 flavor_catalog/website/src/layouts/BaseLayout.astro               |  133 ++++
 flavor_catalog/website/src/lib/notation.ts                        |  381 ++++++++++
 flavor_catalog/website/src/lib/prose.ts                           |  504 +++++++++++++
 flavor_catalog/website/src/pages/browse.astro                     |  262 +++++++
 flavor_catalog/website/src/pages/entries/[id].astro               |  356 +++++++++
 flavor_catalog/website/src/pages/families/[family].astro          |  140 ++++
 flavor_catalog/website/src/pages/index.astro                      |   53 ++
 flavor_catalog/website/src/pages/methodology.astro                |  263 +++++++
 flavor_catalog/website/src/styles/global.css                      | 1216 +++++++++++++++++++++++++++++++
 flavor_catalog/website/tsconfig.json                              |    5 +
 336 files changed, 70231 insertions(+)
```

### 5. `a809cc3` — chore(deploy): trigger Cloudflare Pages deploy on main
- SHA: `a809cc3510cdaac5a64c24d0470b476e03b06eb9`
- Message: chore(deploy): trigger Cloudflare Pages deploy on main
- Physical/numerical summary: Touched deployment metadata to trigger a Cloudflare Pages deploy from the consolidated main branch.
- Files touched (`git show --stat`):
  - `git show --stat` reported no file-level changes.

### 6. `cd20f3a` — docs: single-trunk state; close INFRA-1 (Cloudflare deploys from main)
- SHA: `cd20f3a233bd7ece82e9f2377a3143e7c3792cbd`
- Message: docs: single-trunk state; close INFRA-1 (Cloudflare deploys from main)
- Physical/numerical summary: Documented the single-trunk state and closed the Cloudflare-main deployment infrastructure issue.
- Files touched (`git show --stat`):
```text
 .orchestration/ISSUES.md | 2 +-
 CLAUDE.md                | 5 ++---
 2 files changed, 3 insertions(+), 4 deletions(-)
```

## Phase-1 (24-unit review tracking)

### 1. `6ccf6d8` — chore(gitignore): track paper docs + phase logs, ignore .claude/ and snapshots/
- SHA: `6ccf6d87d0c9dd09b7f18670fda5ad99c38188f0`
- Message: chore(gitignore): track paper docs + phase logs, ignore .claude/ and snapshots/
- Physical/numerical summary: Started tracking paper docs/phase logs and ignored local agent/snapshot outputs; no physical numerics changed.
- Files touched (`git show --stat`):
```text
 .gitignore | 13 +++++++++++++
 1 file changed, 13 insertions(+)
```

### 2. `109ec02` — feat(qcd): add PDG 2024 MS-bar quark mass running
- SHA: `109ec02a72f83a594038c4330a7befc24615ab7d`
- Message: feat(qcd): add PDG 2024 MS-bar quark mass running
- Physical/numerical summary: Added PDG-2024 MS-bar quark-mass running, establishing the mass-evolution numerics used by later quark-scan gates.
- Files touched (`git show --stat`):
```text
 qcd/__init__.py                      |   5 +++-
 qcd/decoupling.py                    | 112 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-----
 qcd/mass_running.py                  | 248 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 quarkConstraints/pdg_quark_masses.py | 225 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 tests/test_mass_running.py           |  83 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 tests/test_pdg_quark_masses.py       |  82 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 6 files changed, 747 insertions(+), 8 deletions(-)
```

### 3. `d4f873c` — feat(quark): integrate PDG targets into scan gates
- SHA: `d4f873c375c078315c991bad08bab7fa79ffa253`
- Message: feat(quark): integrate PDG targets into scan gates
- Physical/numerical summary: Integrated PDG target masses/CKM constraints into scan gates, making physical acceptance depend on updated quark inputs.
- Files touched (`git show --stat`):
```text
 quarkConstraints/benchmarks.py        |  94 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++--------------
 quarkConstraints/diagnostics.py       | 149 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 quarkConstraints/modern/inputs.py     |  62 ++++++++++++++++++++++++++++++++++++++++++++++++------------
 quarkConstraints/scales.py            |  20 +++++++++++++++++++-
 quarkConstraints/scan.py              | 160 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++---
 scripts/calibrate_phase0.py           | 181 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 tests/test_diagnostics.py             | 113 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 tests/test_modern_input_registry.py   |  33 ++++++++++++++++++++++++++++----
 tests/test_quark_benchmarks.py        |  49 ++++++++++++++++++++++++++++++++++++++++-------
 tests/test_quark_fit.py               |  53 +++++++++++++++++++++++++++++++++++++++++++++++----
 tests/test_quark_scan.py              |   3 +++
 tests/test_quark_target_regression.py |  75 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 12 files changed, 947 insertions(+), 45 deletions(-)
```

### 4. `d0a9103` — feat(scan): add RS-anarchy driver and dispatch scripts
- SHA: `d0a91039091233f8cc3c83ea9ead54ba5d9e0031`
- Message: feat(scan): add RS-anarchy driver and dispatch scripts
- Physical/numerical summary: Added the RS-anarchy scan driver and dispatch scripts, creating the numerical ensemble-generation machinery.
- Files touched (`git show --stat`):
```text
 scripts/run_rs_anarchy.py          | 860 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 scripts/run_rs_anarchy.sbatch      |  67 ++++++++++++++
 scripts/run_rs_anarchy_run3.sbatch |  85 ++++++++++++++++++
 scripts/run_rs_anarchy_runA.sbatch |  52 +++++++++++
 scripts/run_rs_anarchy_runB.sbatch |  75 ++++++++++++++++
 scripts/run_rs_anarchy_runC.sbatch |  60 +++++++++++++
 tests/test_rs_anarchy_priors.py    | 106 ++++++++++++++++++++++
 7 files changed, 1305 insertions(+)
```

### 5. `bcca1df` — feat(analysis): add quark-scan follow-up utilities
- SHA: `bcca1df9c4815b64ab495f789efbbcc97935e13c`
- Message: feat(analysis): add quark-scan follow-up utilities
- Physical/numerical summary: Added follow-up analysis utilities for scan outputs, enabling derived numerical summaries and plots.
- Files touched (`git show --stat`):
```text
 scripts/_audit_followup_crossings.py                | 114 +++++++++++++++++++++++++++++++
 scripts/compute_per_quark_residuals.py              | 106 +++++++++++++++++++++++++++++
 scripts/export_accepted_quark_scan.py               | 185 ++++++++++++++++++++++++++++++++++++++++++++++++++
 scripts/export_accepted_quark_scan_with_yukawas.py  | 600 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 scripts/plot_publication_figures.py                 |  33 +++++++--
 scripts/plot_rs_anarchy_summary.py                  | 161 +++++++++++++++++++++++++++++++++++++++++++
 scripts/rs_anarchy_cfw_comparison.py                | 153 +++++++++++++++++++++++++++++++++++++++++
 scripts/rs_anarchy_gate_sensitivity.py              | 176 +++++++++++++++++++++++++++++++++++++++++++++++
 scripts/rs_anarchy_mkk_min_hist.py                  | 218 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 scripts/rs_anarchy_mkk_min_hist_by_cvals.py         | 184 +++++++++++++++++++++++++++++++++++++++++++++++++
 scripts/rs_anarchy_mkk_min_hist_by_pdg_tightness.py | 184 +++++++++++++++++++++++++++++++++++++++++++++++++
 scripts/rs_anarchy_mkk_min_hist_by_yprior.py        | 167 +++++++++++++++++++++++++++++++++++++++++++++
 scripts/rs_anarchy_mkk_min_hist_by_yprior_gsstar.py | 160 +++++++++++++++++++++++++++++++++++++++++++
 scripts/rs_anarchy_mkk_min_hist_gsstar.py           | 160 +++++++++++++++++++++++++++++++++++++++++++
 scripts/run_followups_README.md                     | 203 ++++++++++++++++++++++++++++++++++++++++++++++++++++++
 scripts/yukawa_envelope_vs_anarchic.py              | 132 +++++++++++++++++++++++++++++++++++
 scripts/yukawa_per_element_anatomy.py               | 166 ++++++++++++++++++++++++++++++++++++++++++++
 17 files changed, 3095 insertions(+), 7 deletions(-)
```

### 6. `d3bcac9` — docs(paper): add quark-scan methodology and notebooks
- SHA: `d3bcac9c54bf68d186f8ca755b75eab904261acc`
- Message: docs(paper): add quark-scan methodology and notebooks
- Physical/numerical summary: Added the quark-scan methodology note and notebooks documenting the scan physics and reproducible analysis path.
- Files touched (`git show --stat`):
```text
 docs/paper_execution_decisions.md                      |  108 ++++++++++++
 docs/paper_execution_roadmap.md                        |  214 +++++++++++++++++++++++
 docs/phase_logs/phase1_impl.md                         |   48 ++++++
 docs/quark_scan_assumptions_compact.pdf                |  Bin 212890 -> 204920 bytes
 docs/quark_scan_assumptions_compact.tex                |   53 +++++-
 docs/quark_scan_methodology_note.pdf                   |  Bin 0 -> 551180 bytes
 docs/quark_scan_methodology_note.tex                   |  904 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 notebooks/dense_scan_2sigma_vs_1sigma_comparison.ipynb |  662 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 notebooks/dense_scan_mkk_constraints_pdg2024.ipynb     |  323 ++++++++++++++++++++++++++++++++++
 notebooks/pdg_quark_target_fix_verification.ipynb      | 1477 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 notebooks/rs_anarchy_analysis.ipynb                    |  675 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 11 files changed, 4458 insertions(+), 6 deletions(-)
```

### 7. `b0675cd` — docs(paper): seal Phase 1 logs (impl + review + signoff)
- SHA: `b0675cdae93ccc3a8341979465e125848d35c9cf`
- Message: docs(paper): seal Phase 1 logs (impl + review + signoff)
- Physical/numerical summary: Recorded reviewer/sign-off disposition; review ledger only.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/phase1_review.md  | 25 +++++++++++++++++++++++++
 docs/phase_logs/phase1_signoff.md | 28 ++++++++++++++++++++++++++++
 2 files changed, 53 insertions(+)
```

### 8. `16e6b36` — test(fixtures): capture pre-fix xfail behavior for spurion-seed tests
- SHA: `16e6b36ebe9bf518cbbb56b3fe50e1c4133aeee5`
- Message: test(fixtures): capture pre-fix xfail behavior for spurion-seed tests
- Physical/numerical summary: Captured the pre-fix xfail state for spurion-seed tests, preserving the failing numerical baseline.
- Files touched (`git show --stat`):
```text
 tests/baselines/pre-fix-xfail-output.txt | 135 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 135 insertions(+)
```

### 9. `2871fce` — fix(benchmarks): re-derive spurion seed under PDG targets
- SHA: `2871fce75a9a9411fdc0a74cdcd49785571b8eaf`
- Message: fix(benchmarks): re-derive spurion seed under PDG targets
- Physical/numerical summary: Re-derived the benchmark spurion seed under PDG targets, changing the pinned benchmark flavor numerics.
- Files touched (`git show --stat`):
```text
 quarkConstraints/benchmarks.py | 36 ++++++++++++++++++++++++++++++------
 tests/test_quark_fit.py        | 20 --------------------
 2 files changed, 30 insertions(+), 26 deletions(-)
```

### 10. `559b851` — chore(pytest): enable strict xfail policy
- SHA: `559b851d147037d0600253b7c3324dadd7393c11`
- Message: chore(pytest): enable strict xfail policy
- Physical/numerical summary: Enabled strict xfail handling so repaired numerical tests cannot silently keep passing as expected failures.
- Files touched (`git show --stat`):
```text
 pyproject.toml | 3 +++
 1 file changed, 3 insertions(+)
```

### 11. `fd83d96` — test(fixtures): capture post-fix passing test output
- SHA: `fd83d961bf63f4e07b3b56036ddfbb57c1067d7f`
- Message: test(fixtures): capture post-fix passing test output
- Physical/numerical summary: Captured post-fix passing fixture output for the benchmark repair, locking the corrected numerical behavior.
- Files touched (`git show --stat`):
```text
 tests/baselines/post-fix-test-output.txt | 10 ++++++++++
 1 file changed, 10 insertions(+)
```

### 12. `7e679ea` — docs(paper): seal Phase 2 hole #4 (xfail repair) impl log
- SHA: `7e679ea29d1640d2c7f952bbb0cfcfbe3b78a67f`
- Message: docs(paper): seal Phase 2 hole #4 (xfail repair) impl log
- Physical/numerical summary: Recorded review/signoff state for the paper or phase logs; bookkeeping over already-reviewed numerical work.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/phase2_h4_impl.md | 77 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 77 insertions(+)
```

### 13. `1fec2c7` — docs(paper): seal Phase 2 hole #4 review + signoff
- SHA: `1fec2c76530560fec88e7df661be5bc271dc7a83`
- Message: docs(paper): seal Phase 2 hole #4 review + signoff
- Physical/numerical summary: Recorded reviewer/sign-off disposition; review ledger only.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/phase2_h4_review.md  | 24 ++++++++++++++++++++++++
 docs/phase_logs/phase2_h4_signoff.md | 19 +++++++++++++++++++
 2 files changed, 43 insertions(+)
```

### 14. `dc9c498` — physics(deltaf2): update kaon inputs to FLAG 2024 and BGS 2020
- SHA: `dc9c498f82444c2866ec512a9931333c16fb4bb8`
- Message: physics(deltaf2): update kaon inputs to FLAG 2024 and BGS 2020
- Physical/numerical summary: Updated kaon DeltaF=2 inputs to FLAG 2024 and BGS 2020 values, changing hadronic constants used in epsilon_K/Delta m_K numerics.
- Files touched (`git show --stat`):
```text
 quarkConstraints/deltaf2.py | 8 ++++----
 tests/test_quark_deltaf2.py | 8 ++++++++
 2 files changed, 12 insertions(+), 4 deletions(-)
```

### 15. `82a96f0` — docs(paper): document hadronic input provenance
- SHA: `82a96f00cccb907dbc88eaf8ce4a2c4da43bf8e3`
- Message: docs(paper): document hadronic input provenance
- Physical/numerical summary: Documented hadronic-input provenance so the updated bag constants and SM inputs have explicit literature anchors.
- Files touched (`git show --stat`):
```text
 docs/audits/bag_param_inventory.md   |  86 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 docs/quark_scan_methodology_note.pdf | Bin 551180 -> 574719 bytes
 docs/quark_scan_methodology_note.tex |  76 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++--
 3 files changed, 160 insertions(+), 2 deletions(-)
```

### 16. `695f35e` — docs(paper): seal Phase 2 hole #5 (bag-param audit) impl log
- SHA: `695f35ea238c6c3dea7c173e554bdcf03fee654e`
- Message: docs(paper): seal Phase 2 hole #5 (bag-param audit) impl log
- Physical/numerical summary: Recorded review/signoff state for the paper or phase logs; bookkeeping over already-reviewed numerical work.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/phase2_h5_impl.md | 73 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 73 insertions(+)
```

### 17. `2521f6b` — docs(paper): seal Phase 2 hole #5 review + signoff + audit decisions
- SHA: `2521f6be7aff8cd749a0de37f78440cb9462c0cd`
- Message: docs(paper): seal Phase 2 hole #5 review + signoff + audit decisions
- Physical/numerical summary: Recorded reviewer/sign-off disposition; review ledger only.
- Files touched (`git show --stat`):
```text
 .gitignore                           |   2 ++
 docs/audits/epsilon_k_sm_decision.md | 148 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 docs/phase_logs/phase2_h5_review.md  | 137 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 docs/phase_logs/phase2_h5_signoff.md | 175 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 4 files changed, 462 insertions(+)
```

### 18. `7f71908` — physics(deltaf2): fix BMU LO Wilson RG running
- SHA: `7f71908a26d7015be20741b798552103e900be7d`
- Message: physics(deltaf2): fix BMU LO Wilson RG running
- Physical/numerical summary: Fixed BMU LO Wilson RG running, altering DeltaF=2 Wilson evolution numerics.
- Files touched (`git show --stat`):
```text
 docs/audits/wilson_rg_inventory.md        | 136 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 docs/audits/wilson_rg_reference_values.md | 118 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 docs/quark_scan_methodology_note.pdf      | Bin 574719 -> 580429 bytes
 docs/quark_scan_methodology_note.tex      |  63 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 quarkConstraints/deltaf2.py               |   9 ++++++---
 quarkConstraints/qcd_running.py           |  72 +++++++++++++++++++++++++++++++++++++++++++++++-------------------------
 scripts/audit_wilson_rg.py                | 167 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 tests/test_qcd_running.py                 |  97 ++++++++++++++++++++++++++++++++++++++++++++++---------------------------------------------------
 tests/test_quark_deltaf2.py               |  19 +++++++++++--------
 tests/test_wilson_rg_audit.py             |  66 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 10 files changed, 660 insertions(+), 87 deletions(-)
```

### 19. `561d848` — docs(paper): seal Phase 2 hole #6 (Wilson-RG audit) impl log
- SHA: `561d848479ab058df9d384fe611e664721325ceb`
- Message: docs(paper): seal Phase 2 hole #6 (Wilson-RG audit) impl log
- Physical/numerical summary: Recorded review/signoff state for the paper or phase logs; bookkeeping over already-reviewed numerical work.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/phase2_h6_impl.md | 121 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 121 insertions(+)
```

### 20. `c80a8e0` — fix(deltaf2): correct BMU-to-scalar-LR sign convention in ADM and matrix elements
- SHA: `c80a8e0a15acdddd984bbd19806eaf2c42161917`
- Message: fix(deltaf2): correct BMU-to-scalar-LR sign convention in ADM and matrix elements
- Physical/numerical summary: Corrected the BMU-to-scalar-LR sign convention in ADM/matrix elements, changing DeltaF=2 amplitude signs consistently.
- Files touched (`git show --stat`):
```text
 docs/audits/wilson_rg_inventory.md        |  84 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++----------------
 docs/audits/wilson_rg_reference_values.md |  15 ++++++++-------
 docs/phase_logs/phase2_h6_impl.md         |  80 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++---------------------
 docs/quark_scan_methodology_note.pdf      | Bin 580429 -> 580512 bytes
 docs/quark_scan_methodology_note.tex      |  17 +++++++++--------
 quarkConstraints/deltaf2.py               |  18 +++++++++++++++---
 quarkConstraints/qcd_running.py           |  16 ++++++++--------
 scripts/audit_wilson_rg.py                |   8 ++++----
 tests/test_epsilon_k_physics.py           |   9 +++++----
 tests/test_modern_scan.py                 |  16 +++++++++-------
 tests/test_qcd_running.py                 |   7 ++++---
 tests/test_quark_deltaf2.py               |   4 +++-
 tests/test_quark_plot_data.py             |  15 ++++++++-------
 tests/test_wilson_rg_audit.py             |   2 +-
 14 files changed, 201 insertions(+), 90 deletions(-)
```

### 21. `c83c16d` — docs(paper): seal Phase 2 hole #6 review + revision + signoff
- SHA: `c83c16db70160fb72c37edfeed4eac21afab95cc`
- Message: docs(paper): seal Phase 2 hole #6 review + revision + signoff
- Physical/numerical summary: Recorded reviewer/sign-off disposition; review ledger only.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/phase2_h6_review.md    |  62 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 docs/phase_logs/phase2_h6_review_v2.md |  29 +++++++++++++++++++++++++++++
 docs/phase_logs/phase2_h6_signoff.md   | 142 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 3 files changed, 233 insertions(+)
```

### 22. `29803ff` — scripts(rerun): snapshot pre-audit figures to quark_pre_audit_constants/
- SHA: `29803fffa70ede4ffab736b4bebc705eeceef0f6`
- Message: scripts(rerun): snapshot pre-audit figures to quark_pre_audit_constants/
- Physical/numerical summary: Added or adjusted an analysis/rerun script supporting the documented numerical workflow.
- Files touched (`git show --stat`):
```text
 results/figures/quark_pre_audit_constants/c_values_2sigma_vs_1sigma.pdf                         | Bin 0 -> 30474 bytes
 results/figures/quark_pre_audit_constants/c_values_2sigma_vs_1sigma.png                         | Bin 0 -> 194674 bytes
 results/figures/quark_pre_audit_constants/fig1_exclusion_boundaries.png                         | Bin 0 -> 145529 bytes
 results/figures/quark_pre_audit_constants/fig1_exclusion_boundaries_1sigma.pdf                  | Bin 0 -> 34232 bytes
 results/figures/quark_pre_audit_constants/fig1_exclusion_boundaries_1sigma.png                  | Bin 0 -> 90140 bytes
 results/figures/quark_pre_audit_constants/fig1_exclusion_boundaries_2sigma.pdf                  | Bin 0 -> 34232 bytes
 results/figures/quark_pre_audit_constants/fig1_exclusion_boundaries_2sigma.png                  | Bin 0 -> 90140 bytes
 results/figures/quark_pre_audit_constants/fig1_exclusion_boundaries_pdg2024.pdf                 | Bin 0 -> 34232 bytes
 results/figures/quark_pre_audit_constants/fig1_exclusion_boundaries_pdg2024.png                 | Bin 0 -> 90140 bytes
 results/figures/quark_pre_audit_constants/fig2_mkk_bound_1sigma.pdf                             | Bin 0 -> 20850 bytes
 results/figures/quark_pre_audit_constants/fig2_mkk_bound_1sigma.png                             | Bin 0 -> 87971 bytes
 results/figures/quark_pre_audit_constants/fig2_mkk_bound_2007_vs_modern.png                     | Bin 0 -> 130167 bytes
 results/figures/quark_pre_audit_constants/fig2_mkk_bound_2007_vs_modern_pdg2024.pdf             | Bin 0 -> 20850 bytes
 results/figures/quark_pre_audit_constants/fig2_mkk_bound_2007_vs_modern_pdg2024.png             | Bin 0 -> 87971 bytes
 results/figures/quark_pre_audit_constants/fig2_mkk_bound_2sigma.pdf                             | Bin 0 -> 20850 bytes
 results/figures/quark_pre_audit_constants/fig2_mkk_bound_2sigma.png                             | Bin 0 -> 87971 bytes
 results/figures/quark_pre_audit_constants/rs_anarchy_cfw_comparison.pdf                         | Bin 0 -> 36074 bytes
 results/figures/quark_pre_audit_constants/rs_anarchy_cfw_comparison.png                         | Bin 0 -> 167675 bytes
 results/figures/quark_pre_audit_constants/rs_anarchy_gate_sensitivity.pdf                       | Bin 0 -> 28152 bytes
 results/figures/quark_pre_audit_constants/rs_anarchy_gate_sensitivity.png                       | Bin 0 -> 204557 bytes
 results/figures/quark_pre_audit_constants/rs_anarchy_max_ratio_vs_mkk.pdf                       | Bin 0 -> 24061 bytes
 results/figures/quark_pre_audit_constants/rs_anarchy_max_ratio_vs_mkk.png                       | Bin 0 -> 159006 bytes
 results/figures/quark_pre_audit_constants/rs_anarchy_mkk_min_hist.pdf                           | Bin 0 -> 41452 bytes
 results/figures/quark_pre_audit_constants/rs_anarchy_mkk_min_hist.png                           | Bin 0 -> 219090 bytes
 results/figures/quark_pre_audit_constants/rs_anarchy_mkk_min_hist_by_cvals.pdf                  | Bin 0 -> 31942 bytes
 results/figures/quark_pre_audit_constants/rs_anarchy_mkk_min_hist_by_cvals.png                  | Bin 0 -> 214766 bytes
 results/figures/quark_pre_audit_constants/rs_anarchy_mkk_min_hist_by_pdg_tightness.pdf          | Bin 0 -> 37118 bytes
 results/figures/quark_pre_audit_constants/rs_anarchy_mkk_min_hist_by_pdg_tightness.png          | Bin 0 -> 227456 bytes
 results/figures/quark_pre_audit_constants/rs_anarchy_mkk_min_hist_by_yprior.pdf                 | Bin 0 -> 41202 bytes
 results/figures/quark_pre_audit_constants/rs_anarchy_mkk_min_hist_by_yprior.png                 | Bin 0 -> 233191 bytes
 results/figures/quark_pre_audit_constants/rs_anarchy_mkk_min_hist_by_yprior_gsstar.pdf          | Bin 0 -> 45132 bytes
 results/figures/quark_pre_audit_constants/rs_anarchy_mkk_min_hist_by_yprior_gsstar.png          | Bin 0 -> 238873 bytes
 results/figures/quark_pre_audit_constants/rs_anarchy_mkk_min_hist_gsstar.pdf                    | Bin 0 -> 28458 bytes
 results/figures/quark_pre_audit_constants/rs_anarchy_mkk_min_hist_gsstar.png                    | Bin 0 -> 137511 bytes
 results/figures/quark_pre_audit_constants/rs_anarchy_pdg_pass_fraction.pdf                      | Bin 0 -> 19017 bytes
 results/figures/quark_pre_audit_constants/rs_anarchy_pdg_pass_fraction.png                      | Bin 0 -> 62009 bytes
 results/figures/quark_pre_audit_constants/rs_anarchy_per_system_pass_rate.pdf                   | Bin 0 -> 25806 bytes
 results/figures/quark_pre_audit_constants/rs_anarchy_per_system_pass_rate.png                   | Bin 0 -> 84964 bytes
 results/figures/quark_pre_audit_constants/runA/rs_anarchy_gate_sensitivity_runA_8Mdraws.pdf     | Bin 0 -> 28152 bytes
 results/figures/quark_pre_audit_constants/runA/rs_anarchy_gate_sensitivity_runA_8Mdraws.png     | Bin 0 -> 204557 bytes
 results/figures/quark_pre_audit_constants/runA/rs_anarchy_max_ratio_vs_mkk_runA_8Mdraws.pdf     | Bin 0 -> 24061 bytes
 results/figures/quark_pre_audit_constants/runA/rs_anarchy_max_ratio_vs_mkk_runA_8Mdraws.png     | Bin 0 -> 159006 bytes
 results/figures/quark_pre_audit_constants/runA/rs_anarchy_mkk_min_hist_runA_8Mdraws.pdf         | Bin 0 -> 41452 bytes
 results/figures/quark_pre_audit_constants/runA/rs_anarchy_mkk_min_hist_runA_8Mdraws.png         | Bin 0 -> 219090 bytes
 results/figures/quark_pre_audit_constants/runA/rs_anarchy_pdg_pass_fraction_runA_8Mdraws.pdf    | Bin 0 -> 19017 bytes
 results/figures/quark_pre_audit_constants/runA/rs_anarchy_pdg_pass_fraction_runA_8Mdraws.png    | Bin 0 -> 62009 bytes
 results/figures/quark_pre_audit_constants/runA/rs_anarchy_per_system_pass_rate_runA_8Mdraws.pdf | Bin 0 -> 25806 bytes
 results/figures/quark_pre_audit_constants/runA/rs_anarchy_per_system_pass_rate_runA_8Mdraws.png | Bin 0 -> 84964 bytes
 results/figures/quark_pre_audit_constants/yukawa_d_2sigma_vs_1sigma.pdf                         | Bin 0 -> 36070 bytes
 results/figures/quark_pre_audit_constants/yukawa_d_2sigma_vs_1sigma.png                         | Bin 0 -> 110714 bytes
 results/figures/quark_pre_audit_constants/yukawa_offdiag_spread.pdf                             | Bin 0 -> 21091 bytes
 results/figures/quark_pre_audit_constants/yukawa_offdiag_spread.png                             | Bin 0 -> 68408 bytes
 results/figures/quark_pre_audit_constants/yukawa_per_element_anatomy.pdf                        | Bin 0 -> 32887 bytes
 results/figures/quark_pre_audit_constants/yukawa_per_element_anatomy.png                        | Bin 0 -> 254000 bytes
 results/figures/quark_pre_audit_constants/yukawa_size_envelope_vs_anarchic.pdf                  | Bin 0 -> 26022 bytes
 results/figures/quark_pre_audit_constants/yukawa_size_envelope_vs_anarchic.png                  | Bin 0 -> 121962 bytes
 results/figures/quark_pre_audit_constants/yukawa_u_2sigma_vs_1sigma.pdf                         | Bin 0 -> 37548 bytes
 results/figures/quark_pre_audit_constants/yukawa_u_2sigma_vs_1sigma.png                         | Bin 0 -> 116629 bytes
 58 files changed, 0 insertions(+), 0 deletions(-)
```

### 23. `db02223` — scan(rerun): regenerate 8 RS-anarchy scans under corrected ΔF=2 constants
- SHA: `db0222366cf53184ed90e7009739005580facec9`
- Message: scan(rerun): regenerate 8 RS-anarchy scans under corrected ΔF=2 constants
- Physical/numerical summary: Regenerated eight RS-anarchy scans under corrected DeltaF=2 constants, replacing stale ensemble outputs.
- Files touched (`git show --stat`):
```text
 scan_outputs/rs_anarchy_run3_baseline_20260515T085324/tile_summary.json        | 478 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 scan_outputs/rs_anarchy_run3_moreIR_20260515T085324/tile_summary.json          | 366 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 scan_outputs/rs_anarchy_run3_moreUV_20260515T085324/tile_summary.json          | 366 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 scan_outputs/rs_anarchy_run3_qtop_shifted_20260515T085324/tile_summary.json    | 478 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 scan_outputs/rs_anarchy_runA_20260515T085316/tile_summary.json                 | 478 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 scan_outputs/rs_anarchy_runB_gaussian_3sigma_20260515T085324/tile_summary.json | 478 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 scan_outputs/rs_anarchy_runB_narrow_uniform_20260515T085324/tile_summary.json  | 478 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 scan_outputs/rs_anarchy_runB_wide_uniform_20260515T085325/tile_summary.json    | 478 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 scan_outputs/rs_anarchy_runC_20260515T085323/tile_summary.json                 | 366 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 9 files changed, 3966 insertions(+)
```

### 24. `87c728f` — figures(rerun): regenerate publication plots from corrected scans
- SHA: `87c728fb1652b90f7d2d1bb37c0f237a0e8d9d82`
- Message: figures(rerun): regenerate publication plots from corrected scans
- Physical/numerical summary: Regenerated publication plots from corrected scans, updating displayed M_KK and acceptance numerics.
- Files touched (`git show --stat`):
```text
 results/figures/quark/rs_anarchy_cfw_comparison.pdf                | Bin 0 -> 36453 bytes
 results/figures/quark/rs_anarchy_cfw_comparison.png                | Bin 0 -> 167105 bytes
 results/figures/quark/rs_anarchy_gate_sensitivity.pdf              | Bin 0 -> 27529 bytes
 results/figures/quark/rs_anarchy_gate_sensitivity.png              | Bin 0 -> 179051 bytes
 results/figures/quark/rs_anarchy_max_ratio_vs_mkk.pdf              | Bin 0 -> 22827 bytes
 results/figures/quark/rs_anarchy_max_ratio_vs_mkk.png              | Bin 0 -> 151372 bytes
 results/figures/quark/rs_anarchy_mkk_min_hist.pdf                  | Bin 0 -> 39613 bytes
 results/figures/quark/rs_anarchy_mkk_min_hist.png                  | Bin 0 -> 218846 bytes
 results/figures/quark/rs_anarchy_mkk_min_hist_by_cvals.pdf         | Bin 0 -> 33789 bytes
 results/figures/quark/rs_anarchy_mkk_min_hist_by_cvals.png         | Bin 0 -> 208340 bytes
 results/figures/quark/rs_anarchy_mkk_min_hist_by_pdg_tightness.pdf | Bin 0 -> 40451 bytes
 results/figures/quark/rs_anarchy_mkk_min_hist_by_pdg_tightness.png | Bin 0 -> 240873 bytes
 results/figures/quark/rs_anarchy_mkk_min_hist_by_yprior.pdf        | Bin 0 -> 43388 bytes
 results/figures/quark/rs_anarchy_mkk_min_hist_by_yprior.png        | Bin 0 -> 238740 bytes
 results/figures/quark/rs_anarchy_mkk_min_hist_by_yprior_gsstar.pdf | Bin 0 -> 44928 bytes
 results/figures/quark/rs_anarchy_mkk_min_hist_by_yprior_gsstar.png | Bin 0 -> 246600 bytes
 results/figures/quark/rs_anarchy_mkk_min_hist_gsstar.pdf           | Bin 0 -> 28205 bytes
 results/figures/quark/rs_anarchy_mkk_min_hist_gsstar.png           | Bin 0 -> 137367 bytes
 results/figures/quark/rs_anarchy_pdg_pass_fraction.pdf             | Bin 0 -> 17056 bytes
 results/figures/quark/rs_anarchy_pdg_pass_fraction.png             | Bin 0 -> 58505 bytes
 results/figures/quark/rs_anarchy_per_system_pass_rate.pdf          | Bin 0 -> 24923 bytes
 results/figures/quark/rs_anarchy_per_system_pass_rate.png          | Bin 0 -> 92259 bytes
 22 files changed, 0 insertions(+), 0 deletions(-)
```

### 25. `38a4586` — docs(paper): update methodology-note headline numbers post-audit
- SHA: `38a4586399e8e873e38f0a74e9719e494b0f7ea4`
- Message: docs(paper): update methodology-note headline numbers post-audit
- Physical/numerical summary: Updated methodology-note headline numbers after the audit rerun, aligning text with corrected scan outputs.
- Files touched (`git show --stat`):
```text
 docs/quark_scan_methodology_note.pdf         | Bin 580512 -> 587168 bytes
 docs/quark_scan_methodology_note.tex         | 187 +++++++++++++++++++++++++++++++++++++++++++++++++++++++---------------------------------
 scan_outputs/followup_crossings_summary.json | 355 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 3 files changed, 472 insertions(+), 70 deletions(-)
```

### 26. `217af80` — docs(paper): record invalidation-gate rerun clearance
- SHA: `217af802507c650451a7f244c3dc75120cc0ff47`
- Message: docs(paper): record invalidation-gate rerun clearance
- Physical/numerical summary: Updated paper/repository documentation; physical impact is explanatory/provenance text rather than executable numerics.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/invalidation_gate_rerun.md | 137 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 137 insertions(+)
```

### 27. `f19f17c` — docs(paper): seal invalidation-gate rerun signoff
- SHA: `f19f17cc788c7e62c13594475d02d43e681b6c4d`
- Message: docs(paper): seal invalidation-gate rerun signoff
- Physical/numerical summary: Recorded reviewer/sign-off disposition; review ledger only.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/invalidation_gate_signoff.md | 125 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 125 insertions(+)
```

### 28. `ffec986` — audit(cfw): extract CFW 0804.1954 conventions and reference values
- SHA: `ffec98692449e1403b8aa5ad1c3453e49860ca46`
- Message: audit(cfw): extract CFW 0804.1954 conventions and reference values
- Physical/numerical summary: Extracted CFW 0804.1954 convention/reference values for later comparison, adding the physical benchmark ledger.
- Files touched (`git show --stat`):
```text
 docs/audits/cfw_convention_extract.md | 154 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 docs/audits/cfw_vs_ours.md            |  17 +++++++++++++++++
 2 files changed, 171 insertions(+)
```

### 29. `06d5d85` — scripts(cfw): add convention-override flags to comparison driver
- SHA: `06d5d85d63993db12c8da2f1aa6b8ecfce014300`
- Message: scripts(cfw): add convention-override flags to comparison driver
- Physical/numerical summary: Added convention-override flags to the CFW comparison driver, allowing matched g_s and normalization comparisons.
- Files touched (`git show --stat`):
```text
 scripts/rs_anarchy_cfw_comparison.py | 455 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-------------------------------
 1 file changed, 374 insertions(+), 81 deletions(-)
```

### 30. `da8f647` — figures(cfw): regenerate CFW comparison with matched conventions
- SHA: `da8f647ec3e910c93489a76bde98193af0497d25`
- Message: figures(cfw): regenerate CFW comparison with matched conventions
- Physical/numerical summary: Regenerated the CFW comparison figure under matched conventions, updating visual numerical reconciliation.
- Files touched (`git show --stat`):
```text
 results/figures/quark/rs_anarchy_cfw_comparison.pdf | Bin 36453 -> 31874 bytes
 results/figures/quark/rs_anarchy_cfw_comparison.png | Bin 167105 -> 172165 bytes
 2 files changed, 0 insertions(+), 0 deletions(-)
```

### 31. `508ae69` — docs(paper): cite quantitative CFW reconciliation in methodology note
- SHA: `508ae69619a219c6f23e1cc4b9714dc900e33131`
- Message: docs(paper): cite quantitative CFW reconciliation in methodology note
- Physical/numerical summary: Updated paper/repository documentation; physical impact is explanatory/provenance text rather than executable numerics.
- Files touched (`git show --stat`):
```text
 docs/audits/cfw_reproduction.md      | 159 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 docs/quark_scan_methodology_note.pdf | Bin 587168 -> 583860 bytes
 docs/quark_scan_methodology_note.tex | 110 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++--------------------------------------------------
 3 files changed, 219 insertions(+), 50 deletions(-)
```

### 32. `b47aa76` — test(cfw): regression test for the comparison driver
- SHA: `b47aa769990febaa7761e3cad56e618394496f02`
- Message: test(cfw): regression test for the comparison driver
- Physical/numerical summary: Added or refreshed regression/fixture tests pinning the numerical behavior described by the surrounding unit.
- Files touched (`git show --stat`):
```text
 tests/test_cfw_comparison.py | 70 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 70 insertions(+)
```

### 33. `330ffa9` — docs(paper): record Phase 2 hole #7 implementation
- SHA: `330ffa9fbcd85048957a98391453e3a24f1f232d`
- Message: docs(paper): record Phase 2 hole #7 implementation
- Physical/numerical summary: Updated paper/repository documentation; physical impact is explanatory/provenance text rather than executable numerics.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/phase2_h7_impl.md | 86 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 86 insertions(+)
```

### 34. `11e0d58` — figures(cfw): show CFW markers at both g_s*=3 and g_s*~6 conventions
- SHA: `11e0d58b98a28d4820920ce03fe082125db67e8f`
- Message: figures(cfw): show CFW markers at both g_s*=3 and g_s*~6 conventions
- Physical/numerical summary: Regenerated or annotated figures from existing scan/audit data, changing rendered numerical presentation.
- Files touched (`git show --stat`):
```text
 results/figures/quark/rs_anarchy_cfw_comparison.pdf | Bin 31874 -> 33436 bytes
 results/figures/quark/rs_anarchy_cfw_comparison.png | Bin 172165 -> 194099 bytes
 scripts/rs_anarchy_cfw_comparison.py                |  47 ++++++++++++++++++++++++++++++++++++++++-------
 3 files changed, 40 insertions(+), 7 deletions(-)
```

### 35. `321fe0f` — docs(paper): honest CFW reconciliation - factor-2.2 stronger, not 11% agreement
- SHA: `321fe0faf7acc0922076aa37b798e9de473614ea`
- Message: docs(paper): honest CFW reconciliation - factor-2.2 stronger, not 11% agreement
- Physical/numerical summary: Corrected the CFW reconciliation claim from near-agreement to a factor-2.2 stronger result at matched conventions.
- Files touched (`git show --stat`):
```text
 docs/quark_scan_methodology_note.pdf | Bin 583860 -> 586942 bytes
 docs/quark_scan_methodology_note.tex |  54 +++++++++++++++++++++++++++++++++++-------------------
 2 files changed, 35 insertions(+), 19 deletions(-)
```

### 36. `7c284a8` — audit(cfw): update reconciliation doc + mapping table with corrected agreement claim
- SHA: `7c284a86e7056716e8331a45608f07bceafbd0dd`
- Message: audit(cfw): update reconciliation doc + mapping table with corrected agreement claim
- Physical/numerical summary: Added or updated an audit inventory/reconciliation document, turning numerical/physics checks into explicit provenance.
- Files touched (`git show --stat`):
```text
 docs/audits/cfw_reproduction.md | 47 ++++++++++++++++++++++++++++++++++-------------
 docs/audits/cfw_vs_ours.md      |  9 ++++++---
 2 files changed, 40 insertions(+), 16 deletions(-)
```

### 37. `2a3a2c0` — test(cfw): assert factor-2.2 result at matched conventions
- SHA: `2a3a2c08a25b177240f00b228ff9a95480b14c46`
- Message: test(cfw): assert factor-2.2 result at matched conventions
- Physical/numerical summary: Added a regression test pinning the factor-2.2 CFW comparison at matched conventions.
- Files touched (`git show --stat`):
```text
 tests/test_cfw_comparison.py | 64 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 64 insertions(+)
```

### 38. `f2c29c8` — docs(h7): log CFW peer-review revision
- SHA: `f2c29c8a5e5b9953899458b89a650195da6893aa`
- Message: docs(h7): log CFW peer-review revision
- Physical/numerical summary: Recorded review/signoff state for the paper or phase logs; bookkeeping over already-reviewed numerical work.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/phase2_h7_impl.md | 57 +++++++++++++++++++++++++++++++++++++++++++++++++++++----
 1 file changed, 53 insertions(+), 4 deletions(-)
```

### 39. `ad2e3ee` — docs(paper): seal Phase 2 hole #7 review + revision + signoff
- SHA: `ad2e3ee8a169e55634f20286584a10859e8854ef`
- Message: docs(paper): seal Phase 2 hole #7 review + revision + signoff
- Physical/numerical summary: Recorded reviewer/sign-off disposition; review ledger only.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/phase2_h7_review.md    | 49 +++++++++++++++++++++++++++++++++++++++++++++++++
 docs/phase_logs/phase2_h7_review_v2.md | 23 +++++++++++++++++++++++
 docs/phase_logs/phase2_h7_signoff.md   | 86 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 3 files changed, 158 insertions(+)
```

### 40. `f34f9df` — physics(stats): add Wilson-score upper-limit helper for zero-pass observations
- SHA: `f34f9df203d28a711c0514398a1c4a02fe6476a6`
- Message: physics(stats): add Wilson-score upper-limit helper for zero-pass observations
- Physical/numerical summary: Added a Wilson-score upper-limit helper for finite zero-pass observations, introducing bounded pass-rate numerics.
- Files touched (`git show --stat`):
```text
 quarkConstraints/finite_stats.py | 26 ++++++++++++++++++++++++++
 tests/test_finite_stats.py       | 20 ++++++++++++++++++++
 2 files changed, 46 insertions(+)
```

### 41. `a9d1058` — audit(zero-pass): inventory + 95% Wilson UL on moreUV/moreIR/runC
- SHA: `a9d10589db9f0d205a1c00e44e30ba826ed55f0c`
- Message: audit(zero-pass): inventory + 95% Wilson UL on moreUV/moreIR/runC
- Physical/numerical summary: Inventoried zero-pass ensembles and computed 95% Wilson upper limits for moreUV/moreIR/runC.
- Files touched (`git show --stat`):
```text
 docs/audits/zero_pass_inventory.md | 54 ++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 54 insertions(+)
```

### 42. `212320f` — figures(zero-pass): annotate plot legends with finite-ensemble upper limits
- SHA: `212320fd399295a46a1649d565ddf38ab8ad5fd5`
- Message: figures(zero-pass): annotate plot legends with finite-ensemble upper limits
- Physical/numerical summary: Annotated zero-pass figures with finite-ensemble upper limits, replacing impossible-zero language in plots.
- Files touched (`git show --stat`):
```text
 results/figures/quark/rs_anarchy_mkk_min_hist_by_cvals.pdf | Bin 33789 -> 35472 bytes
 results/figures/quark/rs_anarchy_mkk_min_hist_by_cvals.png | Bin 208340 -> 216619 bytes
 scripts/rs_anarchy_cfw_comparison.py                       |  54 +++++++++++++++++++++++++++++++++++++++++++++++++++++-
 scripts/rs_anarchy_mkk_min_hist_by_cvals.py                |  51 +++++++++++++++++++++++++++++++++++++++++++++++++--
 4 files changed, 102 insertions(+), 3 deletions(-)
```

### 43. `1d9e17e` — docs(paper): replace zero-pass impossibility wording with binomial bounds
- SHA: `1d9e17ea1269ed5c552198fd30acfab7fffa39f5`
- Message: docs(paper): replace zero-pass impossibility wording with binomial bounds
- Physical/numerical summary: Replaced zero-pass impossibility wording with binomial-bound language in the paper text.
- Files touched (`git show --stat`):
```text
 docs/audits/zero_pass_inventory.md           |  25 +++++++++++++------------
 docs/phase_logs/invalidation_gate_rerun.md   |   6 +++---
 docs/phase_logs/invalidation_gate_signoff.md |  12 ++++++++----
 docs/phase_logs/phase2_h8_impl.md            |  50 ++++++++++++++++++++++++++++++++++++++++++++++++++
 docs/quark_scan_methodology_note.pdf         | Bin 586942 -> 590665 bytes
 docs/quark_scan_methodology_note.tex         |  49 +++++++++++++++++++++++++++++++++++++------------
 6 files changed, 111 insertions(+), 31 deletions(-)
```

### 44. `715b679` — docs(paper): seal Phase 2 hole #8 review + signoff
- SHA: `715b67922ddc16e5ada22fb4e220467c19607bcc`
- Message: docs(paper): seal Phase 2 hole #8 review + signoff
- Physical/numerical summary: Recorded reviewer/sign-off disposition; review ledger only.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/phase2_h8_review.md  |  25 +++++++++++++++++++++++++
 docs/phase_logs/phase2_h8_signoff.md | 105 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 2 files changed, 130 insertions(+)
```

### 45. `f911a13` — audit(scope): inventory scope-wording statements across the methodology note
- SHA: `f911a1378540e1940148bda955b535a492abd1a0`
- Message: audit(scope): inventory scope-wording statements across the methodology note
- Physical/numerical summary: Added or updated an audit inventory/reconciliation document, turning numerical/physics checks into explicit provenance.
- Files touched (`git show --stat`):
```text
 docs/audits/scope_wording_inventory.md | 11 +++++++++++
 1 file changed, 11 insertions(+)
```

### 46. `bd8ded1` — docs(paper): add Scope and approximations appendix subsection
- SHA: `bd8ded15dfea030806f0d3446fa8d5c425fd2b4d`
- Message: docs(paper): add Scope and approximations appendix subsection
- Physical/numerical summary: Added a Scope and approximations appendix subsection clarifying the physical domain of the scan claims.
- Files touched (`git show --stat`):
```text
 docs/quark_scan_methodology_note.tex | 53 +++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 53 insertions(+)
```

### 47. `e359faa` — docs(paper): align inline body wording with scope appendix
- SHA: `e359faa0897aa769a4a7380610685622c109a4cc`
- Message: docs(paper): align inline body wording with scope appendix
- Physical/numerical summary: Updated paper/repository documentation; physical impact is explanatory/provenance text rather than executable numerics.
- Files touched (`git show --stat`):
```text
 docs/quark_scan_methodology_note.pdf | Bin 590665 -> 595742 bytes
 docs/quark_scan_methodology_note.tex |  59 +++++++++++++++++++++++++++++++----------------------------
 2 files changed, 31 insertions(+), 28 deletions(-)
```

### 48. `1a107c8` — docs(repo): update README + CLAUDE.md with paper-branch pointer
- SHA: `1a107c8a71d8f047648dcb0f731751d0177f13b1`
- Message: docs(repo): update README + CLAUDE.md with paper-branch pointer
- Physical/numerical summary: Updated paper/repository documentation; physical impact is explanatory/provenance text rather than executable numerics.
- Files touched (`git show --stat`):
```text
 CLAUDE.md                         |  5 +++++
 README.md                         | 11 ++++++++---
 docs/phase_logs/phase2_h9_impl.md | 54 ++++++++++++++++++++++++++++++++++++++++++++++++++++++
 3 files changed, 67 insertions(+), 3 deletions(-)
```

### 49. `760f23d` — docs(paper): seal Phase 2 hole #9 review + signoff
- SHA: `760f23daa6c0db8e126bc4a2cf8b97c53ca1d03b`
- Message: docs(paper): seal Phase 2 hole #9 review + signoff
- Physical/numerical summary: Recorded reviewer/sign-off disposition; review ledger only.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/phase2_h9_review.md  | 26 ++++++++++++++++++++++++++
 docs/phase_logs/phase2_h9_signoff.md | 50 ++++++++++++++++++++++++++++++++++++++++++++++++++
 2 files changed, 76 insertions(+)
```

### 50. `80641cf` — chore(gitignore): track artifacts/ manifest and checksums
- SHA: `80641cf653a3a304afb40b65ffa0287222eea945`
- Message: chore(gitignore): track artifacts/ manifest and checksums
- Physical/numerical summary: Repository/configuration bookkeeping; no direct physical numerical change unless noted by touched generated artifacts.
- Files touched (`git show --stat`):
```text
 .gitignore | 8 ++++++++
 1 file changed, 8 insertions(+)
```

### 51. `cb27631` — docs(artifacts): seal canonical artifact manifest for paper
- SHA: `cb27631ae90471feb97c51b5309ac71581c175c0`
- Message: docs(artifacts): seal canonical artifact manifest for paper
- Physical/numerical summary: Sealed the canonical artifact manifest for the paper, anchoring figures/data/checksums used by the final narrative.
- Files touched (`git show --stat`):
```text
 artifacts/README.md                         |  17 +++++++++
 artifacts/checksums.sha256                  |  20 +++++++++++
 artifacts/quarkscan_paper_rc1_manifest.json | 329 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 docs/artifact_manifest.md                   |  67 +++++++++++++++++++++++++++++++++++
 docs/phase_logs/phase2_h2_impl.md           |  48 +++++++++++++++++++++++++
 5 files changed, 481 insertions(+)
```

### 52. `b86e66d` — docs(paper): cite artifact manifest in methodology note
- SHA: `b86e66d9e99db66238301a64442e62e72292440a`
- Message: docs(paper): cite artifact manifest in methodology note
- Physical/numerical summary: Updated paper/repository documentation; physical impact is explanatory/provenance text rather than executable numerics.
- Files touched (`git show --stat`):
```text
 docs/quark_scan_methodology_note.pdf | Bin 595742 -> 596248 bytes
 docs/quark_scan_methodology_note.tex |   7 +++++++
 2 files changed, 7 insertions(+)
```

### 53. `2729659` — chore(gitignore): track manifest bundle tarball + sidecar
- SHA: `27296592c9174f5e4ef22469ccb266906a9320de`
- Message: chore(gitignore): track manifest bundle tarball + sidecar
- Physical/numerical summary: Repository/configuration bookkeeping; no direct physical numerical change unless noted by touched generated artifacts.
- Files touched (`git show --stat`):
```text
 .gitignore                                                  |   2 ++
 artifacts/quarkscan_paper_rc1_manifest_bundle.tar.gz        | Bin 0 -> 6067 bytes
 artifacts/quarkscan_paper_rc1_manifest_bundle.tar.gz.sha256 |   1 +
 3 files changed, 3 insertions(+)
```

### 54. `b946abf` — docs(artifacts): add CFW mapping rows to canonical manifest
- SHA: `b946abf2421f91190f9cb50e4865f5c58dbf1cf0`
- Message: docs(artifacts): add CFW mapping rows to canonical manifest
- Physical/numerical summary: Updated paper/repository documentation; physical impact is explanatory/provenance text rather than executable numerics.
- Files touched (`git show --stat`):
```text
 docs/artifact_manifest.md | 2 ++
 1 file changed, 2 insertions(+)
```

### 55. `5e60bc0` — docs(paper): add canonical-artifact-manifest anchor phrase in methodology note
- SHA: `5e60bc0e613e7d483bef4c919eac1af86f0aa237`
- Message: docs(paper): add canonical-artifact-manifest anchor phrase in methodology note
- Physical/numerical summary: Updated paper/repository documentation; physical impact is explanatory/provenance text rather than executable numerics.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/phase2_h2_impl.md    |   6 ++++++
 docs/quark_scan_methodology_note.pdf | Bin 596248 -> 596242 bytes
 docs/quark_scan_methodology_note.tex |   2 +-
 3 files changed, 7 insertions(+), 1 deletion(-)
```

### 56. `c9d9cb5` — docs(paper): seal Phase 2 hole #2 review + revision + signoff
- SHA: `c9d9cb55db77e82b788113701321906af52f9166`
- Message: docs(paper): seal Phase 2 hole #2 review + revision + signoff
- Physical/numerical summary: Recorded reviewer/sign-off disposition; review ledger only.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/phase2_h2_review.md    |  23 +++++++++++++++++++++++
 docs/phase_logs/phase2_h2_review_v2.md |  14 ++++++++++++++
 docs/phase_logs/phase2_h2_signoff.md   | 105 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 3 files changed, 142 insertions(+)
```

### 57. `e7d824d` — chore(figures): move unreferenced/exploratory figures aside
- SHA: `e7d824d4151ec22027baec5bcb2a61e050cfcfa4`
- Message: chore(figures): move unreferenced/exploratory figures aside
- Physical/numerical summary: Moved unreferenced or exploratory figures aside, pruning non-load-bearing visual artifacts from the paper surface.
- Files touched (`git show --stat`):
```text
 .gitignore                                                                           |  13 ++++++++++---
 README.md                                                                            |   7 +++----
 results/figures/quark/{ => exploratory}/fig1_exclusion_boundaries.png                | Bin
 results/figures/quark/{ => exploratory}/fig2_mkk_bound_2007_vs_modern.png            | Bin
 results/figures/quark/{ => exploratory}/rs_anarchy_cfw_comparison.png                | Bin
 results/figures/quark/{ => exploratory}/rs_anarchy_gate_sensitivity.png              | Bin
 results/figures/quark/{ => exploratory}/rs_anarchy_max_ratio_vs_mkk.png              | Bin
 results/figures/quark/{ => exploratory}/rs_anarchy_mkk_min_hist.pdf                  | Bin
 results/figures/quark/{ => exploratory}/rs_anarchy_mkk_min_hist.png                  | Bin
 results/figures/quark/{ => exploratory}/rs_anarchy_mkk_min_hist_by_cvals.png         | Bin
 results/figures/quark/{ => exploratory}/rs_anarchy_mkk_min_hist_by_pdg_tightness.png | Bin
 results/figures/quark/{ => exploratory}/rs_anarchy_mkk_min_hist_by_yprior.pdf        | Bin
 results/figures/quark/{ => exploratory}/rs_anarchy_mkk_min_hist_by_yprior.png        | Bin
 results/figures/quark/{ => exploratory}/rs_anarchy_mkk_min_hist_by_yprior_gsstar.png | Bin
 results/figures/quark/{ => exploratory}/rs_anarchy_mkk_min_hist_gsstar.png           | Bin
 results/figures/quark/{ => exploratory}/rs_anarchy_pdg_pass_fraction.pdf             | Bin
 results/figures/quark/{ => exploratory}/rs_anarchy_pdg_pass_fraction.png             | Bin
 results/figures/quark/{ => exploratory}/rs_anarchy_per_system_pass_rate.png          | Bin
 results/figures/quark/yukawa_size_envelope_vs_anarchic.pdf                           | Bin 0 -> 26022 bytes
 19 files changed, 13 insertions(+), 7 deletions(-)
```

### 58. `3951f41` — audit(figures): inventory figure references and dispositions
- SHA: `3951f41afdf155994d84d6fe0bafc7181accd2c1`
- Message: audit(figures): inventory figure references and dispositions
- Physical/numerical summary: Added or updated an audit inventory/reconciliation document, turning numerical/physics checks into explicit provenance.
- Files touched (`git show --stat`):
```text
 docs/audits/figure_prune_inventory.md | 111 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 111 insertions(+)
```

### 59. `bf6186c` — docs(paper): final PDF rebuild after figure prune
- SHA: `bf6186c7a0a7e933178d36e06c3e3b4f3eb85eb9`
- Message: docs(paper): final PDF rebuild after figure prune
- Physical/numerical summary: Rebuilt the final PDF after figure pruning, refreshing the rendered paper artifact without changing physics code.
- Files touched (`git show --stat`):
```text
 docs/quark_scan_methodology_note.pdf | Bin 596242 -> 596242 bytes
 1 file changed, 0 insertions(+), 0 deletions(-)
```

### 60. `00b222d` — docs(phase): report Phase 3 hole 10 figure hygiene
- SHA: `00b222d23c85ff26e9c4f98ee4d6b8d23f4a95a3`
- Message: docs(phase): report Phase 3 hole 10 figure hygiene
- Physical/numerical summary: Updated paper/repository documentation; physical impact is explanatory/provenance text rather than executable numerics.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/phase3_h10_impl.md | 53 +++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 53 insertions(+)
```

### 61. `e3c0c79` — docs(paper): seal Phase 3 hole #10 review + signoff
- SHA: `e3c0c79e545ded087e08036fb3ce69675997b47d`
- Message: docs(paper): seal Phase 3 hole #10 review + signoff
- Physical/numerical summary: Recorded reviewer/sign-off disposition; review ledger only.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/phase3_h10_review.md  | 27 +++++++++++++++++++++++++++
 docs/phase_logs/phase3_h10_signoff.md | 28 ++++++++++++++++++++++++++++
 2 files changed, 55 insertions(+)
```

### 62. `1faba23` — docs(paper): final acceptance log for quarkscan-paper-rc1
- SHA: `1faba236ffa76df2a92784f6be199a7805b68cdf`
- Message: docs(paper): final acceptance log for quarkscan-paper-rc1
- Physical/numerical summary: Recorded review/signoff state for the paper or phase logs; bookkeeping over already-reviewed numerical work.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/phase3_final_acceptance.md |  62 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 docs/quark_scan_assumptions_compact.pdf    | Bin 204920 -> 204920 bytes
 docs/quark_scan_methodology_note.pdf       | Bin 596242 -> 596242 bytes
 3 files changed, 62 insertions(+)
```

### 63. `dce4d18` — docs(paper): final summary for quarkscan-paper-rc1
- SHA: `dce4d1869f53bf61f58852d2e2f5a6d0f2573343`
- Message: docs(paper): final summary for quarkscan-paper-rc1
- Physical/numerical summary: Recorded review/signoff state for the paper or phase logs; bookkeeping over already-reviewed numerical work.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/phase3_final_summary.md | 26 ++++++++++++++++++++++++++
 1 file changed, 26 insertions(+)
```

### 64. `4be26fd` — docs(orchestrator): post-compaction briefing for rc1 close-out tasks
- SHA: `4be26fd4a772d84e7f175cea32f48681e86d75bb`
- Message: docs(orchestrator): post-compaction briefing for rc1 close-out tasks
- Physical/numerical summary: Updated paper/repository documentation; physical impact is explanatory/provenance text rather than executable numerics.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/POST_COMPACTION_BRIEFING.md | 243 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 243 insertions(+)
```

### 65. `e3a0f1e` — docs(paper): seal Phase 3 final Opus end-to-end physics signoff
- SHA: `e3a0f1e74f7d899aae5e24e6720b68d24e069b7b`
- Message: docs(paper): seal Phase 3 final Opus end-to-end physics signoff
- Physical/numerical summary: Recorded reviewer/sign-off disposition; review ledger only.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/phase3_final_opus_signoff.md | 175 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 175 insertions(+)
```

### 66. `e0b24e8` — docs(paper): rc1.1 text fixes — WARNING-1/2/3 from Opus end-to-end review
- SHA: `e0b24e866b8ee8f8ed85819a5bace1d0c7e87783`
- Message: docs(paper): rc1.1 text fixes — WARNING-1/2/3 from Opus end-to-end review
- Physical/numerical summary: Applied rc1.1 text fixes from final review, including warning-driven wording/numerical phrasing corrections.
- Files touched (`git show --stat`):
```text
 docs/quark_scan_methodology_note.pdf | Bin 596242 -> 620832 bytes
 docs/quark_scan_methodology_note.tex |  51 +++++++++++++++++++++++++++++++++------------------
 2 files changed, 33 insertions(+), 18 deletions(-)
```

### 67. `ebe8a3c` — docs(flavor-catalog): planner v0 — design for discovery-mode catalog subdir
- SHA: `ebe8a3c9e5719918c6dcc1228afa86de892df966`
- Message: docs(flavor-catalog): planner v0 — design for discovery-mode catalog subdir
- Physical/numerical summary: Updated paper/repository documentation; physical impact is explanatory/provenance text rather than executable numerics.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/flavor_catalog_plan_v0.md | 584 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 584 insertions(+)
```

### 68. `8deb291` — docs(artifacts): refresh PDF checksum after rc1.1 text fixes
- SHA: `8deb29176442f44721ab417e9b1e5f389ace5c80`
- Message: docs(artifacts): refresh PDF checksum after rc1.1 text fixes
- Physical/numerical summary: Refreshed the PDF checksum after rc1.1 text fixes, updating artifact provenance only.
- Files touched (`git show --stat`):
```text
 artifacts/checksums.sha256                   |  1 +
 docs/phase_logs/phase3_rc1p1_textfix_impl.md | 59 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 2 files changed, 60 insertions(+)
```

### 69. `e2b8349` — notebooks: re-execute dense_scan_2sigma_vs_1sigma_comparison.ipynb under post-audit ΔF=2 constants
- SHA: `e2b83492df97579fa85393dbd4c36ae745bfa8ec`
- Message: notebooks: re-execute dense_scan_2sigma_vs_1sigma_comparison.ipynb under post-audit ΔF=2 constants
- Physical/numerical summary: Re-executed the named notebook under post-audit DeltaF=2 constants, refreshing stored notebook outputs and checksums.
- Files touched (`git show --stat`):
```text
 notebooks/dense_scan_2sigma_vs_1sigma_comparison.ipynb | 72 +++++++++++++++++++++++++++++++++++++++---------------------------------
 1 file changed, 39 insertions(+), 33 deletions(-)
```

### 70. `5471ccd` — notebooks: re-execute dense_scan_mkk_constraints_pdg2024.ipynb under post-audit ΔF=2 constants
- SHA: `5471ccd9257629250b89abaceb8e37f3cf75078c`
- Message: notebooks: re-execute dense_scan_mkk_constraints_pdg2024.ipynb under post-audit ΔF=2 constants
- Physical/numerical summary: Re-executed the named notebook under post-audit DeltaF=2 constants, refreshing stored notebook outputs and checksums.
- Files touched (`git show --stat`):
```text
 notebooks/dense_scan_mkk_constraints_pdg2024.ipynb | 40 ++++++++++++++++++++--------------------
 1 file changed, 20 insertions(+), 20 deletions(-)
```

### 71. `a4ee3be` — notebooks: re-execute pdg_quark_target_fix_verification.ipynb under post-audit ΔF=2 constants
- SHA: `a4ee3bea34fcfe91c94c7cbddd77673a78f63195`
- Message: notebooks: re-execute pdg_quark_target_fix_verification.ipynb under post-audit ΔF=2 constants
- Physical/numerical summary: Re-executed the named notebook under post-audit DeltaF=2 constants, refreshing stored notebook outputs and checksums.
- Files touched (`git show --stat`):
```text
 notebooks/pdg_quark_target_fix_verification.ipynb | 132 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++---------------------------------------------------------------------
 1 file changed, 63 insertions(+), 69 deletions(-)
```

### 72. `8981648` — docs(flavor-catalog): peer review of plan v0
- SHA: `8981648a5884337addbd75d8566badb706568d74`
- Message: docs(flavor-catalog): peer review of plan v0
- Physical/numerical summary: Recorded review/signoff state for the paper or phase logs; bookkeeping over already-reviewed numerical work.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/flavor_catalog_plan_v0_review.md | 56 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 56 insertions(+)
```

### 73. `8541a2b` — docs(paper): peer review of rc1.1 text fixes
- SHA: `8541a2bf929589a0ae225d76843c38bc1664c6f5`
- Message: docs(paper): peer review of rc1.1 text fixes
- Physical/numerical summary: Recorded review/signoff state for the paper or phase logs; bookkeeping over already-reviewed numerical work.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/phase3_rc1p1_textfix_review.md | 30 ++++++++++++++++++++++++++++++
 1 file changed, 30 insertions(+)
```

### 74. `2121e4e` — notebooks: re-execute rs_anarchy_analysis.ipynb under post-audit ΔF=2 constants
- SHA: `2121e4e545827f45998df4d9f0a91ce4440378db`
- Message: notebooks: re-execute rs_anarchy_analysis.ipynb under post-audit ΔF=2 constants
- Physical/numerical summary: Re-executed the named notebook under post-audit DeltaF=2 constants, refreshing stored notebook outputs and checksums.
- Files touched (`git show --stat`):
```text
 notebooks/rs_anarchy_analysis.ipynb | 235 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++---------------------------------------------------------------------------------------
 1 file changed, 119 insertions(+), 116 deletions(-)
```

### 75. `82ff25e` — docs(paper): Opus re-review sign-off on rc1.1 text fixes — PASS
- SHA: `82ff25e9b517ee5805b36cfe7038f4261ddebaee`
- Message: docs(paper): Opus re-review sign-off on rc1.1 text fixes — PASS
- Physical/numerical summary: Recorded reviewer/sign-off disposition; review ledger only.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/phase3_rc1p1_textfix_opus_signoff.md | 113 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 113 insertions(+)
```

### 76. `73d347d` — docs(paper): Task B notebook re-execution impl report
- SHA: `73d347da30bac2cc8f6faecb4a3aa92c299dbd08`
- Message: docs(paper): Task B notebook re-execution impl report
- Physical/numerical summary: Recorded or reviewed notebook re-execution under corrected constants; numerical output provenance only.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/phase3_notebook_rerun_impl.md | 34 ++++++++++++++++++++++++++++++++++
 1 file changed, 34 insertions(+)
```

### 77. `744e70c` — docs(flavor-catalog): planner v1 — address BLOCKERs from v0 peer review
- SHA: `744e70cd4a90fe5c64665bb94f903968f50ba870`
- Message: docs(flavor-catalog): planner v1 — address BLOCKERs from v0 peer review
- Physical/numerical summary: Recorded review/signoff state for the paper or phase logs; bookkeeping over already-reviewed numerical work.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/flavor_catalog_plan_v1.md | 698 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 698 insertions(+)
```

### 78. `33367d9` — docs(flavor-catalog): peer review of plan v1
- SHA: `33367d99b2878dd614de7fa1fdcf7f52b6edc271`
- Message: docs(flavor-catalog): peer review of plan v1
- Physical/numerical summary: Recorded review/signoff state for the paper or phase logs; bookkeeping over already-reviewed numerical work.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/flavor_catalog_plan_v1_review.md | 30 ++++++++++++++++++++++++++++++
 1 file changed, 30 insertions(+)
```

### 79. `fb37473` — docs(paper): peer review of Task B notebook re-execution
- SHA: `fb374736a4915d73f816c6b6f384f27b6c566ad1`
- Message: docs(paper): peer review of Task B notebook re-execution
- Physical/numerical summary: Recorded or reviewed notebook re-execution under corrected constants; numerical output provenance only.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/phase3_notebook_rerun_review.md | 34 ++++++++++++++++++++++++++++++++++
 1 file changed, 34 insertions(+)
```

### 80. `95f9f63` — docs(flavor-catalog): Opus final approval of plan v1 — APPROVE-FOR-EXECUTION
- SHA: `95f9f63a2697b11c4544e0af1c24491e3baeb139`
- Message: docs(flavor-catalog): Opus final approval of plan v1 — APPROVE-FOR-EXECUTION
- Physical/numerical summary: Updated paper/repository documentation; physical impact is explanatory/provenance text rather than executable numerics.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/flavor_catalog_plan_v1_opus_signoff.md | 103 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 103 insertions(+)
```

### 81. `9a13d16` — docs(paper): Opus sign-off on Task B notebook re-execution — PASS
- SHA: `9a13d1640ccc8c9f8270c8186211e86f29f683ea`
- Message: docs(paper): Opus sign-off on Task B notebook re-execution — PASS
- Physical/numerical summary: Recorded reviewer/sign-off disposition; review ledger only.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/phase3_notebook_rerun_opus_signoff.md | 100 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 100 insertions(+)
```

### 82. `c78fd73` — docs(flavor-catalog): orchestrator decisions on Section I PI questions
- SHA: `c78fd73ef145380cf5f283f9c24c749f60089230`
- Message: docs(flavor-catalog): orchestrator decisions on Section I PI questions
- Physical/numerical summary: Updated paper/repository documentation; physical impact is explanatory/provenance text rather than executable numerics.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/flavor_catalog_orchestrator_decisions.md | 142 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 142 insertions(+)
```

### 83. `83c0178` — flavor-catalog: scaffold discovery-mode subdir per plan v1
- SHA: `83c0178d80adf890ea6e0bf76128c54847972a15`
- Message: flavor-catalog: scaffold discovery-mode subdir per plan v1
- Physical/numerical summary: Applied the change described by the subject; based on stat surface, this is documentation/catalog/orchestration work rather than a new physics-code calculation.
- Files touched (`git show --stat`):
```text
 flavor_catalog/README.md                          | 54 ++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/catalog_index.tex                  | 24 ++++++++++++++++++++++++
 flavor_catalog/catalog_index.yaml                 |  7 +++++++
 flavor_catalog/catalog_master.tex                 | 51 +++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/latex/macros.tex                   | 26 ++++++++++++++++++++++++++
 flavor_catalog/latex/process_template.tex         | 36 ++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/beauty/.gitkeep          |  0
 flavor_catalog/processes/beauty/index.tex         |  1 +
 flavor_catalog/processes/charged_lepton/.gitkeep  |  0
 flavor_catalog/processes/charged_lepton/index.tex |  1 +
 flavor_catalog/processes/charm/.gitkeep           |  0
 flavor_catalog/processes/charm/index.tex          |  1 +
 flavor_catalog/processes/edm_neutrino/.gitkeep    |  0
 flavor_catalog/processes/edm_neutrino/index.tex   |  1 +
 flavor_catalog/processes/kaon/.gitkeep            |  0
 flavor_catalog/processes/kaon/index.tex           |  1 +
 flavor_catalog/processes/top_higgs_ew/.gitkeep    |  0
 flavor_catalog/processes/top_higgs_ew/index.tex   |  1 +
 flavor_catalog/references/.gitkeep                |  0
 flavor_catalog/signoff/by_process/.gitkeep        |  0
 flavor_catalog/signoff/round_index/.gitkeep       |  0
 flavor_catalog/worklogs/checker/.gitkeep          |  0
 flavor_catalog/worklogs/discovery/.gitkeep        |  0
 flavor_catalog/worklogs/pka/.gitkeep              |  0
 flavor_catalog/worklogs/writer/.gitkeep           |  0
 25 files changed, 204 insertions(+)
```

### 84. `a516dc5` — flavor-catalog: record scaffold implementation report
- SHA: `a516dc529a025432fc240f5f3d5de528e3c267d9`
- Message: flavor-catalog: record scaffold implementation report
- Physical/numerical summary: Applied the change described by the subject; based on stat surface, this is documentation/catalog/orchestration work rather than a new physics-code calculation.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/flavor_catalog_scaffold_impl.md | 73 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 73 insertions(+)
```

### 85. `e45eadc` — flavor-catalog(kaon): PKA draft for K006 KL to mumu
- SHA: `e45eadc6dab0d7137cfa02e6c021895dd0e51f9c`
- Message: flavor-catalog(kaon): PKA draft for K006 KL to mumu
- Physical/numerical summary: Added an initial catalog process draft for K006 KL to mumu; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/kaon/K006.tex                                    | 85 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/kaon/K006.yaml                                   | 93 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/K006/ambrose_e871_2000_pubmed.txt               | 22 ++++++++++++++++++++++
 flavor_catalog/references/K006/chao_christ_2024_arxiv.txt                 | 23 +++++++++++++++++++++++
 flavor_catalog/references/K006/dambrosio_kitahara_2017_arxiv.txt          | 23 +++++++++++++++++++++++
 flavor_catalog/references/K006/dery_ghosh_grossman_schacht_2021_arxiv.txt | 21 +++++++++++++++++++++
 flavor_catalog/references/K006/isidori_unterdorfer_2004_arxiv.txt         | 21 +++++++++++++++++++++
 flavor_catalog/references/K006/pdg2024_kl_decay_modes.txt                 | 15 +++++++++++++++
 flavor_catalog/references/K006/source_manifest.yaml                       | 58 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/K006.md                                       | 11 +++++++++++
 10 files changed, 372 insertions(+)
```

### 86. `6fb9d78` — flavor-catalog(kaon): PKA draft for K005 KL pi0 nunu
- SHA: `6fb9d7820a594b3addac3e4e23f475b4f3a5e204`
- Message: flavor-catalog(kaon): PKA draft for K005 KL pi0 nunu
- Physical/numerical summary: Added an initial catalog process draft for K005 KL pi0 nunu; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/kaon/K005.tex                                     |  95 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/kaon/K005.yaml                                    | 124 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/K005/brod_gorbahn_stamou2021_arxiv2105_02868.txt |  11 +++++++++++
 flavor_catalog/references/K005/buchalla_buras1996_arxiv_hep-ph_9607447.txt |  13 +++++++++++++
 flavor_catalog/references/K005/buras_venturini2022_arxiv2203_10099.txt     |  14 ++++++++++++++
 flavor_catalog/references/K005/koto2025_prl134081802_arxiv2411_11237.txt   |  16 ++++++++++++++++
 flavor_catalog/references/K005/na62_2026_arxiv2604_12649.txt               |  18 ++++++++++++++++++
 flavor_catalog/references/K005/pdg2026_pdgLive_S013R40.txt                 |  21 +++++++++++++++++++++
 flavor_catalog/references/K005/source_manifest.yaml                        |  59 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/K005.md                                        |  35 +++++++++++++++++++++++++++++++++++
 10 files changed, 406 insertions(+)
```

### 87. `a36f67e` — flavor-catalog(beauty): PKA draft for B011 B to Xs gamma
- SHA: `a36f67e3c4522dfcb85fc44208155bd2da1f197b`
- Message: flavor-catalog(beauty): PKA draft for B011 B to Xs gamma
- Physical/numerical summary: Added an initial catalog process draft for B011 B to Xs gamma; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B011.tex                                | 100 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/beauty/B011.yaml                               | 108 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B011/belle_2016_arxiv.txt                     |  28 ++++++++++++++++++++++++++++
 flavor_catalog/references/B011/belleii_2022_arxiv.txt                   |  30 ++++++++++++++++++++++++++++++
 flavor_catalog/references/B011/hflav_rare_decays_dec2024_table67.txt    |  34 ++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B011/hflav_results_2026-05-16.txt             |  23 +++++++++++++++++++++++
 flavor_catalog/references/B011/misiak_et_al_2015_arxiv.txt              |  41 +++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B011/misiak_rehman_steinhauser_2020_arxiv.txt |  37 +++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B011/pdg2024_b_meson_prod_decay_excerpt.txt   |  33 +++++++++++++++++++++++++++++++++
 flavor_catalog/references/B011/sha256sums.txt                           |   7 +++++++
 flavor_catalog/references/B011/source_manifest.yaml                     |  67 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/B011.md                                     |  34 ++++++++++++++++++++++++++++++++++
 12 files changed, 542 insertions(+)
```

### 88. `6fae247` — flavor-catalog(kaon): PKA draft for K004 K+ to pi+ nu nubar
- SHA: `6fae24759db24572c6d03e40fd500310bb535cf7`
- Message: flavor-catalog(kaon): PKA draft for K004 K+ to pi+ nu nubar
- Physical/numerical summary: Added an initial catalog process draft for K004 K+ to pi+ nu nubar; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/kaon/K004.tex                                      |  91 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/kaon/K004.yaml                                     | 119 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/K004/brod_gorbahn_stamou_2021_arxiv2105_02868.txt |  23 +++++++++++++++++++++++
 flavor_catalog/references/K004/buras_venturini_2022_arxiv2203_10099.txt     |  23 +++++++++++++++++++++++
 flavor_catalog/references/K004/cfw_2008_arxiv0804_1954.txt                  |  24 ++++++++++++++++++++++++
 flavor_catalog/references/K004/na62_2025_arxiv2412_12015.txt                |  26 ++++++++++++++++++++++++++
 flavor_catalog/references/K004/na62_2026_arxiv2604_12649.txt                |  27 +++++++++++++++++++++++++++
 flavor_catalog/references/K004/sha256sums.txt                               |   5 +++++
 flavor_catalog/references/K004/source_manifest.yaml                         |  51 +++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/K004.md                                         |  38 ++++++++++++++++++++++++++++++++++++++
 10 files changed, 427 insertions(+)
```

### 89. `0a1ce41` — flavor-catalog(kaon): PKA draft for K003 epsilon prime over epsilon
- SHA: `0a1ce41cd019a192fbfd14ac5f28678a1ff5e60e`
- Message: flavor-catalog(kaon): PKA draft for K003 epsilon prime over epsilon
- Physical/numerical summary: Added an initial catalog process draft for K003 epsilon prime over epsilon; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/kaon/K003.tex                                      |  97 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/kaon/K003.yaml                                     | 141 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/K003/arxiv_0804.1954_cfw_2008.txt                 |  18 ++++++++++++++++++
 flavor_catalog/references/K003/arxiv_1011.0127_ktev_2011.txt                |  16 ++++++++++++++++
 flavor_catalog/references/K003/arxiv_1807.02520_bsm_master_formula_2019.txt |  25 ++++++++++++++++++++++++
 flavor_catalog/references/K003/arxiv_1808.00466_bsm_anatomy_2019.txt        |  19 +++++++++++++++++++
 flavor_catalog/references/K003/arxiv_2004.09440_rbcukqcd_2020.txt           |  20 ++++++++++++++++++++
 flavor_catalog/references/K003/arxiv_2005.05978_aebischer_buras_2020.txt    |  20 ++++++++++++++++++++
 flavor_catalog/references/K003/arxiv_hep-ex_0208009_na48_2002.txt           |  17 +++++++++++++++++
 flavor_catalog/references/K003/pdg_live_S013EPS_datablock_20260516.txt      |  24 +++++++++++++++++++++++
 flavor_catalog/references/K003/pdg_live_S013_20260516.txt                   |  15 +++++++++++++++
 flavor_catalog/references/K003/source_manifest.yaml                         |  86 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/K003.md                                         |  39 ++++++++++++++++++++++++++++++++++++++
 13 files changed, 537 insertions(+)
```

### 90. `be62498` — flavor-catalog(top_higgs_ew): PKA draft for T001 t to cZ
- SHA: `be624981c1c87b132e90a7805afd295c1d68d08a`
- Message: flavor-catalog(top_higgs_ew): PKA draft for T001 t to cZ
- Physical/numerical summary: Added an initial catalog process draft for T001 t to cZ; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T001.tex                                |  85 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/top_higgs_ew/T001.yaml                               | 131 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/T001/aguilar_saavedra_hepph0409342.txt              |  34 ++++++++++++++++++++++++++++++++++
 flavor_catalog/references/T001/airen_franceschini_2601_14966_arxiv_abs.txt    |  29 +++++++++++++++++++++++++++++
 flavor_catalog/references/T001/atlas_2301_11605_arxiv_abs.txt                 |  42 ++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/T001/cms_1702_01404_arxiv_abs.txt                   |  34 ++++++++++++++++++++++++++++++++++
 flavor_catalog/references/T001/csaki_falkowski_weiler_0804_1954_arxiv_abs.txt |  24 ++++++++++++++++++++++++
 flavor_catalog/references/T001/pdg_2025_top_fcnc_tzq.txt                      |  35 +++++++++++++++++++++++++++++++++++
 flavor_catalog/references/T001/sha256sums.txt                                 |   6 ++++++
 flavor_catalog/references/T001/source_manifest.yaml                           |  58 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/T001.md                                           |  35 +++++++++++++++++++++++++++++++++++
 11 files changed, 513 insertions(+)
```

### 91. `18fe48f` — flavor-catalog(top_higgs_ew): PKA draft for T010 Zbb pole observables
- SHA: `18fe48f2c8c72241b0d512db0c077fb715cd41a8`
- Message: flavor-catalog(top_higgs_ew): PKA draft for T010 Zbb pole observables
- Physical/numerical summary: Added an initial catalog process draft for T010 Zbb pole observables; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T010.tex                |  80 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/top_higgs_ew/T010.yaml               | 105 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/T010/casagrande_2008_rs_ewpt.txt    |  17 +++++++++++++++++
 flavor_catalog/references/T010/cfw_2008_rs_flavor.txt         |  18 ++++++++++++++++++
 flavor_catalog/references/T010/fcc_ee_2025_zbb_projection.txt |  17 +++++++++++++++++
 flavor_catalog/references/T010/freitas_2014_z_widths.txt      |  17 +++++++++++++++++
 flavor_catalog/references/T010/lepslc_2006_z_resonance.txt    |  19 +++++++++++++++++++
 flavor_catalog/references/T010/pdg_2025_z_boson.txt           |  22 ++++++++++++++++++++++
 flavor_catalog/references/T010/source_manifest.yaml           |  58 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/T010.md                           |  34 ++++++++++++++++++++++++++++++++++
 10 files changed, 387 insertions(+)
```

### 92. `0155d7c` — flavor-catalog(kaon): PKA draft for K013 KL pi0 gamma gamma
- SHA: `0155d7c9fad2b556d2da978432903d7bbd918f82`
- Message: flavor-catalog(kaon): PKA draft for K013 KL pi0 gamma gamma
- Physical/numerical summary: Added an initial catalog process draft for K013 KL pi0 gamma gamma; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/kaon/K013.tex                                    | 39 +++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/kaon/K013.yaml                                   | 89 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/K013/cappiello_cata_dambrosio_2018_springer.txt | 24 ++++++++++++++++++++++++
 flavor_catalog/references/K013/gabbiani_valencia_2001_arxiv.txt           | 17 +++++++++++++++++
 flavor_catalog/references/K013/ktev_abouzaid_2008_arxiv.txt               | 21 +++++++++++++++++++++
 flavor_catalog/references/K013/na48_2_batley_2014_arxiv.txt               | 19 +++++++++++++++++++
 flavor_catalog/references/K013/na48_lai_2002_repo.txt                     | 21 +++++++++++++++++++++
 flavor_catalog/references/K013/pdg2025_kl_pi0gammagamma.txt               | 27 +++++++++++++++++++++++++++
 flavor_catalog/references/K013/sha256sums.txt                             |  6 ++++++
 flavor_catalog/references/K013/source_manifest.yaml                       | 57 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/K013.md                                       | 11 +++++++++++
 11 files changed, 331 insertions(+)
```

### 93. `47e4339` — flavor-catalog(beauty): PKA draft for B009 B+ to tau nu
- SHA: `47e4339f464e0a32b28752da145d36ac98bf977d`
- Message: flavor-catalog(beauty): PKA draft for B009 B+ to tau nu
- Physical/numerical summary: Added an initial catalog process draft for B009 B+ to tau nu; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B009.tex                       | 101 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/beauty/B009.yaml                      | 146 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B009/babar2010_arxiv0912_2453.txt    |  15 +++++++++++++++
 flavor_catalog/references/B009/babar2013_arxiv1207_0698.txt    |  14 ++++++++++++++
 flavor_catalog/references/B009/belle2013_arxiv1208_4678.txt    |  14 ++++++++++++++
 flavor_catalog/references/B009/belle2015_arxiv1503_05613.txt   |  14 ++++++++++++++
 flavor_catalog/references/B009/belleii2025_arxiv2502_04885.txt |  15 +++++++++++++++
 flavor_catalog/references/B009/cfw2008_arxiv0804_1954.txt      |  14 ++++++++++++++
 flavor_catalog/references/B009/hflav_dec2025_btaunu.txt        |  18 ++++++++++++++++++
 flavor_catalog/references/B009/pdg2025_btaunu_api.txt          |  23 +++++++++++++++++++++++
 flavor_catalog/references/B009/source_manifest.yaml            |  85 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B009/utfit_summer2024_btaunu.txt     |  18 ++++++++++++++++++
 flavor_catalog/worklogs/pka/B009.md                            |  40 ++++++++++++++++++++++++++++++++++++++++
 13 files changed, 517 insertions(+)
```

### 94. `c4cea8a` — flavor-catalog(beauty): PKA draft for B015 inclusive b to s ll
- SHA: `c4cea8a7cc0c2cf67d60a43428ed050e9d31eb18`
- Message: flavor-catalog(beauty): PKA draft for B015 inclusive b to s ll
- Physical/numerical summary: Added an initial catalog process draft for B015 inclusive b to s ll; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B015.tex                               |  82 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/beauty/B015.yaml                              | 111 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B015/babar_1312_5364_arxiv.txt               |  23 +++++++++++++++++++++++
 flavor_catalog/references/B015/belle_hepex_0503044_arxiv.txt           |  17 +++++++++++++++++
 flavor_catalog/references/B015/cfw_0804_1954_arxiv.txt                 |  21 +++++++++++++++++++++
 flavor_catalog/references/B015/hflav_dec2024_B_to_Xsll.txt             |  21 +++++++++++++++++++++
 flavor_catalog/references/B015/huber_hurth_lunghi_1503_04849_arxiv.txt |  33 +++++++++++++++++++++++++++++++++
 flavor_catalog/references/B015/miscited_prompt_arxiv_ids.txt           |  22 ++++++++++++++++++++++
 flavor_catalog/references/B015/sha256sums.txt                          |   6 ++++++
 flavor_catalog/references/B015/source_manifest.yaml                    |  64 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/B015.md                                    |  11 +++++++++++
 11 files changed, 411 insertions(+)
```

### 95. `1cf18c6` — flavor-catalog(top_higgs_ew): WA batch wa_wave1_top_higgs_ew — polish PKA drafts for T001 T010
- SHA: `1cf18c61d7fbf170ffb27ea75e76435c54371df4`
- Message: flavor-catalog(top_higgs_ew): WA batch wa_wave1_top_higgs_ew — polish PKA drafts for T001 T010
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T001.tex          | 57 ++++++++++++++++++++++++++++++++-------------------------
 flavor_catalog/processes/top_higgs_ew/T001.yaml         | 10 +++++++++-
 flavor_catalog/processes/top_higgs_ew/T010.tex          | 32 ++++++++++++++++++--------------
 flavor_catalog/processes/top_higgs_ew/T010.yaml         | 10 +++++++++-
 flavor_catalog/worklogs/writer/wa_wave1_top_higgs_ew.md | 32 ++++++++++++++++++++++++++++++++
 5 files changed, 100 insertions(+), 41 deletions(-)
```

### 96. `73b138a` — flavor-catalog(beauty): WA batch wa_wave1_beauty — polish PKA drafts for B009 B011 B015
- SHA: `73b138a42cb52ff8f7c210f6fcdaac9c64e7e702`
- Message: flavor-catalog(beauty): WA batch wa_wave1_beauty — polish PKA drafts for B009 B011 B015
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B009.tex          | 25 +++++++++++--------------
 flavor_catalog/processes/beauty/B009.yaml         | 10 +++++++++-
 flavor_catalog/processes/beauty/B011.tex          | 29 +++++++++++++++--------------
 flavor_catalog/processes/beauty/B011.yaml         | 10 +++++++++-
 flavor_catalog/processes/beauty/B015.tex          | 47 +++++++++++++++++++++++++----------------------
 flavor_catalog/processes/beauty/B015.yaml         | 10 +++++++++-
 flavor_catalog/worklogs/writer/wa_wave1_beauty.md | 42 ++++++++++++++++++++++++++++++++++++++++++
 7 files changed, 120 insertions(+), 53 deletions(-)
```

### 97. `50311f6` — flavor-catalog(kaon): WA batch wa_wave1_kaon — polish PKA drafts for K003 K004 K005 K006 K013
- SHA: `50311f6b81af6b6c61ba76b501d90c2c9d76ac29`
- Message: flavor-catalog(kaon): WA batch wa_wave1_kaon — polish PKA drafts for K003 K004 K005 K006 K013
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/kaon/K003.tex          | 31 +++++++++++++++----------------
 flavor_catalog/processes/kaon/K003.yaml         | 11 +++++++++--
 flavor_catalog/processes/kaon/K004.tex          | 30 +++++++++++++++---------------
 flavor_catalog/processes/kaon/K004.yaml         | 10 +++++++++-
 flavor_catalog/processes/kaon/K005.tex          | 14 +++++++-------
 flavor_catalog/processes/kaon/K005.yaml         | 10 +++++++++-
 flavor_catalog/processes/kaon/K006.tex          | 18 ++++++++----------
 flavor_catalog/processes/kaon/K006.yaml         | 10 +++++++++-
 flavor_catalog/processes/kaon/K013.tex          | 56 +++++++++++++++++++++++++++++++++++++++++++++++++-------
 flavor_catalog/processes/kaon/K013.yaml         | 10 +++++++++-
 flavor_catalog/worklogs/writer/wa_wave1_kaon.md | 83 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 11 files changed, 222 insertions(+), 61 deletions(-)
```

### 98. `9bb58af` — flavor-catalog(top_higgs_ew): PKA draft for T002 t -> u Z
- SHA: `9bb58afaac3f933b4d9687a35b087caf9fc5d524`
- Message: flavor-catalog(top_higgs_ew): PKA draft for T002 t -> u Z
- Physical/numerical summary: Added an initial catalog process draft for T002 t -> u Z; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T002.tex                                |  88 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/top_higgs_ew/T002.yaml                               | 131 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/T002/aguilar_saavedra_hepph0409342.txt              |  33 +++++++++++++++++++++++++++++++++
 flavor_catalog/references/T002/atlas_2301_11605_arxiv_abs.txt                 |  42 ++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/T002/casagrande_0807_4937_rs_top_fcnc.txt           |  40 ++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/T002/cms_1702_01404_arxiv_abs.txt                   |  34 ++++++++++++++++++++++++++++++++++
 flavor_catalog/references/T002/csaki_falkowski_weiler_0804_1954_arxiv_abs.txt |  25 +++++++++++++++++++++++++
 flavor_catalog/references/T002/pdg_2025_top_fcnc_tzq.txt                      |  35 +++++++++++++++++++++++++++++++++++
 flavor_catalog/references/T002/source_manifest.yaml                           |  58 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/T002.md                                           |  38 ++++++++++++++++++++++++++++++++++++++
 10 files changed, 524 insertions(+)
```

### 99. `fcbe870` — flavor-catalog(beauty): PKA draft for B002 S_psiK_S
- SHA: `fcbe870e90a92b66b59b3c4717f6680d40662977`
- Message: flavor-catalog(beauty): PKA draft for B002 S_psiK_S
- Physical/numerical summary: Added an initial catalog process draft for B002 S_psiK_S; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B002.tex                                   | 100 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/beauty/B002.yaml                                  | 132 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B002/belleii_2024_sin2phi1_arxiv.txt             |  19 +++++++++++++++++++
 flavor_catalog/references/B002/csaki_falkowski_weiler_2008_arxiv.txt       |  17 +++++++++++++++++
 flavor_catalog/references/B002/faller_fleischer_jung_mannel_2008_arxiv.txt |  17 +++++++++++++++++
 flavor_catalog/references/B002/frings_nierste_wiebusch_2015_arxiv.txt      |  19 +++++++++++++++++++
 flavor_catalog/references/B002/hflav_triangle_summer2025_sin2beta.txt      |  33 +++++++++++++++++++++++++++++++++
 flavor_catalog/references/B002/lhcb_2024_spsiks_arxiv.txt                  |  20 ++++++++++++++++++++
 flavor_catalog/references/B002/source_manifest.yaml                        |  58 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/B002.md                                        |  37 +++++++++++++++++++++++++++++++++++++
 10 files changed, 452 insertions(+)
```

### 100. `2f443a0` — flavor-catalog(edm_neutrino): PKA draft for E001 electron EDM
- SHA: `2f443a0f0856149fd217a58510b3dbc20f1900f7`
- Message: flavor-catalog(edm_neutrino): PKA draft for E001 electron EDM
- Physical/numerical summary: Added an initial catalog process draft for E001 electron EDM; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/edm_neutrino/E001.tex                    | 90 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/edm_neutrino/E001.yaml                   | 94 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/E001/andreev2018_caltechauthors.txt     | 25 +++++++++++++++++++++++++
 flavor_catalog/references/E001/cfw2008_arxiv0804_1954.txt         | 24 ++++++++++++++++++++++++
 flavor_catalog/references/E001/panico2019_arxiv1810_09413.txt     | 19 +++++++++++++++++++
 flavor_catalog/references/E001/pdg2026_electron_edm_datablock.txt | 33 +++++++++++++++++++++++++++++++++
 flavor_catalog/references/E001/roussy2023_arxiv2212_11841.txt     | 27 +++++++++++++++++++++++++++
 flavor_catalog/references/E001/source_manifest.yaml               | 55 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/E001.md                               | 31 +++++++++++++++++++++++++++++++
 9 files changed, 398 insertions(+)
```

### 101. `21661a4` — flavor-catalog(kaon): PKA draft for K001 epsilon_K
- SHA: `21661a4e9048776b496dde7375af962a955f8882`
- Message: flavor-catalog(kaon): PKA draft for K001 epsilon_K
- Physical/numerical summary: Added an initial catalog process draft for K001 epsilon_K; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/kaon/K001.tex                               |  91 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/kaon/K001.yaml                              | 115 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/K001/bgs2020_epsilon_k_arxiv1911_06822.txt |  18 ++++++++++++++++++
 flavor_catalog/references/K001/cfw2008_rs_flavor_arxiv0804_1954.txt  |  18 ++++++++++++++++++
 flavor_catalog/references/K001/flag2024_bk_arxiv2411_04268.txt       |  18 ++++++++++++++++++
 flavor_catalog/references/K001/pdg2026_epsilon_k.txt                 |  17 +++++++++++++++++
 flavor_catalog/references/K001/source_manifest.yaml                  |  45 +++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/K001.md                                  |  31 +++++++++++++++++++++++++++++++
 8 files changed, 353 insertions(+)
```

### 102. `9a92781` — flavor-catalog(beauty): PKA draft for B005 Bs to mu mu
- SHA: `9a92781374bddd04f33f51e229c89b46b7914d1f`
- Message: flavor-catalog(beauty): PKA draft for B005 Bs to mu mu
- Physical/numerical summary: Added an initial catalog process draft for B005 Bs to mu mu; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B005.tex                             |  97 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/beauty/B005.yaml                            | 106 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B005/atlas_1812_03017_bsmumu.txt           |  15 +++++++++++++++
 flavor_catalog/references/B005/cfw_0804_1954_rs_flavor.txt           |  15 +++++++++++++++
 flavor_catalog/references/B005/cms_2212_10311_bsmumu.txt             |  22 ++++++++++++++++++++++
 flavor_catalog/references/B005/czaja_misiak_2407_03810_bsmumu_sm.txt |  17 +++++++++++++++++
 flavor_catalog/references/B005/hflav_2023_bsmumu.txt                 |  19 +++++++++++++++++++
 flavor_catalog/references/B005/lhcb_2108_09283_bsmumu.txt            |  22 ++++++++++++++++++++++
 flavor_catalog/references/B005/pdg_2026_bsmumu.txt                   |  33 +++++++++++++++++++++++++++++++++
 flavor_catalog/references/B005/source_manifest.yaml                  |  76 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/B005.md                                  |  34 ++++++++++++++++++++++++++++++++++
 11 files changed, 456 insertions(+)
```

### 103. `6e60d02` — flavor-catalog(charm): PKA draft for C001 D0 mixing
- SHA: `6e60d020443d8f645c2688a8baaf056e3ff7ebe1`
- Message: flavor-catalog(charm): PKA draft for C001 D0 mixing
- Physical/numerical summary: Added an initial catalog process draft for C001 D0 mixing; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charm/C001.tex                             |  92 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/charm/C001.yaml                            | 125 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/C001/cfw2008_arxiv0804_1954.txt           |  21 +++++++++++++++++++++
 flavor_catalog/references/C001/hflav_ckm25_delta_y_input.txt        |  23 +++++++++++++++++++++++
 flavor_catalog/references/C001/hflav_ckm25_dmixing_global_fit.txt   |  27 +++++++++++++++++++++++++++
 flavor_catalog/references/C001/lhcb2021_arxiv2106_03744.txt         |  18 ++++++++++++++++++
 flavor_catalog/references/C001/pdg2025_deltam_d_pdgLive.txt         |  21 +++++++++++++++++++++
 flavor_catalog/references/C001/pdg2025_dmixing_review_table70_7.txt |  22 ++++++++++++++++++++++
 flavor_catalog/references/C001/sha256sums.txt                       |   6 ++++++
 flavor_catalog/references/C001/source_manifest.yaml                 |  64 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/C001.md                                 |  14 ++++++++++++++
 11 files changed, 433 insertions(+)
```

### 104. `9068956` — flavor-catalog(charged_lepton): PKA draft for L001 mu -> e gamma
- SHA: `90689565e04f01e4662c57169104a7eca12a3c73`
- Message: flavor-catalog(charged_lepton): PKA draft for L001 mu -> e gamma
- Physical/numerical summary: Added an initial catalog process draft for L001 mu -> e gamma; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charged_lepton/L001.tex                     |  84 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/charged_lepton/L001.yaml                    | 151 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/L001/meg2016_arxiv1605_05081.txt           |  23 ++++++++++++++++++++++
 flavor_catalog/references/L001/megii2024_arxiv2310_12614.txt         |  29 ++++++++++++++++++++++++++++
 flavor_catalog/references/L001/megii2025_arxiv2504_15711.txt         |  33 ++++++++++++++++++++++++++++++++
 flavor_catalog/references/L001/pdg2025_muon_listing_s004.txt         |  32 +++++++++++++++++++++++++++++++
 flavor_catalog/references/L001/perez_randall_2008_arxiv0805_4652.txt |  38 ++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/L001/sha256sums.txt                        |   5 +++++
 flavor_catalog/references/L001/source_manifest.yaml                  |  55 ++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/L001.md                                  |  34 ++++++++++++++++++++++++++++++++
 10 files changed, 484 insertions(+)
```

### 105. `6027db0` — flavor-catalog(kaon): PKA draft for K002 Delta m_K
- SHA: `6027db0b328911f793c1d1128b2dfc4643de4227`
- Message: flavor-catalog(kaon): PKA draft for K002 Delta m_K
- Physical/numerical summary: Added an initial catalog process draft for K002 Delta m_K; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/kaon/K002.tex                                          |  89 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/kaon/K002.yaml                                         | 145 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/K002/arxiv_0804.1954_cfw_2008.txt                     |  17 ++++++++++++++++
 flavor_catalog/references/K002/arxiv_1011.0127_ktev_2011.txt                    |  19 ++++++++++++++++++
 flavor_catalog/references/K002/arxiv_1212.5931_christ_2013_ld_delta_mk.txt      |  22 ++++++++++++++++++++
 flavor_catalog/references/K002/arxiv_1406.0916_bai_2014_delta_mk.txt            |  21 +++++++++++++++++++
 flavor_catalog/references/K002/arxiv_2001.06374_wang_2020_physical_delta_mk.txt |  18 +++++++++++++++++
 flavor_catalog/references/K002/arxiv_2301.01387_wang_2023_delta_mk.txt          |  19 ++++++++++++++++++
 flavor_catalog/references/K002/arxiv_2404.02297_boyle_2024_bsm_kaon_mixing.txt  |  28 ++++++++++++++++++++++++++
 flavor_catalog/references/K002/arxiv_2411.04268_flag2024.txt                    |  20 +++++++++++++++++++
 flavor_catalog/references/K002/pdg_live_S013D_20260516.txt                      |  27 +++++++++++++++++++++++++
 flavor_catalog/references/K002/source_manifest.yaml                             |  85 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/K002.md                                             |  39 ++++++++++++++++++++++++++++++++++++
 13 files changed, 549 insertions(+)
```

### 106. `9b49fae` — flavor-catalog(beauty): CA batch ca_wave1_beauty — verify WA polish for B009 B011 B015
- SHA: `9b49fae1619be13f886be7a7b8544c591c6d3cae`
- Message: flavor-catalog(beauty): CA batch ca_wave1_beauty — verify WA polish for B009 B011 B015
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B009.yaml          | 14 ++++++++++++--
 flavor_catalog/processes/beauty/B011.yaml          | 14 +++++++++++---
 flavor_catalog/processes/beauty/B015.yaml          | 15 +++++++++++++--
 flavor_catalog/worklogs/checker/ca_wave1_beauty.md | 40 ++++++++++++++++++++++++++++++++++++++++
 4 files changed, 76 insertions(+), 7 deletions(-)
```

### 107. `4e6b7a9` — flavor-catalog(kaon): CA batch ca_wave1_kaon — verify WA polish for K003 K004 K005 K006 K013
- SHA: `4e6b7a9b04ec32613d3de8900da564d33da591b7`
- Message: flavor-catalog(kaon): CA batch ca_wave1_kaon — verify WA polish for K003 K004 K005 K006 K013
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/kaon/K003.yaml          | 14 +++++++++++---
 flavor_catalog/processes/kaon/K004.yaml          | 14 +++++++++++---
 flavor_catalog/processes/kaon/K005.yaml          | 13 +++++++++++--
 flavor_catalog/processes/kaon/K006.yaml          | 14 +++++++++++---
 flavor_catalog/processes/kaon/K013.yaml          | 14 +++++++++++---
 flavor_catalog/worklogs/checker/ca_wave1_kaon.md | 61 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 6 files changed, 116 insertions(+), 14 deletions(-)
```

### 108. `333042f` — flavor-catalog(top_higgs_ew): CA batch ca_wave1_top_higgs_ew — verify WA polish for T001 T010
- SHA: `333042f736a2d1aca5b46f65817e5b559405bba6`
- Message: flavor-catalog(top_higgs_ew): CA batch ca_wave1_top_higgs_ew — verify WA polish for T001 T010
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T001.yaml          | 14 ++++++++++++--
 flavor_catalog/processes/top_higgs_ew/T010.yaml          | 14 ++++++++++++--
 flavor_catalog/worklogs/checker/ca_wave1_top_higgs_ew.md | 39 +++++++++++++++++++++++++++++++++++++++
 3 files changed, 63 insertions(+), 4 deletions(-)
```

### 109. `e5ba712` — flavor-catalog(kaon): WA v2 rework wa_wave1_kaon_v2 — address CA findings
- SHA: `e5ba712bf86b05f63a6ccbc6d52e7f62467519ff`
- Message: flavor-catalog(kaon): WA v2 rework wa_wave1_kaon_v2 — address CA findings
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/kaon/K005.tex             | 10 +++++-----
 flavor_catalog/processes/kaon/K005.yaml            | 15 ++++++++++++---
 flavor_catalog/worklogs/writer/wa_wave1_kaon_v2.md | 31 +++++++++++++++++++++++++++++++
 3 files changed, 48 insertions(+), 8 deletions(-)
```

### 110. `fb9a170` — flavor-catalog(top_higgs_ew): WA v2 rework wa_wave1_top_higgs_ew_v2 — address CA findings
- SHA: `fb9a170efde052b630e9ef396a8799172fbf2e3d`
- Message: flavor-catalog(top_higgs_ew): WA v2 rework wa_wave1_top_higgs_ew_v2 — address CA findings
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T001.yaml            | 134 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-----------------------------------------------------------
 flavor_catalog/processes/top_higgs_ew/T010.yaml            |  12 ++++++++++--
 flavor_catalog/worklogs/writer/wa_wave1_top_higgs_ew_v2.md |  34 ++++++++++++++++++++++++++++++++++
 3 files changed, 119 insertions(+), 61 deletions(-)
```

### 111. `25afeb0` — flavor-catalog(beauty): WA v2 rework wa_wave1_beauty_v2 — address CA findings
- SHA: `25afeb0db219c1a833e6c7143070adf09344ba2a`
- Message: flavor-catalog(beauty): WA v2 rework wa_wave1_beauty_v2 — address CA findings
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B009.yaml            | 21 ++++++++++++++++++++-
 flavor_catalog/processes/beauty/B015.yaml            | 55 ++++++++++++++++++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/worklogs/writer/wa_wave1_beauty_v2.md | 52 ++++++++++++++++++++++++++++++++++++++++++++++++++++
 3 files changed, 126 insertions(+), 2 deletions(-)
```

### 112. `f902c33` — flavor-catalog(charged_lepton): PKA draft for L007 tau -> mu gamma
- SHA: `f902c33b9a55342191a6dd8e1e11c387242b9515`
- Message: flavor-catalog(charged_lepton): PKA draft for L007 tau -> mu gamma
- Physical/numerical summary: Added an initial catalog process draft for L007 tau -> mu gamma; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charged_lepton/L007.tex                           |  83 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/charged_lepton/L007.yaml                          | 136 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/L007/babar2010_tau_lgamma_arxiv0908_2381.txt     |  29 +++++++++++++++++++++++++++++
 flavor_catalog/references/L007/belle2021_tau_lgamma_arxiv2103_12994.txt    |  28 ++++++++++++++++++++++++++++
 flavor_catalog/references/L007/belleii_physics_book_1808_10567_tau_lfv.txt |  31 +++++++++++++++++++++++++++++++
 flavor_catalog/references/L007/cfw2008_arxiv0804_1954.txt                  |  26 ++++++++++++++++++++++++++
 flavor_catalog/references/L007/pdg2025_tau_mu_gamma_pdglive.txt            |  27 +++++++++++++++++++++++++++
 flavor_catalog/references/L007/perez_randall2008_arxiv0805_4652.txt        |  26 ++++++++++++++++++++++++++
 flavor_catalog/references/L007/source_manifest.yaml                        |  64 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/L007.md                                        |  31 +++++++++++++++++++++++++++++++
 10 files changed, 481 insertions(+)
```

### 113. `56a0a7f` — flavor-catalog(top_higgs_ew): PKA draft for T007 t to c h
- SHA: `56a0a7f0c390942df1f0b7e2722e318c5efe0fc3`
- Message: flavor-catalog(top_higgs_ew): PKA draft for T007 t to c h
- Physical/numerical summary: Added an initial catalog process draft for T007 t to c h; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T007.tex                             |  92 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/top_higgs_ew/T007.yaml                            | 122 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/T007/aguilar_saavedra_hepph0409342_arxiv_abs.txt |  18 ++++++++++++++++++
 flavor_catalog/references/T007/atlas_2404_02123_arxiv_abs.txt              |  28 ++++++++++++++++++++++++++++
 flavor_catalog/references/T007/cfw_0804_1954_arxiv_abs.txt                 |  18 ++++++++++++++++++
 flavor_catalog/references/T007/cms_2407_15172_arxiv_abs.txt                |  25 +++++++++++++++++++++++++
 flavor_catalog/references/T007/pdg2026_top_hc_datablock.txt                |  40 ++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/T007/source_manifest.yaml                        |  49 +++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/T007.md                                        |  33 +++++++++++++++++++++++++++++++++
 9 files changed, 425 insertions(+)
```

### 114. `282fb9c` — flavor-catalog(beauty): PKA draft for B032 Bbar to pi Kbar
- SHA: `282fb9cac4028b710a91f0b341b8363e1197ae31`
- Message: flavor-catalog(beauty): PKA draft for B032 Bbar to pi Kbar
- Physical/numerical summary: Added an initial catalog process draft for B032 Bbar to pi Kbar; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B032.tex                               | 105 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/beauty/B032.yaml                              | 130 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B032/belleii2023_b0_ks_pi0_cp_arxiv.txt      |  18 ++++++++++++++++++
 flavor_catalog/references/B032/belleii2024_btokpi_arxiv.txt            |  19 +++++++++++++++++++
 flavor_catalog/references/B032/cfw2008_rs_flavor_arxiv.txt             |  20 ++++++++++++++++++++
 flavor_catalog/references/B032/hflav_dec2025_btopik_values.txt         |  42 ++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B032/lhcb2021_bplus_kplus_pi0_acp_arxiv.txt  |  20 ++++++++++++++++++++
 flavor_catalog/references/B032/pdg2025_k0pi0_time_dependent_cp_api.txt |  42 ++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B032/source_manifest.yaml                    |  60 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/B032.md                                    |  33 +++++++++++++++++++++++++++++++++
 10 files changed, 489 insertions(+)
```

### 115. `7b2c2c0` — flavor-catalog(charged_lepton): PKA draft for L002 mu to 3e
- SHA: `7b2c2c0f4322606294ef52c15613277f19d98e1e`
- Message: flavor-catalog(charged_lepton): PKA draft for L002 mu to 3e
- Physical/numerical summary: Added an initial catalog process draft for L002 mu to 3e; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charged_lepton/L002.tex                          |  86 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/charged_lepton/L002.yaml                         | 103 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/L002/agashe2006_rs_lfv_arxiv_hepph_0606021.txt  |  28 ++++++++++++++++++++++++++++
 flavor_catalog/references/L002/cfw2008_arxiv0804_1954.txt                 |  28 ++++++++++++++++++++++++++++
 flavor_catalog/references/L002/crivellin2017_mu_e_eft_arxiv1702_03020.txt |  28 ++++++++++++++++++++++++++++
 flavor_catalog/references/L002/mu3e_status_2025_arxiv2501_14667.txt       |  27 +++++++++++++++++++++++++++
 flavor_catalog/references/L002/mu3e_tdr_arxiv2009_11690.txt               |  25 +++++++++++++++++++++++++
 flavor_catalog/references/L002/pdg2026_muon_listing_s004r4.txt            |  38 ++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/L002/sindrum1988_inspire_abstract.txt           |  27 +++++++++++++++++++++++++++
 flavor_catalog/references/L002/source_manifest.yaml                       |  75 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/L002.md                                       |  35 +++++++++++++++++++++++++++++++++++
 11 files changed, 500 insertions(+)
```

### 116. `93c5a85` — flavor-catalog(beauty): PKA draft for B025 R_D
- SHA: `93c5a8513aa115c92e2347e0694ad8818fe94e89`
- Message: flavor-catalog(beauty): PKA draft for B025 R_D
- Physical/numerical summary: Added an initial catalog process draft for B025 R_D; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B025.tex                       |   92 +++++
 flavor_catalog/processes/beauty/B025.yaml                      |  200 +++++++++
 flavor_catalog/references/B025/babar2012_arxiv1205_5442.txt    |   28 ++
 flavor_catalog/references/B025/babar2013_arxiv1303_0571.txt    |   34 ++
 flavor_catalog/references/B025/belle2015_arxiv1507_03233.txt   |  637 +++++++++++++++++++++++++++++
 flavor_catalog/references/B025/belle2019_arxiv1910_05864.txt   |  598 +++++++++++++++++++++++++++
 flavor_catalog/references/B025/belleii2025_arxiv2504_11220.txt | 1383 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B025/cfw2008_arxiv0804_1954.txt      |   34 ++
 flavor_catalog/references/B025/flag2024_arxiv2411_04268.txt    |  118 ++++++
 flavor_catalog/references/B025/hflav_ckm2025_logfile.txt       | 2303 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B025/hflav_ckm2025_rdrds.txt         |  524 ++++++++++++++++++++++++
 flavor_catalog/references/B025/lhcb2023_arxiv2302_02886.txt    |   28 ++
 flavor_catalog/references/B025/lhcb2025_arxiv2406_03387.txt    | 3293 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B025/sha256sums.txt                  |   11 +
 flavor_catalog/references/B025/source_manifest.yaml            |  114 ++++++
 flavor_catalog/worklogs/pka/B025.md                            |   33 ++
 16 files changed, 9430 insertions(+)
```

### 117. `69c33e4` — flavor-catalog(kaon): CA batch ca_wave1_kaon_v2 — verify WA polish for K005
- SHA: `69c33e4dd5432ae33040cb00cbac97fec632019d`
- Message: flavor-catalog(kaon): CA batch ca_wave1_kaon_v2 — verify WA polish for K005
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/kaon/K005.yaml             | 12 ++++++++++--
 flavor_catalog/worklogs/checker/ca_wave1_kaon_v2.md | 32 ++++++++++++++++++++++++++++++++
 2 files changed, 42 insertions(+), 2 deletions(-)
```

### 118. `27de127` — flavor-catalog(top_higgs_ew): PKA draft for T018 h to mu tau
- SHA: `27de127f7ec1d5848a9071930b0014fe14d971d2`
- Message: flavor-catalog(top_higgs_ew): PKA draft for T018 h to mu tau
- Physical/numerical summary: Added an initial catalog process draft for T018 h to mu tau; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T018.tex                          |  90 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/top_higgs_ew/T018.yaml                         | 158 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/T018/atlas2023_higg201911_arxiv_abs.txt       |  32 +++++++++++++++++++++++++++++
 flavor_catalog/references/T018/cms2015_hig14005_arxiv_abs.txt           |  32 +++++++++++++++++++++++++++++
 flavor_catalog/references/T018/cms2021_hig20009_arxiv_abs.txt           |  27 ++++++++++++++++++++++++
 flavor_catalog/references/T018/csaki_falkowski_weiler2008_arxiv_abs.txt |  23 +++++++++++++++++++++
 flavor_catalog/references/T018/harnik_kopp_zupan2012_arxiv_abs.txt      |  23 +++++++++++++++++++++
 flavor_catalog/references/T018/pdg2025_higgs_lfv_review.txt             |  25 ++++++++++++++++++++++
 flavor_catalog/references/T018/perez_randall2008_arxiv_abs.txt          |  21 +++++++++++++++++++
 flavor_catalog/references/T018/source_manifest.yaml                     |  67 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/T018.md                                     |  39 +++++++++++++++++++++++++++++++++++
 11 files changed, 537 insertions(+)
```

### 119. `494309d` — flavor-catalog(beauty): PKA draft for B033 B to phi K_S
- SHA: `494309d7c902651a4b6805f9cbe266b3839e5345`
- Message: flavor-catalog(beauty): PKA draft for B033 B to phi K_S
- Physical/numerical summary: Added an initial catalog process draft for B033 B to phi K_S; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B033.tex                             |  93 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/beauty/B033.yaml                            | 170 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B033/babar_2012_b0_kkks_arxiv.txt          |  30 ++++++++++++++++++++++++++
 flavor_catalog/references/B033/belle_2010_b0_kkks_arxiv.txt          |  33 ++++++++++++++++++++++++++++
 flavor_catalog/references/B033/belleii_2023_phik0_arxiv.txt          |  31 ++++++++++++++++++++++++++
 flavor_catalog/references/B033/csaki_falkowski_weiler_2008_arxiv.txt |  28 ++++++++++++++++++++++++
 flavor_catalog/references/B033/hflav_triangle_summer2025_phiK0.txt   |  50 ++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B033/source_manifest.yaml                  |  49 +++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/B033.md                                  |  40 ++++++++++++++++++++++++++++++++++
 9 files changed, 524 insertions(+)
```

### 120. `00393cd` — flavor-catalog(top_higgs_ew): CA batch ca_wave1_top_higgs_ew_v2 — verify WA polish for T001 T010
- SHA: `00393cddfaacf14f6349e3f4ab63198bb5faae58`
- Message: flavor-catalog(top_higgs_ew): CA batch ca_wave1_top_higgs_ew_v2 — verify WA polish for T001 T010
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T001.yaml             | 12 ++++++++++--
 flavor_catalog/processes/top_higgs_ew/T010.yaml             | 12 ++++++++++--
 flavor_catalog/worklogs/checker/ca_wave1_top_higgs_ew_v2.md | 41 +++++++++++++++++++++++++++++++++++++++++
 3 files changed, 61 insertions(+), 4 deletions(-)
```

### 121. `6ffcdb2` — flavor-catalog(beauty): CA batch ca_wave1_beauty_v2 — verify WA polish for B009 B015
- SHA: `6ffcdb2780f8f39330d3ee84ee53c294ef9652dc`
- Message: flavor-catalog(beauty): CA batch ca_wave1_beauty_v2 — verify WA polish for B009 B015
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B009.yaml             | 13 +++++++++++--
 flavor_catalog/processes/beauty/B015.yaml             | 13 +++++++++++--
 flavor_catalog/worklogs/checker/ca_wave1_beauty_v2.md | 50 ++++++++++++++++++++++++++++++++++++++++++++++++++
 3 files changed, 72 insertions(+), 4 deletions(-)
```

### 122. `f93ed80` — flavor-catalog(beauty): PKA draft for B018 R_K
- SHA: `f93ed807b49c552f4d23795b4fce48cf44f73c20`
- Message: flavor-catalog(beauty): PKA draft for B018 R_K
- Physical/numerical summary: Added an initial catalog process draft for B018 R_K; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B018.tex                                              |   88 ++++
 flavor_catalog/processes/beauty/B018.yaml                                             |  106 +++++
 flavor_catalog/references/B018/bordone_isidori_pattori_2016_rk_sm_arxiv1605_07633.txt |   34 ++
 flavor_catalog/references/B018/cfw_0804_1954_arxiv.txt                                |   34 ++
 flavor_catalog/references/B018/hflav_dec2025_rk_centralq2.txt                         |   23 +
 flavor_catalog/references/B018/hflav_dec2025_rk_fullq2_belle.txt                      |   21 +
 flavor_catalog/references/B018/hflav_dec2025_rk_lowq2.txt                             |   21 +
 flavor_catalog/references/B018/lhcb_2021_rk_arxiv2105_10303.txt                       |   25 +
 flavor_catalog/references/B018/lhcb_2023_rk_rkstar_arxiv2212_09153.txt                | 3217 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B018/source_manifest.yaml                                   |   73 +++
 flavor_catalog/worklogs/pka/B018.md                                                   |   31 ++
 11 files changed, 3673 insertions(+)
```

### 123. `a2c7802` — flavor-catalog(top_higgs_ew): WA batch wa_w23_top_higgs_ew — polish PKA drafts for T002 T007 T018
- SHA: `a2c7802058da410cfdc916fb41fe8292ece64328`
- Message: flavor-catalog(top_higgs_ew): WA batch wa_w23_top_higgs_ew — polish PKA drafts for T002 T007 T018
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T002.tex        | 56 ++++++++++++++++++++++++++++++++------------------------
 flavor_catalog/processes/top_higgs_ew/T002.yaml       | 10 +++++++++-
 flavor_catalog/processes/top_higgs_ew/T007.tex        | 30 +++++++++++++++---------------
 flavor_catalog/processes/top_higgs_ew/T007.yaml       | 10 +++++++++-
 flavor_catalog/processes/top_higgs_ew/T018.tex        | 28 +++++++++++++++-------------
 flavor_catalog/processes/top_higgs_ew/T018.yaml       | 10 +++++++++-
 flavor_catalog/worklogs/writer/wa_w23_top_higgs_ew.md | 62 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 7 files changed, 151 insertions(+), 55 deletions(-)
```

### 124. `19a8a23` — flavor-catalog(charged_lepton): WA batch wa_w23_charged_lepton — polish PKA drafts for L001 L002 L007
- SHA: `19a8a23b80c7a78324445546020c1c15dc1681f4`
- Message: flavor-catalog(charged_lepton): WA batch wa_w23_charged_lepton — polish PKA drafts for L001 L002 L007
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charged_lepton/L001.tex        | 78 ++++++++++++++++++++++++++++++++++++++++++++----------------------------------
 flavor_catalog/processes/charged_lepton/L001.yaml       | 12 ++++++++++--
 flavor_catalog/processes/charged_lepton/L002.tex        | 67 +++++++++++++++++++++++++++++++++++++------------------------------
 flavor_catalog/processes/charged_lepton/L002.yaml       | 12 ++++++++++--
 flavor_catalog/processes/charged_lepton/L007.tex        | 53 +++++++++++++++++++++++++++++------------------------
 flavor_catalog/processes/charged_lepton/L007.yaml       | 12 ++++++++++--
 flavor_catalog/worklogs/writer/wa_w23_charged_lepton.md | 42 ++++++++++++++++++++++++++++++++++++++++++
 7 files changed, 182 insertions(+), 94 deletions(-)
```

### 125. `0fe9e31` — flavor-catalog(mixed): WA batch wa_w23_kaon_charm_edm — polish PKA drafts for K001 K002 C001 E001
- SHA: `0fe9e312636c07a15b07a1b84d7ead49cf7d25bd`
- Message: flavor-catalog(mixed): WA batch wa_w23_kaon_charm_edm — polish PKA drafts for K001 K002 C001 E001
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charm/C001.tex                 | 33 ++++++++++++++++++---------------
 flavor_catalog/processes/charm/C001.yaml                | 10 +++++++++-
 flavor_catalog/processes/edm_neutrino/E001.tex          | 32 +++++++++++++++++++-------------
 flavor_catalog/processes/edm_neutrino/E001.yaml         | 10 +++++++++-
 flavor_catalog/processes/kaon/K001.tex                  | 57 +++++++++++++++++++++++++++++----------------------------
 flavor_catalog/processes/kaon/K001.yaml                 | 10 +++++++++-
 flavor_catalog/processes/kaon/K002.tex                  | 47 +++++++++++++++++++++++++----------------------
 flavor_catalog/processes/kaon/K002.yaml                 | 10 +++++++++-
 flavor_catalog/worklogs/writer/wa_w23_kaon_charm_edm.md | 80 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 9 files changed, 207 insertions(+), 82 deletions(-)
```

### 126. `5031e06` — flavor-catalog(beauty): WA batch wa_w23_beauty — polish PKA drafts for B002 B005 B017 B018 B025 B032 B033
- SHA: `5031e06ec0a937cfb617ea992afb4342346240f2`
- Message: flavor-catalog(beauty): WA batch wa_w23_beauty — polish PKA drafts for B002 B005 B017 B018 B025 B032 B033
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B002.tex        |  39 ++++++++++++++++++++++-----------------
 flavor_catalog/processes/beauty/B002.yaml       |  10 +++++++++-
 flavor_catalog/processes/beauty/B005.tex        |  19 +++++++++++--------
 flavor_catalog/processes/beauty/B005.yaml       |  10 +++++++++-
 flavor_catalog/processes/beauty/B017.tex        |  95 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/beauty/B017.yaml       | 132 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/beauty/B018.tex        |  25 +++++++++++++++----------
 flavor_catalog/processes/beauty/B018.yaml       |  10 +++++++++-
 flavor_catalog/processes/beauty/B025.tex        |   7 ++++---
 flavor_catalog/processes/beauty/B025.yaml       |  10 +++++++++-
 flavor_catalog/processes/beauty/B032.tex        |  29 +++++++++++++++++------------
 flavor_catalog/processes/beauty/B032.yaml       |  10 +++++++++-
 flavor_catalog/processes/beauty/B033.tex        |  24 +++++++++++++++---------
 flavor_catalog/processes/beauty/B033.yaml       |  10 +++++++++-
 flavor_catalog/worklogs/writer/wa_w23_beauty.md | 101 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 15 files changed, 466 insertions(+), 65 deletions(-)
```

### 127. `8795419` — flavor-catalog(top_higgs_ew): CA batch ca_w23_top_higgs_ew — verify WA polish for T002 T007 T018
- SHA: `8795419b6d64fa09a0af2ff7f6fe558784215ee6`
- Message: flavor-catalog(top_higgs_ew): CA batch ca_w23_top_higgs_ew — verify WA polish for T002 T007 T018
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T002.yaml        | 11 ++++++++++-
 flavor_catalog/processes/top_higgs_ew/T007.yaml        | 14 +++++++++++---
 flavor_catalog/processes/top_higgs_ew/T018.yaml        | 11 ++++++++++-
 flavor_catalog/worklogs/checker/ca_w23_top_higgs_ew.md | 49 +++++++++++++++++++++++++++++++++++++++++++++++++
 4 files changed, 80 insertions(+), 5 deletions(-)
```

### 128. `76b1fd6` — flavor-catalog(charged_lepton): CA batch ca_w23_charged_lepton — verify WA polish for L001 L002 L007
- SHA: `76b1fd6879cd3acc79b01d5f75e10426bfbc31e2`
- Message: flavor-catalog(charged_lepton): CA batch ca_w23_charged_lepton — verify WA polish for L001 L002 L007
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charged_lepton/L001.yaml        | 16 ++++++++++++++--
 flavor_catalog/processes/charged_lepton/L002.yaml        | 16 ++++++++++++++--
 flavor_catalog/processes/charged_lepton/L007.yaml        | 15 ++++++++++++---
 flavor_catalog/worklogs/checker/ca_w23_charged_lepton.md | 42 ++++++++++++++++++++++++++++++++++++++++++
 4 files changed, 82 insertions(+), 7 deletions(-)
```

### 129. `bcb34b3` — flavor-catalog(mixed): CA batch ca_w23_kaon_charm_edm — verify WA polish for K001 K002 C001 E001
- SHA: `bcb34b39ae040ecd49c317820fbe1729912f5992`
- Message: flavor-catalog(mixed): CA batch ca_w23_kaon_charm_edm — verify WA polish for K001 K002 C001 E001
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charm/C001.yaml                 | 15 +++++++++++++--
 flavor_catalog/processes/edm_neutrino/E001.yaml          | 15 +++++++++++++--
 flavor_catalog/processes/kaon/K001.yaml                  | 15 ++++++++++++---
 flavor_catalog/processes/kaon/K002.yaml                  | 15 +++++++++++++--
 flavor_catalog/worklogs/checker/ca_w23_kaon_charm_edm.md | 49 +++++++++++++++++++++++++++++++++++++++++++++++++
 5 files changed, 100 insertions(+), 9 deletions(-)
```

### 130. `8b6388c` — flavor-catalog(beauty): CA batch ca_w23_beauty — verify WA polish for B002 B005 B017 B018 B025 B032 B033
- SHA: `8b6388cd45cf4e46891b8473c8f53423fa1e27c0`
- Message: flavor-catalog(beauty): CA batch ca_w23_beauty — verify WA polish for B002 B005 B017 B018 B025 B032 B033
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B002.yaml        | 14 +++++++++++---
 flavor_catalog/processes/beauty/B005.yaml        | 14 +++++++++++---
 flavor_catalog/processes/beauty/B017.yaml        | 14 ++++++++++++--
 flavor_catalog/processes/beauty/B018.yaml        | 14 ++++++++++++--
 flavor_catalog/processes/beauty/B025.yaml        | 14 +++++++++++---
 flavor_catalog/processes/beauty/B032.yaml        | 14 ++++++++++++--
 flavor_catalog/processes/beauty/B033.yaml        | 14 +++++++++++---
 flavor_catalog/worklogs/checker/ca_w23_beauty.md | 69 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 8 files changed, 149 insertions(+), 18 deletions(-)
```

### 131. `df697bc` — flavor-catalog: DA-1 discovery worklog round_001_full_scope
- SHA: `df697bce50886bbf8e9ff53413fbc763a8d0fa45`
- Message: flavor-catalog: DA-1 discovery worklog round_001_full_scope
- Physical/numerical summary: Applied the change described by the subject; based on stat surface, this is documentation/catalog/orchestration work rather than a new physics-code calculation.
- Files touched (`git show --stat`):
```text
 flavor_catalog/worklogs/discovery/round_001_full_scope.md | 66 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 66 insertions(+)
```

### 132. `a45dfed` — flavor-catalog(beauty): WA v2 rework wa_w23_beauty_v2 — address CA findings
- SHA: `a45dfed7768554c050c822c672cdf6dc56355103`
- Message: flavor-catalog(beauty): WA v2 rework wa_w23_beauty_v2 — address CA findings
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B017.yaml                          | 11 ++++++++++-
 flavor_catalog/processes/beauty/B018.yaml                          | 32 +++++++++++++++++++++++++++++++-
 flavor_catalog/processes/beauty/B032.yaml                          | 34 +++++++++++++++++++++++++++++++++-
 flavor_catalog/references/B017/cfw_0804_1954_rs_flavor.txt         | 21 +++++++++++++++++++++
 flavor_catalog/references/B017/hflav_dec2025_B0_to_Kst0_ll.txt     | 18 ++++++++++++++++++
 flavor_catalog/references/B017/hflav_dec2025_B_to_Xsll.txt         | 21 +++++++++++++++++++++
 flavor_catalog/references/B017/hflav_dec2025_Bplus_to_Kplus_ll.txt | 18 ++++++++++++++++++
 flavor_catalog/references/B017/hflav_dec2025_overview.txt          | 19 +++++++++++++++++++
 flavor_catalog/references/B017/lhcb_2003_04831_Kst_angular.txt     | 25 +++++++++++++++++++++++++
 flavor_catalog/references/B017/lhcb_2212_09153_RK_RKst.txt         | 37 +++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B017/source_manifest.yaml                | 75 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B017/wen_xu_2305_19038_global_fit.txt    | 26 ++++++++++++++++++++++++++
 flavor_catalog/worklogs/writer/wa_w23_beauty_v2.md                 | 34 ++++++++++++++++++++++++++++++++++
 13 files changed, 368 insertions(+), 3 deletions(-)
```

### 133. `37d83ca` — flavor-catalog(kaon_charm_edm): WA v2 rework wa_w23_kaon_charm_edm_v2 — address CA findings
- SHA: `37d83cafa790de4b753f35d3bda76e489f40fffb`
- Message: flavor-catalog(kaon_charm_edm): WA v2 rework wa_w23_kaon_charm_edm_v2 — address CA findings
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charm/C001.yaml                   | 40 +++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/processes/edm_neutrino/E001.yaml            | 47 ++++++++++++++++++++++++++++++++++++++++++++---
 flavor_catalog/processes/kaon/K002.yaml                    | 99 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-------------------------------------
 flavor_catalog/worklogs/writer/wa_w23_kaon_charm_edm_v2.md | 55 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
 4 files changed, 200 insertions(+), 41 deletions(-)
```

### 134. `48bffd7` — flavor-catalog(charged_lepton): WA v2 rework wa_w23_charged_lepton_v2 — address CA findings
- SHA: `48bffd72676f4fe32b2d5d8b0cf5eefe36db0c24`
- Message: flavor-catalog(charged_lepton): WA v2 rework wa_w23_charged_lepton_v2 — address CA findings
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charged_lepton/L001.yaml          | 75 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-----
 flavor_catalog/processes/charged_lepton/L002.yaml          | 56 +++++++++++++++++++++++++++++++++++++++++++++++++++-----
 flavor_catalog/worklogs/writer/wa_w23_charged_lepton_v2.md | 42 ++++++++++++++++++++++++++++++++++++++++++
 3 files changed, 163 insertions(+), 10 deletions(-)
```

### 135. `7bf7c21` — flavor-catalog(top_higgs_ew): WA v2 rework wa_w23_top_higgs_ew_v2 — address CA findings
- SHA: `7bf7c2199d97d78613cd1898d5d9cc61eb25daa5`
- Message: flavor-catalog(top_higgs_ew): WA v2 rework wa_w23_top_higgs_ew_v2 — address CA findings
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T002.yaml          | 142 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++----------------------------------------------
 flavor_catalog/processes/top_higgs_ew/T018.yaml          | 200 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-------------------------------------------------------------
 flavor_catalog/worklogs/writer/wa_w23_top_higgs_ew_v2.md |  46 ++++++++++++++++++++++++++++++++++++
 3 files changed, 251 insertions(+), 137 deletions(-)
```

### 136. `03e2dba` — flavor-catalog(top_higgs_ew): CA batch ca_w23_top_higgs_ew_v2 — verify WA polish for T002 T018
- SHA: `03e2dba1d1481a07877b84ca1429e4e36ee4f2c0`
- Message: flavor-catalog(top_higgs_ew): CA batch ca_w23_top_higgs_ew_v2 — verify WA polish for T002 T018
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T002.yaml           | 15 ++++++++++++---
 flavor_catalog/processes/top_higgs_ew/T018.yaml           | 15 ++++++++++++---
 flavor_catalog/worklogs/checker/ca_w23_top_higgs_ew_v2.md | 41 +++++++++++++++++++++++++++++++++++++++++
 3 files changed, 65 insertions(+), 6 deletions(-)
```

### 137. `96446a1` — flavor-catalog(charged_lepton): CA batch ca_w23_charged_lepton_v2 — verify WA polish for L001 L002
- SHA: `96446a1a91bafbbceb106c80546167b1cab9f2ef`
- Message: flavor-catalog(charged_lepton): CA batch ca_w23_charged_lepton_v2 — verify WA polish for L001 L002
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charged_lepton/L001.yaml           | 13 ++++++++++++-
 flavor_catalog/processes/charged_lepton/L002.yaml           | 13 +++++++++++--
 flavor_catalog/worklogs/checker/ca_w23_charged_lepton_v2.md | 38 ++++++++++++++++++++++++++++++++++++++
 3 files changed, 61 insertions(+), 3 deletions(-)
```

### 138. `30465ce` — flavor-catalog(kaon_charm_edm): CA batch ca_w23_kaon_charm_edm_v2 — verify WA polish for K002 C001 E001
- SHA: `30465cee99557595f957946259ae234be3960dfb`
- Message: flavor-catalog(kaon_charm_edm): CA batch ca_w23_kaon_charm_edm_v2 — verify WA polish for K002 C001 E001
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charm/C001.yaml                    | 13 +++++++++++--
 flavor_catalog/processes/edm_neutrino/E001.yaml             | 13 +++++++++++--
 flavor_catalog/processes/kaon/K002.yaml                     | 13 +++++++++++--
 flavor_catalog/worklogs/checker/ca_w23_kaon_charm_edm_v2.md | 43 +++++++++++++++++++++++++++++++++++++++++++
 4 files changed, 76 insertions(+), 6 deletions(-)
```

### 139. `3d8c18d` — flavor-catalog(beauty): CA batch ca_w23_beauty_v2 — verify WA polish for B017 B018 B032
- SHA: `3d8c18df2ad80495d8746dbdc2682671e96d7a3b`
- Message: flavor-catalog(beauty): CA batch ca_w23_beauty_v2 — verify WA polish for B017 B018 B032
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B017.yaml           | 12 ++++++++++--
 flavor_catalog/processes/beauty/B018.yaml           | 12 ++++++++++--
 flavor_catalog/processes/beauty/B032.yaml           | 12 +++++++++++-
 flavor_catalog/worklogs/checker/ca_w23_beauty_v2.md | 45 +++++++++++++++++++++++++++++++++++++++++++++
 4 files changed, 76 insertions(+), 5 deletions(-)
```

### 140. `f048312` — flavor-catalog(beauty): PKA draft for B004 phi_s
- SHA: `f048312c8e9a063c9dbb15948af61eaadc6882bd`
- Message: flavor-catalog(beauty): PKA draft for B004 phi_s
- Physical/numerical summary: Added an initial catalog process draft for B004 phi_s; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B004.tex                       | 101 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/beauty/B004.yaml                      | 106 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B004/cfw_2008_rs_flavor.txt          |  22 ++++++++++++++++++++++
 flavor_catalog/references/B004/hflav_pdg2025_phis_inputs.txt   |  32 ++++++++++++++++++++++++++++++++
 flavor_catalog/references/B004/lhcb_2024_jpsikk_phis_arxiv.txt |  45 +++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B004/sha256sums.txt                  |   3 +++
 flavor_catalog/references/B004/source_manifest.yaml            |  37 +++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/B004.md                            |  31 +++++++++++++++++++++++++++++++
 8 files changed, 377 insertions(+)
```

### 141. `7d3da08` — flavor-catalog(top_higgs_ew): PKA draft for EW001 oblique parameters
- SHA: `7d3da085e8a0da5dd283e6991a916d9ee36bb8f6`
- Message: flavor-catalog(top_higgs_ew): PKA draft for EW001 oblique parameters
- Physical/numerical summary: Added an initial catalog process draft for EW001 oblique parameters; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/edm_neutrino/E004.tex                              |  94 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/edm_neutrino/E004.yaml                             | 142 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/top_higgs_ew/EW001.tex                             |  87 +++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/top_higgs_ew/EW001.yaml                            | 225 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/E004/abel2020_arxiv2001_11966.txt                 |  26 ++++++++++++++++
 flavor_catalog/references/E004/bhattacharya2022_arxiv2203_03746.txt         |  27 +++++++++++++++++
 flavor_catalog/references/E004/cfw2008_arxiv0804_1954.txt                   |  27 +++++++++++++++++
 flavor_catalog/references/E004/koenig_neubert_straub2014_arxiv1403_2756.txt |  27 +++++++++++++++++
 flavor_catalog/references/E004/pdg2026_neutron_edm_datablock.txt            |  35 ++++++++++++++++++++++
 flavor_catalog/references/E004/source_manifest.yaml                         |  57 +++++++++++++++++++++++++++++++++++
 flavor_catalog/references/EW001/cfw_2008_rs_flavor.txt                      |  19 ++++++++++++
 flavor_catalog/references/EW001/gfitter_oblique_parameters.txt              |  26 ++++++++++++++++
 flavor_catalog/references/EW001/hepfit_deblas_2022_cdf_stu.txt              |  45 ++++++++++++++++++++++++++++
 flavor_catalog/references/EW001/hepfit_precision_ew_page.txt                |  25 ++++++++++++++++
 flavor_catalog/references/EW001/pdg_2025_electroweak_stu.txt                |  43 ++++++++++++++++++++++++++
 flavor_catalog/references/EW001/peskin_takeuchi_1992_oblique.txt            |  14 +++++++++
 flavor_catalog/references/EW001/source_manifest.yaml                        |  58 +++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/E004.md                                         |  34 +++++++++++++++++++++
 flavor_catalog/worklogs/pka/EW001.md                                        |  36 ++++++++++++++++++++++
 19 files changed, 1047 insertions(+)
```

### 142. `6a00c06` — flavor-catalog(charged_lepton): PKA draft for L003 mu-e conversion in Al
- SHA: `6a00c063d5662bc00ad4509c3fdae1bb144a5c6f`
- Message: flavor-catalog(charged_lepton): PKA draft for L003 mu-e conversion in Al
- Physical/numerical summary: Added an initial catalog process draft for L003 mu-e conversion in Al; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charged_lepton/L003.tex                       |  87 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/charged_lepton/L003.yaml                      | 154 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/L003/cfw2008_arxiv0804_1954.txt              |  21 ++++++++++++++++++++
 flavor_catalog/references/L003/comet_phase_i_tdr_arxiv1812_09018.txt   |  26 ++++++++++++++++++++++++
 flavor_catalog/references/L003/comet_strategy_2018_arxiv1812_07824.txt |  28 ++++++++++++++++++++++++++
 flavor_catalog/references/L003/mu2e_2016_arxiv1606_05559.txt           |  26 ++++++++++++++++++++++++
 flavor_catalog/references/L003/mu2e_run1_2022_arxiv2210_11380.txt      |  27 +++++++++++++++++++++++++
 flavor_catalog/references/L003/mu2eii_2022_arxiv2203_07569.txt         |  25 +++++++++++++++++++++++
 flavor_catalog/references/L003/sindrumii_2006_gold_juser.txt           |  22 +++++++++++++++++++++
 flavor_catalog/references/L003/source_manifest.yaml                    |  75 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/L003.md                                    |  33 +++++++++++++++++++++++++++++++
 11 files changed, 524 insertions(+)
```

### 143. `cfbee57` — flavor-catalog(charged_lepton): PKA draft for L008 tau e gamma
- SHA: `cfbee57f42f1c6e029bf8671192a19e025360832`
- Message: flavor-catalog(charged_lepton): PKA draft for L008 tau e gamma
- Physical/numerical summary: Added an initial catalog process draft for L008 tau e gamma; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charged_lepton/L008.tex                           |  91 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/charged_lepton/L008.yaml                          | 131 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/L008/babar2010_tau_lgamma_arxiv0908_2381.txt     |  29 +++++++++++++++++++++++++++++
 flavor_catalog/references/L008/belle2021_tau_lgamma_arxiv2103_12994.txt    |  27 +++++++++++++++++++++++++++
 flavor_catalog/references/L008/belleii_physics_book_1808_10567_tau_lfv.txt |  31 +++++++++++++++++++++++++++++++
 flavor_catalog/references/L008/cfw2008_arxiv0804_1954.txt                  |  26 ++++++++++++++++++++++++++
 flavor_catalog/references/L008/pdg2025_tau_e_gamma_pdglive.txt             |  35 +++++++++++++++++++++++++++++++++++
 flavor_catalog/references/L008/perez_randall2008_arxiv0805_4652.txt        |  26 ++++++++++++++++++++++++++
 flavor_catalog/references/L008/source_manifest.yaml                        |  64 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/L008.md                                        |  32 ++++++++++++++++++++++++++++++++
 10 files changed, 492 insertions(+)
```

### 144. `d5f2a59` — flavor-catalog(charm): PKA draft for C003 Delta A_CP
- SHA: `d5f2a592adb98e6c1b9705a3f0e96e6811d2e183`
- Message: flavor-catalog(charm): PKA draft for C003 Delta A_CP
- Physical/numerical summary: Added an initial catalog process draft for C003 Delta A_CP; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charm/C003.tex                          |  83 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/charm/C003.yaml                         | 155 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/C003/cfw2008_arxiv0804_1954.txt        |  24 +++++++++++++++++++++++
 flavor_catalog/references/C003/hflav2025_direct_indirect_cpv.txt |  36 ++++++++++++++++++++++++++++++++++
 flavor_catalog/references/C003/hflav_ckm25_dcpv_inputs.txt       |  37 +++++++++++++++++++++++++++++++++++
 flavor_catalog/references/C003/lhcb2019_arxiv1903_08726.txt      |  29 ++++++++++++++++++++++++++++
 flavor_catalog/references/C003/sha256sums.txt                    |   4 ++++
 flavor_catalog/references/C003/source_manifest.yaml              |  46 ++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/C003.md                              |  36 ++++++++++++++++++++++++++++++++++
 9 files changed, 450 insertions(+)
```

### 145. `e0b7a8f` — flavor-catalog(charm): PKA draft for C004 D0 to mu+ mu-
- SHA: `e0b7a8f1154a36c2a7875fb4a220b5943c0ca5d4`
- Message: flavor-catalog(charm): PKA draft for C004 D0 to mu+ mu-
- Physical/numerical summary: Added an initial catalog process draft for C004 D0 to mu+ mu-; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charm/C004.tex                                       |  97 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/charm/C004.yaml                                      | 127 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/C004/burdman2002_arxiv_hepph0112235.txt             |  29 +++++++++++++++++++++++++++++
 flavor_catalog/references/C004/cfw2008_arxiv0804_1954.txt                     |  20 ++++++++++++++++++++
 flavor_catalog/references/C004/cms2025_d0_mumu_public.txt                     |  22 ++++++++++++++++++++++
 flavor_catalog/references/C004/gisbert_hiller_suelmann2024_rare_charm_eft.txt |  20 ++++++++++++++++++++
 flavor_catalog/references/C004/lhcb2023_arxiv2212_11203.txt                   |  23 +++++++++++++++++++++++
 flavor_catalog/references/C004/pdg2026_d0_mumu_pdgLive.txt                    |  39 +++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/C004/source_manifest.yaml                           |  71 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/C004.md                                           |  35 +++++++++++++++++++++++++++++++++++
 10 files changed, 483 insertions(+)
```

### 146. `36602e4` — flavor-catalog(charged_lepton): PKA draft for L009 tau to 3 mu
- SHA: `36602e43ad4eec46bd9bea385f09cd8e9febc3c6`
- Message: flavor-catalog(charged_lepton): PKA draft for L009 tau to 3 mu
- Physical/numerical summary: Added an initial catalog process draft for L009 tau to 3 mu; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charged_lepton/L009.tex                 |  82 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/charged_lepton/L009.yaml                | 145 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/L009/babar_2010_arxiv_1002_4550.txt    |  29 +++++++++++++++++++++++++++++
 flavor_catalog/references/L009/belle_2010_arxiv_1001_3221.txt    |  34 ++++++++++++++++++++++++++++++++++
 flavor_catalog/references/L009/belleii_2024_arxiv_2405_07386.txt |  29 +++++++++++++++++++++++++++++
 flavor_catalog/references/L009/cfw_2008_arxiv_0804_1954.txt      |  23 +++++++++++++++++++++++
 flavor_catalog/references/L009/lhcb_2015_arxiv_1409_8548.txt     |  25 +++++++++++++++++++++++++
 flavor_catalog/references/L009/lhcb_2026_arxiv_2601_20785.txt    |  33 +++++++++++++++++++++++++++++++++
 flavor_catalog/references/L009/pdg_2025_tau_listing.txt          |  41 +++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/L009/source_manifest.yaml              |  74 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/L009.md                              |  31 +++++++++++++++++++++++++++++++
 11 files changed, 546 insertions(+)
```

### 147. `e5f45ab` — flavor-catalog(charged_lepton): WA v2 rework wa_w23_charged_lepton_v3 — address CA findings
- SHA: `e5f45ab0714943e3986e80771adee4ee29a5e18a`
- Message: flavor-catalog(charged_lepton): WA v2 rework wa_w23_charged_lepton_v3 — address CA findings
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charged_lepton/L001.tex           |  4 ++--
 flavor_catalog/processes/charged_lepton/L001.yaml          | 12 +++++++++++-
 flavor_catalog/worklogs/writer/wa_w23_charged_lepton_v3.md | 25 +++++++++++++++++++++++++
 3 files changed, 38 insertions(+), 3 deletions(-)
```

### 148. `4ec4e95` — flavor-catalog(beauty): WA v2 rework wa_w23_beauty_v3 — address CA findings
- SHA: `4ec4e958229d9de98f1ea9a913b8d9073e4f7785`
- Message: flavor-catalog(beauty): WA v2 rework wa_w23_beauty_v3 — address CA findings
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B032.yaml          | 53 ++++++++++++++++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/worklogs/writer/wa_w23_beauty_v3.md | 40 ++++++++++++++++++++++++++++++++++++++++
 2 files changed, 92 insertions(+), 1 deletion(-)
```

### 149. `741f7c0` — flavor-catalog(edm_neutrino): PKA draft for E006 mercury EDM
- SHA: `741f7c0e2438571201953598aa87662025ddc3f5`
- Message: flavor-catalog(edm_neutrino): PKA draft for E006 mercury EDM
- Physical/numerical summary: Added an initial catalog process draft for E006 mercury EDM; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/edm_neutrino/E006.tex                     | 101 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/edm_neutrino/E006.yaml                    | 135 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/E006/cfw2008_arxiv0804_1954.txt          |  22 ++++++++++++++++++++++
 flavor_catalog/references/E006/graner2016_arxiv1601_04339.txt      |  31 +++++++++++++++++++++++++++++++
 flavor_catalog/references/E006/pdg2026_neutron_edm_hg_crossref.txt |  31 +++++++++++++++++++++++++++++++
 flavor_catalog/references/E006/sahoo2017_arxiv1612_09371.txt       |  28 ++++++++++++++++++++++++++++
 flavor_catalog/references/E006/source_manifest.yaml                |  46 ++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/E006.md                                |  35 +++++++++++++++++++++++++++++++++++
 8 files changed, 429 insertions(+)
```

### 150. `de78fa1` — flavor-catalog(beauty): PKA draft for B006 B_d to mu mu
- SHA: `de78fa12749a602b3b1a910dc2101c35283a5ce7`
- Message: flavor-catalog(beauty): PKA draft for B006 B_d to mu mu
- Physical/numerical summary: Added an initial catalog process draft for B006 B_d to mu mu; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B006.tex                      |  88 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/beauty/B006.yaml                     | 107 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B006/atlas_1812_03017_bdmumu.txt    |  21 +++++++++++++++++++++
 flavor_catalog/references/B006/bobeth_1311_0903_bsdll_sm.txt  |  18 ++++++++++++++++++
 flavor_catalog/references/B006/cfw_0804_1954_rs_flavor.txt    |  16 ++++++++++++++++
 flavor_catalog/references/B006/cms_2212_10311_bdmumu.txt      |  23 +++++++++++++++++++++++
 flavor_catalog/references/B006/hflav_2023_bd_over_bs_mumu.txt |  21 +++++++++++++++++++++
 flavor_catalog/references/B006/lhcb_2108_09283_bdmumu.txt     |  22 ++++++++++++++++++++++
 flavor_catalog/references/B006/pdg_2026_bdmumu.txt            |  39 +++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B006/source_manifest.yaml           |  76 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/B006.md                           |  35 +++++++++++++++++++++++++++++++++++
 11 files changed, 466 insertions(+)
```

### 151. `f734e4e` — flavor-catalog(top_higgs_ew): PKA draft for EW002 first-row CKM unitarity
- SHA: `f734e4e01e4b5f64cc78a0f736f2479207edd1dd`
- Message: flavor-catalog(top_higgs_ew): PKA draft for EW002 first-row CKM unitarity
- Physical/numerical summary: Added an initial catalog process draft for EW002 first-row CKM unitarity; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/EW002.tex                                       |  93 +++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/top_higgs_ew/EW002.yaml                                      | 231 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/EW002/bryman_cirigliano_crivellin_inguglia_2022_lfu_ckm.txt |  24 +++++++++++++
 flavor_catalog/references/EW002/cfw_2008_rs_flavor.txt                                |  25 ++++++++++++++
 flavor_catalog/references/EW002/crivellin_kirk_kitahara_mescia_2022_vlq_caa.txt       |  28 ++++++++++++++++
 flavor_catalog/references/EW002/flag_2024_vus_ckm_unitarity.txt                       |  34 +++++++++++++++++++
 flavor_catalog/references/EW002/hardy_towner_2020_superallowed_vud.txt                |  23 +++++++++++++
 flavor_catalog/references/EW002/pdg_2025_vud_vus_ckm_unitarity.txt                    |  29 ++++++++++++++++
 flavor_catalog/references/EW002/source_manifest.yaml                                  |  58 ++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/EW002.md                                                  |  33 ++++++++++++++++++
 10 files changed, 578 insertions(+)
```

### 152. `fd49402` — flavor-catalog(beauty): PKA draft for B019 R_Kstar
- SHA: `fd49402dbfd1dc4585462aeef85ff81059616150`
- Message: flavor-catalog(beauty): PKA draft for B019 R_Kstar
- Physical/numerical summary: Added an initial catalog process draft for B019 R_Kstar; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B019.tex                                                  |  94 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/beauty/B019.yaml                                                 | 144 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B019/bordone_isidori_pattori_2016_rkstar_sm_arxiv1605_07633.txt |  25 ++++++++++++++++++++++
 flavor_catalog/references/B019/cfw_0804_1954_rs_flavor.txt                                |  15 +++++++++++++
 flavor_catalog/references/B019/hflav_dec2025_rkstar_centralq2.txt                         |  23 ++++++++++++++++++++
 flavor_catalog/references/B019/hflav_dec2025_rkstar_lowq2_lhcb.txt                        |  21 ++++++++++++++++++
 flavor_catalog/references/B019/lhcb_2017_rkstar_arxiv1705_05802.txt                       |  20 +++++++++++++++++
 flavor_catalog/references/B019/lhcb_2023_rk_rkstar_arxiv2212_09153.txt                    |  32 +++++++++++++++++++++++++++
 flavor_catalog/references/B019/source_manifest.yaml                                       |  63 +++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/B019.md                                                       |  29 +++++++++++++++++++++++++
 10 files changed, 466 insertions(+)
```

### 153. `cf21e2b` — flavor-catalog(charm): PKA draft for C007 D+ to pi+ mu+ mu-
- SHA: `cf21e2b0ff3986cd0a6f081b976fad98f7f9ddfa`
- Message: flavor-catalog(charm): PKA draft for C007 D+ to pi+ mu+ mu-
- Physical/numerical summary: Added an initial catalog process draft for C007 D+ to pi+ mu+ mu-; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charm/C007.tex                              |  95 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/charm/C007.yaml                             | 115 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/C007/cfw2008_arxiv0804_1954.txt            |  25 +++++++++++++++++++++++++
 flavor_catalog/references/C007/deboer_hiller2015_arxiv1510_00311.txt |  24 ++++++++++++++++++++++++
 flavor_catalog/references/C007/lhcb2013_arxiv1304_6365.txt           |  35 +++++++++++++++++++++++++++++++++++
 flavor_catalog/references/C007/lhcb2021_arxiv2011_00217.txt          |  45 +++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/C007/pdg2026_dplus_piplus_mumu_pdgLive.txt |  55 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/C007/sha256sums.txt                        |   5 +++++
 flavor_catalog/references/C007/source_manifest.yaml                  |  60 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/C007.md                                  |  35 +++++++++++++++++++++++++++++++++++
 10 files changed, 494 insertions(+)
```

### 154. `519aec5` — flavor-catalog(charged_lepton): PKA draft for L004 mu-e conversion in Au
- SHA: `519aec532e3ebf4dd1414d13350ab0d681529e19`
- Message: flavor-catalog(charged_lepton): PKA draft for L004 mu-e conversion in Au
- Physical/numerical summary: Added an initial catalog process draft for L004 mu-e conversion in Au; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charged_lepton/L004.tex                           |  82 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/charged_lepton/L004.yaml                          | 126 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/L004/cfw2008_arxiv0804_1954.txt                  |  21 +++++++++++++++++++++
 flavor_catalog/references/L004/comet_phase_i_tdr_arxiv1812_09018.txt       |  27 +++++++++++++++++++++++++++
 flavor_catalog/references/L004/mu2e_tdr_arxiv1501_05241.txt                |  25 +++++++++++++++++++++++++
 flavor_catalog/references/L004/pdg2026_muon_listing_s004_au_conversion.txt |  25 +++++++++++++++++++++++++
 flavor_catalog/references/L004/sha256sums.txt                              |   5 +++++
 flavor_catalog/references/L004/sindrumii_2006_gold_crossref.txt            |  24 ++++++++++++++++++++++++
 flavor_catalog/references/L004/source_manifest.yaml                        |  59 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/L004.md                                        |  33 +++++++++++++++++++++++++++++++++
 10 files changed, 427 insertions(+)
```

### 155. `8eec2a0` — flavor-catalog(edm_neutrino): PKA draft for E008 quark chromo-EDM bounds
- SHA: `8eec2a0ebab676d7b768e8c4549ebff13afa46ec`
- Message: flavor-catalog(edm_neutrino): PKA draft for E008 quark chromo-EDM bounds
- Physical/numerical summary: Added an initial catalog process draft for E008 quark chromo-EDM bounds; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/edm_neutrino/E008.tex                                                 |  99 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/edm_neutrino/E008.yaml                                                | 189 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/E008/bhattacharya2024_isovector_qcedm_lattice_arxiv2404_13516.txt    |  25 ++++++++++++++++
 flavor_catalog/references/E008/cfw2008_arxiv0804_1954.txt                                      |  24 +++++++++++++++
 flavor_catalog/references/E008/graner2016_hg199_arxiv1601_04339.txt                            |  24 +++++++++++++++
 flavor_catalog/references/E008/koenig_neubert_straub2014_composite_dipoles_arxiv1403_2756.txt  |  38 ++++++++++++++++++++++++
 flavor_catalog/references/E008/olive_pospelov_ritz_santoso2005_hg_qcedm_arxiv_hepph0506106.txt |  34 +++++++++++++++++++++
 flavor_catalog/references/E008/pdg2026_neutron_edm_datablock.txt                               |  30 +++++++++++++++++++
 flavor_catalog/references/E008/pospelov_ritz2000_qcedm_neutron_arxiv_hepph0010037.txt          |  25 ++++++++++++++++
 flavor_catalog/references/E008/source_manifest.yaml                                            |  79 +++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/E008.md                                                            |  34 +++++++++++++++++++++
 11 files changed, 601 insertions(+)
```

### 156. `d64d457` — flavor-catalog(charm): PKA draft for C002 D mixing CP
- SHA: `d64d4576a5326fe9385d23fd1580f72f62043739`
- Message: flavor-catalog(charm): PKA draft for C002 D mixing CP
- Physical/numerical summary: Added an initial catalog process draft for C002 D mixing CP; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charm/C002.tex                              |  86 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/charm/C002.yaml                             | 121 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/C002/belle_belleii2025_arxiv2410_22961.txt |  23 +++++++++++++++++++++++
 flavor_catalog/references/C002/cfw2008_arxiv0804_1954.txt            |  22 ++++++++++++++++++++++
 flavor_catalog/references/C002/hflav_ckm25_dmixing_cpv_fit.txt       |  35 +++++++++++++++++++++++++++++++++++
 flavor_catalog/references/C002/lhcb2023_arxiv2208_06512.txt          |  24 ++++++++++++++++++++++++
 flavor_catalog/references/C002/pdg2025_dmixing_review_cpv.txt        |  24 ++++++++++++++++++++++++
 flavor_catalog/references/C002/source_manifest.yaml                  |  54 ++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/C002.md                                  |  33 +++++++++++++++++++++++++++++++++
 9 files changed, 422 insertions(+)
```

### 157. `66ea613` — flavor-catalog(charged_lepton): PKA draft for L010 tau to 3e
- SHA: `66ea613306d70f303bd3b43d555e235568d849b5`
- Message: flavor-catalog(charged_lepton): PKA draft for L010 tau to 3e
- Physical/numerical summary: Added an initial catalog process draft for L010 tau to 3e; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charged_lepton/L010.tex                           |  81 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/charged_lepton/L010.yaml                          | 154 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/L010/agashe2006_rs_lfv_arxiv_hepph0606021.txt    |  28 +++++++++++++++++++++++++
 flavor_catalog/references/L010/babar2010_tau3leptons_arxiv1002_4550.txt    |  33 ++++++++++++++++++++++++++++++
 flavor_catalog/references/L010/belle2010_tau3leptons_arxiv1001_3221.txt    |  36 ++++++++++++++++++++++++++++++++
 flavor_catalog/references/L010/belleii_physics_book_1808_10567_tau_lfv.txt |  30 +++++++++++++++++++++++++++
 flavor_catalog/references/L010/cfw2008_arxiv0804_1954.txt                  |  27 ++++++++++++++++++++++++
 flavor_catalog/references/L010/lhcb2026_tau3mu_context_arxiv2601_20785.txt |  29 ++++++++++++++++++++++++++
 flavor_catalog/references/L010/pdg2025_tau3e_pdglive_s035r58.txt           |  39 +++++++++++++++++++++++++++++++++++
 flavor_catalog/references/L010/source_manifest.yaml                        |  74 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/L010.md                                        |  35 +++++++++++++++++++++++++++++++
 11 files changed, 566 insertions(+)
```

### 158. `fbbc3ff` — flavor-catalog(kaon): PKA draft for K017 K to l nu RK
- SHA: `fbbc3fff046d831e36ebee760968f12fdf2c81e0`
- Message: flavor-catalog(kaon): PKA draft for K017 K to l nu RK
- Physical/numerical summary: Added an initial catalog process draft for K017 K to l nu RK; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/kaon/K017.tex                                   |  90 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/kaon/K017.yaml                                  | 113 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/K017/cfw_2008_arxiv0804_1954.txt               |  15 +++++++++++++++
 flavor_catalog/references/K017/cirigliano_rosell_2007_arxiv0707_4464.txt |  17 +++++++++++++++++
 flavor_catalog/references/K017/flavianet_2010_arxiv1005_2323.txt         |  20 ++++++++++++++++++++
 flavor_catalog/references/K017/kloe_2009_arxiv0907_3594.txt              |  18 ++++++++++++++++++
 flavor_catalog/references/K017/na62_2013_arxiv1212_4012.txt              |  20 ++++++++++++++++++++
 flavor_catalog/references/K017/pdg2025_rk_api_20260516.txt               |  57 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/K017/source_manifest.yaml                      |  69 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/K017.md                                      |  40 ++++++++++++++++++++++++++++++++++++++++
 10 files changed, 459 insertions(+)
```

### 159. `a91f599` — flavor-catalog(top_higgs_ew): PKA draft for EW003 Vcb Vub tensions
- SHA: `a91f599f800700a81adc55dd001d2721f835b07a`
- Message: flavor-catalog(top_higgs_ew): PKA draft for EW003 Vcb Vub tensions
- Physical/numerical summary: Added an initial catalog process draft for EW003 Vcb Vub tensions; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/EW003.tex                      |  88 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/top_higgs_ew/EW003.yaml                     | 179 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/EW003/bhatta_2021_bu_tau_np.txt            |  21 +++++++++++++++++
 flavor_catalog/references/EW003/cfw_2008_rs_flavor_baseline.txt      |  10 ++++++++
 flavor_catalog/references/EW003/flag_2024_b_semileptonic_lattice.txt |  17 ++++++++++++++
 flavor_catalog/references/EW003/hflav_2023_semileptonic_report.txt   |  20 ++++++++++++++++
 flavor_catalog/references/EW003/iguro_2024_bc_tau_global_fit.txt     |  15 ++++++++++++
 flavor_catalog/references/EW003/pdg_2024_vcb_vub.txt                 |  18 +++++++++++++++
 flavor_catalog/references/EW003/source_manifest.yaml                 |  58 +++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/EW003.md                                 |  37 ++++++++++++++++++++++++++++++
 10 files changed, 463 insertions(+)
```

### 160. `0c395f2` — flavor-catalog(beauty): PKA draft for B026 R_Dstar
- SHA: `0c395f20cfc970c4150cd2506d89539bf1fb40fd`
- Message: flavor-catalog(beauty): PKA draft for B026 R_Dstar
- Physical/numerical summary: Added an initial catalog process draft for B026 R_Dstar; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B026.tex                       |  105 +++++++
 flavor_catalog/processes/beauty/B026.yaml                      |  192 +++++++++++++
 flavor_catalog/references/B026/babar2012_arxiv1205_5442.txt    |   13 +
 flavor_catalog/references/B026/babar2013_arxiv1303_0571.txt    |   13 +
 flavor_catalog/references/B026/belle2015_arxiv1507_03233.txt   |   13 +
 flavor_catalog/references/B026/belle2017_arxiv1612_00529.txt   |   13 +
 flavor_catalog/references/B026/belle2018_arxiv1709_00129.txt   |   13 +
 flavor_catalog/references/B026/belle2019_arxiv1910_05864.txt   |   13 +
 flavor_catalog/references/B026/belleii2024_arxiv2401_02840.txt |   12 +
 flavor_catalog/references/B026/belleii2025_arxiv2504_11220.txt |   13 +
 flavor_catalog/references/B026/cfw2008_arxiv0804_1954.txt      |   13 +
 flavor_catalog/references/B026/flag2024_arxiv2411_04268.txt    |   13 +
 flavor_catalog/references/B026/hflav_ckm2025_logfile.txt       | 2303 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B026/hflav_ckm2025_rdrds.txt         |  524 ++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B026/lhcb2023_arxiv2302_02886.txt    |   28 ++
 flavor_catalog/references/B026/lhcb2023_arxiv2305_01463.txt    |   13 +
 flavor_catalog/references/B026/lhcb2025_arxiv2406_03387.txt    |   13 +
 flavor_catalog/references/B026/sha256sums.txt                  |   15 +
 flavor_catalog/references/B026/source_manifest.yaml            |  157 +++++++++++
 flavor_catalog/worklogs/pka/B026.md                            |   34 +++
 20 files changed, 3513 insertions(+)
```

### 161. `49c57c2` — flavor-catalog(beauty): CA batch ca_w23_beauty_v3 — verify WA polish for B032
- SHA: `49c57c2b652a7970bba807e1d1b06d5eeb5a0632`
- Message: flavor-catalog(beauty): CA batch ca_w23_beauty_v3 — verify WA polish for B032
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B032.yaml           | 13 +++++++++++--
 flavor_catalog/worklogs/checker/ca_w23_beauty_v3.md | 38 ++++++++++++++++++++++++++++++++++++++
 2 files changed, 49 insertions(+), 2 deletions(-)
```

### 162. `65e3ca5` — flavor-catalog(charged_lepton): CA batch ca_w23_charged_lepton_v3 — verify WA polish for L001
- SHA: `65e3ca5d1715eff7ea72a87da5818926db61a7e5`
- Message: flavor-catalog(charged_lepton): CA batch ca_w23_charged_lepton_v3 — verify WA polish for L001
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charged_lepton/L001.yaml           | 14 +++++++++++++-
 flavor_catalog/worklogs/checker/ca_w23_charged_lepton_v3.md | 33 +++++++++++++++++++++++++++++++++
 2 files changed, 46 insertions(+), 1 deletion(-)
```

### 163. `1268def` — flavor-catalog(beauty): WA batch wa_w4_beauty — polish PKA drafts for B004 B006 B019 B026
- SHA: `1268deff8504452762c83499c2dedfa8009ab242`
- Message: flavor-catalog(beauty): WA batch wa_w4_beauty — polish PKA drafts for B004 B006 B019 B026
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B004.tex       | 41 ++++++++++++++++++++++-------------------
 flavor_catalog/processes/beauty/B004.yaml      | 10 +++++++++-
 flavor_catalog/processes/beauty/B006.tex       | 40 +++++++++++++++++++++++-----------------
 flavor_catalog/processes/beauty/B006.yaml      | 10 +++++++++-
 flavor_catalog/processes/beauty/B019.tex       | 40 +++++++++++++++++++++++-----------------
 flavor_catalog/processes/beauty/B019.yaml      | 10 +++++++++-
 flavor_catalog/processes/beauty/B026.tex       | 66 +++++++++++++++++++++++++++++++++++-------------------------------
 flavor_catalog/processes/beauty/B026.yaml      | 10 +++++++++-
 flavor_catalog/worklogs/writer/wa_w4_beauty.md | 80 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 9 files changed, 219 insertions(+), 88 deletions(-)
```

### 164. `7741e1c` — flavor-catalog(kaon_edm): WA batch wa_w4_kaon_edm — polish PKA drafts for K017 E004 E006 E008
- SHA: `7741e1ceb1ed568cd29f6ecb359b30efe71f3f4b`
- Message: flavor-catalog(kaon_edm): WA batch wa_w4_kaon_edm — polish PKA drafts for K017 E004 E006 E008
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/edm_neutrino/E004.tex   | 26 ++++++++++++++------------
 flavor_catalog/processes/edm_neutrino/E004.yaml  | 10 +++++++++-
 flavor_catalog/processes/edm_neutrino/E006.tex   | 26 ++++++++++++--------------
 flavor_catalog/processes/edm_neutrino/E006.yaml  | 10 +++++++++-
 flavor_catalog/processes/edm_neutrino/E008.tex   | 35 ++++++++++++++++++++++-------------
 flavor_catalog/processes/edm_neutrino/E008.yaml  | 10 +++++++++-
 flavor_catalog/processes/kaon/K017.tex           | 52 +++++++++++++++++++++++++++++-----------------------
 flavor_catalog/processes/kaon/K017.yaml          | 10 +++++++++-
 flavor_catalog/worklogs/writer/wa_w4_kaon_edm.md | 84 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 9 files changed, 197 insertions(+), 66 deletions(-)
```

### 165. `eae8d05` — flavor-catalog(charged_lepton): WA batch wa_w4_charged_lepton — polish PKA drafts for L003 L004 L008 L009 L010
- SHA: `eae8d05b880a49e134bd7295c0359b1021454fa8`
- Message: flavor-catalog(charged_lepton): WA batch wa_w4_charged_lepton — polish PKA drafts for L003 L004 L008 L009 L010
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charged_lepton/L003.tex       | 71 ++++++++++++++++++++++++++++++++++++++---------------------------------
 flavor_catalog/processes/charged_lepton/L003.yaml      | 11 ++++++++++-
 flavor_catalog/processes/charged_lepton/L004.tex       | 31 ++++++++++++++++++++-----------
 flavor_catalog/processes/charged_lepton/L004.yaml      | 11 ++++++++++-
 flavor_catalog/processes/charged_lepton/L008.tex       | 47 ++++++++++++++++++++++++++++-------------------
 flavor_catalog/processes/charged_lepton/L008.yaml      | 11 ++++++++++-
 flavor_catalog/processes/charged_lepton/L009.tex       | 61 +++++++++++++++++++++++++++++++++++++------------------------
 flavor_catalog/processes/charged_lepton/L009.yaml      | 11 ++++++++++-
 flavor_catalog/processes/charged_lepton/L010.tex       | 55 +++++++++++++++++++++++++++++++++----------------------
 flavor_catalog/processes/charged_lepton/L010.yaml      | 11 ++++++++++-
 flavor_catalog/worklogs/writer/wa_w4_charged_lepton.md | 84 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 11 files changed, 290 insertions(+), 114 deletions(-)
```

### 166. `470d498` — flavor-catalog(charged_lepton): CA batch ca_w4_charged_lepton — verify WA polish for L003 L004 L008 L009 L010
- SHA: `470d4985d6f07a0f4854fd269e904b106e35c4d8`
- Message: flavor-catalog(charged_lepton): CA batch ca_w4_charged_lepton — verify WA polish for L003 L004 L008 L009 L010
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charged_lepton/L003.yaml       | 15 ++++++++++++---
 flavor_catalog/processes/charged_lepton/L004.yaml       | 15 ++++++++++++---
 flavor_catalog/processes/charged_lepton/L008.yaml       | 15 ++++++++++++---
 flavor_catalog/processes/charged_lepton/L009.yaml       | 15 ++++++++++++---
 flavor_catalog/processes/charged_lepton/L010.yaml       | 15 ++++++++++++---
 flavor_catalog/worklogs/checker/ca_w4_charged_lepton.md | 59 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 6 files changed, 119 insertions(+), 15 deletions(-)
```

### 167. `0f5972b` — flavor-catalog(top_higgs_ew): CA batch ca_w4_ew — verify WA polish for EW001 EW002 EW003
- SHA: `0f5972bf9f3f7209f07315e63abd9aa1aed81865`
- Message: flavor-catalog(top_higgs_ew): CA batch ca_w4_ew — verify WA polish for EW001 EW002 EW003
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/EW001.yaml | 12 +++++++++++-
 flavor_catalog/processes/top_higgs_ew/EW002.yaml | 12 +++++++++++-
 flavor_catalog/processes/top_higgs_ew/EW003.yaml | 13 ++++++++++++-
 flavor_catalog/worklogs/checker/ca_w4_ew.md      | 54 ++++++++++++++++++++++++++++++++++++++++++++++++++++++
 4 files changed, 88 insertions(+), 3 deletions(-)
```

### 168. `bcd8907` — flavor-catalog(arbitration): Opus signoff on L001 -- APPROVE-OVERRIDE
- SHA: `bcd89070c9c8e08b4931235f483725d55b54e1fb`
- Message: flavor-catalog(arbitration): Opus signoff on L001 -- APPROVE-OVERRIDE
- Physical/numerical summary: Recorded reviewer/sign-off disposition; review ledger only.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charged_lepton/L001.yaml | 14 ++++++++++++--
 flavor_catalog/processes/edm_neutrino/E004.yaml   | 12 ++++++++++--
 flavor_catalog/processes/edm_neutrino/E006.yaml   | 12 +++++++++++-
 flavor_catalog/processes/edm_neutrino/E008.yaml   | 12 ++++++++++--
 flavor_catalog/processes/kaon/K017.yaml           | 12 +++++++++++-
 flavor_catalog/signoff/by_process/L001.md         | 19 +++++++++++++++++++
 flavor_catalog/worklogs/checker/ca_w4_kaon_edm.md | 46 ++++++++++++++++++++++++++++++++++++++++++++++
 7 files changed, 119 insertions(+), 8 deletions(-)
```

### 169. `aab1862` — flavor-catalog(mixed): CA batch ca_w4_kaon_edm — verify WA polish for K017 E004 E006 E008
- SHA: `aab1862846f6f3627a0234d151640217571f79ac`
- Message: flavor-catalog(mixed): CA batch ca_w4_kaon_edm — verify WA polish for K017 E004 E006 E008
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
  - `git show --stat` reported no file-level changes.

### 170. `a5e94cb` — flavor-catalog(charm): CA batch ca_w4_charm — verify WA polish for C002 C003 C004 C007
- SHA: `a5e94cb23411f0faf201b0e96a29170156503d9d`
- Message: flavor-catalog(charm): CA batch ca_w4_charm — verify WA polish for C002 C003 C004 C007
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charm/C002.yaml       | 11 ++++++++++-
 flavor_catalog/processes/charm/C003.yaml       | 11 ++++++++++-
 flavor_catalog/processes/charm/C004.yaml       | 12 +++++++++++-
 flavor_catalog/processes/charm/C007.yaml       | 13 ++++++++++++-
 flavor_catalog/worklogs/checker/ca_w4_charm.md | 52 ++++++++++++++++++++++++++++++++++++++++++++++++++++
 5 files changed, 95 insertions(+), 4 deletions(-)
```

### 171. `c369747` — flavor-catalog(beauty): CA batch ca_w4_beauty — verify WA polish for B004 B006 B019 B026
- SHA: `c3697476ba5313b1c470765b78eb297975d3e52f`
- Message: flavor-catalog(beauty): CA batch ca_w4_beauty — verify WA polish for B004 B006 B019 B026
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B004.yaml       | 12 ++++++++++--
 flavor_catalog/processes/beauty/B006.yaml       | 11 ++++++++++-
 flavor_catalog/processes/beauty/B019.yaml       | 12 ++++++++++--
 flavor_catalog/processes/beauty/B026.yaml       | 12 ++++++++++--
 flavor_catalog/worklogs/checker/ca_w4_beauty.md | 45 +++++++++++++++++++++++++++++++++++++++++++++
 5 files changed, 85 insertions(+), 7 deletions(-)
```

### 172. `4907544` — flavor-catalog: DA-2 discovery worklog round_002_followup
- SHA: `4907544270d8ba4ed1a4ad6e7c0bae0e3536c4f6`
- Message: flavor-catalog: DA-2 discovery worklog round_002_followup
- Physical/numerical summary: Applied the change described by the subject; based on stat surface, this is documentation/catalog/orchestration work rather than a new physics-code calculation.
- Files touched (`git show --stat`):
```text
 flavor_catalog/worklogs/discovery/round_002_followup.md | 51 +++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 51 insertions(+)
```

### 173. `a23b929` — flavor-catalog(kaon_edm): WA v2 rework wa_w4_kaon_edm_v2 — address CA findings
- SHA: `a23b929430fe730b80d07443e6839a215e6276e4`
- Message: flavor-catalog(kaon_edm): WA v2 rework wa_w4_kaon_edm_v2 — address CA findings
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/edm_neutrino/E006.yaml     | 26 ++++++++++++++++++++++++--
 flavor_catalog/processes/kaon/K017.yaml             | 17 +++++++++++++++--
 flavor_catalog/worklogs/writer/wa_w4_kaon_edm_v2.md | 41 +++++++++++++++++++++++++++++++++++++++++
 3 files changed, 80 insertions(+), 4 deletions(-)
```

### 174. `ff366bc` — flavor-catalog(beauty): WA v2 rework wa_w4_beauty_v2 — address CA findings
- SHA: `ff366bc875a35f1e76991916a86f605c9bd09389`
- Message: flavor-catalog(beauty): WA v2 rework wa_w4_beauty_v2 — address CA findings
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B006.yaml         | 35 ++++++++++++++++++++++++++++++++++-
 flavor_catalog/worklogs/writer/wa_w4_beauty_v2.md | 36 ++++++++++++++++++++++++++++++++++++
 2 files changed, 70 insertions(+), 1 deletion(-)
```

### 175. `52acd5e` — flavor-catalog(top_higgs_ew): WA v2 rework wa_w4_ew_v2 — address CA findings
- SHA: `52acd5e5ccb529b6b679063a98c6430da44f21a0`
- Message: flavor-catalog(top_higgs_ew): WA v2 rework wa_w4_ew_v2 — address CA findings
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/EW001.yaml | 12 +++++++++++-
 flavor_catalog/processes/top_higgs_ew/EW002.yaml | 12 +++++++++++-
 flavor_catalog/processes/top_higgs_ew/EW003.yaml | 34 +++++++++++++++++++++++++++++++---
 flavor_catalog/worklogs/writer/wa_w4_ew_v2.md    | 30 ++++++++++++++++++++++++++++++
 4 files changed, 83 insertions(+), 5 deletions(-)
```

### 176. `f863013` — flavor-catalog(charm): WA v2 rework wa_w4_charm_v2 — address CA findings
- SHA: `f863013c38b245829e86fa512d8ef24d3ec7d637`
- Message: flavor-catalog(charm): WA v2 rework wa_w4_charm_v2 — address CA findings
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charm/C002.yaml         | 11 ++++++++++-
 flavor_catalog/processes/charm/C003.yaml         | 11 ++++++++++-
 flavor_catalog/processes/charm/C004.tex          |  9 ++++-----
 flavor_catalog/processes/charm/C004.yaml         | 17 +++++++++++++----
 flavor_catalog/processes/charm/C007.yaml         | 38 +++++++++++++++++++++++++++++++++++++-
 flavor_catalog/worklogs/writer/wa_w4_charm_v2.md | 29 +++++++++++++++++++++++++++++
 6 files changed, 103 insertions(+), 12 deletions(-)
```

### 177. `1397db5` — flavor-catalog(charged_lepton): PKA draft for L005 mu-e conversion in Ti
- SHA: `1397db50cde7bb8050acd1bd7c6ad05dc81ff68f`
- Message: flavor-catalog(charged_lepton): PKA draft for L005 mu-e conversion in Ti
- Physical/numerical summary: Added an initial catalog process draft for L005 mu-e conversion in Ti; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charged_lepton/L005.tex                           |  79 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/charged_lepton/L005.yaml                          | 130 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/L005/cfw2008_arxiv0804_1954.txt                  |  21 +++++++++++++++++++++
 flavor_catalog/references/L005/comet_phase_i_tdr_arxiv1812_09018.txt       |  26 ++++++++++++++++++++++++++
 flavor_catalog/references/L005/mu2e_tdr_arxiv1501_05241.txt                |  26 ++++++++++++++++++++++++++
 flavor_catalog/references/L005/pdg2026_muon_listing_s004_ti_conversion.txt |  31 +++++++++++++++++++++++++++++++
 flavor_catalog/references/L005/sha256sums.txt                              |   5 +++++
 flavor_catalog/references/L005/sindrumii_1993_titanium_crossref.txt        |  25 +++++++++++++++++++++++++
 flavor_catalog/references/L005/source_manifest.yaml                        |  59 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/L005.md                                        |  32 ++++++++++++++++++++++++++++++++
 10 files changed, 434 insertions(+)
```

### 178. `e1ff427` — flavor-catalog(top_higgs_ew): PKA draft for T019 h to e tau
- SHA: `e1ff4275297d3c2b47070553a5c43c804ad37c7e`
- Message: flavor-catalog(top_higgs_ew): PKA draft for T019 h to e tau
- Physical/numerical summary: Added an initial catalog process draft for T019 h to e tau; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T019.tex                          |  85 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/top_higgs_ew/T019.yaml                         | 181 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/T019/atlas2019_higg201906_arxiv_abs.txt       |  25 ++++++++++++++++++++
 flavor_catalog/references/T019/atlas2023_higg201911_arxiv_abs.txt       |  32 +++++++++++++++++++++++++
 flavor_catalog/references/T019/cms2021_hig20009_arxiv_abs.txt           |  25 ++++++++++++++++++++
 flavor_catalog/references/T019/csaki_falkowski_weiler2008_arxiv_abs.txt |  20 ++++++++++++++++
 flavor_catalog/references/T019/harnik_kopp_zupan2012_arxiv_abs.txt      |  19 +++++++++++++++
 flavor_catalog/references/T019/pdg2025_higgs_lfv_review.txt             |  27 +++++++++++++++++++++
 flavor_catalog/references/T019/perez_randall2008_arxiv_abs.txt          |  19 +++++++++++++++
 flavor_catalog/references/T019/source_manifest.yaml                     |  68 +++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/T019.md                                     |  39 ++++++++++++++++++++++++++++++
 11 files changed, 540 insertions(+)
```

### 179. `08bd183` — flavor-catalog(charm): PKA draft for C005 D0 to e+e-
- SHA: `08bd183151f7b92413c37fe928e0410638036ee3`
- Message: flavor-catalog(charm): PKA draft for C005 D0 to e+e-
- Physical/numerical summary: Added an initial catalog process draft for C005 D0 to e+e-; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charm/C005.tex                           |  94 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/charm/C005.yaml                          | 109 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/C005/babar2012_arxiv1206_5419.txt       |  27 +++++++++++++++++++++++++++
 flavor_catalog/references/C005/belle2010_arxiv1003_2345.txt       |  27 +++++++++++++++++++++++++++
 flavor_catalog/references/C005/burdman2002_arxiv_hepph0112235.txt |  23 +++++++++++++++++++++++
 flavor_catalog/references/C005/cfw2008_arxiv0804_1954.txt         |  28 ++++++++++++++++++++++++++++
 flavor_catalog/references/C005/pdg2026_d0_ee_pdgLive.txt          |  67 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/C005/source_manifest.yaml               |  60 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/C005.md                               |  34 ++++++++++++++++++++++++++++++++++
 9 files changed, 469 insertions(+)
```

### 180. `d1e077e` — flavor-catalog(charged_lepton): PKA draft for L006 muonium conversion
- SHA: `d1e077ebb160d02675078741dbf0389e5c55e7ac`
- Message: flavor-catalog(charged_lepton): PKA draft for L006 muonium conversion
- Physical/numerical summary: Added an initial catalog process draft for L006 muonium conversion; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charged_lepton/L006.tex                         |  87 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/charged_lepton/L006.yaml                        | 139 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/L006/agashe2006_rs_lfv_arxiv_hepph0606021.txt  |  26 ++++++++++++++++++++++++++
 flavor_catalog/references/L006/bai2022_snowmass_mace_arxiv2203_11406.txt |  30 ++++++++++++++++++++++++++++++
 flavor_catalog/references/L006/bai2024_mace_cdr_arxiv2410_18817.txt      |  31 +++++++++++++++++++++++++++++++
 flavor_catalog/references/L006/cfw2008_arxiv0804_1954.txt                |  21 +++++++++++++++++++++
 flavor_catalog/references/L006/heeck2024_delta_l_two_arxiv2401_09580.txt |  26 ++++++++++++++++++++++++++
 flavor_catalog/references/L006/pdg2026_muonium_antimuonium_s004mc.txt    |  38 ++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/L006/source_manifest.yaml                      |  75 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/L006/willmann1999_macs_arxiv_hepex9807011.txt  |  35 +++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/L006.md                                      |  34 ++++++++++++++++++++++++++++++++++
 11 files changed, 542 insertions(+)
```

### 181. `ddc55cf` — flavor-catalog(edm_neutrino): PKA draft for E002 muon EDM
- SHA: `ddc55cfa23e38dc7313930d7e07169ec7d3a0648`
- Message: flavor-catalog(edm_neutrino): PKA draft for E002 muon EDM
- Physical/numerical summary: Added an initial catalog process draft for E002 muon EDM; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/edm_neutrino/E002.tex                  |  88 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/edm_neutrino/E002.yaml                 | 160 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/E002/adelmann2021_arxiv2102_08838.txt |  25 +++++++++++++++++++++++
 flavor_catalog/references/E002/bennett2009_arxiv0811_1207.txt   |  25 +++++++++++++++++++++++
 flavor_catalog/references/E002/cfw2008_arxiv0804_1954.txt       |  23 ++++++++++++++++++++++
 flavor_catalog/references/E002/pdg2026_muon_edm_datablock.txt   |  34 ++++++++++++++++++++++++++++++++
 flavor_catalog/references/E002/renga2024_arxiv2409_20050.txt    |  26 ++++++++++++++++++++++++
 flavor_catalog/references/E002/source_manifest.yaml             |  55 +++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/E002.md                             |  27 +++++++++++++++++++++++++
 9 files changed, 463 insertions(+)
```

### 182. `febb3e8` — flavor-catalog(beauty): PKA draft for B021 Lambda_b to Lambda ll
- SHA: `febb3e83e9c6def4cc2c68fc7558e7f3a1ce019e`
- Message: flavor-catalog(beauty): PKA draft for B021 Lambda_b to Lambda ll
- Physical/numerical summary: Added an initial catalog process draft for B021 Lambda_b to Lambda ll; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B021.tex                                  |  86 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/beauty/B021.yaml                                 | 121 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B021/cdf_1107_3753_lambdab_observation.txt      |  22 ++++++++++++++++++++++
 flavor_catalog/references/B021/cfw_0804_1954_rs_flavor.txt                |  23 +++++++++++++++++++++++
 flavor_catalog/references/B021/detmold_1212_4827_lattice_form_factors.txt |  23 +++++++++++++++++++++++
 flavor_catalog/references/B021/lhcb_1503_07138_lambdab_lambdamumu.txt     |  27 +++++++++++++++++++++++++++
 flavor_catalog/references/B021/pdg2025_lambdab_lambdamumu.txt             |  42 ++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B021/sha256sums.txt                             |   5 +++++
 flavor_catalog/references/B021/source_manifest.yaml                       |  58 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/B021.md                                       |  37 +++++++++++++++++++++++++++++++++++++
 10 files changed, 444 insertions(+)
```

### 183. `4ba7e48` — flavor-catalog(beauty): PKA draft for B023 B to Kstar nunu
- SHA: `4ba7e48224525f96ff20148b0e4a4707fc055c2b`
- Message: flavor-catalog(beauty): PKA draft for B023 B to Kstar nunu
- Physical/numerical summary: Added an initial catalog process draft for B023 B to Kstar nunu; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B023.tex                      | 100 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/beauty/B023.yaml                     | 127 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B023/babar2013_arxiv1303_7465.txt   |  17 +++++++++++++
 flavor_catalog/references/B023/belle2013_arxiv1303_3719.txt   |  19 ++++++++++++++
 flavor_catalog/references/B023/belle2017_arxiv1702_03224.txt  |  24 ++++++++++++++++++
 flavor_catalog/references/B023/buras2015_arxiv1409_4557.txt   |  22 ++++++++++++++++
 flavor_catalog/references/B023/cfw2008_arxiv0804_1954.txt     |  15 +++++++++++
 flavor_catalog/references/B023/hflav_dec2025_b0_kst0_nunu.txt |  24 ++++++++++++++++++
 flavor_catalog/references/B023/hflav_dec2025_bp_kstp_nunu.txt |  24 ++++++++++++++++++
 flavor_catalog/references/B023/pdg2025_b0_kst0_nunu_api.txt   | 209 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B023/pdg2025_bp_kstp_nunu_api.txt   | 177 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B023/source_manifest.yaml           |  93 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/B023.md                           |  33 ++++++++++++++++++++++++
 13 files changed, 884 insertions(+)
```

### 184. `68e117c` — flavor-catalog(top_higgs_ew): PKA draft for T015 Z to e mu
- SHA: `68e117c0a3dbd8faf567365f974eed85905be484`
- Message: flavor-catalog(top_higgs_ew): PKA draft for T015 Z to e mu
- Physical/numerical summary: Added an initial catalog process draft for T015 Z to e mu; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B034.tex                                |  95 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/beauty/B034.yaml                               | 183 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/top_higgs_ew/T005.tex                          |  96 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/top_higgs_ew/T005.yaml                         | 166 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/top_higgs_ew/T015.tex                          |  91 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/top_higgs_ew/T015.yaml                         | 161 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B034/cfw2008_arxiv0804_1954.txt               |  25 +++++++++++++++++++
 flavor_catalog/references/B034/hflav2024_bs_phiphi_branching.txt        |  19 +++++++++++++++
 flavor_catalog/references/B034/lhcb2014_bs_phiphi_cp_arxiv.txt          |  23 ++++++++++++++++++
 flavor_catalog/references/B034/lhcb2019_bs_phiphi_cp_cds.txt            |  30 +++++++++++++++++++++++
 flavor_catalog/references/B034/lhcb2023_bs_phiphi_cp_cds.txt            |  33 ++++++++++++++++++++++++++
 flavor_catalog/references/B034/pdg2025_bs_phiphi_cp.txt                 |  23 ++++++++++++++++++
 flavor_catalog/references/B034/sha256sums.txt                           |   6 +++++
 flavor_catalog/references/B034/source_manifest.yaml                     |  70 ++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/T005/aguilar_saavedra_hepph0409342_tcg_sm.txt |  31 ++++++++++++++++++++++++
 flavor_catalog/references/T005/atlas_2112_01302_arxiv_abs.txt           |  41 ++++++++++++++++++++++++++++++++
 flavor_catalog/references/T005/cfw_0804_1954_rs_flavor.txt              |  31 ++++++++++++++++++++++++
 flavor_catalog/references/T005/cms_1610_03545_arxiv_abs.txt             |  32 +++++++++++++++++++++++++
 flavor_catalog/references/T005/pdg_2025_top_tug_tcg.txt                 |  44 ++++++++++++++++++++++++++++++++++
 flavor_catalog/references/T005/sha256sums.txt                           |   5 ++++
 flavor_catalog/references/T005/source_manifest.yaml                     |  49 ++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/T015/atlas2014_exot201302_arxiv_abs.txt       |  22 +++++++++++++++++
 flavor_catalog/references/T015/atlas2023_exot201835_arxiv_abs.txt       |  21 ++++++++++++++++
 flavor_catalog/references/T015/calibbi_marcano_roy2021_arxiv_abs.txt    |  23 ++++++++++++++++++
 flavor_catalog/references/T015/cms2025_smp23003_arxiv_abs.txt           |  27 +++++++++++++++++++++
 flavor_catalog/references/T015/csaki_falkowski_weiler2008_arxiv_abs.txt |  20 ++++++++++++++++
 flavor_catalog/references/T015/pdg2025_z_lfv_listing.txt                |  24 +++++++++++++++++++
 flavor_catalog/references/T015/perez_randall2008_arxiv_abs.txt          |  18 ++++++++++++++
 flavor_catalog/references/T015/sha256sums.txt                           |   7 ++++++
 flavor_catalog/references/T015/source_manifest.yaml                     |  67 +++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/B034.md                                     |  33 ++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/T005.md                                     |  37 +++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/T015.md                                     |  38 +++++++++++++++++++++++++++++
 33 files changed, 1591 insertions(+)
```

### 185. `945ede5` — flavor-catalog(beauty): PKA draft for B022 B to K nu nubar
- SHA: `945ede5a23ab0d527e98ce83f538afa140b5b5ee`
- Message: flavor-catalog(beauty): PKA draft for B022 B to K nu nubar
- Physical/numerical summary: Added an initial catalog process draft for B022 B to K nu nubar; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B022.tex                              |  93 +++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/beauty/B022.yaml                             | 133 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B022/babar2013_btokstarnunu_arxiv.txt       |  28 +++++++++++++
 flavor_catalog/references/B022/belleii2024_bplus_kplus_nunu_pubdb.txt | 305 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B022/belleii2025_likelihood_pubdb.txt       | 276 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B022/cfw2008_rs_flavor_arxiv.txt            |  34 +++++++++++++++
 flavor_catalog/references/B022/hflav_dec2025_bplus_kplus_nunu.txt     |  30 ++++++++++++++
 flavor_catalog/references/B022/hpqcd2023_btoknunu_sm_arxiv.txt        |  87 +++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B022/pdg2026_bplus_kplus_nunu_api.txt       | 321 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B022/source_manifest.yaml                   |  73 +++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/B022.md                                   |  35 ++++++++++++++++
 11 files changed, 1415 insertions(+)
```

### 186. `3aba841` — flavor-catalog(top_higgs_ew): PKA draft for T005 t to c g
- SHA: `3aba8415723e4871482bab23cf1bf5a7d64d0366`
- Message: flavor-catalog(top_higgs_ew): PKA draft for T005 t to c g
- Physical/numerical summary: Added an initial catalog process draft for T005 t to c g; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
  - `git show --stat` reported no file-level changes.

### 187. `c8dd354` — flavor-catalog(beauty): PKA draft for B034 Bs to phi phi
- SHA: `c8dd354fa5219a85eb2e20ee4a05354c30a63021`
- Message: flavor-catalog(beauty): PKA draft for B034 Bs to phi phi
- Physical/numerical summary: Added an initial catalog process draft for B034 Bs to phi phi; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B034.yaml                                                           |   9 ++++++-
 flavor_catalog/processes/kaon/K008.tex                                                              |  89 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/kaon/K008.yaml                                                             | 148 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/K008/aebischer_buras_kumar2022_rare_kaon_arxiv.txt                        |  19 +++++++++++++++
 flavor_catalog/references/K008/cfw2008_rs_flavor_arxiv.txt                                          |  24 +++++++++++++++++++
 flavor_catalog/references/K008/christ_feng_juettner_lawson_portelli_sachrajda2016_lattice_arxiv.txt |  26 ++++++++++++++++++++
 flavor_catalog/references/K008/isidori_smith_unterdorfer2004_kl_pi0ll_arxiv.txt                     |  29 ++++++++++++++++++++++
 flavor_catalog/references/K008/ktev2004_kl_pi0ee_arxiv.txt                                          |  27 +++++++++++++++++++++
 flavor_catalog/references/K008/na48_ks_pi0ee_pdgLive.txt                                            |  34 ++++++++++++++++++++++++++
 flavor_catalog/references/K008/pdg2026_kl_pi0ee_pdgLive.txt                                         |  46 +++++++++++++++++++++++++++++++++++
 flavor_catalog/references/K008/source_manifest.yaml                                                 |  74 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/K008.md                                                                 |  39 ++++++++++++++++++++++++++++++
 12 files changed, 563 insertions(+), 1 deletion(-)
```

### 188. `fabad04` — flavor-catalog(kaon_edm): CA batch ca_w4_kaon_edm_v2 — verify WA polish for K017 E006
- SHA: `fabad04a966086fb79c84b245a6143e1c883f23b`
- Message: flavor-catalog(kaon_edm): CA batch ca_w4_kaon_edm_v2 — verify WA polish for K017 E006
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/edm_neutrino/E006.yaml      | 15 ++++++++++++---
 flavor_catalog/processes/kaon/K017.yaml              | 15 ++++++++++++---
 flavor_catalog/worklogs/checker/ca_w4_kaon_edm_v2.md | 35 +++++++++++++++++++++++++++++++++++
 3 files changed, 59 insertions(+), 6 deletions(-)
```

### 189. `049fa43` — flavor-catalog(beauty): CA batch ca_w4_beauty_v2 — verify WA polish for B006
- SHA: `049fa43d8fc66260e4071dbd01bc3eb84592b2b2`
- Message: flavor-catalog(beauty): CA batch ca_w4_beauty_v2 — verify WA polish for B006
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B006.yaml          | 15 ++++++++++++---
 flavor_catalog/worklogs/checker/ca_w4_beauty_v2.md | 33 +++++++++++++++++++++++++++++++++
 2 files changed, 45 insertions(+), 3 deletions(-)
```

### 190. `d6e78b6` — flavor-catalog(top_higgs_ew): CA batch ca_w4_ew_v2 — verify WA polish for EW001 EW002 EW003
- SHA: `d6e78b6162ed814edc708b71a35e7584be0f35ce`
- Message: flavor-catalog(top_higgs_ew): CA batch ca_w4_ew_v2 — verify WA polish for EW001 EW002 EW003
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/EW001.yaml | 16 +++++++++++++---
 flavor_catalog/processes/top_higgs_ew/EW002.yaml | 16 +++++++++++++---
 flavor_catalog/processes/top_higgs_ew/EW003.yaml | 16 +++++++++++++---
 flavor_catalog/worklogs/checker/ca_w4_ew_v2.md   | 51 +++++++++++++++++++++++++++++++++++++++++++++++++++
 4 files changed, 90 insertions(+), 9 deletions(-)
```

### 191. `8a0d6fa` — flavor-catalog(charm): CA batch ca_w4_charm_v2 — verify WA polish for C002 C003 C004 C007
- SHA: `8a0d6fa5ba6cce4954967b758edeaa23aae4682d`
- Message: flavor-catalog(charm): CA batch ca_w4_charm_v2 — verify WA polish for C002 C003 C004 C007
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charm/C002.yaml          | 15 ++++++++++++---
 flavor_catalog/processes/charm/C003.yaml          | 15 ++++++++++++---
 flavor_catalog/processes/charm/C004.yaml          | 15 ++++++++++++---
 flavor_catalog/processes/charm/C007.yaml          | 15 ++++++++++++---
 flavor_catalog/worklogs/checker/ca_w4_charm_v2.md | 50 ++++++++++++++++++++++++++++++++++++++++++++++++++
 5 files changed, 98 insertions(+), 12 deletions(-)
```

### 192. `8aa7b96` — flavor-catalog(top_higgs_ew): WA batch wa_w5a_top_higgs_ew — polish PKA drafts for T005 T015 T019
- SHA: `8aa7b9607fe7fad6269a7a42857052d7aadc8dda`
- Message: flavor-catalog(top_higgs_ew): WA batch wa_w5a_top_higgs_ew — polish PKA drafts for T005 T015 T019
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T005.tex        | 37 +++++++++++++++++++++++--------------
 flavor_catalog/processes/top_higgs_ew/T005.yaml       | 11 ++++++++++-
 flavor_catalog/processes/top_higgs_ew/T015.tex        | 31 ++++++++++++++++++-------------
 flavor_catalog/processes/top_higgs_ew/T015.yaml       | 11 ++++++++++-
 flavor_catalog/processes/top_higgs_ew/T019.tex        | 34 +++++++++++++++++++++-------------
 flavor_catalog/processes/top_higgs_ew/T019.yaml       | 11 ++++++++++-
 flavor_catalog/worklogs/writer/wa_w5a_top_higgs_ew.md | 50 ++++++++++++++++++++++++++++++++++++++++++++++++++
 7 files changed, 142 insertions(+), 43 deletions(-)
```

### 193. `fd3c97e` — flavor-catalog(kaon): PKA draft for K010 K_S pi0 ee
- SHA: `fd3c97ed4e383d091faa927fb878a2351f85d3af`
- Message: flavor-catalog(kaon): PKA draft for K010 K_S pi0 ee
- Physical/numerical summary: Added an initial catalog process draft for K010 K_S pi0 ee; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/kaon/K010.tex                                  |  94 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/kaon/K010.yaml                                 | 103 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/K010/aebischer_buras_kumar_2022_arxiv.txt     |  18 ++++++++++++++++++
 flavor_catalog/references/K010/cfw2008_rs_flavor_arxiv.txt              |  18 ++++++++++++++++++
 flavor_catalog/references/K010/isidori_smith_unterdorfer_2004_arxiv.txt |  25 +++++++++++++++++++++++++
 flavor_catalog/references/K010/na48_2003_ks_pi0ee_arxiv.txt             |  22 ++++++++++++++++++++++
 flavor_catalog/references/K010/pdg2026_ks_pi0ee_pdgLive.txt             |  52 ++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/K010/sha256sums.txt                           |   5 +++++
 flavor_catalog/references/K010/source_manifest.yaml                     |  60 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/K010.md                                     |  31 +++++++++++++++++++++++++++++++
 10 files changed, 428 insertions(+)
```

### 194. `121559e` — flavor-catalog(charged_lepton_edm): WA batch wa_w5a_charged_lepton_edm - polish PKA drafts for L005 L006 E002
- SHA: `121559e0893526df6ad14084bda0729d401d0ed6`
- Message: flavor-catalog(charged_lepton_edm): WA batch wa_w5a_charged_lepton_edm - polish PKA drafts for L005 L006 E002
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charged_lepton/L005.tex            | 19 +++++++++++--------
 flavor_catalog/processes/charged_lepton/L005.yaml           | 10 +++++++++-
 flavor_catalog/processes/charged_lepton/L006.tex            | 24 ++++++++++++++----------
 flavor_catalog/processes/charged_lepton/L006.yaml           | 10 +++++++++-
 flavor_catalog/processes/edm_neutrino/E002.tex              | 19 ++++++++++++-------
 flavor_catalog/processes/edm_neutrino/E002.yaml             | 10 +++++++++-
 flavor_catalog/worklogs/writer/wa_w5a_charged_lepton_edm.md | 51 +++++++++++++++++++++++++++++++++++++++++++++++++++
 7 files changed, 115 insertions(+), 28 deletions(-)
```

### 195. `74ff5bc` — flavor-catalog(beauty): WA batch wa_w5a_beauty — polish PKA drafts for B021 B022 B023 B034
- SHA: `74ff5bc10775eb20087143173cd07c73d106adb7`
- Message: flavor-catalog(beauty): WA batch wa_w5a_beauty — polish PKA drafts for B021 B022 B023 B034
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B021.tex        | 28 ++++++++++++++++------------
 flavor_catalog/processes/beauty/B021.yaml       | 10 +++++++++-
 flavor_catalog/processes/beauty/B022.tex        | 30 +++++++++++++++++-------------
 flavor_catalog/processes/beauty/B022.yaml       | 10 +++++++++-
 flavor_catalog/processes/beauty/B023.tex        | 42 +++++++++++++++++++++---------------------
 flavor_catalog/processes/beauty/B023.yaml       | 10 +++++++++-
 flavor_catalog/processes/beauty/B034.tex        | 28 ++++++++++++++++------------
 flavor_catalog/processes/beauty/B034.yaml       | 10 +++++++++-
 flavor_catalog/worklogs/writer/wa_w5a_beauty.md | 72 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 9 files changed, 178 insertions(+), 62 deletions(-)
```

### 196. `140ab97` — flavor-catalog(charm): PKA draft for C006 D0 e mu LFV
- SHA: `140ab97459135760a840ccd9d2ce1d6ebd91f36f`
- Message: flavor-catalog(charm): PKA draft for C006 D0 e mu LFV
- Physical/numerical summary: Added an initial catalog process draft for C006 D0 e mu LFV; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charm/C006.tex                     |  95 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/charm/C006.yaml                    | 115 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/C006/babar2012_arxiv1206_5419.txt |  23 +++++++++++++++++++++++
 flavor_catalog/references/C006/belle2010_arxiv1003_2345.txt |  22 ++++++++++++++++++++++
 flavor_catalog/references/C006/cfw2008_arxiv0804_1954.txt   |  24 ++++++++++++++++++++++++
 flavor_catalog/references/C006/lhcb2016_arxiv1512_00322.txt |  24 ++++++++++++++++++++++++
 flavor_catalog/references/C006/pdg2026_d0_emu_pdgLive.txt   |  50 ++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/C006/source_manifest.yaml         |  60 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/C006.md                         |  33 +++++++++++++++++++++++++++++++++
 9 files changed, 446 insertions(+)
```

### 197. `92f3286` — flavor-catalog(kaon_charm): WA batch wa_w5a_kaon_charm — polish PKA drafts for K008 C005
- SHA: `92f32869d7a665c0baf3b0cc74c9b2ad298d0cb0`
- Message: flavor-catalog(kaon_charm): WA batch wa_w5a_kaon_charm — polish PKA drafts for K008 C005
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charm/C005.tex             | 76 ++++++++++++++++++++++++++++++++++++++++++----------------------------------
 flavor_catalog/processes/charm/C005.yaml            | 11 ++++++++++-
 flavor_catalog/processes/kaon/K008.tex              | 83 ++++++++++++++++++++++++++++++++++++++++++++++++-----------------------------------
 flavor_catalog/processes/kaon/K008.yaml             | 11 ++++++++++-
 flavor_catalog/worklogs/writer/wa_w5a_kaon_charm.md | 45 +++++++++++++++++++++++++++++++++++++++++++++
 5 files changed, 155 insertions(+), 71 deletions(-)
```

### 198. `de608a6` — flavor-catalog(charm): PKA draft for C008 D+ pi e mu
- SHA: `de608a6d370f89deb8cb1493bc658e27747bcd19`
- Message: flavor-catalog(charm): PKA draft for C008 D+ pi e mu
- Physical/numerical summary: Added an initial catalog process draft for C008 D+ pi e mu; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charm/C008.tex                                |  90 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/charm/C008.yaml                               | 114 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/C008/babar2011_arxiv1107_4465.txt            |  37 +++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/C008/cfw2008_arxiv0804_1954.txt              |  20 ++++++++++++++++++++
 flavor_catalog/references/C008/lhcb2021_arxiv2011_00217.txt            |  45 +++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/C008/pdg2026_dplus_piplus_em_mup_pdgLive.txt |  39 +++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/C008/pdg2026_dplus_piplus_ep_mum_pdgLive.txt |  39 +++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/C008/sha256sums.txt                          |   6 ++++++
 flavor_catalog/references/C008/source_manifest.yaml                    |  60 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/C008.md                                    |  32 ++++++++++++++++++++++++++++++++
 10 files changed, 482 insertions(+)
```

### 199. `a4519a6` — flavor-catalog(top_higgs_ew): PKA draft for T006 t to u g
- SHA: `a4519a6c837aa68d21b7f2a0ec241adf32f416e9`
- Message: flavor-catalog(top_higgs_ew): PKA draft for T006 t to u g
- Physical/numerical summary: Added an initial catalog process draft for T006 t to u g; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T006.tex                          |  94 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/top_higgs_ew/T006.yaml                         | 184 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/T006/aguilar_saavedra_hepph0409342_tug_sm.txt |  27 +++++++++++++++++++++
 flavor_catalog/references/T006/atlas_2112_01302_arxiv_abs.txt           |  37 ++++++++++++++++++++++++++++
 flavor_catalog/references/T006/cfw_0804_1954_rs_flavor.txt              |  28 ++++++++++++++++++++++
 flavor_catalog/references/T006/cms_1610_03545_arxiv_abs.txt             |  30 +++++++++++++++++++++++
 flavor_catalog/references/T006/pdg_2025_top_tug.txt                     |  44 ++++++++++++++++++++++++++++++++++
 flavor_catalog/references/T006/source_manifest.yaml                     |  49 ++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/T006.md                                     |  40 +++++++++++++++++++++++++++++++
 9 files changed, 533 insertions(+)
```

### 200. `d5dc2a6` — flavor-catalog(edm_neutrino): PKA draft for E007 radium xenon EDMs
- SHA: `d5dc2a6556f846bd7ffef8cb30534d0f419a34cb`
- Message: flavor-catalog(edm_neutrino): PKA draft for E007 radium xenon EDMs
- Physical/numerical summary: Added an initial catalog process draft for E007 radium xenon EDMs; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/edm_neutrino/E007.tex                      | 116 +++++++++++++++++++++++
 flavor_catalog/processes/edm_neutrino/E007.yaml                     | 179 +++++++++++++++++++++++++++++++++++
 flavor_catalog/references/E007/allmendinger2019_arxiv1904_12295.txt | 705 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/E007/argonne_ra225_edm_page.txt           | 231 +++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/E007/bishof2016_arxiv1606_04931.txt       | 703 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/E007/cfw2008_arxiv0804_1954.txt           | 723 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/E007/dfg_heil_schmidt_129xe_upgrade.txt   | 334 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/E007/kuleuven_raf_acf_project.txt         | 535 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/E007/parker2015_arxiv1504_07477.txt       | 685 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/E007/sachdeva2019_arxiv1909_12800.txt     | 735 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/E007/source_manifest.yaml                 |  86 +++++++++++++++++
 flavor_catalog/worklogs/pka/E007.md                                 |  34 +++++++
 12 files changed, 5066 insertions(+)
```

### 201. `74c6aa6` — flavor-catalog(top_higgs_ew): PKA draft for T017 Z -> mu tau
- SHA: `74c6aa67cbc6a7193983669dcd4d014a799d911a`
- Message: flavor-catalog(top_higgs_ew): PKA draft for T017 Z -> mu tau
- Physical/numerical summary: Added an initial catalog process draft for T017 Z -> mu tau; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T017.tex                          |  92 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/top_higgs_ew/T017.yaml                         | 157 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/T017/atlas2021_z_lfv_tau_arxiv_abs.txt        |  30 +++++++++++++++++++++++++++
 flavor_catalog/references/T017/calibbi_marcano_roy2021_arxiv_abs.txt    |  24 ++++++++++++++++++++++
 flavor_catalog/references/T017/cms2025_smp23003_arxiv_abs.txt           |  30 +++++++++++++++++++++++++++
 flavor_catalog/references/T017/csaki_falkowski_weiler2008_arxiv_abs.txt |  27 ++++++++++++++++++++++++
 flavor_catalog/references/T017/pdg2025_z_lfv_listing.txt                |  24 ++++++++++++++++++++++
 flavor_catalog/references/T017/perez_randall2008_arxiv_abs.txt          |  26 ++++++++++++++++++++++++
 flavor_catalog/references/T017/source_manifest.yaml                     |  58 ++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/T017.md                                     |  33 ++++++++++++++++++++++++++++++
 10 files changed, 501 insertions(+)
```

### 202. `a873420` — flavor-catalog(top_higgs_ew): PKA draft for T016 Z e tau
- SHA: `a873420c7740c4a7aed9a9a3cd7657a65580e7f2`
- Message: flavor-catalog(top_higgs_ew): PKA draft for T016 Z e tau
- Physical/numerical summary: Added an initial catalog process draft for T016 Z e tau; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T016.tex                          |  90 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/top_higgs_ew/T016.yaml                         | 180 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/T016/atlas2020_natphys_z_ltau_arxiv_abs.txt   |  24 +++++++++++++++++++
 flavor_catalog/references/T016/atlas2021_prl_z_ltau_arxiv_abs.txt       |  29 +++++++++++++++++++++++
 flavor_catalog/references/T016/calibbi_marcano_roy2021_arxiv_abs.txt    |  23 ++++++++++++++++++
 flavor_catalog/references/T016/cms2025_smp23003_arxiv_abs.txt           |  26 +++++++++++++++++++++
 flavor_catalog/references/T016/csaki_falkowski_weiler2008_arxiv_abs.txt |  23 ++++++++++++++++++
 flavor_catalog/references/T016/pdg2025_z_lfv_listing.txt                |  24 +++++++++++++++++++
 flavor_catalog/references/T016/perez_randall2008_arxiv_abs.txt          |  22 +++++++++++++++++
 flavor_catalog/references/T016/sha256sums.txt                           |   7 ++++++
 flavor_catalog/references/T016/source_manifest.yaml                     |  67 ++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/T016.md                                     |  34 +++++++++++++++++++++++++++
 12 files changed, 549 insertions(+)
```

### 203. `6a978ea` — flavor-catalog(kaon): PKA draft for K009 KL pi0 mu mu
- SHA: `6a978ea97872bdc0d1de3095318be54e344f76ed`
- Message: flavor-catalog(kaon): PKA draft for K009 KL pi0 mu mu
- Physical/numerical summary: Added an initial catalog process draft for K009 KL pi0 mu mu; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/kaon/K009.tex                                                              |  90 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/kaon/K009.yaml                                                             | 131 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/K009/aebischer_buras_kumar2022_rare_kaon_arxiv.txt                        |  21 ++++++++++++++++++
 flavor_catalog/references/K009/cfw2008_rs_flavor_arxiv.txt                                          |  23 ++++++++++++++++++++
 flavor_catalog/references/K009/christ_feng_juettner_lawson_portelli_sachrajda2016_lattice_arxiv.txt |  25 ++++++++++++++++++++++
 flavor_catalog/references/K009/isidori_smith_unterdorfer2004_kl_pi0mumu_arxiv.txt                   |  40 ++++++++++++++++++++++++++++++++++
 flavor_catalog/references/K009/ktev2000_kl_pi0mumu_arxiv.txt                                        |  21 ++++++++++++++++++
 flavor_catalog/references/K009/na48_ks_pi0mumu_arxiv.txt                                            |  21 ++++++++++++++++++
 flavor_catalog/references/K009/na48_ks_pi0mumu_pdgLive.txt                                          |  29 +++++++++++++++++++++++++
 flavor_catalog/references/K009/pdg2026_kl_pi0mumu_pdgLive.txt                                       |  29 +++++++++++++++++++++++++
 flavor_catalog/references/K009/source_manifest.yaml                                                 |  84 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/K009.md                                                                 |  37 ++++++++++++++++++++++++++++++++
 12 files changed, 551 insertions(+)
```

### 204. `a1d969e` — flavor-catalog(charged_lepton_edm): CA batch ca_w5a_charged_lepton_edm — verify WA polish for L005 L006 E002
- SHA: `a1d969e53a5ce86660ce43973c3dfa8e464986b3`
- Message: flavor-catalog(charged_lepton_edm): CA batch ca_w5a_charged_lepton_edm — verify WA polish for L005 L006 E002
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charged_lepton/L005.yaml            | 12 +++++++++++-
 flavor_catalog/processes/charged_lepton/L006.yaml            | 12 +++++++++++-
 flavor_catalog/processes/edm_neutrino/E002.yaml              | 12 +++++++++++-
 flavor_catalog/worklogs/checker/ca_w5a_charged_lepton_edm.md | 38 ++++++++++++++++++++++++++++++++++++++
 4 files changed, 71 insertions(+), 3 deletions(-)
```

### 205. `82780c6` — flavor-catalog(beauty): CA batch ca_w5a_beauty — verify WA polish for B021 B022 B023 B034
- SHA: `82780c69932920f47cabadf54add8a1e3ba523dd`
- Message: flavor-catalog(beauty): CA batch ca_w5a_beauty — verify WA polish for B021 B022 B023 B034
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B021.yaml        | 14 +++++++++++++-
 flavor_catalog/processes/beauty/B022.yaml        | 12 +++++++++++-
 flavor_catalog/processes/beauty/B023.yaml        | 12 +++++++++++-
 flavor_catalog/processes/beauty/B034.yaml        | 14 +++++++++++---
 flavor_catalog/worklogs/checker/ca_w5a_beauty.md | 70 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 5 files changed, 116 insertions(+), 6 deletions(-)
```

### 206. `ca886e0` — flavor-catalog(top_higgs_ew): CA batch ca_w5a_top_higgs_ew — verify WA polish for T005 T015 T019
- SHA: `ca886e02d9cf4c41648a8483341d31889bc1e5d9`
- Message: flavor-catalog(top_higgs_ew): CA batch ca_w5a_top_higgs_ew — verify WA polish for T005 T015 T019
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T005.yaml        | 14 ++++++++++++--
 flavor_catalog/processes/top_higgs_ew/T015.yaml        | 15 ++++++++++++---
 flavor_catalog/processes/top_higgs_ew/T019.yaml        | 14 ++++++++++++--
 flavor_catalog/worklogs/checker/ca_w5a_top_higgs_ew.md | 53 +++++++++++++++++++++++++++++++++++++++++++++++++++++
 4 files changed, 89 insertions(+), 7 deletions(-)
```

### 207. `8d02741` — flavor-catalog(kaon_charm): CA batch ca_w5a_kaon_charm — verify WA polish for K008 C005
- SHA: `8d02741bba170c66d829f156eff0e49b77aa1356`
- Message: flavor-catalog(kaon_charm): CA batch ca_w5a_kaon_charm — verify WA polish for K008 C005
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charm/C005.yaml             | 15 ++++++++++++---
 flavor_catalog/processes/kaon/K008.yaml              | 12 +++++++++++-
 flavor_catalog/worklogs/checker/ca_w5a_kaon_charm.md | 44 ++++++++++++++++++++++++++++++++++++++++++++
 3 files changed, 67 insertions(+), 4 deletions(-)
```

### 208. `7849dcc` — flavor-catalog(kaon): WA batch wa_w5b_kaon — polish PKA drafts for K009 K010
- SHA: `7849dccee4d6ad6cc990ca33c2209d732979b1a6`
- Message: flavor-catalog(kaon): WA batch wa_w5b_kaon — polish PKA drafts for K009 K010
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/kaon/K009.tex        | 49 ++++++++++++++++++++++++++++---------------------
 flavor_catalog/processes/kaon/K009.yaml       | 10 +++++++++-
 flavor_catalog/processes/kaon/K010.tex        | 23 +++++++++++------------
 flavor_catalog/processes/kaon/K010.yaml       | 10 +++++++++-
 flavor_catalog/worklogs/writer/wa_w5b_kaon.md | 26 ++++++++++++++++++++++++++
 5 files changed, 83 insertions(+), 35 deletions(-)
```

### 209. `d288119` — flavor-catalog(top_higgs_ew): WA batch wa_w5b_top_higgs_ew — polish PKA drafts for T006 T016 T017
- SHA: `d2881197592c996538f2a1efc15c9b93f56ad498`
- Message: flavor-catalog(top_higgs_ew): WA batch wa_w5b_top_higgs_ew — polish PKA drafts for T006 T016 T017
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T006.tex        | 53 ++++++++++++++++++++++++++++++++---------------------
 flavor_catalog/processes/top_higgs_ew/T006.yaml       | 10 +++++++++-
 flavor_catalog/processes/top_higgs_ew/T016.tex        | 12 +++++++-----
 flavor_catalog/processes/top_higgs_ew/T016.yaml       | 10 +++++++++-
 flavor_catalog/processes/top_higgs_ew/T017.tex        | 11 +++++------
 flavor_catalog/processes/top_higgs_ew/T017.yaml       | 10 +++++++++-
 flavor_catalog/worklogs/writer/wa_w5b_top_higgs_ew.md | 43 +++++++++++++++++++++++++++++++++++++++++++
 7 files changed, 114 insertions(+), 35 deletions(-)
```

### 210. `a9dd56c` — docs(flavor-catalog): codex quota pause + resumption plan
- SHA: `a9dd56c90965325a16bb17622838a1a3d2ca9234`
- Message: docs(flavor-catalog): codex quota pause + resumption plan
- Physical/numerical summary: Updated paper/repository documentation; physical impact is explanatory/provenance text rather than executable numerics.
- Files touched (`git show --stat`):
```text
 docs/phase_logs/flavor_catalog_codex_quota_pause.md | 89 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 89 insertions(+)
```

### 211. `4eca813` — flavor-catalog(charm_edm): WA batch wa_w5b_charm_edm — polish PKA drafts for C006 C008 E007
- SHA: `4eca813efc95a29834d4f0cde21d7c947663d6c3`
- Message: flavor-catalog(charm_edm): WA batch wa_w5b_charm_edm — polish PKA drafts for C006 C008 E007
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charm/C006.tex            | 49 ++++++++++++++++++++++++++++---------------------
 flavor_catalog/processes/charm/C006.yaml           | 11 ++++++++++-
 flavor_catalog/processes/charm/C008.tex            | 37 +++++++++++++++++++++----------------
 flavor_catalog/processes/charm/C008.yaml           | 11 ++++++++++-
 flavor_catalog/processes/edm_neutrino/E007.tex     | 53 +++++++++++++++++++++++++++++++----------------------
 flavor_catalog/processes/edm_neutrino/E007.yaml    | 11 ++++++++++-
 flavor_catalog/worklogs/writer/wa_w5b_charm_edm.md | 51 +++++++++++++++++++++++++++++++++++++++++++++++++++
 7 files changed, 161 insertions(+), 62 deletions(-)
```

### 212. `bb942a1` — flavor-catalog(signoff): Opus round 1 sign-off on 50 processes
- SHA: `bb942a1ffcba075e305c54b0df83ba95756ca341`
- Message: flavor-catalog(signoff): Opus round 1 sign-off on 50 processes
- Physical/numerical summary: Recorded reviewer/sign-off disposition (50 processes); review ledger only.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B002.yaml         |  13 +++++++++++--
 flavor_catalog/processes/beauty/B004.yaml         |  13 +++++++++++--
 flavor_catalog/processes/beauty/B005.yaml         |  13 +++++++++++--
 flavor_catalog/processes/beauty/B006.yaml         |  13 +++++++++++--
 flavor_catalog/processes/beauty/B009.yaml         |  13 +++++++++++--
 flavor_catalog/processes/beauty/B011.yaml         |  13 +++++++++++--
 flavor_catalog/processes/beauty/B015.yaml         |  13 +++++++++++--
 flavor_catalog/processes/beauty/B017.yaml         |  13 +++++++++++--
 flavor_catalog/processes/beauty/B018.yaml         |  13 +++++++++++--
 flavor_catalog/processes/beauty/B019.yaml         |  13 +++++++++++--
 flavor_catalog/processes/beauty/B025.yaml         |  13 +++++++++++--
 flavor_catalog/processes/beauty/B026.yaml         |  13 +++++++++++--
 flavor_catalog/processes/beauty/B032.yaml         |  13 +++++++++++--
 flavor_catalog/processes/beauty/B033.yaml         |  13 +++++++++++--
 flavor_catalog/processes/beauty/B034.yaml         |  13 +++++++++++--
 flavor_catalog/processes/charged_lepton/L001.yaml |  13 +++++++++++--
 flavor_catalog/processes/charged_lepton/L002.yaml |  13 +++++++++++--
 flavor_catalog/processes/charged_lepton/L003.yaml |  13 +++++++++++--
 flavor_catalog/processes/charged_lepton/L004.yaml |  13 +++++++++++--
 flavor_catalog/processes/charged_lepton/L007.yaml |  13 +++++++++++--
 flavor_catalog/processes/charged_lepton/L008.yaml |  13 +++++++++++--
 flavor_catalog/processes/charged_lepton/L009.yaml |  13 +++++++++++--
 flavor_catalog/processes/charged_lepton/L010.yaml |  13 +++++++++++--
 flavor_catalog/processes/charm/C001.yaml          |  13 +++++++++++--
 flavor_catalog/processes/charm/C002.yaml          |  13 +++++++++++--
 flavor_catalog/processes/charm/C003.yaml          |  13 +++++++++++--
 flavor_catalog/processes/charm/C004.yaml          |  13 +++++++++++--
 flavor_catalog/processes/charm/C005.yaml          |  13 +++++++++++--
 flavor_catalog/processes/charm/C007.yaml          |  13 +++++++++++--
 flavor_catalog/processes/edm_neutrino/E001.yaml   |  13 +++++++++++--
 flavor_catalog/processes/edm_neutrino/E004.yaml   |  13 +++++++++++--
 flavor_catalog/processes/edm_neutrino/E006.yaml   |  13 +++++++++++--
 flavor_catalog/processes/edm_neutrino/E008.yaml   |  13 +++++++++++--
 flavor_catalog/processes/kaon/K001.yaml           |  13 +++++++++++--
 flavor_catalog/processes/kaon/K002.yaml           |  13 +++++++++++--
 flavor_catalog/processes/kaon/K003.yaml           |  13 +++++++++++--
 flavor_catalog/processes/kaon/K004.yaml           |  13 +++++++++++--
 flavor_catalog/processes/kaon/K005.yaml           |  13 +++++++++++--
 flavor_catalog/processes/kaon/K006.yaml           |  13 +++++++++++--
 flavor_catalog/processes/kaon/K013.yaml           |  13 +++++++++++--
 flavor_catalog/processes/kaon/K017.yaml           |  13 +++++++++++--
 flavor_catalog/processes/top_higgs_ew/EW001.yaml  |  13 +++++++++++--
 flavor_catalog/processes/top_higgs_ew/EW002.yaml  |  13 +++++++++++--
 flavor_catalog/processes/top_higgs_ew/EW003.yaml  |  13 +++++++++++--
 flavor_catalog/processes/top_higgs_ew/T001.yaml   |  13 +++++++++++--
 flavor_catalog/processes/top_higgs_ew/T002.yaml   |  13 +++++++++++--
 flavor_catalog/processes/top_higgs_ew/T007.yaml   |  13 +++++++++++--
 flavor_catalog/processes/top_higgs_ew/T010.yaml   |  13 +++++++++++--
 flavor_catalog/processes/top_higgs_ew/T015.yaml   |  13 +++++++++++--
 flavor_catalog/processes/top_higgs_ew/T018.yaml   |  13 +++++++++++--
 flavor_catalog/signoff/round_001_index.md         | 155 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 51 files changed, 705 insertions(+), 100 deletions(-)
```

### 213. `7618cd0` — flavor-catalog: DA-3 discovery worklog round_003_final_sweep
- SHA: `7618cd0bb9b43e4383ac00b0b874377ade9e712b`
- Message: flavor-catalog: DA-3 discovery worklog round_003_final_sweep
- Physical/numerical summary: Applied the change described by the subject; based on stat surface, this is documentation/catalog/orchestration work rather than a new physics-code calculation.
- Files touched (`git show --stat`):
```text
 flavor_catalog/worklogs/discovery/round_003_final_sweep.md | 41 +++++++++++++++++++++++++++++++++++++++++
 1 file changed, 41 insertions(+)
```

### 214. `8ce0306` — flavor-catalog(beauty): WA v2 rework wa_w5a_beauty_v2 — address CA findings
- SHA: `8ce0306c99e0e4cc72c1b6b118476d10cb6bf76e`
- Message: flavor-catalog(beauty): WA v2 rework wa_w5a_beauty_v2 — address CA findings
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B021.yaml           | 75 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/processes/beauty/B022.yaml           | 24 +++++++++++++++++++++++-
 flavor_catalog/processes/beauty/B023.yaml           | 35 ++++++++++++++++++++++++++++++++++-
 flavor_catalog/references/B021/source_manifest.yaml | 10 +++++-----
 flavor_catalog/worklogs/writer/wa_w5a_beauty_v2.md  | 46 ++++++++++++++++++++++++++++++++++++++++++++++
 5 files changed, 182 insertions(+), 8 deletions(-)
```

### 215. `058c7f1` — flavor-catalog(kaon): CA batch ca_w5b_kaon — verify WA polish for K009 K010
- SHA: `058c7f15f81fd03d854d48b0e82ed68c30f7830a`
- Message: flavor-catalog(kaon): CA batch ca_w5b_kaon — verify WA polish for K009 K010
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/kaon/K009.yaml        | 12 +++++++++++-
 flavor_catalog/processes/kaon/K010.yaml        | 12 +++++++++++-
 flavor_catalog/worklogs/checker/ca_w5b_kaon.md | 38 ++++++++++++++++++++++++++++++++++++++
 3 files changed, 60 insertions(+), 2 deletions(-)
```

### 216. `45f051a` — flavor-catalog(kaon_charm): WA v2 rework wa_w5a_kaon_charm_v2 — address CA findings
- SHA: `45f051adab99d132c289cec0190d354fe42053dd`
- Message: flavor-catalog(kaon_charm): WA v2 rework wa_w5a_kaon_charm_v2 — address CA findings
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/kaon/K008.yaml                | 244 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/worklogs/writer/wa_w5a_kaon_charm_v2.md |  40 ++++++++++++++++++++++++++
 2 files changed, 283 insertions(+), 1 deletion(-)
```

### 217. `90b60ea` — flavor-catalog(top_higgs_ew): CA batch ca_w5b_top_higgs_ew — verify WA polish for T006 T016 T017
- SHA: `90b60ea350dafc0bf1bd0407e78e9f766149b672`
- Message: flavor-catalog(top_higgs_ew): CA batch ca_w5b_top_higgs_ew — verify WA polish for T006 T016 T017
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T006.yaml        | 14 ++++++++++++--
 flavor_catalog/processes/top_higgs_ew/T016.yaml        | 15 ++++++++++++---
 flavor_catalog/processes/top_higgs_ew/T017.yaml        | 15 ++++++++++++---
 flavor_catalog/worklogs/checker/ca_w5b_top_higgs_ew.md | 48 ++++++++++++++++++++++++++++++++++++++++++++++++
 4 files changed, 84 insertions(+), 8 deletions(-)
```

### 218. `b55bed6` — flavor-catalog(charm_edm): CA batch ca_w5b_charm_edm - verify WA polish for C006 C008 E007
- SHA: `b55bed6bc0e9f834ef2045377dc9f04c705700ea`
- Message: flavor-catalog(charm_edm): CA batch ca_w5b_charm_edm - verify WA polish for C006 C008 E007
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charm/C006.yaml            | 15 ++++++++++++---
 flavor_catalog/processes/charm/C008.yaml            | 15 ++++++++++++---
 flavor_catalog/processes/edm_neutrino/E007.yaml     | 14 ++++++++++++--
 flavor_catalog/worklogs/checker/ca_w5b_charm_edm.md | 35 +++++++++++++++++++++++++++++++++++
 4 files changed, 71 insertions(+), 8 deletions(-)
```

### 219. `a496085` — flavor-catalog(charged_lepton_edm): WA v2 rework wa_w5a_charged_lepton_edm_v2 — address CA findings
- SHA: `a496085ae75978ead7a63578bf49edda5f14f374`
- Message: flavor-catalog(charged_lepton_edm): WA v2 rework wa_w5a_charged_lepton_edm_v2 — address CA findings
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charged_lepton/L005.yaml              | 47 ++++++++++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/processes/charged_lepton/L006.yaml              | 36 +++++++++++++++++++++++++++++++++++-
 flavor_catalog/processes/edm_neutrino/E002.yaml                | 90 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/worklogs/writer/wa_w5a_charged_lepton_edm_v2.md | 40 ++++++++++++++++++++++++++++++++++++++++
 4 files changed, 210 insertions(+), 3 deletions(-)
```

### 220. `16ef205` — flavor-catalog(top_higgs_ew): WA v2 rework wa_w5a_top_higgs_ew_v2 — address CA findings
- SHA: `16ef205e7136a360495d203d78061819d5b9b46f`
- Message: flavor-catalog(top_higgs_ew): WA v2 rework wa_w5a_top_higgs_ew_v2 — address CA findings
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T005.yaml          | 49 +++++++++++++++++++++++++++++++++++++++++++++++--
 flavor_catalog/processes/top_higgs_ew/T019.yaml          | 29 +++++++++++++++++++++++++++--
 flavor_catalog/worklogs/writer/wa_w5a_top_higgs_ew_v2.md | 41 +++++++++++++++++++++++++++++++++++++++++
 3 files changed, 115 insertions(+), 4 deletions(-)
```

### 221. `51fe52f` — flavor-catalog(charm_edm): WA v2 rework wa_w5b_charm_edm_v2 — address CA findings
- SHA: `51fe52f8e32f3a20adc999a92ee375c04af6f958`
- Message: flavor-catalog(charm_edm): WA v2 rework wa_w5b_charm_edm_v2 — address CA findings
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/edm_neutrino/E007.tex        | 21 +++++++++------------
 flavor_catalog/processes/edm_neutrino/E007.yaml       | 12 +++++++++++-
 flavor_catalog/worklogs/writer/wa_w5b_charm_edm_v2.md | 24 ++++++++++++++++++++++++
 3 files changed, 44 insertions(+), 13 deletions(-)
```

### 222. `0fcb364` — flavor-catalog(top_higgs_ew): WA v2 rework wa_w5b_top_higgs_ew_v2 — address CA findings
- SHA: `0fcb3643260510d912805c698ae52aab2a3198c7`
- Message: flavor-catalog(top_higgs_ew): WA v2 rework wa_w5b_top_higgs_ew_v2 — address CA findings
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T006.yaml          | 31 +++++++++++++++++++++----------
 flavor_catalog/worklogs/writer/wa_w5b_top_higgs_ew_v2.md | 40 ++++++++++++++++++++++++++++++++++++++++
 2 files changed, 61 insertions(+), 10 deletions(-)
```

### 223. `3b76b20` — flavor-catalog(kaon_charm): CA batch ca_w5a_kaon_charm_v2 — verify WA polish for K008
- SHA: `3b76b20bbc77141d7341b7abb32d54ba7807ba8a`
- Message: flavor-catalog(kaon_charm): CA batch ca_w5a_kaon_charm_v2 — verify WA polish for K008
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/kaon/K008.yaml                 | 13 +++++++++++--
 flavor_catalog/worklogs/checker/ca_w5a_kaon_charm_v2.md | 38 ++++++++++++++++++++++++++++++++++++++
 2 files changed, 49 insertions(+), 2 deletions(-)
```

### 224. `4f07abf` — flavor-catalog(charged_lepton_edm): CA batch ca_w5a_charged_lepton_edm_v2 — verify WA polish for L005 L006 E002
- SHA: `4f07abfee08c950539090c4c08b368e0cbc6cc5e`
- Message: flavor-catalog(charged_lepton_edm): CA batch ca_w5a_charged_lepton_edm_v2 — verify WA polish for L005 L006 E002
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charged_lepton/L005.yaml               | 13 +++++++++++--
 flavor_catalog/processes/charged_lepton/L006.yaml               | 13 +++++++++++--
 flavor_catalog/processes/edm_neutrino/E002.yaml                 | 13 +++++++++++--
 flavor_catalog/worklogs/checker/ca_w5a_charged_lepton_edm_v2.md | 38 ++++++++++++++++++++++++++++++++++++++
 4 files changed, 71 insertions(+), 6 deletions(-)
```

### 225. `923b8cd` — flavor-catalog(kaon): WA v2 rework wa_w5b_kaon_v2 — address CA findings
- SHA: `923b8cd93870fd08da208572c9b457e2cf13450e`
- Message: flavor-catalog(kaon): WA v2 rework wa_w5b_kaon_v2 — address CA findings
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/kaon/K009.yaml          | 185 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/processes/kaon/K010.yaml          |  80 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/worklogs/writer/wa_w5b_kaon_v2.md |  27 ++++++++++++++++++++++++
 3 files changed, 290 insertions(+), 2 deletions(-)
```

### 226. `64b5043` — flavor-catalog(beauty): PKA draft for B001 Delta m_d
- SHA: `64b50438a19697dbf19660136f1e18adfd5971dc`
- Message: flavor-catalog(beauty): PKA draft for B001 Delta m_d
- Physical/numerical summary: Added an initial catalog process draft for B001 Delta m_d; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B001.tex                        |  91 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/beauty/B001.yaml                       | 136 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B001/belleii_2023_bd_mixing_arxiv.txt |  24 ++++++++++++++++++++++++
 flavor_catalog/references/B001/cfw_2008_rs_flavor_arxiv.txt     |  22 ++++++++++++++++++++++
 flavor_catalog/references/B001/flag2024_b_mixing_arxiv.txt      |  22 ++++++++++++++++++++++
 flavor_catalog/references/B001/hflav_pdg2025_bd_mixing.txt      |  26 ++++++++++++++++++++++++++
 flavor_catalog/references/B001/source_manifest.yaml             |  49 +++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/B001.md                             |  28 ++++++++++++++++++++++++++++
 8 files changed, 398 insertions(+)
```

### 227. `d596b3e` — flavor-catalog(beauty): CA batch ca_w5a_beauty_v2 — verify WA polish for B021 B022 B023
- SHA: `d596b3ec3361df89490bada087c8fd6f6fc29435`
- Message: flavor-catalog(beauty): CA batch ca_w5a_beauty_v2 — verify WA polish for B021 B022 B023
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B021.yaml           | 12 +++++++++++-
 flavor_catalog/processes/beauty/B022.yaml           | 12 ++++++++++--
 flavor_catalog/processes/beauty/B023.yaml           | 12 +++++++++++-
 flavor_catalog/worklogs/checker/ca_w5a_beauty_v2.md | 65 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 4 files changed, 97 insertions(+), 4 deletions(-)
```

### 228. `17ab228` — flavor-catalog(beauty): PKA draft for B016 B to K ll
- SHA: `17ab228e5f0e4ae7b8fcde2660155adc38411b6b`
- Message: flavor-catalog(beauty): PKA draft for B016 B to K ll
- Physical/numerical summary: Added an initial catalog process draft for B016 B to K ll; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B016.tex                           | 88 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/beauty/B016.yaml                          | 84 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B016/babar_0807_4119_BtoKll.txt          | 12 ++++++++++++
 flavor_catalog/references/B016/belle_1908_01848_BtoKll.txt         | 12 ++++++++++++
 flavor_catalog/references/B016/cfw_0804_1954_rs_flavor.txt         | 12 ++++++++++++
 flavor_catalog/references/B016/hflav_dec2025_B0_to_K0_ll.txt       | 23 +++++++++++++++++++++++
 flavor_catalog/references/B016/hflav_dec2025_Bplus_to_Kplus_ll.txt | 18 ++++++++++++++++++
 flavor_catalog/references/B016/lhcb_1403_8044_BtoKmumu.txt         | 12 ++++++++++++
 flavor_catalog/references/B016/source_manifest.yaml                | 65 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/B016.md                                | 37 +++++++++++++++++++++++++++++++++++++
 10 files changed, 363 insertions(+)
```

### 229. `bc4f253` — flavor-catalog(top_higgs_ew): CA batch ca_w5a_top_higgs_ew_v2 — verify WA polish for T005 T019
- SHA: `bc4f25344363cdd678d3f92d3938fba996f92a4e`
- Message: flavor-catalog(top_higgs_ew): CA batch ca_w5a_top_higgs_ew_v2 — verify WA polish for T005 T019
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T005.yaml           | 14 ++++++++++++--
 flavor_catalog/processes/top_higgs_ew/T019.yaml           | 14 ++++++++++++--
 flavor_catalog/worklogs/checker/ca_w5a_top_higgs_ew_v2.md | 45 +++++++++++++++++++++++++++++++++++++++++++++
 3 files changed, 69 insertions(+), 4 deletions(-)
```

### 230. `6a16bf4` — flavor-catalog(kaon): PKA draft for K012 K_S to mu mu
- SHA: `6a16bf498c7f86efaf8e3d1b1d30da2fb5e4befd`
- Message: flavor-catalog(kaon): PKA draft for K012 K_S to mu mu
- Physical/numerical summary: Added an initial catalog process draft for K012 K_S to mu mu; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/kaon/K012.tex                                          |  77 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/kaon/K012.yaml                                         | 106 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/K012/chobanova_2018_ksmumu_susy_arxiv.txt             |  21 +++++++++++++++++++++
 flavor_catalog/references/K012/dery_ghosh_grossman_schacht_2021_kmumu_arxiv.txt |  20 ++++++++++++++++++++
 flavor_catalog/references/K012/lhcb_2017_ksmumu_arxiv.txt                       |  20 ++++++++++++++++++++
 flavor_catalog/references/K012/lhcb_2020_ksmumu_arxiv.txt                       |  20 ++++++++++++++++++++
 flavor_catalog/references/K012/pdg2025_ks_decay_modes.txt                       |  43 +++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/K012/source_manifest.yaml                             |  54 ++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/K012.md                                             |  31 +++++++++++++++++++++++++++++++
 9 files changed, 392 insertions(+)
```

### 231. `3d2f833` — flavor-catalog(beauty): PKA draft for B003 Delta m_s
- SHA: `3d2f833cf7a556b05a05f56ca2ab4f0bdaa173f2`
- Message: flavor-catalog(beauty): PKA draft for B003 Delta m_s
- Physical/numerical summary: Added an initial catalog process draft for B003 Delta m_s; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B003.tex                   | 106 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/beauty/B003.yaml                  | 125 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B003/cfw_2008_rs_flavor.txt      |  26 ++++++++++++++++++++++++++
 flavor_catalog/references/B003/flag_2024_bmixing.txt       |  40 ++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B003/hflav_2024_dms.txt          |  33 +++++++++++++++++++++++++++++++++
 flavor_catalog/references/B003/hflav_pdg2025_dms.txt       |  32 ++++++++++++++++++++++++++++++++
 flavor_catalog/references/B003/lhcb_2021_deltams_arxiv.txt |  36 ++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B003/sha256sums.txt              |   5 +++++
 flavor_catalog/references/B003/source_manifest.yaml        |  59 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/B003.md                        |  30 ++++++++++++++++++++++++++++++
 10 files changed, 492 insertions(+)
```

### 232. `98b203a` — flavor-catalog(kaon): PKA draft for K018 Kl3 Vus
- SHA: `98b203ae177cd89df5028140a940cee8603ba18b`
- Message: flavor-catalog(kaon): PKA draft for K018 Kl3 Vus
- Physical/numerical summary: Added an initial catalog process draft for K018 Kl3 Vus; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/kaon/K018.tex                        |  85 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/kaon/K018.yaml                       | 220 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/K018/cfw2008_arxiv0804_1954.txt     |  25 +++++++++++++++++
 flavor_catalog/references/K018/flag2024_kaon_semileptonic.txt |  42 +++++++++++++++++++++++++++++
 flavor_catalog/references/K018/flavianet2010_kl3.txt          |  51 +++++++++++++++++++++++++++++++++++
 flavor_catalog/references/K018/pdg2025_vud_vus_kl3.txt        |  49 ++++++++++++++++++++++++++++++++++
 flavor_catalog/references/K018/sha256sums.txt                 |   4 +++
 flavor_catalog/references/K018/source_manifest.yaml           |  46 ++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/K018.md                           |  30 +++++++++++++++++++++
 9 files changed, 552 insertions(+)
```

### 233. `fea0961` — flavor-catalog(top_higgs_ew): PKA draft for T020 h to e mu
- SHA: `fea0961431dd55b47071332887c312e36b928f1c`
- Message: flavor-catalog(top_higgs_ew): PKA draft for T020 h to e mu
- Physical/numerical summary: Added an initial catalog process draft for T020 h to e mu; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T020.tex                          |  94 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/top_higgs_ew/T020.yaml                         | 170 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/T020/atlas2020_higgs_ee_emu_arxiv_abs.txt     |  25 +++++++++++++++++++++
 flavor_catalog/references/T020/cms2023_hig22002_public_page.txt         |  30 +++++++++++++++++++++++++
 flavor_catalog/references/T020/csaki_falkowski_weiler2008_arxiv_abs.txt |  23 +++++++++++++++++++
 flavor_catalog/references/T020/harnik_kopp_zupan2012_arxiv_abs.txt      |  24 ++++++++++++++++++++
 flavor_catalog/references/T020/pdg2025_higgs_lfv_review.txt             |  32 +++++++++++++++++++++++++++
 flavor_catalog/references/T020/perez_randall2008_arxiv_abs.txt          |  22 ++++++++++++++++++
 flavor_catalog/references/T020/sha256sums.txt                           |   6 +++++
 flavor_catalog/references/T020/source_manifest.yaml                     |  59 +++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/T020.md                                     |  32 +++++++++++++++++++++++++++
 11 files changed, 517 insertions(+)
```

### 234. `bbc91de` — flavor-catalog(edm_neutrino): PKA draft for E009 Weinberg three-gluon operator
- SHA: `bbc91de773db515b4a0a6a02087f09ce9dac406f`
- Message: flavor-catalog(edm_neutrino): PKA draft for E009 Weinberg three-gluon operator
- Physical/numerical summary: Added an initial catalog process draft for E009 Weinberg three-gluon operator; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/edm_neutrino/E009.tex                                                |  91 +++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/edm_neutrino/E009.yaml                                               | 210 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/E009/abel2020_neutron_edm_arxiv2001_11966.txt                       |  23 +++++++++++++
 flavor_catalog/references/E009/bhattacharya2022_weinberg_lattice_arxiv2203_03746.txt          |  31 ++++++++++++++++++
 flavor_catalog/references/E009/cfw2008_arxiv0804_1954.txt                                     |  16 +++++++++
 flavor_catalog/references/E009/haisch_hala2019_weinberg_sum_rules_arxiv1909_08955.txt         |  33 +++++++++++++++++++
 flavor_catalog/references/E009/koenig_neubert_straub2014_composite_dipoles_arxiv1403_2756.txt |  35 ++++++++++++++++++++
 flavor_catalog/references/E009/pdg2026_neutron_edm_datablock.txt                              |  29 +++++++++++++++++
 flavor_catalog/references/E009/pospelov_ritz2005_edm_review_arxiv_hepph0504231.txt            |  27 ++++++++++++++++
 flavor_catalog/references/E009/ramsey_musolf2014_edm_global_analysis_arxiv1407_1064.txt       |  30 +++++++++++++++++
 flavor_catalog/references/E009/source_manifest.yaml                                           |  99 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/E009/weinberg1989_inspire_metadata.txt                              |  27 ++++++++++++++++
 flavor_catalog/worklogs/pka/E009.md                                                           |  36 +++++++++++++++++++++
 13 files changed, 687 insertions(+)
```

### 235. `821fe14` — flavor-catalog(charged_lepton): PKA draft for L023 neutrino trident
- SHA: `821fe14713dcaf7f3ec8b87369a3cd63b52dffeb`
- Message: flavor-catalog(charged_lepton): PKA draft for L023 neutrino trident
- Physical/numerical summary: Added an initial catalog process draft for L023 neutrino trident; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charged_lepton/L023.tex                   |  85 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/charged_lepton/L023.yaml                  | 154 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/L023/altmannshofer2014_trident_arxiv.txt |  25 ++++++++++++++++++++++++
 flavor_catalog/references/L023/ccfr1991_prl_snapshot.txt           |  25 ++++++++++++++++++++++++
 flavor_catalog/references/L023/cfw2008_arxiv0804_1954.txt          |  18 +++++++++++++++++
 flavor_catalog/references/L023/charmii1990_cds_snapshot.txt        |  26 +++++++++++++++++++++++++
 flavor_catalog/references/L023/dune2019_tridents_arxiv.txt         |  29 ++++++++++++++++++++++++++++
 flavor_catalog/references/L023/kaneta2017_belleii_arxiv.txt        |  27 ++++++++++++++++++++++++++
 flavor_catalog/references/L023/nutev1998_arxiv_snapshot.txt        |  24 +++++++++++++++++++++++
 flavor_catalog/references/L023/source_manifest.yaml                |  68 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/L023.md                                |  33 +++++++++++++++++++++++++++++++
 11 files changed, 514 insertions(+)
```

### 236. `f6ab87a` — flavor-catalog(arbitration): Opus signoff on B021 + B023 cycle-3 cap
- SHA: `f6ab87a4df6ff678ebac25a2be2793c15cea0a34`
- Message: flavor-catalog(arbitration): Opus signoff on B021 + B023 cycle-3 cap
- Physical/numerical summary: Recorded reviewer/sign-off disposition; review ledger only.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B021.yaml      |  15 ++++++++++++---
 flavor_catalog/processes/beauty/B023.tex       |   4 ++--
 flavor_catalog/processes/beauty/B023.yaml      |  15 ++++++++++++---
 flavor_catalog/signoff/by_process/B021_B023.md | 108 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 4 files changed, 134 insertions(+), 8 deletions(-)
```

### 237. `7f87327` — flavor-catalog(charm_edm): CA batch ca_w5b_charm_edm_v2 - verify WA polish for E007
- SHA: `7f87327b7b26fde56e3c8d8901d193d20cfacce1`
- Message: flavor-catalog(charm_edm): CA batch ca_w5b_charm_edm_v2 - verify WA polish for E007
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/edm_neutrino/E007.yaml        | 13 +++++++++++--
 flavor_catalog/worklogs/checker/ca_w5b_charm_edm_v2.md | 30 ++++++++++++++++++++++++++++++
 2 files changed, 41 insertions(+), 2 deletions(-)
```

### 238. `b709912` — flavor-catalog(kaon): CA batch ca_w5b_kaon_v2 — verify WA polish for K009 K010
- SHA: `b7099122f0a7dcb76e75e366181a01c225d6be62`
- Message: flavor-catalog(kaon): CA batch ca_w5b_kaon_v2 — verify WA polish for K009 K010
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/kaon/K009.yaml           | 13 +++++++++++--
 flavor_catalog/processes/kaon/K010.yaml           | 13 +++++++++++--
 flavor_catalog/worklogs/checker/ca_w5b_kaon_v2.md | 37 +++++++++++++++++++++++++++++++++++++
 3 files changed, 59 insertions(+), 4 deletions(-)
```

### 239. `398b939` — flavor-catalog(edm_neutrino): WA batch wa_w6_edm — polish PKA drafts for E009
- SHA: `398b93937667f99443c4e2ab6344cdd528258afc`
- Message: flavor-catalog(edm_neutrino): WA batch wa_w6_edm — polish PKA drafts for E009
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/edm_neutrino/E009.tex  | 49 ++++++++++++++++++++++++++++---------------------
 flavor_catalog/processes/edm_neutrino/E009.yaml | 11 ++++++++++-
 flavor_catalog/worklogs/writer/wa_w6_edm.md     | 30 ++++++++++++++++++++++++++++++
 3 files changed, 68 insertions(+), 22 deletions(-)
```

### 240. `50793c3` — flavor-catalog(charged_lepton_top): WA batch wa_w6_charged_lepton_top — polish PKA drafts for L023 T020
- SHA: `50793c39cccbc266d4c09e91f23dc973499609d0`
- Message: flavor-catalog(charged_lepton_top): WA batch wa_w6_charged_lepton_top — polish PKA drafts for L023 T020
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charged_lepton/L023.tex           | 60 +++++++++++++++++++++++++++++++++---------------------------
 flavor_catalog/processes/charged_lepton/L023.yaml          | 11 ++++++++++-
 flavor_catalog/processes/top_higgs_ew/T020.tex             | 37 ++++++++++++++++++++-----------------
 flavor_catalog/processes/top_higgs_ew/T020.yaml            | 11 ++++++++++-
 flavor_catalog/worklogs/writer/wa_w6_charged_lepton_top.md | 33 +++++++++++++++++++++++++++++++++
 5 files changed, 106 insertions(+), 46 deletions(-)
```

### 241. `a28aa73` — flavor-catalog(beauty): WA batch wa_w6_beauty — polish PKA drafts for B001 B003 B016
- SHA: `a28aa736e056c9059b457a668357a9cba3253526`
- Message: flavor-catalog(beauty): WA batch wa_w6_beauty — polish PKA drafts for B001 B003 B016
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B001.tex       | 34 ++++++++++++++++++----------------
 flavor_catalog/processes/beauty/B001.yaml      | 10 +++++++++-
 flavor_catalog/processes/beauty/B003.tex       | 24 +++++++++++++-----------
 flavor_catalog/processes/beauty/B003.yaml      | 10 +++++++++-
 flavor_catalog/processes/beauty/B016.tex       | 44 ++++++++++++++++++++++++--------------------
 flavor_catalog/processes/beauty/B016.yaml      | 10 +++++++++-
 flavor_catalog/worklogs/writer/wa_w6_beauty.md | 65 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 7 files changed, 147 insertions(+), 50 deletions(-)
```

### 242. `297b820` — flavor-catalog(kaon): WA batch wa_w6_kaon — polish PKA drafts for K012 K018
- SHA: `297b820dbd616e2b8134f17da5cd8b5faa7920c6`
- Message: flavor-catalog(kaon): WA batch wa_w6_kaon — polish PKA drafts for K012 K018
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/kaon/K012.tex       | 52 +++++++++++++++++++++++++++-------------------------
 flavor_catalog/processes/kaon/K012.yaml      | 11 ++++++++++-
 flavor_catalog/processes/kaon/K018.tex       | 57 +++++++++++++++++++++++++++++++--------------------------
 flavor_catalog/processes/kaon/K018.yaml      | 11 ++++++++++-
 flavor_catalog/worklogs/writer/wa_w6_kaon.md | 52 ++++++++++++++++++++++++++++++++++++++++++++++++++++
 5 files changed, 130 insertions(+), 53 deletions(-)
```

### 243. `f115228` — flavor-catalog(top_higgs_ew): CA batch ca_w5b_top_higgs_ew_v2 — verify WA polish for T006
- SHA: `f11522863b429c5008635da7ae804c453921e53d`
- Message: flavor-catalog(top_higgs_ew): CA batch ca_w5b_top_higgs_ew_v2 — verify WA polish for T006
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T006.yaml           | 13 +++++++++++--
 flavor_catalog/worklogs/checker/ca_w5b_top_higgs_ew_v2.md | 35 +++++++++++++++++++++++++++++++++++
 2 files changed, 46 insertions(+), 2 deletions(-)
```

### 244. `aaa22d6` — flavor-catalog(edm_neutrino): CA batch ca_w6_edm — verify WA polish for E009
- SHA: `aaa22d63ef4b9f3afdf8a8c1ab06006a1dc8a7a8`
- Message: flavor-catalog(edm_neutrino): CA batch ca_w6_edm — verify WA polish for E009
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/edm_neutrino/E009.yaml | 12 +++++++++++-
 flavor_catalog/worklogs/checker/ca_w6_edm.md    | 32 ++++++++++++++++++++++++++++++++
 2 files changed, 43 insertions(+), 1 deletion(-)
```

### 245. `a15edea` — flavor-catalog(kaon): CA batch ca_w6_kaon — verify WA polish for K012 K018
- SHA: `a15edea39615a642ad48db49fc3ac1a5712830a4`
- Message: flavor-catalog(kaon): CA batch ca_w6_kaon — verify WA polish for K012 K018
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charged_lepton/L023.yaml           | 12 +++++++++++-
 flavor_catalog/processes/kaon/K012.yaml                     | 11 ++++++++++-
 flavor_catalog/processes/kaon/K018.yaml                     | 11 ++++++++++-
 flavor_catalog/processes/top_higgs_ew/T020.yaml             | 12 +++++++++++-
 flavor_catalog/worklogs/checker/ca_w6_charged_lepton_top.md | 38 ++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/checker/ca_w6_kaon.md               | 47 +++++++++++++++++++++++++++++++++++++++++++++++
 6 files changed, 127 insertions(+), 4 deletions(-)
```

### 246. `bdb8371` — flavor-catalog(beauty): CA batch ca_w6_beauty - verify WA polish for B001 B003 B016
- SHA: `bdb837124781e0b49e52cd816cd10fa619b9c034`
- Message: flavor-catalog(beauty): CA batch ca_w6_beauty - verify WA polish for B001 B003 B016
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B001.yaml       | 14 +++++++++++++-
 flavor_catalog/processes/beauty/B003.yaml       | 15 ++++++++++++++-
 flavor_catalog/processes/beauty/B016.yaml       | 13 ++++++++++++-
 flavor_catalog/worklogs/checker/ca_w6_beauty.md | 47 +++++++++++++++++++++++++++++++++++++++++++++++
 4 files changed, 86 insertions(+), 3 deletions(-)
```

### 247. `c6ec949` — flavor-catalog(charged_lepton_top): CA batch ca_w6_charged_lepton_top — verify WA polish for L023 T020
- SHA: `c6ec949a617e280db496ce8de0f9b27c8f409b62`
- Message: flavor-catalog(charged_lepton_top): CA batch ca_w6_charged_lepton_top — verify WA polish for L023 T020
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
  - `git show --stat` reported no file-level changes.

### 248. `e821611` — flavor-catalog: DA-4 convergence check; recommend lock at 75 + DEFERRED-SCOPE tail
- SHA: `e8216113e6e0aef8f1d01af05c54143afcee14d3`
- Message: flavor-catalog: DA-4 convergence check; recommend lock at 75 + DEFERRED-SCOPE tail
- Physical/numerical summary: Applied the change described by the subject; based on stat surface, this is documentation/catalog/orchestration work rather than a new physics-code calculation.
- Files touched (`git show --stat`):
```text
 flavor_catalog/worklogs/discovery/round_004_convergence.md | 44 ++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 44 insertions(+)
```

### 249. `56ee949` — flavor-catalog(charged_lepton_top): WA v2 rework wa_w6_charged_lepton_top_v2 — address CA findings
- SHA: `56ee949b211472496e2a8ef4c6b6ba4304b38cb4`
- Message: flavor-catalog(charged_lepton_top): WA v2 rework wa_w6_charged_lepton_top_v2 — address CA findings
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charged_lepton/L023.yaml             | 90 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/processes/top_higgs_ew/T020.yaml               | 22 +++++++++++++++++++++-
 flavor_catalog/worklogs/writer/wa_w6_charged_lepton_top_v2.md | 38 ++++++++++++++++++++++++++++++++++++++
 3 files changed, 148 insertions(+), 2 deletions(-)
```

### 250. `e9d13b7` — flavor-catalog(edm_neutrino): WA v2 rework wa_w6_edm_v2 — address CA findings
- SHA: `e9d13b7c6cc218af30ba2cd26ca9c11e2d0d322a`
- Message: flavor-catalog(edm_neutrino): WA v2 rework wa_w6_edm_v2 — address CA findings
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/edm_neutrino/E009.yaml | 141 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++----------------------------------------
 flavor_catalog/worklogs/writer/wa_w6_edm_v2.md  |  33 +++++++++++++++++++++++++++++++++
 2 files changed, 134 insertions(+), 40 deletions(-)
```

### 251. `4666517` — flavor-catalog(kaon): WA v2 rework wa_w6_kaon_v2 — address CA findings
- SHA: `4666517bf2220ae87666f9e06f6984366335b47f`
- Message: flavor-catalog(kaon): WA v2 rework wa_w6_kaon_v2 — address CA findings
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/kaon/K012.tex          | 14 +++++++-------
 flavor_catalog/processes/kaon/K012.yaml         | 73 +++++++++++++++++++++++++++++++++++++++++++++++++------------------------
 flavor_catalog/processes/kaon/K018.tex          |  8 ++++----
 flavor_catalog/processes/kaon/K018.yaml         | 16 ++++++++++++++--
 flavor_catalog/worklogs/writer/wa_w6_kaon_v2.md | 33 +++++++++++++++++++++++++++++++++
 5 files changed, 107 insertions(+), 37 deletions(-)
```

### 252. `ee5c95b` — flavor-catalog(beauty): WA v2 rework wa_w6_beauty_v2 — address CA findings
- SHA: `ee5c95bd80b9b2f7de51acf7de0de44ab8262128`
- Message: flavor-catalog(beauty): WA v2 rework wa_w6_beauty_v2 — address CA findings
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B001.yaml         |  48 ++++++++++++++++++++++++++++++++++++++++++++----
 flavor_catalog/processes/beauty/B003.yaml         | 105 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++----
 flavor_catalog/processes/beauty/B016.yaml         |  16 +++++++++++++++-
 flavor_catalog/worklogs/writer/wa_w6_beauty_v2.md |  56 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 4 files changed, 216 insertions(+), 9 deletions(-)
```

### 253. `67912c0` — flavor-catalog(edm_neutrino): CA batch ca_w6_edm_v2 — verify WA polish for E009
- SHA: `67912c093142eac476668aba7f38ab771a85333f`
- Message: flavor-catalog(edm_neutrino): CA batch ca_w6_edm_v2 — verify WA polish for E009
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/edm_neutrino/E009.yaml | 16 +++++++++++++---
 flavor_catalog/worklogs/checker/ca_w6_edm_v2.md | 32 ++++++++++++++++++++++++++++++++
 2 files changed, 45 insertions(+), 3 deletions(-)
```

### 254. `cd8816d` — flavor-catalog(beauty): CA batch ca_w6_beauty_v2 — verify WA polish for B001 B003 B016
- SHA: `cd8816d3c1a4b1ad49a99c26a1f5268e6c6ccc80`
- Message: flavor-catalog(beauty): CA batch ca_w6_beauty_v2 — verify WA polish for B001 B003 B016
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B001.yaml          | 14 +++++++++++++-
 flavor_catalog/processes/beauty/B003.yaml          | 13 ++++++++++++-
 flavor_catalog/processes/beauty/B016.yaml          | 15 ++++++++++++---
 flavor_catalog/worklogs/checker/ca_w6_beauty_v2.md | 45 +++++++++++++++++++++++++++++++++++++++++++++
 4 files changed, 82 insertions(+), 5 deletions(-)
```

### 255. `7500794` — flavor-catalog(charged_lepton_top): CA batch ca_w6_charged_lepton_top_v2 — verify WA polish for L023 T020
- SHA: `7500794538191067d8a0e8f837b752b828be8a66`
- Message: flavor-catalog(charged_lepton_top): CA batch ca_w6_charged_lepton_top_v2 — verify WA polish for L023 T020
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/charged_lepton/L023.yaml              | 15 ++++++++++++---
 flavor_catalog/processes/top_higgs_ew/T020.yaml                | 15 ++++++++++++---
 flavor_catalog/worklogs/checker/ca_w6_charged_lepton_top_v2.md | 39 +++++++++++++++++++++++++++++++++++++++
 3 files changed, 63 insertions(+), 6 deletions(-)
```

### 256. `6d97db8` — flavor-catalog(kaon): CA batch ca_w6_kaon_v2 — verify WA polish for K012 K018
- SHA: `6d97db8848f392dc1631f9f57f6292f482c5544c`
- Message: flavor-catalog(kaon): CA batch ca_w6_kaon_v2 — verify WA polish for K012 K018
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/kaon/K012.yaml          | 15 ++++++++++++---
 flavor_catalog/processes/kaon/K018.yaml          | 15 ++++++++++++---
 flavor_catalog/worklogs/checker/ca_w6_kaon_v2.md | 46 ++++++++++++++++++++++++++++++++++++++++++++++
 3 files changed, 70 insertions(+), 6 deletions(-)
```

### 257. `7e1b80b` — flavor-catalog(arbitration): Opus signoff on B001 + B003 (already-in-code) cycle-2 cap
- SHA: `7e1b80bc95471a591f7b594d3aecfc6758f8672c`
- Message: flavor-catalog(arbitration): Opus signoff on B001 + B003 (already-in-code) cycle-2 cap
- Physical/numerical summary: Recorded reviewer/sign-off disposition; review ledger only.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B001.yaml      |  16 +++++++++++++---
 flavor_catalog/processes/beauty/B003.yaml      |  16 +++++++++++++---
 flavor_catalog/signoff/by_process/B001_B003.md | 137 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 3 files changed, 163 insertions(+), 6 deletions(-)
```

### 258. `aa8e3e8` — flavor-catalog(signoff): Opus round 2 sign-off on 21 processes
- SHA: `aa8e3e8eb9e7c2dbd2cd5064b09e29894634ec64`
- Message: flavor-catalog(signoff): Opus round 2 sign-off on 21 processes
- Physical/numerical summary: Recorded reviewer/sign-off disposition (21 processes); review ledger only.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B016.yaml         |  13 +++++++++++--
 flavor_catalog/processes/beauty/B022.yaml         |  13 +++++++++++--
 flavor_catalog/processes/charged_lepton/L005.yaml |  13 +++++++++++--
 flavor_catalog/processes/charged_lepton/L006.yaml |  13 +++++++++++--
 flavor_catalog/processes/charged_lepton/L023.yaml |  13 +++++++++++--
 flavor_catalog/processes/charm/C006.yaml          |  13 +++++++++++--
 flavor_catalog/processes/charm/C008.yaml          |  13 +++++++++++--
 flavor_catalog/processes/edm_neutrino/E002.yaml   |  13 +++++++++++--
 flavor_catalog/processes/edm_neutrino/E007.yaml   |  13 +++++++++++--
 flavor_catalog/processes/edm_neutrino/E009.yaml   |  13 +++++++++++--
 flavor_catalog/processes/kaon/K008.yaml           |  13 +++++++++++--
 flavor_catalog/processes/kaon/K009.yaml           |  13 +++++++++++--
 flavor_catalog/processes/kaon/K010.yaml           |  13 +++++++++++--
 flavor_catalog/processes/kaon/K012.yaml           |  13 +++++++++++--
 flavor_catalog/processes/kaon/K018.yaml           |  13 +++++++++++--
 flavor_catalog/processes/top_higgs_ew/T005.yaml   |  13 +++++++++++--
 flavor_catalog/processes/top_higgs_ew/T006.yaml   |  13 +++++++++++--
 flavor_catalog/processes/top_higgs_ew/T016.yaml   |  13 +++++++++++--
 flavor_catalog/processes/top_higgs_ew/T017.yaml   |  13 +++++++++++--
 flavor_catalog/processes/top_higgs_ew/T019.yaml   |  13 +++++++++++--
 flavor_catalog/processes/top_higgs_ew/T020.yaml   |  13 +++++++++++--
 flavor_catalog/signoff/round_002_index.md         | 122 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 22 files changed, 353 insertions(+), 42 deletions(-)
```

### 259. `9f65578` — flavor-catalog: master compile - 75 OPUS-APPROVED processes, 4 DA rounds, 2 Opus rounds + 5 arbitrations
- SHA: `9f655782cdd183c380d9c9be2b27e92e476b4d74`
- Message: flavor-catalog: master compile - 75 OPUS-APPROVED processes, 4 DA rounds, 2 Opus rounds + 5 arbitrations
- Physical/numerical summary: Updated the flavor-catalog master compile and headline counts (75, 4, 2, 5 cited in subject); catalog aggregation only.
- Files touched (`git show --stat`):
```text
 flavor_catalog/README.md                          |   9 +++++++++
 flavor_catalog/catalog_master.pdf                 | Bin 0 -> 699241 bytes
 flavor_catalog/master_compile_report.md           |  71 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/beauty/index.tex         |  23 ++++++++++++++++++++++-
 flavor_catalog/processes/charged_lepton/index.tex |  13 ++++++++++++-
 flavor_catalog/processes/charm/index.tex          |  10 +++++++++-
 flavor_catalog/processes/edm_neutrino/index.tex   |   9 ++++++++-
 flavor_catalog/processes/kaon/index.tex           |  15 ++++++++++++++-
 flavor_catalog/processes/top_higgs_ew/index.tex   |  17 ++++++++++++++++-
 9 files changed, 161 insertions(+), 6 deletions(-)
```

### 260. `a4245e2` — flavor-catalog(audits): initialize factcheck_status.md master checklist for v0.1
- SHA: `a4245e2f680aaec24bced60f84914a93e323970a`
- Message: flavor-catalog(audits): initialize factcheck_status.md master checklist for v0.1
- Physical/numerical summary: Recorded fact-check results; verification bookkeeping, not executable physics code.
- Files touched (`git show --stat`):
```text
 flavor_catalog/audits/factcheck_status.md | 37 +++++++++++++++++++++++++++++++++++++
 1 file changed, 37 insertions(+)
```

### 261. `9fb5763` — flavor-catalog(audits): factcheck edm_neutrino family - 6 verified / 0 mismatch / 1 unresolvable
- SHA: `9fb576397d5f68cd67c707812c2a46e32b68392d`
- Message: flavor-catalog(audits): factcheck edm_neutrino family - 6 verified / 0 mismatch / 1 unresolvable
- Physical/numerical summary: Recorded fact-check results (6 verified; 0 mismatch; 1 unresolvable); verification bookkeeping, not executable physics code.
- Files touched (`git show --stat`):
```text
 flavor_catalog/audits/factcheck_edm_neutrino.md | 128 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/edm_neutrino/E001.yaml |  14 +++++++++++++-
 flavor_catalog/processes/edm_neutrino/E002.yaml |  14 +++++++++++++-
 flavor_catalog/processes/edm_neutrino/E004.yaml |  14 +++++++++++++-
 flavor_catalog/processes/edm_neutrino/E006.yaml |  14 +++++++++++++-
 flavor_catalog/processes/edm_neutrino/E007.yaml |  14 +++++++++++++-
 flavor_catalog/processes/edm_neutrino/E008.yaml |  14 +++++++++++++-
 flavor_catalog/processes/edm_neutrino/E009.yaml |  14 +++++++++++++-
 8 files changed, 219 insertions(+), 7 deletions(-)
```

### 262. `a04a101` — flavor-catalog(audits): factcheck kaon family - 13 verified / 0 mismatch / 0 unresolvable
- SHA: `a04a101a42a30643b5c7940892f508eccfd2a005`
- Message: flavor-catalog(audits): factcheck kaon family - 13 verified / 0 mismatch / 0 unresolvable
- Physical/numerical summary: Recorded fact-check results (13 verified; 0 mismatch; 0 unresolvable); verification bookkeeping, not executable physics code.
- Files touched (`git show --stat`):
```text
 flavor_catalog/audits/factcheck_kaon.md | 250 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/kaon/K001.yaml |  16 ++++++++++-
 flavor_catalog/processes/kaon/K002.yaml |  16 ++++++++++-
 flavor_catalog/processes/kaon/K003.yaml |  16 ++++++++++-
 flavor_catalog/processes/kaon/K004.yaml |  16 ++++++++++-
 flavor_catalog/processes/kaon/K005.yaml |  16 ++++++++++-
 flavor_catalog/processes/kaon/K006.yaml |  16 ++++++++++-
 flavor_catalog/processes/kaon/K008.yaml |  16 ++++++++++-
 flavor_catalog/processes/kaon/K009.yaml |  16 ++++++++++-
 flavor_catalog/processes/kaon/K010.yaml |  16 ++++++++++-
 flavor_catalog/processes/kaon/K012.yaml |  16 ++++++++++-
 flavor_catalog/processes/kaon/K013.yaml |  16 ++++++++++-
 flavor_catalog/processes/kaon/K017.yaml |  16 ++++++++++-
 flavor_catalog/processes/kaon/K018.yaml |  16 ++++++++++-
 14 files changed, 445 insertions(+), 13 deletions(-)
```

### 263. `7f79bab` — flavor-catalog(audits): factcheck charged_lepton family - 11 verified / 0 mismatch / 0 unresolvable
- SHA: `7f79bab9c607e8c1a84979be385d9e211d4d1412`
- Message: flavor-catalog(audits): factcheck charged_lepton family - 11 verified / 0 mismatch / 0 unresolvable
- Physical/numerical summary: Recorded fact-check results (11 verified; 0 mismatch; 0 unresolvable); verification bookkeeping, not executable physics code.
- Files touched (`git show --stat`):
```text
 flavor_catalog/audits/factcheck_charged_lepton.md | 181 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/charged_lepton/L001.yaml |  10 ++++++++-
 flavor_catalog/processes/charged_lepton/L002.yaml |  10 ++++++++-
 flavor_catalog/processes/charged_lepton/L003.yaml |  10 ++++++++-
 flavor_catalog/processes/charged_lepton/L004.yaml |  10 ++++++++-
 flavor_catalog/processes/charged_lepton/L005.yaml |  10 ++++++++-
 flavor_catalog/processes/charged_lepton/L006.yaml |  10 ++++++++-
 flavor_catalog/processes/charged_lepton/L007.yaml |  10 ++++++++-
 flavor_catalog/processes/charged_lepton/L008.yaml |  10 ++++++++-
 flavor_catalog/processes/charged_lepton/L009.yaml |  10 ++++++++-
 flavor_catalog/processes/charged_lepton/L010.yaml |  10 ++++++++-
 flavor_catalog/processes/charged_lepton/L023.yaml |  10 ++++++++-
 12 files changed, 280 insertions(+), 11 deletions(-)
```

### 264. `0b07505` — flavor-catalog(audits): factcheck top_higgs_ew family — 14 verified / 1 mismatch / 0 unresolvable
- SHA: `0b075052ba002dd983a926f91d637e4af8901198`
- Message: flavor-catalog(audits): factcheck top_higgs_ew family — 14 verified / 1 mismatch / 0 unresolvable
- Physical/numerical summary: Recorded fact-check results (14 verified; 1 mismatch; 0 unresolvable); verification bookkeeping, not executable physics code.
- Files touched (`git show --stat`):
```text
 flavor_catalog/audits/factcheck_top_higgs_ew.md  | 293 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/top_higgs_ew/EW001.yaml |  14 +++++++-
 flavor_catalog/processes/top_higgs_ew/EW002.yaml |  14 +++++++-
 flavor_catalog/processes/top_higgs_ew/EW003.yaml |  14 +++++++-
 flavor_catalog/processes/top_higgs_ew/T001.yaml  |  14 +++++++-
 flavor_catalog/processes/top_higgs_ew/T002.yaml  |  14 +++++++-
 flavor_catalog/processes/top_higgs_ew/T005.yaml  |  14 +++++++-
 flavor_catalog/processes/top_higgs_ew/T006.yaml  |  14 +++++++-
 flavor_catalog/processes/top_higgs_ew/T007.yaml  |  14 +++++++-
 flavor_catalog/processes/top_higgs_ew/T010.yaml  |  14 +++++++-
 flavor_catalog/processes/top_higgs_ew/T015.yaml  |  14 +++++++-
 flavor_catalog/processes/top_higgs_ew/T016.yaml  |  14 +++++++-
 flavor_catalog/processes/top_higgs_ew/T017.yaml  |  14 +++++++-
 flavor_catalog/processes/top_higgs_ew/T018.yaml  |  14 +++++++-
 flavor_catalog/processes/top_higgs_ew/T019.yaml  |  14 +++++++-
 flavor_catalog/processes/top_higgs_ew/T020.yaml  |  14 +++++++-
 16 files changed, 488 insertions(+), 15 deletions(-)
```

### 265. `04b1ab1` — flavor-catalog(audits): factcheck charm family - 8 verified / 0 mismatch / 0 unresolvable
- SHA: `04b1ab119c8db7c14dc3b3fe46aa0262b51693f1`
- Message: flavor-catalog(audits): factcheck charm family - 8 verified / 0 mismatch / 0 unresolvable
- Physical/numerical summary: Recorded fact-check results (8 verified; 0 mismatch; 0 unresolvable); verification bookkeeping, not executable physics code.
- Files touched (`git show --stat`):
```text
 flavor_catalog/audits/factcheck_charm.md | 135 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/charm/C001.yaml |  14 +++++++++++++-
 flavor_catalog/processes/charm/C002.yaml |  14 +++++++++++++-
 flavor_catalog/processes/charm/C003.yaml |  14 +++++++++++++-
 flavor_catalog/processes/charm/C004.yaml |  14 +++++++++++++-
 flavor_catalog/processes/charm/C005.yaml |  14 +++++++++++++-
 flavor_catalog/processes/charm/C006.yaml |  14 +++++++++++++-
 flavor_catalog/processes/charm/C007.yaml |  14 +++++++++++++-
 flavor_catalog/processes/charm/C008.yaml |  14 +++++++++++++-
 9 files changed, 239 insertions(+), 8 deletions(-)
```

### 266. `994e346` — flavor-catalog(audits): factcheck beauty family - 21 verified / 0 mismatch / 0 unresolvable
- SHA: `994e346676c09b38f0085edbf8d2be4de74696f4`
- Message: flavor-catalog(audits): factcheck beauty family - 21 verified / 0 mismatch / 0 unresolvable
- Physical/numerical summary: Recorded fact-check results (21 verified; 0 mismatch; 0 unresolvable); verification bookkeeping, not executable physics code.
- Files touched (`git show --stat`):
```text
 flavor_catalog/audits/factcheck_beauty.md | 223 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/beauty/B001.yaml |  10 +++++++-
 flavor_catalog/processes/beauty/B002.yaml |  10 +++++++-
 flavor_catalog/processes/beauty/B003.yaml |  10 +++++++-
 flavor_catalog/processes/beauty/B004.yaml |  10 +++++++-
 flavor_catalog/processes/beauty/B005.yaml |  10 +++++++-
 flavor_catalog/processes/beauty/B006.yaml |  10 +++++++-
 flavor_catalog/processes/beauty/B009.yaml |  10 +++++++-
 flavor_catalog/processes/beauty/B011.yaml |  10 +++++++-
 flavor_catalog/processes/beauty/B015.yaml |  10 +++++++-
 flavor_catalog/processes/beauty/B016.yaml |  10 +++++++-
 flavor_catalog/processes/beauty/B017.yaml |  10 +++++++-
 flavor_catalog/processes/beauty/B018.yaml |  10 +++++++-
 flavor_catalog/processes/beauty/B019.yaml |  10 +++++++-
 flavor_catalog/processes/beauty/B021.yaml |  10 +++++++-
 flavor_catalog/processes/beauty/B022.yaml |  10 +++++++-
 flavor_catalog/processes/beauty/B023.yaml |  10 +++++++-
 flavor_catalog/processes/beauty/B025.yaml |  10 +++++++-
 flavor_catalog/processes/beauty/B026.yaml |  10 +++++++-
 flavor_catalog/processes/beauty/B032.yaml |  10 +++++++-
 flavor_catalog/processes/beauty/B033.yaml |  10 +++++++-
 flavor_catalog/processes/beauty/B034.yaml |  10 +++++++-
 22 files changed, 412 insertions(+), 21 deletions(-)
```

### 267. `6498fad` — flavor-catalog(audits): fix T020 ATLAS h to e mu limit digit (6.1e-5 -> 6.2e-5)
- SHA: `6498fadd98e75066f42d2ab627ecd5da2402c5b7`
- Message: flavor-catalog(audits): fix T020 ATLAS h to e mu limit digit (6.1e-5 -> 6.2e-5)
- Physical/numerical summary: Updated catalog audit/fact-check ledgers; physical impact is verification/provenance, not new calculations.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T020.tex  |  9 +++++----
 flavor_catalog/processes/top_higgs_ew/T020.yaml | 19 ++++++++++++-------
 2 files changed, 17 insertions(+), 11 deletions(-)
```

### 268. `42ac647` — flavor-catalog(audits): consolidated factcheck_status.md (75 processes; 73-75 VERIFIED)
- SHA: `42ac647dcbaafc070279afca403da9f642aa6e60`
- Message: flavor-catalog(audits): consolidated factcheck_status.md (75 processes; 73-75 VERIFIED)
- Physical/numerical summary: Recorded fact-check results (75 VERIFIED); verification bookkeeping, not executable physics code.
- Files touched (`git show --stat`):
```text
 flavor_catalog/audits/factcheck_status.md | 124 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++---
 1 file changed, 121 insertions(+), 3 deletions(-)
```

### 269. `022a20c` — flavor-catalog: import external GPT Deep Research PDFs (may15, may16) on RS flavor constraints
- SHA: `022a20c3c37d46d9480bf3b472ed5fe3420e07a5`
- Message: flavor-catalog: import external GPT Deep Research PDFs (may15, may16) on RS flavor constraints
- Physical/numerical summary: Imported external Deep Research PDFs/texts for RS flavor constraints as literature/context artifacts.
- Files touched (`git show --stat`):
```text
 flavor_catalog/external_research/deepresearch_may15.pdf |  Bin 0 -> 203806 bytes
 flavor_catalog/external_research/deepresearch_may15.txt | 1109 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/external_research/deepresearch_may16.pdf |  Bin 0 -> 343307 bytes
 flavor_catalog/external_research/deepresearch_may16.txt | 2092 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 4 files changed, 3201 insertions(+)
```

### 270. `8fb5f91` — flavor-catalog(external-review): deepresearch_may16 comparison against catalog v0.1
- SHA: `8fb5f91be63952176f3bf20253b43dc407704cb3`
- Message: flavor-catalog(external-review): deepresearch_may16 comparison against catalog v0.1
- Physical/numerical summary: Imported or compared external RS flavor research artifacts, adding reference material for catalog coverage decisions.
- Files touched (`git show --stat`):
```text
 flavor_catalog/external_research/deepresearch_may16_review.md | 61 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 61 insertions(+)
```

### 271. `2c00d84` — flavor-catalog(external-review): deepresearch_may15 comparison against catalog v0.1
- SHA: `2c00d8446990dd9ddaadd2652ef3b27ab65cfcdb`
- Message: flavor-catalog(external-review): deepresearch_may15 comparison against catalog v0.1
- Physical/numerical summary: Imported or compared external RS flavor research artifacts, adding reference material for catalog coverage decisions.
- Files touched (`git show --stat`):
```text
 flavor_catalog/external_research/deepresearch_may15_review.md | 208 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 208 insertions(+)
```

### 272. `dbf2b96` — flavor-catalog(wave7): DA-4 addendum closing plan-v1 Section C bookkeeping for deferred-scope items
- SHA: `dbf2b960d3854804385afadb53bef2796a48e7ff`
- Message: flavor-catalog(wave7): DA-4 addendum closing plan-v1 Section C bookkeeping for deferred-scope items
- Physical/numerical summary: Applied the change described by the subject; based on stat surface, this is documentation/catalog/orchestration work rather than a new physics-code calculation.
- Files touched (`git show --stat`):
```text
 flavor_catalog/worklogs/discovery/round_004_addendum_deferred_scope.md | 74 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 74 insertions(+)
```

### 273. `837e709` — flavor-catalog(wave7): thread 6 cross-cutting subtleties from external review into existing process entries
- SHA: `837e709be1e3764eaca02a7655b49e9b81b026e6`
- Message: flavor-catalog(wave7): thread 6 cross-cutting subtleties from external review into existing process entries
- Physical/numerical summary: Applied the change described by the subject; based on stat surface, this is documentation/catalog/orchestration work rather than a new physics-code calculation.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B011.tex           |   4 ++++
 flavor_catalog/processes/beauty/B011.yaml          |   4 ++++
 flavor_catalog/processes/beauty/B015.tex           |   4 ++++
 flavor_catalog/processes/beauty/B015.yaml          |   4 ++++
 flavor_catalog/processes/beauty/B017.tex           |   6 ++++++
 flavor_catalog/processes/beauty/B017.yaml          |   4 ++++
 flavor_catalog/processes/beauty/B018.tex           |   6 ++++++
 flavor_catalog/processes/beauty/B018.yaml          |   4 ++++
 flavor_catalog/processes/beauty/B019.tex           |   6 ++++++
 flavor_catalog/processes/beauty/B019.yaml          |   4 ++++
 flavor_catalog/processes/charm/C001.tex            |   3 +++
 flavor_catalog/processes/charm/C001.yaml           |   4 ++++
 flavor_catalog/processes/charm/C002.tex            |   3 +++
 flavor_catalog/processes/charm/C002.yaml           |   4 ++++
 flavor_catalog/processes/charm/C003.tex            |   3 +++
 flavor_catalog/processes/charm/C003.yaml           |   4 ++++
 flavor_catalog/processes/charm/C004.tex            |   3 +++
 flavor_catalog/processes/charm/C004.yaml           |   4 ++++
 flavor_catalog/processes/edm_neutrino/E001.tex     |   4 ++++
 flavor_catalog/processes/edm_neutrino/E001.yaml    |   4 ++++
 flavor_catalog/processes/edm_neutrino/E004.tex     |   4 ++++
 flavor_catalog/processes/edm_neutrino/E004.yaml    |   4 ++++
 flavor_catalog/processes/edm_neutrino/E006.tex     |   4 ++++
 flavor_catalog/processes/edm_neutrino/E006.yaml    |   4 ++++
 flavor_catalog/processes/edm_neutrino/E008.tex     |   4 ++++
 flavor_catalog/processes/edm_neutrino/E008.yaml    |   4 ++++
 flavor_catalog/processes/edm_neutrino/E009.tex     |   4 ++++
 flavor_catalog/processes/edm_neutrino/E009.yaml    |   4 ++++
 flavor_catalog/processes/top_higgs_ew/EW001.tex    |   4 ++++
 flavor_catalog/processes/top_higgs_ew/EW001.yaml   |   4 ++++
 flavor_catalog/processes/top_higgs_ew/EW002.tex    |   5 +++++
 flavor_catalog/processes/top_higgs_ew/EW002.yaml   |   4 ++++
 flavor_catalog/processes/top_higgs_ew/T001.tex     |   3 +++
 flavor_catalog/processes/top_higgs_ew/T001.yaml    |   4 ++++
 flavor_catalog/processes/top_higgs_ew/T002.tex     |   3 +++
 flavor_catalog/processes/top_higgs_ew/T002.yaml    |   4 ++++
 flavor_catalog/processes/top_higgs_ew/T005.tex     |   3 +++
 flavor_catalog/processes/top_higgs_ew/T005.yaml    |   4 ++++
 flavor_catalog/processes/top_higgs_ew/T006.tex     |   3 +++
 flavor_catalog/processes/top_higgs_ew/T006.yaml    |   4 ++++
 flavor_catalog/processes/top_higgs_ew/T007.tex     |   5 +++++
 flavor_catalog/processes/top_higgs_ew/T007.yaml    |   4 ++++
 flavor_catalog/processes/top_higgs_ew/T010.tex     |   5 +++++
 flavor_catalog/processes/top_higgs_ew/T010.yaml    |   4 ++++
 flavor_catalog/worklogs/writer/wave7_subtleties.md | 100 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 45 files changed, 277 insertions(+)
```

### 274. `e3ecf8b` — flavor-catalog(top_higgs_ew): PKA draft for T012 Z to ccbar
- SHA: `e3ecf8bcebd97b91a46ccc4c54e156543081086a`
- Message: flavor-catalog(top_higgs_ew): PKA draft for T012 Z to ccbar
- Physical/numerical summary: Added an initial catalog process draft for T012 Z to ccbar; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T012.tex                   |  91 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/top_higgs_ew/T012.yaml                  | 108 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/T012/alcaraz_2021_z_lineshape_fcc.txt  |  18 ++++++++++++++++++
 flavor_catalog/references/T012/casagrande_2008_rs_ewpt.txt       |  19 +++++++++++++++++++
 flavor_catalog/references/T012/cfw_2008_rs_flavor.txt            |  20 ++++++++++++++++++++
 flavor_catalog/references/T012/freitas_2014_z_widths.txt         |  17 +++++++++++++++++
 flavor_catalog/references/T012/lepslc_2006_z_resonance_charm.txt |  20 ++++++++++++++++++++
 flavor_catalog/references/T012/pdg_2025_z_boson_charm.txt        |  23 +++++++++++++++++++++++
 flavor_catalog/references/T012/source_manifest.yaml              |  64 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/T012.md                              |  34 ++++++++++++++++++++++++++++++++++
 10 files changed, 414 insertions(+)
```

### 275. `55e2299` — flavor-catalog(top_higgs_ew): PKA draft for T008 t to H u
- SHA: `55e229918d13016f57219c1fa1768aebf2b8d654`
- Message: flavor-catalog(top_higgs_ew): PKA draft for T008 t to H u
- Physical/numerical summary: Added an initial catalog process draft for T008 t to H u; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T008.tex                |  82 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/top_higgs_ew/T008.yaml               | 140 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/T008/atlas_2404_02123_arxiv_abs.txt |  27 +++++++++++++++++++++++++++
 flavor_catalog/references/T008/cfw_0804_1954_arxiv_abs.txt    |  19 +++++++++++++++++++
 flavor_catalog/references/T008/cms_2111_02219_arxiv_abs.txt   |  23 +++++++++++++++++++++++
 flavor_catalog/references/T008/cms_2407_15172_arxiv_abs.txt   |  25 +++++++++++++++++++++++++
 flavor_catalog/references/T008/pdg2026_top_hu_datablock.txt   |  55 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/T008/source_manifest.yaml           |  49 +++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/T008.md                           |  33 +++++++++++++++++++++++++++++++++
 9 files changed, 453 insertions(+)
```

### 276. `64d3f6a` — flavor-catalog(top_higgs_ew): PKA draft for T004 t to u gamma
- SHA: `64d3f6a3404bd6a3e44bc0b75bc5dcccdcec89f4`
- Message: flavor-catalog(top_higgs_ew): PKA draft for T004 t to u gamma
- Physical/numerical summary: Added an initial catalog process draft for T004 t to u gamma; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T004.tex                                |  86 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/top_higgs_ew/T004.yaml                               | 124 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/T004/aguilar_saavedra_1709_03975_arxiv_abs.txt      |  20 ++++++++++++++++++++
 flavor_catalog/references/T004/atlas_2205_02537_arxiv_abs.txt                 |  22 ++++++++++++++++++++++
 flavor_catalog/references/T004/cms_2312_08229_public_page.txt                 |  20 ++++++++++++++++++++
 flavor_catalog/references/T004/csaki_falkowski_weiler_0804_1954_arxiv_abs.txt |  19 +++++++++++++++++++
 flavor_catalog/references/T004/pdg_2026_top_t_gammaq.txt                      |  15 +++++++++++++++
 flavor_catalog/references/T004/source_manifest.yaml                           |  49 +++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/T004.md                                           |  33 +++++++++++++++++++++++++++++++++
 9 files changed, 388 insertions(+)
```

### 277. `1f9e16a` — flavor-catalog(top_higgs_ew): PKA draft for T003 t to c gamma
- SHA: `1f9e16a09f773a87530882908d59da691bfc8b89`
- Message: flavor-catalog(top_higgs_ew): PKA draft for T003 t to c gamma
- Physical/numerical summary: Added an initial catalog process draft for T003 t to c gamma; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T003.tex                              |  91 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/top_higgs_ew/T003.yaml                             | 152 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/T003/aguilar_saavedra_hepph0409342_tcgamma_sm.txt |  35 ++++++++++++++++++++++++++++++++
 flavor_catalog/references/T003/atlas_2023_2205_02537_tqgamma.txt            |  35 ++++++++++++++++++++++++++++++++
 flavor_catalog/references/T003/cfw_0804_1954_rs_flavor.txt                  |  27 ++++++++++++++++++++++++
 flavor_catalog/references/T003/cms_2024_top21013_tqgamma.txt                |  34 +++++++++++++++++++++++++++++++
 flavor_catalog/references/T003/pdg_2026_top_gammaq_summary.txt              |  23 +++++++++++++++++++++
 flavor_catalog/references/T003/source_manifest.yaml                         |  49 ++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/T003.md                                         |  32 +++++++++++++++++++++++++++++
 9 files changed, 478 insertions(+)
```

### 278. `60e85b9` — flavor-catalog(beauty): PKA draft for B012 B to Kstar gamma
- SHA: `60e85b94c529606bceb5c923f441ac2550c34243`
- Message: flavor-catalog(beauty): PKA draft for B012 B to Kstar gamma
- Physical/numerical summary: Added an initial catalog process draft for B012 B to Kstar gamma; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B012.tex                                    |  91 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/beauty/B012.yaml                                   | 139 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B012/arxiv_0804_1954_cfw_rs_flavor.txt            |  17 +++++++++++++++++
 flavor_catalog/references/B012/arxiv_0906_2177_babar_btokstargamma.txt      |  20 ++++++++++++++++++++
 flavor_catalog/references/B012/arxiv_1402_6852_lhcb_photon_polarization.txt |  19 +++++++++++++++++++
 flavor_catalog/references/B012/arxiv_1707_00394_belle_isospin_cp.txt        |  19 +++++++++++++++++++
 flavor_catalog/references/B012/arxiv_2411_10127_belleii_btokstargamma.txt   |  21 +++++++++++++++++++++
 flavor_catalog/references/B012/hflav_rare_decays_dec2024_btokstargamma.txt  |  53 ++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B012/pdglive_2025_kstargamma_cp.txt               |  26 ++++++++++++++++++++++++++
 flavor_catalog/references/B012/sha256sums.txt                               |   7 +++++++
 flavor_catalog/references/B012/source_manifest.yaml                         |  67 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/B012.md                                         |  35 ++++++++++++++++++++++++++++++++++
 12 files changed, 514 insertions(+)
```

### 279. `99003f3` — flavor-catalog(wave7): peer review of DA-4 deferred-scope addendum
- SHA: `99003f3e167920cd3ee82ef0eea816b8da8f62c0`
- Message: flavor-catalog(wave7): peer review of DA-4 deferred-scope addendum
- Physical/numerical summary: Applied the change described by the subject; based on stat surface, this is documentation/catalog/orchestration work rather than a new physics-code calculation.
- Files touched (`git show --stat`):
```text
 flavor_catalog/audits/wave7_deferred_scope_review.md | 25 +++++++++++++++++++++++++
 1 file changed, 25 insertions(+)
```

### 280. `40248b8` — flavor-catalog(wave7): peer review of subtlety writer (6 phrasings on ~22 .tex)
- SHA: `40248b80dbc1333eaafe20d971d649b8c3330bba`
- Message: flavor-catalog(wave7): peer review of subtlety writer (6 phrasings on ~22 .tex)
- Physical/numerical summary: Applied the change described by the subject; based on stat surface, this is documentation/catalog/orchestration work rather than a new physics-code calculation.
- Files touched (`git show --stat`):
```text
 flavor_catalog/audits/wave7_subtleties_review.md | 70 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 70 insertions(+)
```

### 281. `5f21f3e` — flavor-catalog(wave7): T007 subtlety micro-fix (125 GeV -> observed)
- SHA: `5f21f3ebfca3856d68db6949f7d4957a14c09c30`
- Message: flavor-catalog(wave7): T007 subtlety micro-fix (125 GeV -> observed)
- Physical/numerical summary: Applied the change described by the subject; based on stat surface, this is documentation/catalog/orchestration work rather than a new physics-code calculation.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T007.tex  | 2 +-
 flavor_catalog/processes/top_higgs_ew/T007.yaml | 4 ++++
 2 files changed, 5 insertions(+), 1 deletion(-)
```

### 282. `52f6660` — flavor-catalog(mixed): WA batch wa_w7_new_processes — polish PKA drafts for T003 T004 T008 T012 B012
- SHA: `52f66603e94ea802f64ab49a9043ccb5722269db`
- Message: flavor-catalog(mixed): WA batch wa_w7_new_processes — polish PKA drafts for T003 T004 T008 T012 B012
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B012.tex              | 39 ++++++++++++++++++++++++---------------
 flavor_catalog/processes/beauty/B012.yaml             | 11 ++++++++++-
 flavor_catalog/processes/top_higgs_ew/T003.tex        | 49 +++++++++++++++++++++++++++----------------------
 flavor_catalog/processes/top_higgs_ew/T003.yaml       | 11 ++++++++++-
 flavor_catalog/processes/top_higgs_ew/T004.tex        | 23 ++++++++++++-----------
 flavor_catalog/processes/top_higgs_ew/T004.yaml       | 11 ++++++++++-
 flavor_catalog/processes/top_higgs_ew/T008.tex        | 31 ++++++++++++++++++-------------
 flavor_catalog/processes/top_higgs_ew/T008.yaml       | 11 ++++++++++-
 flavor_catalog/processes/top_higgs_ew/T012.tex        | 44 +++++++++++++++++++++++++-------------------
 flavor_catalog/processes/top_higgs_ew/T012.yaml       | 11 ++++++++++-
 flavor_catalog/worklogs/writer/wa_w7_new_processes.md | 72 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 11 files changed, 228 insertions(+), 85 deletions(-)
```

### 283. `3c21e30` — flavor-catalog(mixed): CA batch ca_w7_new_processes - verify WA polish for T003 T004 T008 T012 B012
- SHA: `3c21e304d0ea79adabbce5ee55394360fc3d5d4c`
- Message: flavor-catalog(mixed): CA batch ca_w7_new_processes - verify WA polish for T003 T004 T008 T012 B012
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B012.yaml              | 14 +++++++++++++-
 flavor_catalog/processes/top_higgs_ew/T003.yaml        | 13 ++++++++++++-
 flavor_catalog/processes/top_higgs_ew/T004.yaml        | 14 +++++++++++++-
 flavor_catalog/processes/top_higgs_ew/T008.yaml        | 13 +++++++++++--
 flavor_catalog/processes/top_higgs_ew/T012.yaml        | 13 +++++++++++--
 flavor_catalog/worklogs/checker/ca_w7_new_processes.md | 54 ++++++++++++++++++++++++++++++++++++++++++++++++++++++
 6 files changed, 114 insertions(+), 7 deletions(-)
```

### 284. `3f46ca0` — flavor-catalog(mixed): WA v2 rework wa_w7_new_v2 — address CA findings
- SHA: `3f46ca0c3d4e77ce58bf7130a97c0672d729b692`
- Message: flavor-catalog(mixed): WA v2 rework wa_w7_new_v2 — address CA findings
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B012.tex        |  8 ++++----
 flavor_catalog/processes/beauty/B012.yaml       | 40 ++++++++++++++++++++++++++++------------
 flavor_catalog/processes/top_higgs_ew/T003.tex  |  8 +++-----
 flavor_catalog/processes/top_higgs_ew/T003.yaml | 12 +++++++++++-
 flavor_catalog/processes/top_higgs_ew/T004.tex  | 11 +++--------
 flavor_catalog/processes/top_higgs_ew/T004.yaml | 12 +++++++++++-
 flavor_catalog/worklogs/writer/wa_w7_new_v2.md  | 43 +++++++++++++++++++++++++++++++++++++++++++
 7 files changed, 103 insertions(+), 31 deletions(-)
```

### 285. `91358f9` — flavor-catalog(mixed): CA batch ca_w7_new_v2 - verify WA polish for T003 T004 B012
- SHA: `91358f9b326ec29be35cb16f380a6995c103a6aa`
- Message: flavor-catalog(mixed): CA batch ca_w7_new_v2 - verify WA polish for T003 T004 B012
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B012.yaml       | 16 +++++++++++++---
 flavor_catalog/processes/top_higgs_ew/T003.yaml | 14 +++++++++++++-
 flavor_catalog/processes/top_higgs_ew/T004.yaml | 16 +++++++++++++---
 flavor_catalog/worklogs/checker/ca_w7_new_v2.md | 45 +++++++++++++++++++++++++++++++++++++++++++++
 4 files changed, 84 insertions(+), 7 deletions(-)
```

### 286. `0cfc3ed` — flavor-catalog(top_higgs_ew): WA v2 rework wa_w7_T003_v3 — address CA findings
- SHA: `0cfc3ed3dbc2966f146dbdc597c5131b15dbfe54`
- Message: flavor-catalog(top_higgs_ew): WA v2 rework wa_w7_T003_v3 — address CA findings
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T003.tex  |  6 +++---
 flavor_catalog/processes/top_higgs_ew/T003.yaml | 12 +++++++++++-
 flavor_catalog/worklogs/writer/wa_w7_T003_v3.md | 25 +++++++++++++++++++++++++
 3 files changed, 39 insertions(+), 4 deletions(-)
```

### 287. `a6ffc9f` — flavor-catalog(top_higgs_ew): CA batch ca_w7_T003_v3 — verify WA polish for T003
- SHA: `a6ffc9fae717e631cfcd4644977c2d7e3b878bc9`
- Message: flavor-catalog(top_higgs_ew): CA batch ca_w7_T003_v3 — verify WA polish for T003
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/top_higgs_ew/T003.yaml  | 16 +++++++++++++---
 flavor_catalog/worklogs/checker/ca_w7_T003_v3.md | 33 +++++++++++++++++++++++++++++++++
 2 files changed, 46 insertions(+), 3 deletions(-)
```

### 288. `ca2beba` — flavor-catalog(audits): factcheck beauty family — 1 verified / 0 mismatch / 0 unresolvable
- SHA: `ca2beba360d8030cdb15dfa1411fe2e9f485d38e`
- Message: flavor-catalog(audits): factcheck beauty family — 1 verified / 0 mismatch / 0 unresolvable
- Physical/numerical summary: Recorded fact-check results (1 verified; 0 mismatch; 0 unresolvable); verification bookkeeping, not executable physics code.
- Files touched (`git show --stat`):
```text
 flavor_catalog/audits/factcheck_beauty.md | 22 +++++++++++++++++++---
 flavor_catalog/processes/beauty/B012.yaml | 14 +++++++++++++-
 2 files changed, 32 insertions(+), 4 deletions(-)
```

### 289. `7a6bd34` — flavor-catalog(audits): factcheck top_higgs_ew family — 4 verified / 0 mismatch / 0 unresolvable
- SHA: `7a6bd342f73cacb131bb8534ca4d9ffe76af4ca7`
- Message: flavor-catalog(audits): factcheck top_higgs_ew family — 4 verified / 0 mismatch / 0 unresolvable
- Physical/numerical summary: Recorded fact-check results (4 verified; 0 mismatch; 0 unresolvable); verification bookkeeping, not executable physics code.
- Files touched (`git show --stat`):
```text
 flavor_catalog/audits/factcheck_top_higgs_ew.md | 94 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/top_higgs_ew/T003.yaml | 14 +++++++++++++-
 flavor_catalog/processes/top_higgs_ew/T004.yaml | 14 +++++++++++++-
 flavor_catalog/processes/top_higgs_ew/T008.yaml | 14 +++++++++++++-
 flavor_catalog/processes/top_higgs_ew/T012.yaml | 14 +++++++++++++-
 5 files changed, 146 insertions(+), 4 deletions(-)
```

### 290. `1aa6b37` — flavor-catalog(signoff): Opus round 3 sign-off on Wave-7 (5 new + subtleties + deferred-scope addendum)
- SHA: `1aa6b379c28e22c23db43398f886761079961f4e`
- Message: flavor-catalog(signoff): Opus round 3 sign-off on Wave-7 (5 new + subtleties + deferred-scope addendum)
- Physical/numerical summary: Recorded reviewer/sign-off disposition (5 new); review ledger only.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/beauty/B012.yaml       | 13 +++++++++++--
 flavor_catalog/processes/top_higgs_ew/T003.yaml | 13 +++++++++++--
 flavor_catalog/processes/top_higgs_ew/T004.yaml | 13 +++++++++++--
 flavor_catalog/processes/top_higgs_ew/T008.yaml | 13 +++++++++++--
 flavor_catalog/processes/top_higgs_ew/T012.yaml | 13 +++++++++++--
 flavor_catalog/signoff/round_003_index.md       | 92 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 6 files changed, 147 insertions(+), 10 deletions(-)
```

### 291. `835cf48` — flavor-catalog: master compile v0.2 - 80 OPUS-APPROVED processes (75 + 5 Wave-7 + subtleties + addendum)
- SHA: `835cf48873e5addf89e58a8aec11cac8d9ee7837`
- Message: flavor-catalog: master compile v0.2 - 80 OPUS-APPROVED processes (75 + 5 Wave-7 + subtleties + addendum)
- Physical/numerical summary: Updated the flavor-catalog master compile and headline counts (2, 80, 75, 5, 7 cited in subject); catalog aggregation only.
- Files touched (`git show --stat`):
```text
 flavor_catalog/README.md                          |   8 +++++---
 flavor_catalog/catalog_master.pdf                 | Bin 699241 -> 733862 bytes
 flavor_catalog/master_compile_v02_report.md       |  50 ++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/beauty/index.tex         |  45 +++++++++++++++++++++++----------------------
 flavor_catalog/processes/charged_lepton/index.tex |  24 ++++++++++++------------
 flavor_catalog/processes/charm/index.tex          |  18 +++++++++---------
 flavor_catalog/processes/edm_neutrino/index.tex   |  16 ++++++++--------
 flavor_catalog/processes/kaon/index.tex           |  28 ++++++++++++++--------------
 flavor_catalog/processes/top_higgs_ew/index.tex   |  36 ++++++++++++++++++++----------------
 9 files changed, 141 insertions(+), 84 deletions(-)
```

### 292. `5086c99` — flavor-catalog: add AGENTIC_WORKFLOW.md (reproducible playbook) + SESSION_NOTES.md (handoff)
- SHA: `5086c9993bfc1a3ff224dac9320091ab55d43adc`
- Message: flavor-catalog: add AGENTIC_WORKFLOW.md (reproducible playbook) + SESSION_NOTES.md (handoff)
- Physical/numerical summary: Updated catalog orchestration/methodology handoff documentation; no numerical data or executable physics changed.
- Files touched (`git show --stat`):
```text
 flavor_catalog/AGENTIC_WORKFLOW.md | 505 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/SESSION_NOTES.md    | 302 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 2 files changed, 807 insertions(+)
```

### 293. `5a6fae0` — flavor-catalog(handoff): add Wave-8 candidate tiering to SESSION_NOTES (top/middle/tail tiers from DA-4 addendum)
- SHA: `5a6fae0d72b449a11f436ce72361e1d376c1b069`
- Message: flavor-catalog(handoff): add Wave-8 candidate tiering to SESSION_NOTES (top/middle/tail tiers from DA-4 addendum)
- Physical/numerical summary: Updated catalog orchestration/methodology handoff documentation; no numerical data or executable physics changed.
- Files touched (`git show --stat`):
```text
 flavor_catalog/SESSION_NOTES.md | 39 +++++++++++++++++++++++++++++++++++----
 1 file changed, 35 insertions(+), 4 deletions(-)
```

### 294. `3be3a7b` — flavor-catalog(handoff): add HANDOFF_PROMPT.md (cold-boot prompt for a new Claude orchestrator)
- SHA: `3be3a7b26f99814b84498fa46a7a7193e93ef900`
- Message: flavor-catalog(handoff): add HANDOFF_PROMPT.md (cold-boot prompt for a new Claude orchestrator)
- Physical/numerical summary: Updated catalog orchestration/methodology handoff documentation; no numerical data or executable physics changed.
- Files touched (`git show --stat`):
```text
 flavor_catalog/HANDOFF_PROMPT.md | 160 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 160 insertions(+)
```

### 295. `24284d7` — flavor-catalog(policy): add PRIORITY_TIERS.md (PRIMARY = Waves 1-7; SECONDARY = Wave-8+, under processes/secondary/)
- SHA: `24284d73512c7d185175fa4f21e35627d6ee313f`
- Message: flavor-catalog(policy): add PRIORITY_TIERS.md (PRIMARY = Waves 1-7; SECONDARY = Wave-8+, under processes/secondary/)
- Physical/numerical summary: Updated catalog orchestration/methodology handoff documentation; no numerical data or executable physics changed.
- Files touched (`git show --stat`):
```text
 flavor_catalog/PRIORITY_TIERS.md | 109 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 109 insertions(+)
```

### 296. `82a2c01` — flavor-catalog(wave8): add wave_008_runbook.md (orchestration state file)
- SHA: `82a2c0167372dc053dd1dfba896855929c1cb8dc`
- Message: flavor-catalog(wave8): add wave_008_runbook.md (orchestration state file)
- Physical/numerical summary: Updated orchestration runbook state for the named wave/stage; workflow bookkeeping only.
- Files touched (`git show --stat`):
```text
 flavor_catalog/worklogs/orchestration/wave_008_runbook.md | 202 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 202 insertions(+)
```

### 297. `bab5bd0` — flavor-catalog(wave8): PKA-K020 initial draft (SECONDARY)
- SHA: `bab5bd04cad2437220d2610ed71a3ff51e70f2da`
- Message: flavor-catalog(wave8): PKA-K020 initial draft (SECONDARY)
- Physical/numerical summary: Added an initial catalog process draft for SECONDARY; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/secondary/kaon/K020.tex                                     | 113 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/secondary/kaon/K020.yaml                                    | 156 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/K020/angelescu_faroughy_sumensari_2020_arxiv2002_05684.txt | 260 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/K020/beneke_moch_rohrwild_2015_arxiv1508_01705.txt         | 260 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/K020/cfw_2008_arxiv0804_1954.txt                           | 260 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/K020/na62_2021_arxiv2105_06759.txt                         | 260 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/K020/pdg2025_kplus_lfv_semileptonic_api.txt                | 233 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/K020/sha256sums.txt                                        |   6 +++
 flavor_catalog/references/K020/sher_2005_arxiv_hep_ex_0502020.txt                    | 260 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/K020/source_manifest.yaml                                  |  68 +++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/K020.md                                                  |  93 ++++++++++++++++++++++++++++++++++++++++++++++
 11 files changed, 1969 insertions(+)
```

### 298. `2630168` — flavor-catalog(wave8): PKA-B008 initial draft (SECONDARY)
- SHA: `2630168f2824bb025c423a583259b8e85ccd585b`
- Message: flavor-catalog(wave8): PKA-B008 initial draft (SECONDARY)
- Physical/numerical summary: Added an initial catalog process draft for SECONDARY; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/secondary/beauty/B008.tex                         | 115 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/secondary/beauty/B008.yaml                        | 200 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B008/arxiv_0804_1954_cfw_rs_flavor.txt           |  20 ++++++++++++++
 flavor_catalog/references/B008/arxiv_1311_0903_bobeth_bqll_sm.txt          |  23 ++++++++++++++++
 flavor_catalog/references/B008/arxiv_1605_09637_babar_bktautau_related.txt |  22 +++++++++++++++
 flavor_catalog/references/B008/arxiv_1703_02508_lhcb_bq_tautau.txt         |  21 +++++++++++++++
 flavor_catalog/references/B008/arxiv_1712_01919_capdevila_bstautau_np.txt  |  20 ++++++++++++++
 flavor_catalog/references/B008/arxiv_2307_07013_bordone_bstautau_np.txt    |  23 ++++++++++++++++
 flavor_catalog/references/B008/arxiv_hep_ex_0511015_babar_bd_tautau.txt    |  16 +++++++++++
 flavor_catalog/references/B008/pdg_2026_bd_tautau.txt                      |  45 +++++++++++++++++++++++++++++++
 flavor_catalog/references/B008/pdg_2026_bs_tautau.txt                      |  27 +++++++++++++++++++
 flavor_catalog/references/B008/sha256sums.txt                              |   9 +++++++
 flavor_catalog/references/B008/source_manifest.yaml                        |  95 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/B008.md                                        |  69 +++++++++++++++++++++++++++++++++++++++++++++++
 14 files changed, 705 insertions(+)
```

### 299. `37beabb` — flavor-catalog(wave8): PKA-K020 initial draft (SECONDARY)
- SHA: `37beabb085ea01ea547f4afa605eb4fa85c54aff`
- Message: flavor-catalog(wave8): PKA-K020 initial draft (SECONDARY)
- Physical/numerical summary: Added an initial catalog process draft for SECONDARY; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/secondary/beauty/B013.tex                        | 113 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/secondary/beauty/B013.yaml                       | 253 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B013/belle_2014_bsphigamma_arxiv.txt            |  22 ++++++++++++
 flavor_catalog/references/B013/cfw_2008_arxiv.txt                         |  19 +++++++++++
 flavor_catalog/references/B013/hflav_rare_decays_dec2024_bs_phi_gamma.txt |  33 ++++++++++++++++++
 flavor_catalog/references/B013/lhcb_2019_bsphigamma_cp_arxiv.txt          |  26 +++++++++++++++
 flavor_catalog/references/B013/lhcb_2024_bskkgamma_arxiv.txt              |  24 +++++++++++++
 flavor_catalog/references/B013/lhcb_2024_bsphiee_arxiv.txt                |  26 +++++++++++++++
 flavor_catalog/references/B013/muheim_xie_zwicky_2008_arxiv.txt           |  23 +++++++++++++
 flavor_catalog/references/B013/pdg2025_bs_phi_gamma.txt                   |  50 ++++++++++++++++++++++++++++
 flavor_catalog/references/B013/sha256sums.txt                             |   8 +++++
 flavor_catalog/references/B013/source_manifest.yaml                       |  83 +++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/B013.md                                       | 100 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
 13 files changed, 780 insertions(+)
```

### 300. `ebd066c` — flavor-catalog(wave8): PKA-B013 initial draft (SECONDARY)
- SHA: `ebd066c6b7942a84399dc9df6ff156f915208d03`
- Message: flavor-catalog(wave8): PKA-B013 initial draft (SECONDARY)
- Physical/numerical summary: Added an initial catalog process draft for SECONDARY; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/worklogs/pka/B013.md | 7 +++++++
 1 file changed, 7 insertions(+)
```

### 301. `e6e1cc3` — flavor-catalog(wave8): PKA-B007 initial draft (SECONDARY)
- SHA: `e6e1cc37fa03826bec53445c00769beed1335d5e`
- Message: flavor-catalog(wave8): PKA-B007 initial draft (SECONDARY)
- Physical/numerical summary: Added an initial catalog process draft for SECONDARY; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/secondary/beauty/B007.tex              | 110 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/secondary/beauty/B007.yaml             | 220 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B007/babar_0712_1516_bdee.txt         |  17 ++++++++++++
 flavor_catalog/references/B007/bobeth_1311_0903_bqll_sm.txt     |  39 +++++++++++++++++++++++++++
 flavor_catalog/references/B007/cdf_0901_3803_bsee_bdee.txt      |  23 ++++++++++++++++
 flavor_catalog/references/B007/cfw_0804_1954_rs_flavor.txt      |  23 ++++++++++++++++
 flavor_catalog/references/B007/fleischer_1703_10160_bqll_np.txt |  24 +++++++++++++++++
 flavor_catalog/references/B007/lhcb_2003_03999_bsee_bdee.txt    |  41 ++++++++++++++++++++++++++++
 flavor_catalog/references/B007/pdg_2026_bsee_bdee.txt           |  70 +++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B007/sha256sums.txt                   |   7 +++++
 flavor_catalog/references/B007/source_manifest.yaml             |  73 +++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/B007.md                             |  72 +++++++++++++++++++++++++++++++++++++++++++++++++
 12 files changed, 719 insertions(+)
```

### 302. `5b68b43` — flavor-catalog(wave8): PKA-K021 initial draft (SECONDARY)
- SHA: `5b68b439b4d4991788b4dcd82c0e2f96cdc5ec54`
- Message: flavor-catalog(wave8): PKA-K021 initial draft (SECONDARY)
- Physical/numerical summary: Added an initial catalog process draft for SECONDARY; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/secondary/kaon/K021.tex                              | 114 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/secondary/kaon/K021.yaml                             | 180 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/K021/angelescu_2020_arxiv2002_05684.txt             |  21 ++++++++++++++++
 flavor_catalog/references/K021/crivellin_2016_arxiv1601_00970.txt             |  22 +++++++++++++++++
 flavor_catalog/references/K021/csaki_falkowski_weiler_2008_arxiv0804_1954.txt |  23 +++++++++++++++++
 flavor_catalog/references/K021/delzanno_2024_arxiv2411_13497.txt              |  22 +++++++++++++++++
 flavor_catalog/references/K021/ktev_abouzaid_2008_arxiv0711_3472.txt          |  23 +++++++++++++++++
 flavor_catalog/references/K021/na62_aliberti_2021_arxiv2105_06759.txt         |  25 +++++++++++++++++++
 flavor_catalog/references/K021/pdg2025_kl_pi0emu.txt                          |  33 +++++++++++++++++++++++++
 flavor_catalog/references/K021/perez_randall_2008_arxiv0805_4652.txt          |  22 +++++++++++++++++
 flavor_catalog/references/K021/roy_valencia_2024_arxiv2410_05859.txt          |  20 +++++++++++++++
 flavor_catalog/references/K021/sha256sums.txt                                 |   9 +++++++
 flavor_catalog/references/K021/source_manifest.yaml                           |  93 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/K021.md                                           |  92 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 14 files changed, 699 insertions(+)
```

### 303. `6b64a18` — flavor-catalog(wave8): PKA-T014 initial draft (SECONDARY)
- SHA: `6b64a18ad7db0dfca3b832c885fd09dfbdcb2837`
- Message: flavor-catalog(wave8): PKA-T014 initial draft (SECONDARY)
- Physical/numerical summary: Added an initial catalog process draft for SECONDARY; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/secondary/top_higgs_ew/T014.tex                      | 111 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/secondary/top_higgs_ew/T014.yaml                     | 218 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/T014/abu_ajamieh_2026_fv_z_quark_extract.txt        |  55 ++++++++++++++++++++++++++++++++++
 flavor_catalog/references/T014/blanke_buras_2009_rare_kb_warped_arxiv_abs.txt |  25 ++++++++++++++++
 flavor_catalog/references/T014/cfw_2008_warped_flavor_arxiv_abs.txt           |  23 +++++++++++++++
 flavor_catalog/references/T014/ecfa2025_fv_z_decays_extract.txt               |  53 +++++++++++++++++++++++++++++++++
 flavor_catalog/references/T014/eilam_2002_rare_z_decays_extract.txt           |  38 ++++++++++++++++++++++++
 flavor_catalog/references/T014/kamenik_2024_fv_higgs_z_extract.txt            |  47 +++++++++++++++++++++++++++++
 flavor_catalog/references/T014/pdg2025_z_hadronic_listing.txt                 |  30 +++++++++++++++++++
 flavor_catalog/references/T014/sha256sums.txt                                 |   7 +++++
 flavor_catalog/references/T014/source_manifest.yaml                           |  66 +++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/T014.md                                           |  63 +++++++++++++++++++++++++++++++++++++++
 12 files changed, 736 insertions(+)
```

### 304. `e348a35` — flavor-catalog(wave8): PKA-K019 initial draft (SECONDARY)
- SHA: `e348a358c3a9f994327f9dfa5afb0ce003d35542`
- Message: flavor-catalog(wave8): PKA-K019 initial draft (SECONDARY)
- Physical/numerical summary: Added an initial catalog process draft for SECONDARY; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/secondary/kaon/K019.tex                              | 108 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/secondary/kaon/K019.yaml                             | 169 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/K019/ambrose_1998_arxiv_hepex9811038.txt            |  21 +++++++++++++++++
 flavor_catalog/references/K019/beneke_moch_rohrwild_2015_arxiv1508_01705.txt  |  25 ++++++++++++++++++++
 flavor_catalog/references/K019/blanke_crivellin_2018_arxiv1801_07256.txt      |  25 ++++++++++++++++++++
 flavor_catalog/references/K019/buras_duling_2010_arxiv1006_5356.txt           |  23 +++++++++++++++++++
 flavor_catalog/references/K019/csaki_falkowski_weiler_2008_arxiv0804_1954.txt |  23 +++++++++++++++++++
 flavor_catalog/references/K019/dambrosio_iyer_2018_arxiv1712_08122.txt        |  25 ++++++++++++++++++++
 flavor_catalog/references/K019/lhcb_strange_2018_arxiv1808_03477.txt          |  25 ++++++++++++++++++++
 flavor_catalog/references/K019/pdg2025_conservation_laws_clfv.txt             |  36 +++++++++++++++++++++++++++++
 flavor_catalog/references/K019/pdg2025_kl_emu.txt                             |  36 +++++++++++++++++++++++++++++
 flavor_catalog/references/K019/perez_randall_2008_arxiv0805_4652.txt          |  25 ++++++++++++++++++++
 flavor_catalog/references/K019/sha256sums.txt                                 |  10 ++++++++
 flavor_catalog/references/K019/source_manifest.yaml                           | 103 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/K019.md                                           | 115 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 15 files changed, 769 insertions(+)
```

### 305. `17599d5` — flavor-catalog(wave8): PKA-B014 initial draft (SECONDARY)
- SHA: `17599d52e331fb245c2565a04f67f385273bbd07`
- Message: flavor-catalog(wave8): PKA-B014 initial draft (SECONDARY)
- Physical/numerical summary: Added an initial catalog process draft for SECONDARY; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/secondary/beauty/B014.tex                         | 139 ++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/secondary/beauty/B014.yaml                        | 383 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/B014/arxiv_0804_1954_cfw_rs_flavor.txt           |  15 ++++++
 flavor_catalog/references/B014/arxiv_0804_4770_belle_bd_gamma.txt          |  20 ++++++++
 flavor_catalog/references/B014/arxiv_0808_1379_babar_rho_omega_gamma.txt   |  19 +++++++
 flavor_catalog/references/B014/arxiv_1005_1224_hurth_nakao_penguins.txt    |  15 ++++++
 flavor_catalog/references/B014/arxiv_1005_4087_babar_xdgamma.txt           |  20 ++++++++
 flavor_catalog/references/B014/arxiv_1104_3342_c7_c7prime.txt              |  15 ++++++
 flavor_catalog/references/B014/arxiv_1503_05534_lcsr_form_factors.txt      |  15 ++++++
 flavor_catalog/references/B014/arxiv_2407_08984_belle_belleii_rhogamma.txt |  18 +++++++
 flavor_catalog/references/B014/arxiv_2507_14401_lhcb_b0rho0gamma.txt       |  21 ++++++++
 flavor_catalog/references/B014/hflav_rare_decays_dec2024_bd_gamma.txt      |  36 +++++++++++++
 flavor_catalog/references/B014/pdg2025_bplus_rhop_gamma.txt                |  21 ++++++++
 flavor_catalog/references/B014/pdg2025_bzero_rho_omega_gamma.txt           |  32 ++++++++++++
 flavor_catalog/references/B014/sha256sums.txt                              |  12 +++++
 flavor_catalog/references/B014/source_manifest.yaml                        | 124 +++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/B014.md                                        | 115 +++++++++++++++++++++++++++++++++++++++++
 17 files changed, 1020 insertions(+)
```

### 306. `0ad4a07` — flavor-catalog(wave8): runbook stage-1 complete (8/8 PKAs landed)
- SHA: `0ad4a074b2ab16fd403ecd11bd993fb269b43d8f`
- Message: flavor-catalog(wave8): runbook stage-1 complete (8/8 PKAs landed)
- Physical/numerical summary: Added an initial catalog process draft for s landed); physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/worklogs/orchestration/wave_008_runbook.md | 26 +++++++++++++++++---------
 1 file changed, 17 insertions(+), 9 deletions(-)
```

### 307. `4a6df17` — flavor-catalog(wave8): runbook - scaffold2 + 4 WAs dispatched
- SHA: `4a6df1788587aeb8aceb2409526a35034ee1cf24`
- Message: flavor-catalog(wave8): runbook - scaffold2 + 4 WAs dispatched
- Physical/numerical summary: Updated orchestration runbook state for the named wave/stage; workflow bookkeeping only.
- Files touched (`git show --stat`):
```text
 flavor_catalog/worklogs/orchestration/wave_008_runbook.md | 12 ++++++------
 1 file changed, 6 insertions(+), 6 deletions(-)
```

### 308. `103833f` — flavor-catalog(wave8): wire secondary/ section into catalog_master.tex
- SHA: `103833f57bff9f0223cd7605fc6ddb0f9ce65a00`
- Message: flavor-catalog(wave8): wire secondary/ section into catalog_master.tex
- Physical/numerical summary: Wired the named catalog family/section into catalog_master.tex so compiled catalog counts include those entries.
- Files touched (`git show --stat`):
```text
 flavor_catalog/catalog_master.tex                         | 17 +++++++++++++++++
 flavor_catalog/processes/secondary/beauty/index.tex       |  5 +++++
 flavor_catalog/processes/secondary/kaon/index.tex         |  4 ++++
 flavor_catalog/processes/secondary/top_higgs_ew/index.tex |  2 ++
 4 files changed, 28 insertions(+)
```

### 309. `d198787` — flavor-catalog(wave8): WA polish WA-w8-T014 (cycle 1)
- SHA: `d198787f41ce9552b3f47fd06a19c40d0aad0c28`
- Message: flavor-catalog(wave8): WA polish WA-w8-T014 (cycle 1)
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/secondary/top_higgs_ew/T014.tex  | 35 ++++++++++++++++-------------------
 flavor_catalog/processes/secondary/top_higgs_ew/T014.yaml | 10 +++++++++-
 flavor_catalog/worklogs/writer/wa_w8_T014.md              | 71 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 3 files changed, 96 insertions(+), 20 deletions(-)
```

### 310. `2b23464` — flavor-catalog(wave8): WA polish WA-w8-B-radiative (cycle 1)
- SHA: `2b2346428b2b39be4257d55f6e790f3a128b3e91`
- Message: flavor-catalog(wave8): WA polish WA-w8-B-radiative (cycle 1)
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/secondary/beauty/B013.tex  | 58 ++++++++++++++++++++++++++++++----------------------------
 flavor_catalog/processes/secondary/beauty/B013.yaml | 38 ++++++++++++++++++++++++++++++++++++--
 flavor_catalog/processes/secondary/beauty/B014.tex  | 61 +++++++++++++++++++++++++++++--------------------------------
 flavor_catalog/processes/secondary/beauty/B014.yaml | 26 ++++++++++++++++++++++----
 flavor_catalog/worklogs/writer/wa_w8_B_radiative.md | 82 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 5 files changed, 199 insertions(+), 66 deletions(-)
```

### 311. `f9c7ff4` — flavor-catalog(wave8): WA polish WA-w8-kaon-LFV (cycle 1)
- SHA: `f9c7ff4def0e7da52312fc3cded7e22aab9d584b`
- Message: flavor-catalog(wave8): WA polish WA-w8-kaon-LFV (cycle 1)
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/secondary/kaon/K019.tex  |  14 +++++++++-----
 flavor_catalog/processes/secondary/kaon/K019.yaml |  46 +++++++++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/processes/secondary/kaon/K020.tex  |  10 ++++------
 flavor_catalog/processes/secondary/kaon/K020.yaml |  20 ++++++++++++++++----
 flavor_catalog/processes/secondary/kaon/K021.tex  |  10 +++-------
 flavor_catalog/processes/secondary/kaon/K021.yaml | 113 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++--
 flavor_catalog/worklogs/writer/wa_w8_kaon_LFV.md  |  50 ++++++++++++++++++++++++++++++++++++++++++++++++++
 7 files changed, 238 insertions(+), 25 deletions(-)
```

### 312. `88f2cde` — flavor-catalog(wave8): WA polish WA-w8-B-rare-leptonic (cycle 1)
- SHA: `88f2cdebe207926239b37e21d7784e723ab9dc7f`
- Message: flavor-catalog(wave8): WA polish WA-w8-B-rare-leptonic (cycle 1)
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/secondary/beauty/B007.tex      |  17 ++++++++-------
 flavor_catalog/processes/secondary/beauty/B007.yaml     | 183 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-----------------------------------------
 flavor_catalog/processes/secondary/beauty/B008.tex      |  38 +++++++++++++++++----------------
 flavor_catalog/processes/secondary/beauty/B008.yaml     |  68 ++++++++++++++++++++++++++++++++--------------------------
 flavor_catalog/worklogs/writer/wa_w8_B_rare_leptonic.md |  76 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 5 files changed, 278 insertions(+), 104 deletions(-)
```

### 313. `cd2f81b` — flavor-catalog(wave8): runbook - stage 2+3 done, stage 4 CAs dispatched
- SHA: `cd2f81b2d93a03d36f3437bf7437f9a406078200`
- Message: flavor-catalog(wave8): runbook - stage 2+3 done, stage 4 CAs dispatched
- Physical/numerical summary: Updated orchestration runbook state for the named wave/stage; workflow bookkeeping only.
- Files touched (`git show --stat`):
```text
 flavor_catalog/worklogs/orchestration/wave_008_runbook.md | 26 ++++++++++++++++----------
 1 file changed, 16 insertions(+), 10 deletions(-)
```

### 314. `b3c35d7` — flavor-catalog(wave8): CA cycle-1 verdict CA-w8-B-radiative
- SHA: `b3c35d75e9848135f0a8b04a1375317f75d97c59`
- Message: flavor-catalog(wave8): CA cycle-1 verdict CA-w8-B-radiative
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/secondary/beauty/B013.yaml  |  9 +++++++++
 flavor_catalog/processes/secondary/beauty/B014.yaml  |  9 +++++++++
 flavor_catalog/worklogs/checker/ca_w8_B_radiative.md | 44 ++++++++++++++++++++++++++++++++++++++++++++
 3 files changed, 62 insertions(+)
```

### 315. `e9fcdf3` — flavor-catalog(wave8): CA cycle-1 verdict CA-w8-B-rare-leptonic
- SHA: `e9fcdf3e4015723c6ad756f33e7ac18415417556`
- Message: flavor-catalog(wave8): CA cycle-1 verdict CA-w8-B-rare-leptonic
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/secondary/beauty/B007.yaml      |  9 +++++++++
 flavor_catalog/processes/secondary/beauty/B008.yaml      |  9 +++++++++
 flavor_catalog/worklogs/checker/ca_w8_B_rare_leptonic.md | 46 ++++++++++++++++++++++++++++++++++++++++++++++
 3 files changed, 64 insertions(+)
```

### 316. `41069e7` — flavor-catalog(wave8): CA cycle-1 verdict CA-w8-kaon-LFV
- SHA: `41069e7b38a840dcbb142fe4ab16b9273b5ca2d6`
- Message: flavor-catalog(wave8): CA cycle-1 verdict CA-w8-kaon-LFV
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/secondary/kaon/K019.yaml | 10 ++++++++++
 flavor_catalog/processes/secondary/kaon/K020.yaml |  9 +++++++++
 flavor_catalog/processes/secondary/kaon/K021.yaml | 10 ++++++++++
 flavor_catalog/worklogs/checker/ca_w8_kaon_LFV.md | 42 ++++++++++++++++++++++++++++++++++++++++++
 4 files changed, 71 insertions(+)
```

### 317. `749d789` — flavor-catalog(wave8): CA cycle-1 verdict CA-w8-T014
- SHA: `749d789d559a1c822f6cae60e5cef9fdc7b7eac9`
- Message: flavor-catalog(wave8): CA cycle-1 verdict CA-w8-T014
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/secondary/top_higgs_ew/T014.yaml |  9 +++++++++
 flavor_catalog/worklogs/checker/ca_w8_T014.md             | 42 ++++++++++++++++++++++++++++++++++++++++++
 2 files changed, 51 insertions(+)
```

### 318. `5c14d2f` — flavor-catalog(wave8): WA-v2 K020 cycle-2 (promote Sher/E865 to pdg_or_equivalent)
- SHA: `5c14d2f4bb491668b76340d5e4149fda00143e54`
- Message: flavor-catalog(wave8): WA-v2 K020 cycle-2 (promote Sher/E865 to pdg_or_equivalent)
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/secondary/kaon/K020.yaml   | 30 +++++++++++++++++++++++++++++-
 flavor_catalog/worklogs/writer/wa_w8_kaon_LFV_v2.md | 42 ++++++++++++++++++++++++++++++++++++++++++
 2 files changed, 71 insertions(+), 1 deletion(-)
```

### 319. `a20d75a` — flavor-catalog(wave8): CA-v2 K020 cycle-2 verdict
- SHA: `a20d75a7ba37dabb97b610bb1687f972c26c0c4b`
- Message: flavor-catalog(wave8): CA-v2 K020 cycle-2 verdict
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/secondary/kaon/K020.yaml    | 13 ++++++++++++-
 flavor_catalog/worklogs/checker/ca_w8_kaon_LFV_v2.md | 33 +++++++++++++++++++++++++++++++++
 2 files changed, 45 insertions(+), 1 deletion(-)
```

### 320. `980b45b` — flavor-catalog(wave8): runbook - stage 4 DONE (8/8 CHECKER-DONE incl. K020 cycle-2), stage 5 dispatched
- SHA: `980b45baf31383fcf716836b903b1cfda1832799`
- Message: flavor-catalog(wave8): runbook - stage 4 DONE (8/8 CHECKER-DONE incl. K020 cycle-2), stage 5 dispatched
- Physical/numerical summary: Updated orchestration runbook state for the named wave/stage; workflow bookkeeping only.
- Files touched (`git show --stat`):
```text
 flavor_catalog/worklogs/orchestration/wave_008_runbook.md | 38 +++++++++++++++++++++++++-------------
 1 file changed, 25 insertions(+), 13 deletions(-)
```

### 321. `118de1c` — flavor-catalog(wave8): fact-check addendum top_higgs_ew (SECONDARY)
- SHA: `118de1c86c3c5eb23190b9e0ac4404ddccd041af`
- Message: flavor-catalog(wave8): fact-check addendum top_higgs_ew (SECONDARY)
- Physical/numerical summary: Recorded fact-check results; verification bookkeeping, not executable physics code.
- Files touched (`git show --stat`):
```text
 flavor_catalog/audits/factcheck_top_higgs_ew.md | 40 ++++++++++++++++++++++++++++++++++++++++
 1 file changed, 40 insertions(+)
```

### 322. `37514e9` — flavor-catalog(wave8): fact-check addendum beauty (SECONDARY)
- SHA: `37514e9ba7bb8b72ba8f04af4817defa86cc84d4`
- Message: flavor-catalog(wave8): fact-check addendum beauty (SECONDARY)
- Physical/numerical summary: Recorded fact-check results; verification bookkeeping, not executable physics code.
- Files touched (`git show --stat`):
```text
 flavor_catalog/audits/factcheck_beauty.md | 87 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 87 insertions(+)
```

### 323. `7e121c3` — flavor-catalog(wave8): fact-check addendum kaon (SECONDARY)
- SHA: `7e121c3a1764553cdb5a76da6546da3705dfb82f`
- Message: flavor-catalog(wave8): fact-check addendum kaon (SECONDARY)
- Physical/numerical summary: Recorded fact-check results; verification bookkeeping, not executable physics code.
- Files touched (`git show --stat`):
```text
 flavor_catalog/audits/factcheck_kaon.md | 61 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 61 insertions(+)
```

### 324. `df7af95` — flavor-catalog(wave8): Opus round-4 sign-off (8/8 APPROVE; 7 VERIFIED + 1 PARTIAL K020 metadata)
- SHA: `df7af952046447e51a3e0116105588570ec5e2ed`
- Message: flavor-catalog(wave8): Opus round-4 sign-off (8/8 APPROVE; 7 VERIFIED + 1 PARTIAL K020 metadata)
- Physical/numerical summary: Recorded reviewer/sign-off disposition (8/8); review ledger only.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/secondary/beauty/B007.yaml       | 21 ++++++++++++++++++++-
 flavor_catalog/processes/secondary/beauty/B008.yaml       | 21 ++++++++++++++++++++-
 flavor_catalog/processes/secondary/beauty/B013.yaml       | 21 ++++++++++++++++++++-
 flavor_catalog/processes/secondary/beauty/B014.yaml       | 21 ++++++++++++++++++++-
 flavor_catalog/processes/secondary/kaon/K019.yaml         | 21 ++++++++++++++++++++-
 flavor_catalog/processes/secondary/kaon/K020.yaml         | 21 ++++++++++++++++++++-
 flavor_catalog/processes/secondary/kaon/K021.yaml         | 21 ++++++++++++++++++++-
 flavor_catalog/processes/secondary/top_higgs_ew/T014.yaml | 21 ++++++++++++++++++++-
 flavor_catalog/signoff/round_004_index.md                 | 58 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 9 files changed, 218 insertions(+), 8 deletions(-)
```

### 325. `ef5a417` — flavor-catalog(wave8): runbook - stage 5 DONE (7V/1P), stage 6 DONE (8/8 APPROVE), stage 7 dispatched
- SHA: `ef5a41737d6e73c3c6ba168d0a0c7e884f838d21`
- Message: flavor-catalog(wave8): runbook - stage 5 DONE (7V/1P), stage 6 DONE (8/8 APPROVE), stage 7 dispatched
- Physical/numerical summary: Updated orchestration runbook state for the named wave/stage; workflow bookkeeping only.
- Files touched (`git show --stat`):
```text
 flavor_catalog/worklogs/orchestration/wave_008_runbook.md | 25 ++++++++++++++++++-------
 1 file changed, 18 insertions(+), 7 deletions(-)
```

### 326. `76c87ae` — flavor-catalog(wave8): master compile v0.3 - PRIMARY 80 + SECONDARY 8 OPUS-APPROVED
- SHA: `76c87ae4051f123768bbe46f5292cc7576ec12d1`
- Message: flavor-catalog(wave8): master compile v0.3 - PRIMARY 80 + SECONDARY 8 OPUS-APPROVED
- Physical/numerical summary: Updated the flavor-catalog master compile and headline counts (3, 80, 8 cited in subject); catalog aggregation only.
- Files touched (`git show --stat`):
```text
 flavor_catalog/catalog_master.pdf           | Bin 733862 -> 796017 bytes
 flavor_catalog/master_compile_v03_report.md |  56 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 2 files changed, 56 insertions(+)
```

### 327. `fff1898` — flavor-catalog(wave8): close-out - SESSION_NOTES + runbook updated, v0.3 tagged
- SHA: `fff18981c85c107362aef79d72dafc9fdea8d35e`
- Message: flavor-catalog(wave8): close-out - SESSION_NOTES + runbook updated, v0.3 tagged
- Physical/numerical summary: Updated orchestration runbook state for the named wave/stage; workflow bookkeeping only.
- Files touched (`git show --stat`):
```text
 flavor_catalog/SESSION_NOTES.md                           | 58 +++++++++++++++++++++++++++++++++++++++++++---------------
 flavor_catalog/worklogs/orchestration/wave_008_runbook.md | 23 +++++++++++++++++++++--
 2 files changed, 64 insertions(+), 17 deletions(-)
```

### 328. `6a2506c` — flavor-catalog(wave8): WA-v3 K020 cycle-3 metadata cleanup (align NA62 2021 author with PDG / K021 convention)
- SHA: `6a2506c5463dae2de27e2fb283e3ae178213c763`
- Message: flavor-catalog(wave8): WA-v3 K020 cycle-3 metadata cleanup (align NA62 2021 author with PDG / K021 convention)
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/secondary/kaon/K020.yaml   | 11 +++++++++--
 flavor_catalog/references/K020/source_manifest.yaml |  2 +-
 flavor_catalog/worklogs/writer/wa_w8_kaon_LFV_v3.md | 11 +++++++++++
 3 files changed, 21 insertions(+), 3 deletions(-)
```

### 329. `b5c2375` — flavor-catalog(wave8): re-fact-check K020 - VERIFIED (NA62 author metadata aligned after WA-v3)
- SHA: `b5c23756cfab30606225483575f21969a30c51af`
- Message: flavor-catalog(wave8): re-fact-check K020 - VERIFIED (NA62 author metadata aligned after WA-v3)
- Physical/numerical summary: Recorded fact-check results; verification bookkeeping, not executable physics code.
- Files touched (`git show --stat`):
```text
 flavor_catalog/audits/factcheck_kaon.md           |  8 +++++---
 flavor_catalog/processes/secondary/kaon/K020.yaml | 15 ++++++++++++---
 2 files changed, 17 insertions(+), 6 deletions(-)
```

### 330. `7ed9117` — flavor-catalog(wave9): scaffold - new collider_rs/ family + wave_009_runbook
- SHA: `7ed91179f19f6773f8c5f2b1dd63acdfdc16a473`
- Message: flavor-catalog(wave9): scaffold - new collider_rs/ family + wave_009_runbook
- Physical/numerical summary: Updated orchestration runbook state for the named wave/stage; workflow bookkeeping only.
- Files touched (`git show --stat`):
```text
 flavor_catalog/worklogs/orchestration/wave_009_runbook.md | 150 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 150 insertions(+)
```

### 331. `108cae4` — flavor-catalog(wave9): runbook - 14 PKAs dispatched
- SHA: `108cae4a01e65691383f2f322a87561cd698de0c`
- Message: flavor-catalog(wave9): runbook - 14 PKAs dispatched
- Physical/numerical summary: Added an initial catalog process draft for s dispatched; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/worklogs/orchestration/wave_009_runbook.md | 34 ++++++++++++++++++++--------------
 1 file changed, 20 insertions(+), 14 deletions(-)
```

### 332. `5911918` — flavor-catalog(wave9): PKA-CR012 initial draft (collider_rs PRIMARY)
- SHA: `591191861b37dfe3e610113a333220965797e0c7`
- Message: flavor-catalog(wave9): PKA-CR012 initial draft (collider_rs PRIMARY)
- Physical/numerical summary: Added an initial catalog process draft for collider_rs PRIMARY; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/collider_rs/CR012.tex                              | 121 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/collider_rs/CR012.yaml                             | 295 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/CR012/atlas_2015_arxiv1506_00962.txt              |  17 ++++++++
 flavor_catalog/references/CR012/atlas_2017_arxiv1708_04445.txt              |  17 ++++++++
 flavor_catalog/references/CR012/atlas_2019_arxiv1906_08589.txt              |  17 ++++++++
 flavor_catalog/references/CR012/cms_2017_combination_arxiv1705_09171.txt    |  17 ++++++++
 flavor_catalog/references/CR012/cms_2018_dijet_arxiv1708_05379.txt          |  17 ++++++++
 flavor_catalog/references/CR012/cms_2019_arxiv1906_05977.txt                |  18 +++++++++
 flavor_catalog/references/CR012/cms_2022_semileptonic_arxiv2109_06055.txt   |  17 ++++++++
 flavor_catalog/references/CR012/cms_2023_b2g_20_009_arxiv2210_00043.txt     |  34 ++++++++++++++++
 flavor_catalog/references/CR012/hepdata_cms_2023_b2g_20_009_record.txt      |  24 +++++++++++
 flavor_catalog/references/CR012/internal_quark_scan_methodology_extract.txt |  21 ++++++++++
 flavor_catalog/references/CR012/pappadopulo_2014_hvt_arxiv1402_4431.txt     |  16 ++++++++
 flavor_catalog/references/CR012/pdg2025_heavy_bosons_wprime_wz.txt          |  28 +++++++++++++
 flavor_catalog/references/CR012/sha256sums.txt                              |  12 ++++++
 flavor_catalog/references/CR012/source_manifest.yaml                        | 123 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/CR012.md                                        | 124 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 17 files changed, 918 insertions(+)
```

### 333. `f195c45` — flavor-catalog(wave9): PKA-CR013 initial draft (collider_rs PRIMARY)
- SHA: `f195c453d7d2fa8cd1b03eee1c8427876b1427ff`
- Message: flavor-catalog(wave9): PKA-CR013 initial draft (collider_rs PRIMARY)
- Physical/numerical summary: Added an initial catalog process draft for collider_rs PRIMARY; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/collider_rs/CR013.tex                                               |  125 +++++++
 flavor_catalog/processes/collider_rs/CR013.yaml                                              |  344 ++++++++++++++++++++
 flavor_catalog/references/CR013/atlas2013_arxiv1210_8389_abstract.xml                        |   28 ++
 flavor_catalog/references/CR013/atlas2015_arxiv1504_05511_abstract.xml                       |   28 ++
 flavor_catalog/references/CR013/atlas2015_arxiv1504_05511_fulltext.txt                       | 1675 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/CR013/atlas2021_arxiv2102_13405_abstract.xml                       |   28 ++
 flavor_catalog/references/CR013/atlas2021_arxiv2102_13405_fulltext.txt                       | 1709 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/CR013/cms2012_arxiv1112_0688_abstract.xml                          |   27 ++
 flavor_catalog/references/CR013/cms2018_arxiv1809_00327_abstract.xml                         |   28 ++
 flavor_catalog/references/CR013/cms2018_arxiv1809_00327_fulltext.txt                         | 2103 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/CR013/cms2024_arxiv2405_09320_abstract.xml                         |   28 ++
 flavor_catalog/references/CR013/cms2024_arxiv2405_09320_fulltext.txt                         | 2107 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/CR013/cms2024_exo_22_024_public_page.html                          | 1298 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/CR013/davoudiasl_hewett_rizzo_1999_arxiv_hepph9909255_abstract.xml |   36 ++
 flavor_catalog/references/CR013/hepdata2024_cms_diphoton_record.json                         |    1 +
 flavor_catalog/references/CR013/pdg2025_extra_dimensions_listing.txt                         | 1360 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/CR013/pdg2025_extra_dimensions_review.txt                          |  905 +++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/CR013/randall_sundrum_1999_arxiv_hepph9905221_abstract.xml         |   32 ++
 flavor_catalog/references/CR013/sha256sums.txt                                               |   16 +
 flavor_catalog/references/CR013/source_manifest.yaml                                         |  163 ++++++++++
 flavor_catalog/worklogs/pka/CR013.md                                                         |   83 +++++
 21 files changed, 12124 insertions(+)
```

### 334. `144cc61` — flavor-catalog(wave9): PKA-CR004 initial draft (collider_rs PRIMARY)
- SHA: `144cc6143531f4c5fec41398e893ea4d0ac26d89`
- Message: flavor-catalog(wave9): PKA-CR004 initial draft (collider_rs PRIMARY)
- Physical/numerical summary: Added an initial catalog process draft for collider_rs PRIMARY; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/collider_rs/CR004.tex                   |  129 +++
 flavor_catalog/processes/collider_rs/CR004.yaml                  |  219 +++++
 flavor_catalog/references/CR004/arxiv_0801_1679.txt              |  292 +++++++
 flavor_catalog/references/CR004/arxiv_1808_02343.txt             |  286 +++++++
 flavor_catalog/references/CR004/arxiv_2008_09835.txt             |  291 +++++++
 flavor_catalog/references/CR004/arxiv_2209_07327.txt             |  290 +++++++
 flavor_catalog/references/CR004/arxiv_2210_15413.txt             |  293 +++++++
 flavor_catalog/references/CR004/arxiv_2212_05263.txt             |  288 +++++++
 flavor_catalog/references/CR004/arxiv_2402_13808.txt             |  287 +++++++
 flavor_catalog/references/CR004/arxiv_2405_17605.txt             |  292 +++++++
 flavor_catalog/references/CR004/cms_b2g_20_014_public.txt        | 2545 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/CR004/cms_exo_23_006_review_public.txt | 6774 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/CR004/hepdata_2024_bvlq_record.txt     |  114 +++
 flavor_catalog/references/CR004/pdg2025_bprime_listing.txt       |  570 +++++++++++++
 flavor_catalog/references/CR004/sha256sums.txt                   |   12 +
 flavor_catalog/references/CR004/source_manifest.yaml             |  123 +++
 flavor_catalog/worklogs/pka/CR004.md                             |  143 ++++
 17 files changed, 12948 insertions(+)
```

### 335. `4054e52` — flavor-catalog(wave9): PKA-CR011 initial draft (collider_rs PRIMARY)
- SHA: `4054e5284b59ecd99ce459c2c97fe018f6e7ec1b`
- Message: flavor-catalog(wave9): PKA-CR011 initial draft (collider_rs PRIMARY)
- Physical/numerical summary: Added an initial catalog process draft for collider_rs PRIMARY; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/collider_rs/CR011.tex                   | 133 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/collider_rs/CR011.yaml                  | 238 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/CR011/atlas_2017_arxiv1611_02428.txt   |  25 ++++++++++++++++
 flavor_catalog/references/CR011/atlas_2020_arxiv2004_10612.txt   |  28 ++++++++++++++++++
 flavor_catalog/references/CR011/atlas_2024_arxiv2312_00420.txt   |  29 ++++++++++++++++++
 flavor_catalog/references/CR011/atlas_2025_arxiv2503_11317.txt   |  29 ++++++++++++++++++
 flavor_catalog/references/CR011/atlas_2026_arxiv2603_18630.txt   |  31 ++++++++++++++++++++
 flavor_catalog/references/CR011/cms_2020_arxiv2009_09429.txt     |  30 +++++++++++++++++++
 flavor_catalog/references/CR011/contino_2011_arxiv1109_1570.txt  |  23 +++++++++++++++
 flavor_catalog/references/CR011/pdg2025_wz_quartic_couplings.txt |  25 ++++++++++++++++
 flavor_catalog/references/CR011/sha256sums.txt                   |   8 +++++
 flavor_catalog/references/CR011/source_manifest.yaml             |  83 +++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/CR011.md                             |  82 +++++++++++++++++++++++++++++++++++++++++++++++++++
 13 files changed, 764 insertions(+)
```

### 336. `7614f6d` — flavor-catalog(wave9): PKA-CR007 initial draft (collider_rs PRIMARY)
- SHA: `7614f6d546b327eadbeafbca415e4d3b4de309d8`
- Message: flavor-catalog(wave9): PKA-CR007 initial draft (collider_rs PRIMARY)
- Physical/numerical summary: Added an initial catalog process draft for collider_rs PRIMARY; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/collider_rs/CR007.tex                                           | 140 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/collider_rs/CR007.yaml                                          | 254 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/CR007/agashe_davoudiasl_perez_soni_2007_arxiv_hepph0701186.txt |  37 ++++++++++++++++++
 flavor_catalog/references/CR007/agashe_et_al_2007_warped_gauge_arxiv0709_0007.txt        |  49 ++++++++++++++++++++++++
 flavor_catalog/references/CR007/atlas_2018_combination_arxiv1808_02380.txt               |  28 ++++++++++++++
 flavor_catalog/references/CR007/atlas_2019_hadronic_diboson_arxiv1906_08589.txt          |  28 ++++++++++++++
 flavor_catalog/references/CR007/atlas_2020_semileptonic_diboson_arxiv2004_14636.txt      |  28 ++++++++++++++
 flavor_catalog/references/CR007/cms_2021_dilepton_arxiv2103_02708.txt                    |  28 ++++++++++++++
 flavor_catalog/references/CR007/cms_2023_alljets_boson_pairs_arxiv2210_00043.txt         |  28 ++++++++++++++
 flavor_catalog/references/CR007/cms_2024_diphoton_arxiv2405_09320.txt                    |  28 ++++++++++++++
 flavor_catalog/references/CR007/pdg2025_extra_dimensions_rsg_extract.txt                 | 172 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/CR007/randall_sundrum_1999_arxiv_hepph9905221.txt              |  32 ++++++++++++++++
 flavor_catalog/references/CR007/sha256sums.txt                                           |  10 +++++
 flavor_catalog/references/CR007/source_manifest.yaml                                     | 103 ++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/CR007.md                                                     | 119 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 15 files changed, 1084 insertions(+)
```

### 337. `6176694` — flavor-catalog(wave9): PKA-CR010 initial draft (collider_rs PRIMARY)
- SHA: `617669476a71b20c68f107437f08673965a1ec6c`
- Message: flavor-catalog(wave9): PKA-CR010 initial draft (collider_rs PRIMARY)
- Physical/numerical summary: Added an initial catalog process draft for collider_rs PRIMARY; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/collider_rs/CR010.tex                            | 136 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/collider_rs/CR010.yaml                           | 289 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/CR010/aguilar_saavedra_2009_arxiv0907_3155.txt  |  24 ++++++++++++
 flavor_catalog/references/CR010/atlas_2018_arxiv1808_02343.txt            |  24 ++++++++++++
 flavor_catalog/references/CR010/atlas_2024_arxiv2401_17165.txt            |  24 ++++++++++++
 flavor_catalog/references/CR010/cms_2017_arxiv1706_03408.txt              |  27 +++++++++++++
 flavor_catalog/references/CR010/cms_2023_arxiv2209_07327.txt              |  28 ++++++++++++++
 flavor_catalog/references/CR010/cms_review_2025_arxiv2405_17605.txt       |  30 +++++++++++++++
 flavor_catalog/references/CR010/garberson_golling_2013_arxiv1301_4454.txt |  24 ++++++++++++
 flavor_catalog/references/CR010/pdg2025_tprime_bprime_vlq_limits.txt      |  66 ++++++++++++++++++++++++++++++++
 flavor_catalog/references/CR010/sha256sums.txt                            |   8 ++++
 flavor_catalog/references/CR010/source_manifest.yaml                      |  93 +++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/CR010.md                                      | 113 ++++++++++++++++++++++++++++++++++++++++++++++++++++++
 13 files changed, 886 insertions(+)
```

### 338. `0a8e20a` — flavor-catalog(wave9): PKA-CR009 initial draft (collider_rs PRIMARY)
- SHA: `0a8e20a40675ce6c2de911c9d52dbcf17a299337`
- Message: flavor-catalog(wave9): PKA-CR009 initial draft (collider_rs PRIMARY)
- Physical/numerical summary: Added an initial catalog process draft for collider_rs PRIMARY; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/collider_rs/CR009.tex                         | 141 +++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/collider_rs/CR009.yaml                        | 375 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/CR009/atlas_2020_arxiv2006_12946.txt         |  30 ++++++++++++
 flavor_catalog/references/CR009/cms_2021_arxiv2103_02708.txt           |  34 +++++++++++++
 flavor_catalog/references/CR009/cms_hepdata_2021_ci_limits_extract.txt |  29 +++++++++++
 flavor_catalog/references/CR009/lhc_search_history_arxiv_abstracts.txt |  52 ++++++++++++++++++++
 flavor_catalog/references/CR009/pdg2025_compositeness_extract.txt      |  31 ++++++++++++
 flavor_catalog/references/CR009/sha256sums.txt                         |   6 +++
 flavor_catalog/references/CR009/smeft_theory_arxiv_extracts.txt        |  46 ++++++++++++++++++
 flavor_catalog/references/CR009/source_manifest.yaml                   | 133 ++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/CR009.md                                   | 111 ++++++++++++++++++++++++++++++++++++++++++
 11 files changed, 988 insertions(+)
```

### 339. `28879ea` — flavor-catalog(wave9): PKA-CR001 initial draft (collider_rs PRIMARY)
- SHA: `28879ea2af3b22377617b725eb80812f21b51bd5`
- Message: flavor-catalog(wave9): PKA-CR001 initial draft (collider_rs PRIMARY)
- Physical/numerical summary: Added an initial catalog process draft for collider_rs PRIMARY; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/collider_rs/CR001.tex                                  | 129 ++++++++++++++++++++++++
 flavor_catalog/processes/collider_rs/CR001.yaml                                 | 321 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/CR001/atlas_2012_arxiv1207_2409.txt                   | 715 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/CR001/atlas_2013_arxiv1305_2756.txt                   | 706 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/CR001/atlas_2018_arxiv1804_10823.txt                  |  28 ++++++
 flavor_catalog/references/CR001/atlas_2020_arxiv2005_05138.txt                  |  28 ++++++
 flavor_catalog/references/CR001/atlas_2025_arxiv2512_17856.txt                  |  25 +++++
 flavor_catalog/references/CR001/cms_2017_arxiv1704_03366.txt                    | 712 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/CR001/cms_2019_arxiv1810_05905.txt                    |  28 ++++++
 flavor_catalog/references/CR001/cms_2026_arxiv2603_23454.txt                    |  25 +++++
 flavor_catalog/references/CR001/hepdata_2026_ins3134005_record_extract.json     |  82 +++++++++++++++
 flavor_catalog/references/CR001/lillie_randall_wang_2007_arxiv_hepph0701166.txt | 725 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/CR001/pdg_live_2026_s071kkg.txt                       |  15 +++
 flavor_catalog/references/CR001/quark_scan_methodology_note_mkk_extract.txt     | 101 +++++++++++++++++++
 flavor_catalog/references/CR001/sha256sums.txt                                  |  12 +++
 flavor_catalog/references/CR001/source_manifest.yaml                            | 123 +++++++++++++++++++++++
 flavor_catalog/worklogs/pka/CR001.md                                            | 116 +++++++++++++++++++++
 17 files changed, 3891 insertions(+)
```

### 340. `d09c6b3` — flavor-catalog(wave9): PKA-CR003 initial draft (collider_rs PRIMARY)
- SHA: `d09c6b3dc57fcd4ab1cd3b27542302566db5a925`
- Message: flavor-catalog(wave9): PKA-CR003 initial draft (collider_rs PRIMARY)
- Physical/numerical summary: Added an initial catalog process draft for collider_rs PRIMARY; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/collider_rs/CR003.tex                              | 120 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/collider_rs/CR003.yaml                             | 296 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/CR003/atlas_2017_arxiv1707_03347.txt              |  32 +++++++++++++++
 flavor_catalog/references/CR003/atlas_2018_hadronic_arxiv1808_01771.txt     |  34 ++++++++++++++++
 flavor_catalog/references/CR003/atlas_2018_multib_arxiv1803_09678.txt       |  34 ++++++++++++++++
 flavor_catalog/references/CR003/atlas_2022_arxiv2210_15413.txt              |  35 ++++++++++++++++
 flavor_catalog/references/CR003/atlas_2024_arxiv2401_17165.txt              |  33 ++++++++++++++++
 flavor_catalog/references/CR003/cms_2023_arxiv2209_07327.txt                |  33 ++++++++++++++++
 flavor_catalog/references/CR003/contino_servant_2008_arxiv0801_1679.txt     |  31 +++++++++++++++
 flavor_catalog/references/CR003/pdg_live_q009_tprime_pair_limits.txt        |  38 ++++++++++++++++++
 flavor_catalog/references/CR003/quark_scan_methodology_note_mkk_extract.txt |  24 +++++++++++
 flavor_catalog/references/CR003/sha256sums.txt                              |   9 +++++
 flavor_catalog/references/CR003/source_manifest.yaml                        |  93 +++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/CR003.md                                        | 114 ++++++++++++++++++++++++++++++++++++++++++++++++++++
 14 files changed, 926 insertions(+)
```

### 341. `835cc83` — flavor-catalog(wave9): PKA-CR002 initial draft (collider_rs PRIMARY)
- SHA: `835cc83a52a92a2a27e61e85f01f2369b365596f`
- Message: flavor-catalog(wave9): PKA-CR002 initial draft (collider_rs PRIMARY)
- Physical/numerical summary: Added an initial catalog process draft for collider_rs PRIMARY; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/collider_rs/CR002.tex                          | 124 +++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/collider_rs/CR002.yaml                         | 410 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/CR002/atlas_2018_arxiv1807_11883.txt          |  39 ++++++++++++++
 flavor_catalog/references/CR002/atlas_2023_arxiv2212_05263.txt          |  39 ++++++++++++++
 flavor_catalog/references/CR002/cms_2014_arxiv1312_2391.txt             |  29 ++++++++++
 flavor_catalog/references/CR002/cms_2017_arxiv1705_10967.txt            |  34 ++++++++++++
 flavor_catalog/references/CR002/cms_2019_arxiv1810_03188.txt            |  36 +++++++++++++
 flavor_catalog/references/CR002/contino_servant_2008_arxiv0801_1679.txt |  37 +++++++++++++
 flavor_catalog/references/CR002/mrazek_wulzer_2009_arxiv0909_3977.txt   |  35 ++++++++++++
 flavor_catalog/references/CR002/pdg_encoder_q009_tprime_5over3_2026.txt |  56 +++++++++++++++++++
 flavor_catalog/references/CR002/sha256sums.txt                          |   8 +++
 flavor_catalog/references/CR002/source_manifest.yaml                    |  83 +++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/CR002.md                                    | 163 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 13 files changed, 1093 insertions(+)
```

### 342. `92b78c7` — flavor-catalog(wave9): PKA-CR008 initial draft (collider_rs PRIMARY)
- SHA: `92b78c7ad1a02231dd2978de56ff262782c122c6`
- Message: flavor-catalog(wave9): PKA-CR008 initial draft (collider_rs PRIMARY)
- Physical/numerical summary: Added an initial catalog process draft for collider_rs PRIMARY; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/collider_rs/CR008.tex                                 | 112 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/collider_rs/CR008.yaml                                | 226 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/CR008/aguilar_saavedra_2009_arxiv0907_3155.txt       |  26 +++++++++++++
 flavor_catalog/references/CR008/aguilar_saavedra_et_al_2013_arxiv1306_0572.txt |  22 +++++++++++
 flavor_catalog/references/CR008/atlas_2015_arxiv1505_04306.txt                 |  32 ++++++++++++++++
 flavor_catalog/references/CR008/atlas_2023_arxiv2212_05263.txt                 |  29 +++++++++++++++
 flavor_catalog/references/CR008/atlas_2024_arxiv2401_17165.txt                 |  26 +++++++++++++
 flavor_catalog/references/CR008/cms_2017_arxiv1706_03408.txt                   |  30 +++++++++++++++
 flavor_catalog/references/CR008/cms_2018_arxiv1710_01539.txt                   |  27 ++++++++++++++
 flavor_catalog/references/CR008/cms_2018_arxiv1805_04758.txt                   |  25 +++++++++++++
 flavor_catalog/references/CR008/cms_2023_arxiv2209_07327.txt                   |  29 +++++++++++++++
 flavor_catalog/references/CR008/lari_et_al_2008_arxiv0801_1800.txt             |  21 +++++++++++
 flavor_catalog/references/CR008/pdg2025_tprime_vlq_limits.txt                  | 271 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/CR008/sha256sums.txt                                 |  11 ++++++
 flavor_catalog/references/CR008/source_manifest.yaml                           | 113 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/CR008.md                                           | 112 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
 16 files changed, 1112 insertions(+)
```

### 343. `c22fb07` — flavor-catalog(wave9): PKA-CR014 initial draft (collider_rs PRIMARY)
- SHA: `c22fb073f84646074b4b4f43285d6692ce74af47`
- Message: flavor-catalog(wave9): PKA-CR014 initial draft (collider_rs PRIMARY)
- Physical/numerical summary: Added an initial catalog process draft for collider_rs PRIMARY; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/collider_rs/CR014.tex                                               |  154 ++++++++++++
 flavor_catalog/processes/collider_rs/CR014.yaml                                              |  356 +++++++++++++++++++++++++++
 flavor_catalog/references/CR014/atlas_2019_arxiv1811_02305.txt                               |   28 +++
 flavor_catalog/references/CR014/atlas_2021_arxiv2106_11683.txt                               |   28 +++
 flavor_catalog/references/CR014/atlas_2023_arxiv2303_15061.txt                               |   28 +++
 flavor_catalog/references/CR014/atlas_top_philic_2024_arxiv2304_01678.txt                    |   28 +++
 flavor_catalog/references/CR014/banelli_salvioni_serra_theil_weiler_2021_arxiv2010_05915.txt |   41 ++++
 flavor_catalog/references/CR014/cao_chen_liu_2016_arxiv1602_01934.txt                        |   34 +++
 flavor_catalog/references/CR014/cms_2014_arxiv1409_7339.txt                                  |   28 +++
 flavor_catalog/references/CR014/cms_2017_arxiv1702_06164.txt                                 |   28 +++
 flavor_catalog/references/CR014/cms_2018_arxiv1710_10614.txt                                 |   28 +++
 flavor_catalog/references/CR014/cms_2020_arxiv1908_06463.txt                                 |   28 +++
 flavor_catalog/references/CR014/cms_2023_arxiv2305_13439.txt                                 |   28 +++
 flavor_catalog/references/CR014/cms_2023_evidence_arxiv2303_03864.txt                        |   28 +++
 flavor_catalog/references/CR014/cms_2026_arxiv2604_14058.txt                                 |   26 ++
 flavor_catalog/references/CR014/cms_b2g_25_005_public_result.txt                             | 1563 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/CR014/darme_fuks_li_maltoni_toucheque_2025_arxiv2507_05334.txt     |   41 ++++
 flavor_catalog/references/CR014/darme_fuks_maltoni_2021_arxiv2104_09512.txt                  |   35 +++
 flavor_catalog/references/CR014/de_blas_eberhardt_krause_2018_arxiv1803_00939.txt            |   35 +++
 flavor_catalog/references/CR014/pdg2025_top_four_top_extract.txt                             |  266 ++++++++++++++++++++
 flavor_catalog/references/CR014/quark_scan_methodology_note_extract.txt                      |  202 ++++++++++++++++
 flavor_catalog/references/CR014/sha256sums.txt                                               |   19 ++
 flavor_catalog/references/CR014/source_manifest.yaml                                         |  193 +++++++++++++++
 flavor_catalog/worklogs/pka/CR014.md                                                         |  142 +++++++++++
 24 files changed, 3387 insertions(+)
```

### 344. `42428d7` — flavor-catalog(wave9): PKA-CR006 initial draft (collider_rs PRIMARY)
- SHA: `42428d7c5f0f5e0fce32599d949c496339727c50`
- Message: flavor-catalog(wave9): PKA-CR006 initial draft (collider_rs PRIMARY)
- Physical/numerical summary: Added an initial catalog process draft for collider_rs PRIMARY; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/collider_rs/CR006.tex                                      | 155 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/collider_rs/CR006.yaml                                     | 303 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/CR006/agashe_2008_warped_charged_gauge_arxiv0810_1497.txt |  22 ++++++++++
 flavor_catalog/references/CR006/agashe_servant_2004_arxiv_hepph0403143.txt          |  21 +++++++++
 flavor_catalog/references/CR006/atlas_2018_wprime_tb_hadronic_arxiv1801_07893.txt   |  22 ++++++++++
 flavor_catalog/references/CR006/atlas_2018_wprime_tb_leptonjets_arxiv1807_10473.txt |  23 ++++++++++
 flavor_catalog/references/CR006/atlas_2019_wprime_lnu_arxiv1906_05609.txt           |  24 +++++++++++
 flavor_catalog/references/CR006/cms_2021_wprime_tb_hadronic_arxiv2104_04831.txt     |  21 +++++++++
 flavor_catalog/references/CR006/cms_2022_lnu_arxiv2202_06075.txt                    |  22 ++++++++++
 flavor_catalog/references/CR006/cms_2024_wprime_tb_leptonic_arxiv2310_19893.txt     |  23 ++++++++++
 flavor_catalog/references/CR006/pdg2025_wprime_searches.txt                         |  39 +++++++++++++++++
 flavor_catalog/references/CR006/sha256sums.txt                                      |   9 ++++
 flavor_catalog/references/CR006/source_manifest.yaml                                |  93 +++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/CR006.md                                                | 117 ++++++++++++++++++++++++++++++++++++++++++++++++++
 14 files changed, 894 insertions(+)
```

### 345. `1cd8b57` — flavor-catalog(wave9): PKA-CR005 initial draft (collider_rs PRIMARY)
- SHA: `1cd8b57fb942f4b63ecb775d5ff57f5084f617dc`
- Message: flavor-catalog(wave9): PKA-CR005 initial draft (collider_rs PRIMARY)
- Physical/numerical summary: Added an initial catalog process draft for collider_rs PRIMARY; physically records the process observable/limit, references, YAML, TeX, and worklog material.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/collider_rs/CR005.tex                                         | 134 ++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/collider_rs/CR005.yaml                                        | 338 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/CR005/agashe_delgado_may_sundrum_2003_arxiv_hepph0308036.txt |  27 ++++++++++
 flavor_catalog/references/CR005/agashe_servant_2004_arxiv_hepph0411254.txt             |  24 +++++++++
 flavor_catalog/references/CR005/atlas_2012_arxiv1209_2535.txt                          |  27 ++++++++++
 flavor_catalog/references/CR005/atlas_2017_arxiv1707_02424.txt                         |  26 ++++++++++
 flavor_catalog/references/CR005/atlas_2019_arxiv1903_06248.txt                         |  34 +++++++++++++
 flavor_catalog/references/CR005/bella_2010_arxiv1004_2432.txt                          |  25 ++++++++++
 flavor_catalog/references/CR005/cms_2018_arxiv1803_06292.txt                           |  26 ++++++++++
 flavor_catalog/references/CR005/cms_2021_arxiv2103_02708.txt                           |  33 +++++++++++++
 flavor_catalog/references/CR005/davoudiasl_hewett_rizzo_1999_arxiv_hepph9911262.txt    |  26 ++++++++++
 flavor_catalog/references/CR005/hepdata_atlas_2019_record.txt                          |  26 ++++++++++
 flavor_catalog/references/CR005/hepdata_cms_2021_record.txt                            |  28 +++++++++++
 flavor_catalog/references/CR005/pdg2024_zprime_searches.txt                            |  28 +++++++++++
 flavor_catalog/references/CR005/sha256sums.txt                                         |  12 +++++
 flavor_catalog/references/CR005/source_manifest.yaml                                   | 123 ++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/pka/CR005.md                                                   | 137 +++++++++++++++++++++++++++++++++++++++++++++++++++
 17 files changed, 1074 insertions(+)
```

### 346. `a13edbb` — flavor-catalog(wave9): clean CR005 worklog whitespace
- SHA: `a13edbba3802a5a02445cb651981c4f4c7159f46`
- Message: flavor-catalog(wave9): clean CR005 worklog whitespace
- Physical/numerical summary: Applied the change described by the subject; based on stat surface, this is documentation/catalog/orchestration work rather than a new physics-code calculation.
- Files touched (`git show --stat`):
```text
 flavor_catalog/worklogs/pka/CR005.md | 8 ++++----
 1 file changed, 4 insertions(+), 4 deletions(-)
```

### 347. `b96036b` — flavor-catalog(wave9): wire collider_rs/ section into catalog_master.tex
- SHA: `b96036b769a24c943ddbcfe765e82576adfc38ee`
- Message: flavor-catalog(wave9): wire collider_rs/ section into catalog_master.tex
- Physical/numerical summary: Wired the named catalog family/section into catalog_master.tex so compiled catalog counts include those entries.
- Files touched (`git show --stat`):
```text
 flavor_catalog/catalog_master.tex              | 12 ++++++++++++
 flavor_catalog/processes/collider_rs/index.tex | 15 +++++++++++++++
 2 files changed, 27 insertions(+)
```

### 348. `1eac13e` — flavor-catalog(wave9): WA polish WA-w9-vlq-4top (cycle 1)
- SHA: `1eac13e86e66e95903c6f00631df0a5f1f5f924f`
- Message: flavor-catalog(wave9): WA polish WA-w9-vlq-4top (cycle 1)
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/collider_rs/CR008.tex   | 15 ++++++++-------
 flavor_catalog/processes/collider_rs/CR008.yaml  | 11 ++++++++++-
 flavor_catalog/processes/collider_rs/CR010.tex   |  9 ++++++---
 flavor_catalog/processes/collider_rs/CR010.yaml  | 11 ++++++++++-
 flavor_catalog/processes/collider_rs/CR014.tex   |  6 +++---
 flavor_catalog/processes/collider_rs/CR014.yaml  | 17 ++++++++++++++++-
 flavor_catalog/worklogs/writer/wa_w9_vlq_4top.md | 93 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 7 files changed, 146 insertions(+), 16 deletions(-)
```

### 349. `c6d55b0` — flavor-catalog(wave9): WA polish WA-w9-kk-gauge (cycle 1)
- SHA: `c6d55b0f38f84ed13900572d4656f5786f21c498`
- Message: flavor-catalog(wave9): WA polish WA-w9-kk-gauge (cycle 1)
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/collider_rs/CR001.tex   |  8 ++++----
 flavor_catalog/processes/collider_rs/CR001.yaml  | 13 +++++++++++--
 flavor_catalog/processes/collider_rs/CR005.tex   | 10 +++++-----
 flavor_catalog/processes/collider_rs/CR005.yaml  | 13 +++++++++++--
 flavor_catalog/processes/collider_rs/CR006.tex   | 16 ++++++++--------
 flavor_catalog/processes/collider_rs/CR006.yaml  | 13 +++++++++++--
 flavor_catalog/processes/collider_rs/CR007.tex   | 10 +++++-----
 flavor_catalog/processes/collider_rs/CR007.yaml  | 15 ++++++++++++---
 flavor_catalog/worklogs/writer/wa_w9_kk_gauge.md | 71 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 9 files changed, 138 insertions(+), 31 deletions(-)
```

### 350. `2286a39` — flavor-catalog(wave9): WA polish WA-w9-ew-tail (cycle 1)
- SHA: `2286a393d73e9157067412f72b749b6599638535`
- Message: flavor-catalog(wave9): WA polish WA-w9-ew-tail (cycle 1)
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/collider_rs/CR009.tex  |  16 +++++++++-------
 flavor_catalog/processes/collider_rs/CR009.yaml |  28 ++++++++++++++--------------
 flavor_catalog/processes/collider_rs/CR011.tex  |  10 ++++++----
 flavor_catalog/processes/collider_rs/CR011.yaml |  19 ++++++++++++++-----
 flavor_catalog/processes/collider_rs/CR012.tex  |  14 ++++++++------
 flavor_catalog/processes/collider_rs/CR012.yaml |  21 +++++++++++++++------
 flavor_catalog/processes/collider_rs/CR013.tex  |  23 +++++++++++++----------
 flavor_catalog/processes/collider_rs/CR013.yaml |  42 +++++++++++++++++++-----------------------
 flavor_catalog/worklogs/writer/wa_w9_ew_tail.md | 106 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 9 files changed, 204 insertions(+), 75 deletions(-)
```

### 351. `cd8a3fe` — flavor-catalog(wave9): WA polish WA-w9-custodial (cycle 1)
- SHA: `cd8a3fe4b1bcf37c63fbcfcd8ab6bd022973871f`
- Message: flavor-catalog(wave9): WA polish WA-w9-custodial (cycle 1)
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/collider_rs/CR002.tex    |  66 +++++++++++++++++++++++++++++++++++-------------------------------
 flavor_catalog/processes/collider_rs/CR002.yaml   |  13 +++++++++++--
 flavor_catalog/processes/collider_rs/CR003.tex    |  66 +++++++++++++++++++++++++++++++++++++-----------------------------
 flavor_catalog/processes/collider_rs/CR003.yaml   |  13 +++++++++++--
 flavor_catalog/processes/collider_rs/CR004.tex    | 111 ++++++++++++++++++++++++++++++++++++++++++++++++++++++---------------------------------------------------------
 flavor_catalog/processes/collider_rs/CR004.yaml   |  20 ++++++++++++++++++--
 flavor_catalog/worklogs/writer/wa_w9_custodial.md | 142 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 7 files changed, 308 insertions(+), 123 deletions(-)
```

### 352. `e8daa00` — flavor-catalog(wave9): CA cycle-1 verdict CA-w9-vlq-4top
- SHA: `e8daa0064b110008b7426d56dcff5408af6dd724`
- Message: flavor-catalog(wave9): CA cycle-1 verdict CA-w9-vlq-4top
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/collider_rs/CR008.yaml   | 11 ++++++++++-
 flavor_catalog/processes/collider_rs/CR010.yaml   | 11 ++++++++++-
 flavor_catalog/processes/collider_rs/CR014.yaml   | 12 +++++++++++-
 flavor_catalog/worklogs/checker/ca_w9_vlq_4top.md | 53 +++++++++++++++++++++++++++++++++++++++++++++++++++++
 4 files changed, 84 insertions(+), 3 deletions(-)
```

### 353. `7a6fa3d` — flavor-catalog(wave9): CA cycle-1 verdict CA-w9-kk-gauge
- SHA: `7a6fa3d5dc4947d403cd6b12aa2e559eda08a869`
- Message: flavor-catalog(wave9): CA cycle-1 verdict CA-w9-kk-gauge
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/collider_rs/CR001.yaml   |  9 +++++++++
 flavor_catalog/processes/collider_rs/CR005.yaml   |  9 +++++++++
 flavor_catalog/processes/collider_rs/CR006.yaml   |  9 +++++++++
 flavor_catalog/processes/collider_rs/CR007.yaml   |  9 +++++++++
 flavor_catalog/worklogs/checker/ca_w9_kk_gauge.md | 54 ++++++++++++++++++++++++++++++++++++++++++++++++++++++
 5 files changed, 90 insertions(+)
```

### 354. `bc35a06` — flavor-catalog(wave9): CA cycle-1 verdict CA-w9-custodial
- SHA: `bc35a068fe21a9191a00405408ab33aee870c63b`
- Message: flavor-catalog(wave9): CA cycle-1 verdict CA-w9-custodial
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/collider_rs/CR002.yaml    |  9 +++++++++
 flavor_catalog/processes/collider_rs/CR003.yaml    |  9 +++++++++
 flavor_catalog/processes/collider_rs/CR004.yaml    |  9 +++++++++
 flavor_catalog/worklogs/checker/ca_w9_custodial.md | 53 +++++++++++++++++++++++++++++++++++++++++++++++++++++
 4 files changed, 80 insertions(+)
```

### 355. `e846f14` — flavor-catalog(wave9): CA cycle-1 verdict CA-w9-ew-tail
- SHA: `e846f142a1bdc80364fc64939ac3922ea37695d8`
- Message: flavor-catalog(wave9): CA cycle-1 verdict CA-w9-ew-tail
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/collider_rs/CR009.yaml  | 11 ++++++++++-
 flavor_catalog/processes/collider_rs/CR011.yaml  | 11 ++++++++++-
 flavor_catalog/processes/collider_rs/CR012.yaml  | 12 +++++++++++-
 flavor_catalog/processes/collider_rs/CR013.yaml  | 13 +++++++++++--
 flavor_catalog/worklogs/checker/ca_w9_ew_tail.md | 56 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 5 files changed, 98 insertions(+), 5 deletions(-)
```

### 356. `a8758ac` — flavor-catalog(wave9): WA-v2 vlq_4top cycle-2 (CR008 CHK-1; CR010 CHK-1+CHK-2)
- SHA: `a8758accdd06d36e6f477cba2a1f7682f20852a1`
- Message: flavor-catalog(wave9): WA-v2 vlq_4top cycle-2 (CR008 CHK-1; CR010 CHK-1+CHK-2)
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/collider_rs/CR008.tex      | 20 ++++++++++----------
 flavor_catalog/processes/collider_rs/CR008.yaml     | 10 +++++++++-
 flavor_catalog/processes/collider_rs/CR010.tex      | 24 ++++++++++++------------
 flavor_catalog/processes/collider_rs/CR010.yaml     | 10 +++++++++-
 flavor_catalog/worklogs/writer/wa_w9_vlq_4top_v2.md | 35 +++++++++++++++++++++++++++++++++++
 5 files changed, 75 insertions(+), 24 deletions(-)
```

### 357. `e1aec33` — flavor-catalog(wave9): WA-v2 ew_tail cycle-2 (CR009 CHK-1; CR011 CHK-2)
- SHA: `e1aec33a67ea49de48e29f0acced5e413ff79eb2`
- Message: flavor-catalog(wave9): WA-v2 ew_tail cycle-2 (CR009 CHK-1; CR011 CHK-2)
- Physical/numerical summary: Applied writer/polish revisions to catalog YAML/TeX for the named processes; adjusts prose/provenance/status fields over existing physical limits.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/collider_rs/CR009.tex     | 21 ++++++++++-----------
 flavor_catalog/processes/collider_rs/CR009.yaml    | 18 +++++++++++++-----
 flavor_catalog/processes/collider_rs/CR011.tex     | 10 +++++-----
 flavor_catalog/processes/collider_rs/CR011.yaml    | 10 +++++++++-
 flavor_catalog/worklogs/writer/wa_w9_ew_tail_v2.md | 34 ++++++++++++++++++++++++++++++++++
 5 files changed, 71 insertions(+), 22 deletions(-)
```

### 358. `950ca36` — flavor-catalog(wave9): CA-v2 vlq_4top cycle-2 verdict
- SHA: `950ca367feb0f91c92acd04f03aaed48d75d17e9`
- Message: flavor-catalog(wave9): CA-v2 vlq_4top cycle-2 verdict
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/collider_rs/CR008.yaml      | 12 +++++++++++-
 flavor_catalog/processes/collider_rs/CR010.yaml      | 12 +++++++++++-
 flavor_catalog/worklogs/checker/ca_w9_vlq_4top_v2.md | 41 +++++++++++++++++++++++++++++++++++++++++
 3 files changed, 63 insertions(+), 2 deletions(-)
```

### 359. `82daa9b` — flavor-catalog(wave9): CA-v2 ew_tail cycle-2 verdict
- SHA: `82daa9b5da93d986af4d8e7f91c11c2151d1f069`
- Message: flavor-catalog(wave9): CA-v2 ew_tail cycle-2 verdict
- Physical/numerical summary: Recorded checker verdict/status updates for the named catalog batch; physics content is validation metadata, not new executable numerics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/collider_rs/CR009.yaml     | 11 ++++++++++-
 flavor_catalog/processes/collider_rs/CR011.yaml     | 11 ++++++++++-
 flavor_catalog/worklogs/checker/ca_w9_ew_tail_v2.md | 41 +++++++++++++++++++++++++++++++++++++++++
 3 files changed, 61 insertions(+), 2 deletions(-)
```

### 360. `0c5dacc` — flavor-catalog(wave9): fact-check report collider_rs family (14 entries)
- SHA: `0c5dacc254ada530158ada38e80a87384c3f16ef`
- Message: flavor-catalog(wave9): fact-check report collider_rs family (14 entries)
- Physical/numerical summary: Recorded fact-check results; verification bookkeeping, not executable physics code.
- Files touched (`git show --stat`):
```text
 flavor_catalog/audits/factcheck_collider_rs.md | 220 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/audits/factcheck_status.md      |  19 +++++++++++++--
 2 files changed, 237 insertions(+), 2 deletions(-)
```

### 361. `cdd0238` — flavor-catalog(wave9): Opus round-5 sign-off (14/14 APPROVE; 14 VERIFIED)
- SHA: `cdd023879568ad106829db85284780cdb77d2fa6`
- Message: flavor-catalog(wave9): Opus round-5 sign-off (14/14 APPROVE; 14 VERIFIED)
- Physical/numerical summary: Recorded reviewer/sign-off disposition (14/14); review ledger only.
- Files touched (`git show --stat`):
```text
 flavor_catalog/processes/collider_rs/CR001.yaml | 21 ++++++++++++++++++++-
 flavor_catalog/processes/collider_rs/CR002.yaml | 21 ++++++++++++++++++++-
 flavor_catalog/processes/collider_rs/CR003.yaml | 21 ++++++++++++++++++++-
 flavor_catalog/processes/collider_rs/CR004.yaml | 21 ++++++++++++++++++++-
 flavor_catalog/processes/collider_rs/CR005.yaml | 21 ++++++++++++++++++++-
 flavor_catalog/processes/collider_rs/CR006.yaml | 21 ++++++++++++++++++++-
 flavor_catalog/processes/collider_rs/CR007.yaml | 21 ++++++++++++++++++++-
 flavor_catalog/processes/collider_rs/CR008.yaml | 21 ++++++++++++++++++++-
 flavor_catalog/processes/collider_rs/CR009.yaml | 21 ++++++++++++++++++++-
 flavor_catalog/processes/collider_rs/CR010.yaml | 21 ++++++++++++++++++++-
 flavor_catalog/processes/collider_rs/CR011.yaml | 21 ++++++++++++++++++++-
 flavor_catalog/processes/collider_rs/CR012.yaml | 21 ++++++++++++++++++++-
 flavor_catalog/processes/collider_rs/CR013.yaml | 24 ++++++++++++++++++++----
 flavor_catalog/processes/collider_rs/CR014.yaml | 21 ++++++++++++++++++++-
 flavor_catalog/signoff/round_005_index.md       | 64 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 15 files changed, 344 insertions(+), 17 deletions(-)
```

### 362. `2ad34b1` — flavor-catalog(wave9): master compile v0.4 - PRIMARY 94 (80 flavor + 14 collider_rs) + SECONDARY 8
- SHA: `2ad34b1f0988932fa2dc0c21f298207a971d7677`
- Message: flavor-catalog(wave9): master compile v0.4 - PRIMARY 94 (80 flavor + 14 collider_rs) + SECONDARY 8
- Physical/numerical summary: Updated the flavor-catalog master compile and headline counts (4, 94, 80, 14, 8 cited in subject); catalog aggregation only.
- Files touched (`git show --stat`):
```text
 flavor_catalog/catalog_master.pdf              | Bin 796017 -> 913614 bytes
 flavor_catalog/master_compile_v04_report.md    |  60 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/collider_rs/index.tex |   2 +-
 3 files changed, 61 insertions(+), 1 deletion(-)
```

### 363. `864cd6d` — flavor-catalog(wave9): close-out - SESSION_NOTES + runbook updated, v0.4 tagged
- SHA: `864cd6d6050169b0c8727880a4a358adab9a83b4`
- Message: flavor-catalog(wave9): close-out - SESSION_NOTES + runbook updated, v0.4 tagged
- Physical/numerical summary: Updated orchestration runbook state for the named wave/stage; workflow bookkeeping only.
- Files touched (`git show --stat`):
```text
 flavor_catalog/SESSION_NOTES.md                           | 48 ++++++++++++++++++++++++++++++++++++++++--------
 flavor_catalog/worklogs/orchestration/wave_009_runbook.md | 37 ++++++++++++++++++++++++++++++++++++-
 2 files changed, 76 insertions(+), 9 deletions(-)
```

### 364. `6f072d5` — flavor-catalog(meta): add CATALOG_METHODOLOGY.tex (one-page collaborator-facing pipeline pitch)
- SHA: `6f072d54f25dc3b08a0aff5df22f46b164fe09ea`
- Message: flavor-catalog(meta): add CATALOG_METHODOLOGY.tex (one-page collaborator-facing pipeline pitch)
- Physical/numerical summary: Updated catalog orchestration/methodology handoff documentation; no numerical data or executable physics changed.
- Files touched (`git show --stat`):
```text
 flavor_catalog/CATALOG_METHODOLOGY.pdf | Bin 0 -> 123868 bytes
 flavor_catalog/CATALOG_METHODOLOGY.tex | 105 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 2 files changed, 105 insertions(+)
```

### 365. `d1bbaf4` — flavor-catalog(meta): CATALOG_METHODOLOGY - replace 'sister artifact' with concrete 'companion quark-scan paper'
- SHA: `d1bbaf4c9d8992bb4e22bd4a6847b9ec80deade4`
- Message: flavor-catalog(meta): CATALOG_METHODOLOGY - replace 'sister artifact' with concrete 'companion quark-scan paper'
- Physical/numerical summary: Updated catalog orchestration/methodology handoff documentation; no numerical data or executable physics changed.
- Files touched (`git show --stat`):
```text
 flavor_catalog/CATALOG_METHODOLOGY.pdf | Bin 123868 -> 123868 bytes
 1 file changed, 0 insertions(+), 0 deletions(-)
```

### 366. `781bc20` — flavor-catalog(meta): CATALOG_METHODOLOGY - replace 'sister artifact' with 'companion quark-scan paper'
- SHA: `781bc201075e68167dbf567d77e141bc1092bfcd`
- Message: flavor-catalog(meta): CATALOG_METHODOLOGY - replace 'sister artifact' with 'companion quark-scan paper'
- Physical/numerical summary: Updated catalog orchestration/methodology handoff documentation; no numerical data or executable physics changed.
- Files touched (`git show --stat`):
```text
 flavor_catalog/CATALOG_METHODOLOGY.pdf | Bin 123868 -> 123897 bytes
 flavor_catalog/CATALOG_METHODOLOGY.tex |   6 +++---
 2 files changed, 3 insertions(+), 3 deletions(-)
```

### 367. `2bda5f1` — flavor-catalog(meta): refresh all orchestration docs for post-Wave-9 v0.4 state
- SHA: `2bda5f10eea1d5150d1f0c61118e4b78e43159ad`
- Message: flavor-catalog(meta): refresh all orchestration docs for post-Wave-9 v0.4 state
- Physical/numerical summary: Updated catalog meta-documentation/policy for orchestration and handoff; no physical numerics changed.
- Files touched (`git show --stat`):
```text
 flavor_catalog/AGENTIC_WORKFLOW.md |  52 ++++++++++++++++++++++++++++++++---
 flavor_catalog/HANDOFF_PROMPT.md   | 264 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++---------------------------------------------------------
 flavor_catalog/PRIORITY_TIERS.md   | 115 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-----------------
 flavor_catalog/README.md           | 158 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++------------------------
 4 files changed, 442 insertions(+), 147 deletions(-)
```

### 368. `e3b37c7` — flavor-catalog(meta): add WEBSITE_BUILD_PROMPT.md (cold-boot prompt for website build)
- SHA: `e3b37c7a9bfe2bc96910f50655dca1c3fb8c4a17`
- Message: flavor-catalog(meta): add WEBSITE_BUILD_PROMPT.md (cold-boot prompt for website build)
- Physical/numerical summary: Updated catalog orchestration/methodology handoff documentation; no numerical data or executable physics changed.
- Files touched (`git show --stat`):
```text
 flavor_catalog/WEBSITE_BUILD_PROMPT.md | 478 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1 file changed, 478 insertions(+)
```

### 369. `6d99b17` — flavor-catalog(meta): WEBSITE_BUILD_PROMPT refactored — orchestrator-only + stack freedom
- SHA: `6d99b17b8cac73f441559a9dac014e413e18c609`
- Message: flavor-catalog(meta): WEBSITE_BUILD_PROMPT refactored — orchestrator-only + stack freedom
- Physical/numerical summary: Updated catalog orchestration/methodology handoff documentation; no numerical data or executable physics changed.
- Files touched (`git show --stat`):
```text
 flavor_catalog/WEBSITE_BUILD_PROMPT.md | 552 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++---------------------------------------------------------------
 1 file changed, 350 insertions(+), 202 deletions(-)
```

### 370. `460ece9` — flavor-catalog-website(phase-1): scaffold Astro site + ingest + K001 template
- SHA: `460ece94b4d1749972e8e221dd579ac0d12056bb`
- Message: flavor-catalog-website(phase-1): scaffold Astro site + ingest + K001 template
- Physical/numerical summary: Scaffolded the Astro catalog website, ingest pipeline, and K001 template; converts catalog YAML/TeX into static web content.
- Files touched (`git show --stat`):
```text
 flavor_catalog/website/.gitignore                        |    5 +
 flavor_catalog/website/WEBSITE_RUNBOOK.md                |   51 ++
 flavor_catalog/website/astro.config.mjs                  |   20 +
 flavor_catalog/website/package-lock.json                 | 5599 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/website/package.json                      |   19 +
 flavor_catalog/website/scripts/ingest_catalog.py         |  573 ++++++++++++++++
 flavor_catalog/website/src/components/Badges.astro       |   64 ++
 flavor_catalog/website/src/content.config.ts             |   97 +++
 flavor_catalog/website/src/content/catalog_index.json    | 1039 +++++++++++++++++++++++++++++
 flavor_catalog/website/src/content/entries/B001.json     |  233 +++++++
 flavor_catalog/website/src/content/entries/B002.json     |  219 ++++++
 flavor_catalog/website/src/content/entries/B003.json     |  318 +++++++++
 flavor_catalog/website/src/content/entries/B004.json     |  177 +++++
 flavor_catalog/website/src/content/entries/B005.json     |  163 +++++
 flavor_catalog/website/src/content/entries/B006.json     |  187 ++++++
 flavor_catalog/website/src/content/entries/B007.json     |  333 ++++++++++
 flavor_catalog/website/src/content/entries/B008.json     |  250 +++++++
 flavor_catalog/website/src/content/entries/B009.json     |  215 ++++++
 flavor_catalog/website/src/content/entries/B011.json     |  194 ++++++
 flavor_catalog/website/src/content/entries/B012.json     |  258 ++++++++
 flavor_catalog/website/src/content/entries/B013.json     |  369 +++++++++++
 flavor_catalog/website/src/content/entries/B014.json     |  449 +++++++++++++
 flavor_catalog/website/src/content/entries/B015.json     |  145 ++++
 flavor_catalog/website/src/content/entries/B016.json     |  143 ++++
 flavor_catalog/website/src/content/entries/B017.json     |  144 ++++
 flavor_catalog/website/src/content/entries/B018.json     |  145 ++++
 flavor_catalog/website/src/content/entries/B019.json     |  121 ++++
 flavor_catalog/website/src/content/entries/B021.json     |  256 +++++++
 flavor_catalog/website/src/content/entries/B022.json     |  157 +++++
 flavor_catalog/website/src/content/entries/B023.json     |  185 ++++++
 flavor_catalog/website/src/content/entries/B025.json     |  167 +++++
 flavor_catalog/website/src/content/entries/B026.json     |  160 +++++
 flavor_catalog/website/src/content/entries/B032.json     |  162 +++++
 flavor_catalog/website/src/content/entries/B033.json     |  121 ++++
 flavor_catalog/website/src/content/entries/B034.json     |  136 ++++
 flavor_catalog/website/src/content/entries/C001.json     |  283 ++++++++
 flavor_catalog/website/src/content/entries/C002.json     |  234 +++++++
 flavor_catalog/website/src/content/entries/C003.json     |  276 ++++++++
 flavor_catalog/website/src/content/entries/C004.json     |  219 ++++++
 flavor_catalog/website/src/content/entries/C005.json     |  191 ++++++
 flavor_catalog/website/src/content/entries/C006.json     |  164 +++++
 flavor_catalog/website/src/content/entries/C007.json     |  231 +++++++
 flavor_catalog/website/src/content/entries/C008.json     |  164 +++++
 flavor_catalog/website/src/content/entries/CR001.json    |  321 +++++++++
 flavor_catalog/website/src/content/entries/CR002.json    |  465 +++++++++++++
 flavor_catalog/website/src/content/entries/CR003.json    |  369 +++++++++++
 flavor_catalog/website/src/content/entries/CR004.json    |  253 +++++++
 flavor_catalog/website/src/content/entries/CR005.json    |  381 +++++++++++
 flavor_catalog/website/src/content/entries/CR006.json    |  396 +++++++++++
 flavor_catalog/website/src/content/entries/CR007.json    |  292 ++++++++
 flavor_catalog/website/src/content/entries/CR008.json    |  264 ++++++++
 flavor_catalog/website/src/content/entries/CR009.json    |  363 ++++++++++
 flavor_catalog/website/src/content/entries/CR010.json    |  318 +++++++++
 flavor_catalog/website/src/content/entries/CR011.json    |  340 ++++++++++
 flavor_catalog/website/src/content/entries/CR012.json    |  365 ++++++++++
 flavor_catalog/website/src/content/entries/CR013.json    |  340 ++++++++++
 flavor_catalog/website/src/content/entries/CR014.json    |  455 +++++++++++++
 flavor_catalog/website/src/content/entries/E001.json     |  218 ++++++
 flavor_catalog/website/src/content/entries/E002.json     |  278 ++++++++
 flavor_catalog/website/src/content/entries/E004.json     |  193 ++++++
 flavor_catalog/website/src/content/entries/E006.json     |  232 +++++++
 flavor_catalog/website/src/content/entries/E007.json     |  232 +++++++
 flavor_catalog/website/src/content/entries/E008.json     |  144 ++++
 flavor_catalog/website/src/content/entries/E009.json     |  329 +++++++++
 flavor_catalog/website/src/content/entries/EW001.json    |  147 +++++
 flavor_catalog/website/src/content/entries/EW002.json    |  156 +++++
 flavor_catalog/website/src/content/entries/EW003.json    |  152 +++++
 flavor_catalog/website/src/content/entries/K001.json     |  171 +++++
 flavor_catalog/website/src/content/entries/K002.json     |  228 +++++++
 flavor_catalog/website/src/content/entries/K003.json     |  148 +++++
 flavor_catalog/website/src/content/entries/K004.json     |  199 ++++++
 flavor_catalog/website/src/content/entries/K005.json     |  158 +++++
 flavor_catalog/website/src/content/entries/K006.json     |  130 ++++
 flavor_catalog/website/src/content/entries/K008.json     |  464 +++++++++++++
 flavor_catalog/website/src/content/entries/K009.json     |  384 +++++++++++
 flavor_catalog/website/src/content/entries/K010.json     |  231 +++++++
 flavor_catalog/website/src/content/entries/K012.json     |  142 ++++
 flavor_catalog/website/src/content/entries/K013.json     |  133 ++++
 flavor_catalog/website/src/content/entries/K017.json     |  229 +++++++
 flavor_catalog/website/src/content/entries/K018.json     |  192 ++++++
 flavor_catalog/website/src/content/entries/K019.json     |  253 +++++++
 flavor_catalog/website/src/content/entries/K020.json     |  280 ++++++++
 flavor_catalog/website/src/content/entries/K021.json     |  292 ++++++++
 flavor_catalog/website/src/content/entries/L001.json     |  264 ++++++++
 flavor_catalog/website/src/content/entries/L002.json     |  187 ++++++
 flavor_catalog/website/src/content/entries/L003.json     |  154 +++++
 flavor_catalog/website/src/content/entries/L004.json     |  153 +++++
 flavor_catalog/website/src/content/entries/L005.json     |  247 +++++++
 flavor_catalog/website/src/content/entries/L006.json     |  185 ++++++
 flavor_catalog/website/src/content/entries/L007.json     |  211 ++++++
 flavor_catalog/website/src/content/entries/L008.json     |  206 ++++++
 flavor_catalog/website/src/content/entries/L009.json     |  188 ++++++
 flavor_catalog/website/src/content/entries/L010.json     |  188 ++++++
 flavor_catalog/website/src/content/entries/L023.json     |  268 ++++++++
 flavor_catalog/website/src/content/entries/T001.json     |  306 +++++++++
 flavor_catalog/website/src/content/entries/T002.json     |  306 +++++++++
 flavor_catalog/website/src/content/entries/T003.json     |  315 +++++++++
 flavor_catalog/website/src/content/entries/T004.json     |  232 +++++++
 flavor_catalog/website/src/content/entries/T005.json     |  346 ++++++++++
 flavor_catalog/website/src/content/entries/T006.json     |  323 +++++++++
 flavor_catalog/website/src/content/entries/T007.json     |  239 +++++++
 flavor_catalog/website/src/content/entries/T008.json     |  268 ++++++++
 flavor_catalog/website/src/content/entries/T010.json     |  154 +++++
 flavor_catalog/website/src/content/entries/T012.json     |  158 +++++
 flavor_catalog/website/src/content/entries/T014.json     |  217 ++++++
 flavor_catalog/website/src/content/entries/T015.json     |  249 +++++++
 flavor_catalog/website/src/content/entries/T016.json     |  264 ++++++++
 flavor_catalog/website/src/content/entries/T017.json     |  215 ++++++
 flavor_catalog/website/src/content/entries/T018.json     |  351 ++++++++++
 flavor_catalog/website/src/content/entries/T019.json     |  291 ++++++++
 flavor_catalog/website/src/content/entries/T020.json     |  281 ++++++++
 flavor_catalog/website/src/content/families.json         |   50 ++
 flavor_catalog/website/src/layouts/BaseLayout.astro      |  124 ++++
 flavor_catalog/website/src/pages/entries/[id].astro      |  214 ++++++
 flavor_catalog/website/src/pages/families/[family].astro |   42 ++
 flavor_catalog/website/src/pages/index.astro             |   74 +++
 flavor_catalog/website/src/styles/global.css             |  512 ++++++++++++++
 flavor_catalog/website/tsconfig.json                     |    5 +
 118 files changed, 33060 insertions(+)
```

### 371. `6ffbb33` — flavor-catalog-website(phase-2): codex citation-anchor resolution for 102 entries
- SHA: `6ffbb338d8bdd0991b2844d1c398d02d9228cf08`
- Message: flavor-catalog-website(phase-2): codex citation-anchor resolution for 102 entries
- Physical/numerical summary: Resolved citation anchors for all 102 catalog entries, adding website-side source-target metadata without changing source physics.
- Files touched (`git show --stat`):
```text
 flavor_catalog/website/WEBSITE_RUNBOOK.md                         |  13 ++-
 flavor_catalog/website/_data/citation_anchors/B001.yaml           |  77 +++++++++++++
 flavor_catalog/website/_data/citation_anchors/B002.yaml           | 127 +++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/B003.yaml           | 211 ++++++++++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/B004.yaml           |  80 +++++++++++++
 flavor_catalog/website/_data/citation_anchors/B005.yaml           | 254 +++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/B006.yaml           | 412 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/B007.yaml           | 488 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/B008.yaml           |  93 +++++++++++++++
 flavor_catalog/website/_data/citation_anchors/B009.yaml           | 157 +++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/B011.yaml           |  61 ++++++++++
 flavor_catalog/website/_data/citation_anchors/B012.yaml           | 154 +++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/B013.yaml           | 284 +++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/B014.yaml           | 338 ++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/B015.yaml           | 197 ++++++++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/B016.yaml           |  42 +++++++
 flavor_catalog/website/_data/citation_anchors/B017.yaml           | 155 +++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/B018.yaml           |  98 ++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/B019.yaml           | 163 ++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/B021.yaml           | 176 ++++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/B022.yaml           | 180 +++++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/B023.yaml           | 176 ++++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/B025.yaml           | 203 +++++++++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/B026.yaml           | 310 ++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/B032.yaml           | 445 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/B033.yaml           | 154 +++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/B034.yaml           | 271 +++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/C001.yaml           | 115 +++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/C002.yaml           |  70 ++++++++++++
 flavor_catalog/website/_data/citation_anchors/C003.yaml           | 119 +++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/C004.yaml           |  79 +++++++++++++
 flavor_catalog/website/_data/citation_anchors/C005.yaml           |  48 ++++++++
 flavor_catalog/website/_data/citation_anchors/C006.yaml           |  36 ++++++
 flavor_catalog/website/_data/citation_anchors/C007.yaml           | 102 +++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/C008.yaml           |  60 ++++++++++
 flavor_catalog/website/_data/citation_anchors/CR001.yaml          | 135 ++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/CR002.yaml          | 175 ++++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/CR003.yaml          |  80 +++++++++++++
 flavor_catalog/website/_data/citation_anchors/CR004.yaml          |  42 +++++++
 flavor_catalog/website/_data/citation_anchors/CR005.yaml          |  97 ++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/CR006.yaml          |  96 ++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/CR007.yaml          |  79 +++++++++++++
 flavor_catalog/website/_data/citation_anchors/CR008.yaml          |  37 ++++++
 flavor_catalog/website/_data/citation_anchors/CR009.yaml          | 156 +++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/CR010.yaml          | 100 ++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/CR011.yaml          |  41 +++++++
 flavor_catalog/website/_data/citation_anchors/CR012.yaml          |  76 +++++++++++++
 flavor_catalog/website/_data/citation_anchors/CR013.yaml          | 154 +++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/CR014.yaml          |  71 ++++++++++++
 flavor_catalog/website/_data/citation_anchors/E001.yaml           |  59 ++++++++++
 flavor_catalog/website/_data/citation_anchors/E002.yaml           | 194 +++++++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/E004.yaml           |  61 ++++++++++
 flavor_catalog/website/_data/citation_anchors/E006.yaml           |  82 +++++++++++++
 flavor_catalog/website/_data/citation_anchors/E007.yaml           | 187 ++++++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/E008.yaml           | 142 +++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/E009.yaml           |  97 ++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/EW001.yaml          | 390 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/EW002.yaml          | 262 ++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/EW003.yaml          | 232 +++++++++++++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/K001.yaml           | 145 +++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/K002.yaml           |  79 +++++++++++++
 flavor_catalog/website/_data/citation_anchors/K003.yaml           |  23 ++++
 flavor_catalog/website/_data/citation_anchors/K004.yaml           |  91 +++++++++++++++
 flavor_catalog/website/_data/citation_anchors/K005.yaml           |  20 ++++
 flavor_catalog/website/_data/citation_anchors/K006.yaml           |  23 ++++
 flavor_catalog/website/_data/citation_anchors/K008.yaml           | 877 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/K009.yaml           | 445 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/K010.yaml           | 348 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/K012.yaml           | 197 ++++++++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/K013.yaml           |  23 ++++
 flavor_catalog/website/_data/citation_anchors/K017.yaml           |  80 +++++++++++++
 flavor_catalog/website/_data/citation_anchors/K018.yaml           | 259 +++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/K019.yaml           |  55 +++++++++
 flavor_catalog/website/_data/citation_anchors/K020.yaml           |  96 ++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/K021.yaml           | 149 ++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/L001.yaml           |  48 ++++++++
 flavor_catalog/website/_data/citation_anchors/L002.yaml           |  33 ++++++
 flavor_catalog/website/_data/citation_anchors/L003.yaml           |  18 +++
 flavor_catalog/website/_data/citation_anchors/L004.yaml           |  18 +++
 flavor_catalog/website/_data/citation_anchors/L005.yaml           |  55 +++++++++
 flavor_catalog/website/_data/citation_anchors/L006.yaml           |  66 +++++++++++
 flavor_catalog/website/_data/citation_anchors/L007.yaml           |  55 +++++++++
 flavor_catalog/website/_data/citation_anchors/L008.yaml           |  55 +++++++++
 flavor_catalog/website/_data/citation_anchors/L009.yaml           |  42 +++++++
 flavor_catalog/website/_data/citation_anchors/L010.yaml           |  51 +++++++++
 flavor_catalog/website/_data/citation_anchors/L023.yaml           | 142 +++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/T001.yaml           | 190 +++++++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/T002.yaml           | 181 +++++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/T003.yaml           | 137 ++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/T004.yaml           |  80 +++++++++++++
 flavor_catalog/website/_data/citation_anchors/T005.yaml           | 210 ++++++++++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/T006.yaml           | 172 ++++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/T007.yaml           | 116 +++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/T008.yaml           | 154 +++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/T010.yaml           |  90 +++++++++++++++
 flavor_catalog/website/_data/citation_anchors/T012.yaml           |  61 ++++++++++
 flavor_catalog/website/_data/citation_anchors/T014.yaml           |  61 ++++++++++
 flavor_catalog/website/_data/citation_anchors/T015.yaml           | 144 +++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/T016.yaml           | 164 ++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/T017.yaml           | 106 +++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/T018.yaml           | 221 +++++++++++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/T019.yaml           | 168 +++++++++++++++++++++++++++
 flavor_catalog/website/_data/citation_anchors/T020.yaml           | 138 ++++++++++++++++++++++
 flavor_catalog/website/scripts/resolve_beauty_citation_anchors.py | 917 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/website/scripts/resolve_citation_anchors.py        | 790 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 105 files changed, 16585 insertions(+), 1 deletion(-)
```

### 372. `d62db33` — flavor-catalog-website(phase-3): citation modal, family pages, search + filters
- SHA: `d62db33c6cf5301cf6e0a84598cfd7a44ec49847`
- Message: flavor-catalog-website(phase-3): citation modal, family pages, search + filters
- Physical/numerical summary: Generalized the website to all entries with family pages, citation modals, search, and filters; presentation-layer only.
- Files touched (`git show --stat`):
```text
 flavor_catalog/website/.gitignore                         |   3 +
 flavor_catalog/website/WEBSITE_RUNBOOK.md                 |   7 ++
 flavor_catalog/website/package-lock.json                  | 120 +++++++++++++++++++++++
 flavor_catalog/website/package.json                       |   5 +-
 flavor_catalog/website/scripts/ingest_catalog.py          |  90 +++++++++++++++++
 flavor_catalog/website/scripts/screenshot_cdp.mjs         |  99 +++++++++++++++++++
 flavor_catalog/website/src/components/CitationModal.astro | 211 ++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/website/src/components/EntryTable.astro    |  99 +++++++++++++++++++
 flavor_catalog/website/src/content.config.ts              |  39 ++++++++
 flavor_catalog/website/src/content/catalog_index.json     | 822 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-------------------
 flavor_catalog/website/src/content/entries/B001.json      |  80 ++++++++++++++-
 flavor_catalog/website/src/content/entries/B002.json      | 136 +++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/B003.json      | 208 ++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/B004.json      |  92 +++++++++++++++++-
 flavor_catalog/website/src/content/entries/B005.json      | 200 +++++++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/B006.json      | 287 +++++++++++++++++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/B007.json      | 328 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/B008.json      |  88 ++++++++++++++++-
 flavor_catalog/website/src/content/entries/B009.json      | 180 +++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/B011.json      |  72 +++++++++++++-
 flavor_catalog/website/src/content/entries/B012.json      | 148 +++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/B013.json      | 288 +++++++++++++++++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/B014.json      | 348 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/B015.json      | 204 +++++++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/B016.json      |  52 +++++++++-
 flavor_catalog/website/src/content/entries/B017.json      | 160 +++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/B018.json      | 100 ++++++++++++++++++-
 flavor_catalog/website/src/content/entries/B019.json      | 164 ++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/B021.json      | 192 +++++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/B022.json      | 148 +++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/B023.json      | 136 +++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/B025.json      | 216 ++++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/B026.json      | 299 +++++++++++++++++++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/B032.json      | 348 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/B033.json      | 150 +++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/B034.json      | 225 +++++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/C001.json      | 152 ++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/C002.json      |  97 ++++++++++++++++++-
 flavor_catalog/website/src/content/entries/C003.json      | 162 ++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/C004.json      | 102 ++++++++++++++++++-
 flavor_catalog/website/src/content/entries/C005.json      |  72 +++++++++++++-
 flavor_catalog/website/src/content/entries/C006.json      |  52 +++++++++-
 flavor_catalog/website/src/content/entries/C007.json      | 130 ++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/C008.json      |  60 +++++++++++-
 flavor_catalog/website/src/content/entries/CR001.json     | 128 +++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/CR002.json     | 192 +++++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/CR003.json     | 102 ++++++++++++++++++-
 flavor_catalog/website/src/content/entries/CR004.json     |  62 +++++++++++-
 flavor_catalog/website/src/content/entries/CR005.json     | 112 ++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/CR006.json     | 112 ++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/CR007.json     |  87 ++++++++++++++++-
 flavor_catalog/website/src/content/entries/CR008.json     |  52 +++++++++-
 flavor_catalog/website/src/content/entries/CR009.json     | 172 ++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/CR010.json     | 132 ++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/CR011.json     |  52 +++++++++-
 flavor_catalog/website/src/content/entries/CR012.json     |  92 +++++++++++++++++-
 flavor_catalog/website/src/content/entries/CR013.json     | 168 +++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/CR014.json     |  86 ++++++++++++++++-
 flavor_catalog/website/src/content/entries/E001.json      |  72 +++++++++++++-
 flavor_catalog/website/src/content/entries/E002.json      | 212 +++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/E004.json      |  72 +++++++++++++-
 flavor_catalog/website/src/content/entries/E006.json      |  92 +++++++++++++++++-
 flavor_catalog/website/src/content/entries/E007.json      | 124 +++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/E008.json      | 152 ++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/E009.json      | 112 ++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/EW001.json     | 379 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/EW002.json     | 297 +++++++++++++++++++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/EW003.json     | 252 ++++++++++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/K001.json      | 146 +++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/K002.json      |  92 +++++++++++++++++-
 flavor_catalog/website/src/content/entries/K003.json      |  32 +++++-
 flavor_catalog/website/src/content/entries/K004.json      | 107 +++++++++++++++++++-
 flavor_catalog/website/src/content/entries/K005.json      |  32 +++++-
 flavor_catalog/website/src/content/entries/K006.json      |  32 +++++-
 flavor_catalog/website/src/content/entries/K008.json      | 632 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/K009.json      | 372 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/K010.json      | 236 +++++++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/K012.json      | 172 ++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/K013.json      |  32 +++++-
 flavor_catalog/website/src/content/entries/K017.json      |  92 +++++++++++++++++-
 flavor_catalog/website/src/content/entries/K018.json      | 264 +++++++++++++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/K019.json      |  72 +++++++++++++-
 flavor_catalog/website/src/content/entries/K020.json      | 102 ++++++++++++++++++-
 flavor_catalog/website/src/content/entries/K021.json      | 148 +++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/L001.json      |  71 +++++++++++++-
 flavor_catalog/website/src/content/entries/L002.json      |  52 +++++++++-
 flavor_catalog/website/src/content/entries/L003.json      |  32 +++++-
 flavor_catalog/website/src/content/entries/L004.json      |  32 +++++-
 flavor_catalog/website/src/content/entries/L005.json      |  80 ++++++++++++++-
 flavor_catalog/website/src/content/entries/L006.json      |  84 +++++++++++++++-
 flavor_catalog/website/src/content/entries/L007.json      |  87 ++++++++++++++++-
 flavor_catalog/website/src/content/entries/L008.json      |  87 ++++++++++++++++-
 flavor_catalog/website/src/content/entries/L009.json      |  67 ++++++++++++-
 flavor_catalog/website/src/content/entries/L010.json      |  76 ++++++++++++++-
 flavor_catalog/website/src/content/entries/L023.json      | 156 +++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/T001.json      | 212 +++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/T002.json      | 207 ++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/T003.json      | 152 ++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/T004.json      |  92 +++++++++++++++++-
 flavor_catalog/website/src/content/entries/T005.json      | 232 +++++++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/T006.json      | 192 +++++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/T007.json      | 142 ++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/T008.json      | 182 +++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/T010.json      | 107 +++++++++++++++++++-
 flavor_catalog/website/src/content/entries/T012.json      |  72 +++++++++++++-
 flavor_catalog/website/src/content/entries/T014.json      |  72 +++++++++++++-
 flavor_catalog/website/src/content/entries/T015.json      | 167 +++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/T016.json      | 186 ++++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/T017.json      | 127 +++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/T018.json      | 257 +++++++++++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/T019.json      | 192 +++++++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/content/entries/T020.json      | 162 ++++++++++++++++++++++++++++++-
 flavor_catalog/website/src/layouts/BaseLayout.astro       |  11 ++-
 flavor_catalog/website/src/pages/browse.astro             | 209 +++++++++++++++++++++++++++++++++++++++
 flavor_catalog/website/src/pages/entries/[id].astro       | 219 ++++++++++++++++++++++++++++++++++-------
 flavor_catalog/website/src/pages/families/[family].astro  | 106 +++++++++++++++-----
 flavor_catalog/website/src/pages/index.astro              |  40 ++++----
 flavor_catalog/website/src/styles/global.css              | 627 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
 118 files changed, 17550 insertions(+), 286 deletions(-)
```

### 373. `5f087fd` — flavor-catalog-website(phase-5): methodology page + Cloudflare config + README
- SHA: `5f087fdb4ca27c3378e2241a5736721a6b733f07`
- Message: flavor-catalog-website(phase-5): methodology page + Cloudflare config + README
- Physical/numerical summary: Added methodology page, Cloudflare config, and website README; deployment/docs layer only.
- Files touched (`git show --stat`):
```text
 flavor_catalog/website/.node-version                |   1 +
 flavor_catalog/website/README.md                    | 227 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/website/WEBSITE_RUNBOOK.md           |   4 +-
 flavor_catalog/website/cloudflare-pages.config.md   | 111 ++++++++++++++++++++++++++++++++
 flavor_catalog/website/public/_headers              |  19 ++++++
 flavor_catalog/website/public/_redirects            |  11 ++++
 flavor_catalog/website/src/layouts/BaseLayout.astro |   1 +
 flavor_catalog/website/src/pages/index.astro        |   1 +
 flavor_catalog/website/src/pages/methodology.astro  | 560 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 9 files changed, 934 insertions(+), 1 deletion(-)
```

### 374. `344e147` — flavor-catalog-website(phase-6): LaTeX delimiters + UI polish + methodology rewrite
- SHA: `344e147ff1f9e00e13451763e32727bcc08cc4c4`
- Message: flavor-catalog-website(phase-6): LaTeX delimiters + UI polish + methodology rewrite
- Physical/numerical summary: Fixed LaTeX delimiter handling and UI polish, improving rendering of existing catalog numerics without source-data changes.
- Files touched (`git show --stat`):
```text
 flavor_catalog/website/WEBSITE_RUNBOOK.md              |   5 +-
 flavor_catalog/website/src/components/Badges.astro     |   2 +-
 flavor_catalog/website/src/components/EntryTable.astro |  63 +++++++++++----
 flavor_catalog/website/src/layouts/BaseLayout.astro    |  29 ++++---
 flavor_catalog/website/src/pages/entries/[id].astro    |   6 +-
 flavor_catalog/website/src/pages/index.astro           |  69 +++++++++--------
 flavor_catalog/website/src/pages/methodology.astro     | 660 +++++++++++++++++++++++++++++++++++++++++++++++++++++--------------------------------------------------------------------------------------------------------
 flavor_catalog/website/src/styles/global.css           |  37 +++++++--
 8 files changed, 372 insertions(+), 499 deletions(-)
```

### 375. `5180988` — flavor-catalog-website(phase-7): shorthand-to-LaTeX normalizer for standard_notation
- SHA: `518098805bf6e3f165eda97bf55c88ecc8122624`
- Message: flavor-catalog-website(phase-7): shorthand-to-LaTeX normalizer for standard_notation
- Physical/numerical summary: Added shorthand-to-LaTeX normalization for standard_notation, changing display rendering but not catalog values.
- Files touched (`git show --stat`):
```text
 flavor_catalog/website/src/components/EntryTable.astro |  13 ++-----------
 flavor_catalog/website/src/lib/notation.ts             | 148 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/website/src/pages/entries/[id].astro    |   3 ++-
 3 files changed, 152 insertions(+), 12 deletions(-)
```

### 376. `eb6b8b7` — flavor-catalog-website(phase-8a): codex constraint-priority ranking + sonnet UI fixes
- SHA: `eb6b8b767115a91619a77f1e429236ec9ebf2951`
- Message: flavor-catalog-website(phase-8a): codex constraint-priority ranking + sonnet UI fixes
- Physical/numerical summary: Added constraint-priority ranking workflow and related UI fixes, classifying 102 entries for website display.
- Files touched (`git show --stat`):
```text
 flavor_catalog/website/_data/priority/B001.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/B002.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/B003.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/B004.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/B005.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/B006.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/B007.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/B008.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/B009.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/B011.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/B012.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/B013.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/B014.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/B015.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/B016.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/B017.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/B018.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/B019.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/B021.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/B022.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/B023.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/B025.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/B026.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/B032.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/B033.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/B034.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/C001.yaml          |  12 ++++++++++++
 flavor_catalog/website/_data/priority/C002.yaml          |  13 +++++++++++++
 flavor_catalog/website/_data/priority/C003.yaml          |  12 ++++++++++++
 flavor_catalog/website/_data/priority/C004.yaml          |  12 ++++++++++++
 flavor_catalog/website/_data/priority/C005.yaml          |  12 ++++++++++++
 flavor_catalog/website/_data/priority/C006.yaml          |  12 ++++++++++++
 flavor_catalog/website/_data/priority/C007.yaml          |  11 +++++++++++
 flavor_catalog/website/_data/priority/C008.yaml          |  12 ++++++++++++
 flavor_catalog/website/_data/priority/CR001.yaml         |  12 ++++++++++++
 flavor_catalog/website/_data/priority/CR002.yaml         |  11 +++++++++++
 flavor_catalog/website/_data/priority/CR003.yaml         |  12 ++++++++++++
 flavor_catalog/website/_data/priority/CR004.yaml         |  12 ++++++++++++
 flavor_catalog/website/_data/priority/CR005.yaml         |  12 ++++++++++++
 flavor_catalog/website/_data/priority/CR006.yaml         |  12 ++++++++++++
 flavor_catalog/website/_data/priority/CR007.yaml         |  12 ++++++++++++
 flavor_catalog/website/_data/priority/CR008.yaml         |  11 +++++++++++
 flavor_catalog/website/_data/priority/CR009.yaml         |  12 ++++++++++++
 flavor_catalog/website/_data/priority/CR010.yaml         |  12 ++++++++++++
 flavor_catalog/website/_data/priority/CR011.yaml         |  11 +++++++++++
 flavor_catalog/website/_data/priority/CR012.yaml         |  11 +++++++++++
 flavor_catalog/website/_data/priority/CR013.yaml         |  11 +++++++++++
 flavor_catalog/website/_data/priority/CR014.yaml         |  12 ++++++++++++
 flavor_catalog/website/_data/priority/E001.yaml          |  12 ++++++++++++
 flavor_catalog/website/_data/priority/E002.yaml          |  11 +++++++++++
 flavor_catalog/website/_data/priority/E004.yaml          |  12 ++++++++++++
 flavor_catalog/website/_data/priority/E006.yaml          |  13 +++++++++++++
 flavor_catalog/website/_data/priority/E007.yaml          |  12 ++++++++++++
 flavor_catalog/website/_data/priority/E008.yaml          |  13 +++++++++++++
 flavor_catalog/website/_data/priority/E009.yaml          |  12 ++++++++++++
 flavor_catalog/website/_data/priority/EW001.yaml         |   8 ++++++++
 flavor_catalog/website/_data/priority/EW002.yaml         |   8 ++++++++
 flavor_catalog/website/_data/priority/EW003.yaml         |   8 ++++++++
 flavor_catalog/website/_data/priority/K001.yaml          |  11 +++++++++++
 flavor_catalog/website/_data/priority/K002.yaml          |  12 ++++++++++++
 flavor_catalog/website/_data/priority/K003.yaml          |  12 ++++++++++++
 flavor_catalog/website/_data/priority/K004.yaml          |  12 ++++++++++++
 flavor_catalog/website/_data/priority/K005.yaml          |  12 ++++++++++++
 flavor_catalog/website/_data/priority/K006.yaml          |  12 ++++++++++++
 flavor_catalog/website/_data/priority/K008.yaml          |  12 ++++++++++++
 flavor_catalog/website/_data/priority/K009.yaml          |  12 ++++++++++++
 flavor_catalog/website/_data/priority/K010.yaml          |  11 +++++++++++
 flavor_catalog/website/_data/priority/K012.yaml          |  12 ++++++++++++
 flavor_catalog/website/_data/priority/K013.yaml          |  12 ++++++++++++
 flavor_catalog/website/_data/priority/K017.yaml          |  12 ++++++++++++
 flavor_catalog/website/_data/priority/K018.yaml          |  11 +++++++++++
 flavor_catalog/website/_data/priority/K019.yaml          |  12 ++++++++++++
 flavor_catalog/website/_data/priority/K020.yaml          |  12 ++++++++++++
 flavor_catalog/website/_data/priority/K021.yaml          |  11 +++++++++++
 flavor_catalog/website/_data/priority/L001.yaml          |  12 ++++++++++++
 flavor_catalog/website/_data/priority/L002.yaml          |  12 ++++++++++++
 flavor_catalog/website/_data/priority/L003.yaml          |  13 +++++++++++++
 flavor_catalog/website/_data/priority/L004.yaml          |  12 ++++++++++++
 flavor_catalog/website/_data/priority/L005.yaml          |  12 ++++++++++++
 flavor_catalog/website/_data/priority/L006.yaml          |  12 ++++++++++++
 flavor_catalog/website/_data/priority/L007.yaml          |  11 +++++++++++
 flavor_catalog/website/_data/priority/L008.yaml          |  11 +++++++++++
 flavor_catalog/website/_data/priority/L009.yaml          |  11 +++++++++++
 flavor_catalog/website/_data/priority/L010.yaml          |  11 +++++++++++
 flavor_catalog/website/_data/priority/L023.yaml          |  12 ++++++++++++
 flavor_catalog/website/_data/priority/T001.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/T002.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/T003.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/T004.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/T005.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/T006.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/T007.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/T008.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/T010.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/T012.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/T014.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/T015.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/T016.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/T017.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/T018.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/T019.yaml          |   8 ++++++++
 flavor_catalog/website/_data/priority/T020.yaml          |   8 ++++++++
 flavor_catalog/website/src/components/EntryTable.astro   |   2 --
 flavor_catalog/website/src/lib/prose.ts                  | 148 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/website/src/pages/browse.astro            |  15 +++------------
 flavor_catalog/website/src/pages/entries/[id].astro      |  13 +++++++------
 flavor_catalog/website/src/pages/families/[family].astro |  16 ++++++++++++++++
 flavor_catalog/website/src/pages/index.astro             |  12 ------------
 108 files changed, 1203 insertions(+), 32 deletions(-)
```

### 377. `0833d0e` — flavor-catalog-website(phase-9): priority UI + site-wide LaTeX sweep + em-dash purge + jargon scrub + methodology rewrite
- SHA: `0833d0e404c24df34fc84a285ca01256130aea74`
- Message: flavor-catalog-website(phase-9): priority UI + site-wide LaTeX sweep + em-dash purge + jargon scrub + methodology rewrite
- Physical/numerical summary: Wired priority UI and swept LaTeX/jargon/em-dashes site-wide; presentation cleanup over existing catalog data.
- Files touched (`git show --stat`):
```text
 flavor_catalog/website/WEBSITE_RUNBOOK.md                 |  12 ++---
 flavor_catalog/website/cloudflare-pages.config.md         |  18 +++----
 flavor_catalog/website/scripts/ingest_catalog.py          |  61 ++++++++++++++++++++---
 flavor_catalog/website/src/components/Badges.astro        |  19 ++++++--
 flavor_catalog/website/src/components/CitationModal.astro |   2 +-
 flavor_catalog/website/src/components/EntryTable.astro    |  30 ++++++++++--
 flavor_catalog/website/src/content.config.ts              |   5 ++
 flavor_catalog/website/src/content/catalog_index.json     | 414 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++--------------------------------------
 flavor_catalog/website/src/content/entries/B001.json      |   5 +-
 flavor_catalog/website/src/content/entries/B002.json      |   5 +-
 flavor_catalog/website/src/content/entries/B003.json      |   5 +-
 flavor_catalog/website/src/content/entries/B004.json      |   5 +-
 flavor_catalog/website/src/content/entries/B005.json      |   5 +-
 flavor_catalog/website/src/content/entries/B006.json      |   5 +-
 flavor_catalog/website/src/content/entries/B007.json      |   5 +-
 flavor_catalog/website/src/content/entries/B008.json      |   5 +-
 flavor_catalog/website/src/content/entries/B009.json      |   5 +-
 flavor_catalog/website/src/content/entries/B011.json      |   5 +-
 flavor_catalog/website/src/content/entries/B012.json      |   5 +-
 flavor_catalog/website/src/content/entries/B013.json      |   5 +-
 flavor_catalog/website/src/content/entries/B014.json      |   5 +-
 flavor_catalog/website/src/content/entries/B015.json      |   5 +-
 flavor_catalog/website/src/content/entries/B016.json      |   5 +-
 flavor_catalog/website/src/content/entries/B017.json      |   5 +-
 flavor_catalog/website/src/content/entries/B018.json      |   5 +-
 flavor_catalog/website/src/content/entries/B019.json      |   5 +-
 flavor_catalog/website/src/content/entries/B021.json      |   5 +-
 flavor_catalog/website/src/content/entries/B022.json      |   5 +-
 flavor_catalog/website/src/content/entries/B023.json      |   5 +-
 flavor_catalog/website/src/content/entries/B025.json      |   5 +-
 flavor_catalog/website/src/content/entries/B026.json      |   5 +-
 flavor_catalog/website/src/content/entries/B032.json      |   5 +-
 flavor_catalog/website/src/content/entries/B033.json      |   5 +-
 flavor_catalog/website/src/content/entries/B034.json      |   5 +-
 flavor_catalog/website/src/content/entries/C001.json      |   5 +-
 flavor_catalog/website/src/content/entries/C002.json      |   5 +-
 flavor_catalog/website/src/content/entries/C003.json      |   5 +-
 flavor_catalog/website/src/content/entries/C004.json      |   5 +-
 flavor_catalog/website/src/content/entries/C005.json      |   5 +-
 flavor_catalog/website/src/content/entries/C006.json      |   5 +-
 flavor_catalog/website/src/content/entries/C007.json      |   5 +-
 flavor_catalog/website/src/content/entries/C008.json      |   5 +-
 flavor_catalog/website/src/content/entries/CR001.json     |   5 +-
 flavor_catalog/website/src/content/entries/CR002.json     |   5 +-
 flavor_catalog/website/src/content/entries/CR003.json     |   5 +-
 flavor_catalog/website/src/content/entries/CR004.json     |   5 +-
 flavor_catalog/website/src/content/entries/CR005.json     |   5 +-
 flavor_catalog/website/src/content/entries/CR006.json     |   5 +-
 flavor_catalog/website/src/content/entries/CR007.json     |   5 +-
 flavor_catalog/website/src/content/entries/CR008.json     |   5 +-
 flavor_catalog/website/src/content/entries/CR009.json     |   5 +-
 flavor_catalog/website/src/content/entries/CR010.json     |   5 +-
 flavor_catalog/website/src/content/entries/CR011.json     |   5 +-
 flavor_catalog/website/src/content/entries/CR012.json     |   5 +-
 flavor_catalog/website/src/content/entries/CR013.json     |   5 +-
 flavor_catalog/website/src/content/entries/CR014.json     |   5 +-
 flavor_catalog/website/src/content/entries/E001.json      |   5 +-
 flavor_catalog/website/src/content/entries/E002.json      |   5 +-
 flavor_catalog/website/src/content/entries/E004.json      |   5 +-
 flavor_catalog/website/src/content/entries/E006.json      |   5 +-
 flavor_catalog/website/src/content/entries/E007.json      |   5 +-
 flavor_catalog/website/src/content/entries/E008.json      |   5 +-
 flavor_catalog/website/src/content/entries/E009.json      |   5 +-
 flavor_catalog/website/src/content/entries/EW001.json     |   5 +-
 flavor_catalog/website/src/content/entries/EW002.json     |   5 +-
 flavor_catalog/website/src/content/entries/EW003.json     |   5 +-
 flavor_catalog/website/src/content/entries/K001.json      |   5 +-
 flavor_catalog/website/src/content/entries/K002.json      |   5 +-
 flavor_catalog/website/src/content/entries/K003.json      |   5 +-
 flavor_catalog/website/src/content/entries/K004.json      |   5 +-
 flavor_catalog/website/src/content/entries/K005.json      |   5 +-
 flavor_catalog/website/src/content/entries/K006.json      |   5 +-
 flavor_catalog/website/src/content/entries/K008.json      |   5 +-
 flavor_catalog/website/src/content/entries/K009.json      |   5 +-
 flavor_catalog/website/src/content/entries/K010.json      |   5 +-
 flavor_catalog/website/src/content/entries/K012.json      |   5 +-
 flavor_catalog/website/src/content/entries/K013.json      |   5 +-
 flavor_catalog/website/src/content/entries/K017.json      |   5 +-
 flavor_catalog/website/src/content/entries/K018.json      |   5 +-
 flavor_catalog/website/src/content/entries/K019.json      |   5 +-
 flavor_catalog/website/src/content/entries/K020.json      |   5 +-
 flavor_catalog/website/src/content/entries/K021.json      |   5 +-
 flavor_catalog/website/src/content/entries/L001.json      |   5 +-
 flavor_catalog/website/src/content/entries/L002.json      |   5 +-
 flavor_catalog/website/src/content/entries/L003.json      |   5 +-
 flavor_catalog/website/src/content/entries/L004.json      |   5 +-
 flavor_catalog/website/src/content/entries/L005.json      |   5 +-
 flavor_catalog/website/src/content/entries/L006.json      |   5 +-
 flavor_catalog/website/src/content/entries/L007.json      |   5 +-
 flavor_catalog/website/src/content/entries/L008.json      |   5 +-
 flavor_catalog/website/src/content/entries/L009.json      |   5 +-
 flavor_catalog/website/src/content/entries/L010.json      |   5 +-
 flavor_catalog/website/src/content/entries/L023.json      |   5 +-
 flavor_catalog/website/src/content/entries/T001.json      |   5 +-
 flavor_catalog/website/src/content/entries/T002.json      |   5 +-
 flavor_catalog/website/src/content/entries/T003.json      |   5 +-
 flavor_catalog/website/src/content/entries/T004.json      |   5 +-
 flavor_catalog/website/src/content/entries/T005.json      |   5 +-
 flavor_catalog/website/src/content/entries/T006.json      |   5 +-
 flavor_catalog/website/src/content/entries/T007.json      |   5 +-
 flavor_catalog/website/src/content/entries/T008.json      |   5 +-
 flavor_catalog/website/src/content/entries/T010.json      |   5 +-
 flavor_catalog/website/src/content/entries/T012.json      |   5 +-
 flavor_catalog/website/src/content/entries/T014.json      |   5 +-
 flavor_catalog/website/src/content/entries/T015.json      |   5 +-
 flavor_catalog/website/src/content/entries/T016.json      |   5 +-
 flavor_catalog/website/src/content/entries/T017.json      |   5 +-
 flavor_catalog/website/src/content/entries/T018.json      |   5 +-
 flavor_catalog/website/src/content/entries/T019.json      |   5 +-
 flavor_catalog/website/src/content/entries/T020.json      |   5 +-
 flavor_catalog/website/src/content/families.json          |  14 +++---
 flavor_catalog/website/src/lib/notation.ts                | 138 ++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/website/src/lib/prose.ts                   | 345 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++---
 flavor_catalog/website/src/pages/browse.astro             |  76 +++++++++++++++++++++++++++--
 flavor_catalog/website/src/pages/entries/[id].astro       |  61 +++++++++++++----------
 flavor_catalog/website/src/pages/families/[family].astro  |  42 +++++++++++-----
 flavor_catalog/website/src/pages/index.astro              |  12 ++++-
 flavor_catalog/website/src/pages/methodology.astro        | 187 +++++++++++++++++++++++++++++++++-------------------------------------
 flavor_catalog/website/src/styles/global.css              |  56 +++++++++++++++++++++
 119 files changed, 1614 insertions(+), 388 deletions(-)
```

### 378. `7053fb7` — Remove tier/anchor UI, fix PDG observable column and unit duplication
- SHA: `7053fb76b6a5da44562ddaa6e47392e9f571e627`
- Message: Remove tier/anchor UI, fix PDG observable column and unit duplication
- Physical/numerical summary: Polished website navigation/rendering and removed secondary/tier/verdict-count UI surfaces; no source physics changed.
- Files touched (`git show --stat`):
```text
 flavor_catalog/website/src/components/EntryTable.astro |  4 +---
 flavor_catalog/website/src/pages/browse.astro          | 14 +++-----------
 flavor_catalog/website/src/pages/entries/[id].astro    | 28 ++--------------------------
 flavor_catalog/website/src/styles/global.css           |  2 --
 4 files changed, 6 insertions(+), 42 deletions(-)
```

### 379. `4eba297` — LaTeX rendering sweep + methodology page rewrite
- SHA: `4eba297d19eed0530fb893f22d061d5f80396b07`
- Message: LaTeX rendering sweep + methodology page rewrite
- Physical/numerical summary: Polished website navigation/rendering and removed secondary/tier/verdict-count UI surfaces; no source physics changed.
- Files touched (`git show --stat`):
```text
 flavor_catalog/website/src/lib/notation.ts         | 133 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-------------
 flavor_catalog/website/src/lib/prose.ts            |  39 ++++++++++++++++++++------
 flavor_catalog/website/src/pages/methodology.astro | 247 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++--------------------------------------------------------------------------
 3 files changed, 278 insertions(+), 141 deletions(-)
```

### 380. `b6864f7` — Update runbook with phases 10-12
- SHA: `b6864f710966d2a9ab31151c7f61df5e608fbae8`
- Message: Update runbook with phases 10-12
- Physical/numerical summary: Updated orchestration runbook state for the named wave/stage; workflow bookkeeping only.
- Files touched (`git show --stat`):
```text
 flavor_catalog/website/WEBSITE_RUNBOOK.md | 3 +++
 1 file changed, 3 insertions(+)
```

### 381. `f485725` — Remove V/P verdict counts and Secondary nav from header
- SHA: `f485725ec821dc6e1aef65cd50160171a70ac77d`
- Message: Remove V/P verdict counts and Secondary nav from header
- Physical/numerical summary: Polished website navigation/rendering and removed secondary/tier/verdict-count UI surfaces; no source physics changed.
- Files touched (`git show --stat`):
```text
 flavor_catalog/website/src/layouts/BaseLayout.astro | 3 ---
 1 file changed, 3 deletions(-)
```

### 382. `5f31f2d` — UI cleanup: remove Secondary, stats strip, methodology nav; polish home
- SHA: `5f31f2d47bf243765aa8b108904a4c36e394bd28`
- Message: UI cleanup: remove Secondary, stats strip, methodology nav; polish home
- Physical/numerical summary: Polished website navigation/rendering and removed secondary/tier/verdict-count UI surfaces; no source physics changed.
- Files touched (`git show --stat`):
```text
 flavor_catalog/website/src/layouts/BaseLayout.astro |   7 -----
 flavor_catalog/website/src/pages/browse.astro       |   2 +-
 flavor_catalog/website/src/pages/index.astro        |  27 ++-----------------
 flavor_catalog/website/src/pages/methodology.astro  | 233 +++++++++++++++++++++++++++++++++++++++++++++++++---------------------------------------------------------------------------------------------------------------
 4 files changed, 74 insertions(+), 195 deletions(-)
```

### 383. `a05832e` — chore(orchestration): commit 24/24 atomic-review reports + ISSUES.md
- SHA: `a05832eea71f2dd8d7af6e20481efbe82c0b517e`
- Message: chore(orchestration): commit 24/24 atomic-review reports + ISSUES.md
- Physical/numerical summary: Committed the 24/24 atomic review reports, review queue, progress ledger, and ISSUES.md issue tracker for R01-R22.
- Files touched (`git show --stat`):
```text
 .orchestration/ISSUES.md          | 484 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 .orchestration/PRE_MERGE_STATE.md |  56 +++++++++++++++++++++
 .orchestration/REVIEW_QUEUE.md    |  30 +++++++++++
 .orchestration/progress.json      | 103 ++++++++++++++++++++++++++++++++++++++
 .orchestration/reviews/R01.md     |  57 +++++++++++++++++++++
 .orchestration/reviews/R02.md     |  66 +++++++++++++++++++++++++
 .orchestration/reviews/R03.md     |  80 ++++++++++++++++++++++++++++++
 .orchestration/reviews/R04.md     | 134 ++++++++++++++++++++++++++++++++++++++++++++++++++
 .orchestration/reviews/R05.md     | 107 ++++++++++++++++++++++++++++++++++++++++
 .orchestration/reviews/R06.md     | 100 +++++++++++++++++++++++++++++++++++++
 .orchestration/reviews/R07.md     |  64 ++++++++++++++++++++++++
 .orchestration/reviews/R08.md     |  75 ++++++++++++++++++++++++++++
 .orchestration/reviews/R09.md     |  91 ++++++++++++++++++++++++++++++++++
 .orchestration/reviews/R10a.md    | 134 ++++++++++++++++++++++++++++++++++++++++++++++++++
 .orchestration/reviews/R10b.md    | 157 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 .orchestration/reviews/R10c.md    | 243 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 .orchestration/reviews/R11.md     |  71 ++++++++++++++++++++++++++
 .orchestration/reviews/R12.md     | 118 ++++++++++++++++++++++++++++++++++++++++++++
 .orchestration/reviews/R13.md     | 115 +++++++++++++++++++++++++++++++++++++++++++
 .orchestration/reviews/R14.md     | 106 +++++++++++++++++++++++++++++++++++++++
 .orchestration/reviews/R15.md     | 115 +++++++++++++++++++++++++++++++++++++++++++
 .orchestration/reviews/R16.md     | 125 ++++++++++++++++++++++++++++++++++++++++++++++
 .orchestration/reviews/R17.md     | 125 ++++++++++++++++++++++++++++++++++++++++++++++
 .orchestration/reviews/R18.md     |  98 ++++++++++++++++++++++++++++++++++++
 .orchestration/reviews/R19.md     | 115 +++++++++++++++++++++++++++++++++++++++++++
 .orchestration/reviews/R20.md     |  88 +++++++++++++++++++++++++++++++++
 .orchestration/reviews/R21.md     |  84 +++++++++++++++++++++++++++++++
 .orchestration/reviews/R22.md     |  83 +++++++++++++++++++++++++++++++
 28 files changed, 3224 insertions(+)
```

## Phase-2 (21-unit cleanup)

### 1. `3ab1f8f` — cleanup(C01a): sync vendored kaon constants in modern/phenomenology.py to post-audit canonical values (R03-I1, R03-I2)
- SHA: `3ab1f8f7ed338be9e35158f13c671b200f6c121d`
- Message: cleanup(C01a): sync vendored kaon constants in modern/phenomenology.py to post-audit canonical values (R03-I1, R03-I2)
- Physical/numerical summary: Updated the modern kaon constants to the post-audit canonical values, tightening epsilon_K ratio behavior without crossing the modern-lane firewall.
- Files touched (`git show --stat`):
```text
 quarkConstraints/modern/phenomenology.py | 16 +++++++++++-----
 tests/test_quark_deltaf2.py              | 38 ++++++++++++++++++++++++++++++++++++++
 2 files changed, 49 insertions(+), 5 deletions(-)
```

### 2. `7899205` — cleanup(C01b): re-export collaborator artifacts post-C01a (R03-I1 followup)
- SHA: `78992057b60bbe090dc05908a96e584b2b6a9eea`
- Message: cleanup(C01b): re-export collaborator artifacts post-C01a (R03-I1 followup)
- Physical/numerical summary: Re-exported collaborator CSV/provenance artifacts after the kaon-constant correction so shared 5/10 TeV point data matches the audited numerics.
- Files touched (`git show --stat`):
```text
 artifacts/collaborator_5tev_points.csv                              | 10 +++++-----
 artifacts/collaborator_5tev_points.provenance.json                  |  4 ++--
 artifacts/collaborator_direct_affine_5_10tev_points.csv             | 10 +++++-----
 artifacts/collaborator_direct_affine_5_10tev_points.provenance.json | 40 ++++++++++++++++++++--------------------
 4 files changed, 32 insertions(+), 32 deletions(-)
```

### 3. `3e1e07c` — cleanup(C02a-code): add --epsilon-k-budget CLI flag + methodology note band paragraph (R03-I3 tasks 1-4)
- SHA: `3e1e07c1a4285a40f63a1294aee8328bc5ed2040`
- Message: cleanup(C02a-code): add --epsilon-k-budget CLI flag + methodology note band paragraph (R03-I3 tasks 1-4)
- Physical/numerical summary: Added the epsilon_K budget-band CLI plumbing and methodology text; central runs remain bit-for-bit identical while low/high scale ratios by 6.7 and 0.2233.
- Files touched (`git show --stat`):
```text
 docs/quark_scan_methodology_note.pdf | Bin 620832 -> 626853 bytes
 docs/quark_scan_methodology_note.tex |  64 +++++++++++++++++++++++++++++++++++++++++++++++
 quarkConstraints/deltaf2.py          |  54 +++++++++++++++++++++++++++++++++++++---
 scripts/run_rs_anarchy.py            |  83 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
 tests/test_run_rs_anarchy.py         | 240 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 5 files changed, 436 insertions(+), 5 deletions(-)
```

### 4. `9dc72d7` — cleanup(C03): Wilson-RG / BMU follow-ups (R04-I1..I4)
- SHA: `9dc72d70cee717ae7383f3f267c9e00b6d9c34a7`
- Message: cleanup(C03): Wilson-RG / BMU follow-ups (R04-I1..I4)
- Physical/numerical summary: Closed Wilson-RG follow-ups with an upper-triangular guard, import cleanup, LR-sign cross-check, and audit-status documentation; numerical audit output unchanged.
- Files touched (`git show --stat`):
```text
 .orchestration/CLEANUP_QUEUE.md       |  29 ++++++++++++++++++++++++
 .orchestration/ISSUES.md              |  81 ++++++++++++++++++++++++++++++++-----------------------------------
 .orchestration/cleanup_progress.json  | 138 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 .orchestration/cleanup_reports/C03.md | 212 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 docs/audits/wilson_rg_inventory.md    |  57 +++++++++++++++++++++++++++++++++++++++++++++++
 quarkConstraints/deltaf2.py           |   7 ++++--
 scripts/audit_wilson_rg.py            |  18 +++++++++++++++
 7 files changed, 497 insertions(+), 45 deletions(-)
```

### 5. `fd26f26` — cleanup(C04): generalize k>0 Wilson-UL test coverage + add website snapshot suite (R07-I2, R22-I1; R02-I2 deferred to C13)
- SHA: `fd26f26da2492bb20f85ab88d5f6a0a00e794a1e`
- Message: cleanup(C04): generalize k>0 Wilson-UL test coverage + add website snapshot suite (R07-I2, R22-I1; R02-I2 deferred to C13)
- Physical/numerical summary: Generalized Wilson upper-limit tests for k>0 and added website notation/prose snapshot coverage; no physics-code shift beyond test coverage.
- Files touched (`git show --stat`):
```text
 .orchestration/CLEANUP_QUEUE.md                                              |   2 +-
 .orchestration/ISSUES.md                                                     |  22 ++++-----
 .orchestration/cleanup_progress.json                                         |   2 +-
 .orchestration/cleanup_reports/C04.md                                        | 144 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/website/package-lock.json                                     | 348 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
 flavor_catalog/website/package.json                                          |   7 ++-
 flavor_catalog/website/src/lib/__tests__/__snapshots__/notation.test.ts.snap |  41 ++++++++++++++++
 flavor_catalog/website/src/lib/__tests__/__snapshots__/prose.test.ts.snap    |  19 ++++++++
 flavor_catalog/website/src/lib/__tests__/notation.test.ts                    |  59 +++++++++++++++++++++++
 flavor_catalog/website/src/lib/__tests__/prose.test.ts                       |  66 ++++++++++++++++++++++++++
 tests/test_finite_stats.py                                                   |  81 ++++++++++++++++++++++++++++++++
 tests/test_quark_fit.py                                                      |  29 +++++++++++-
 12 files changed, 802 insertions(+), 18 deletions(-)
```

### 6. `12c4bd8` — cleanup(C05): add MANIFEST.md + sha256sums.txt to flavor_catalog/external_research/ (R17-I1)
- SHA: `12c4bd8153485eff07e526b9136c6467dd5cc6b2`
- Message: cleanup(C05): add MANIFEST.md + sha256sums.txt to flavor_catalog/external_research/ (R17-I1)
- Physical/numerical summary: Added external-research manifest and six-file SHA256 ledger, fixing provenance without changing catalog physics.
- Files touched (`git show --stat`):
```text
 .orchestration/CLEANUP_QUEUE.md                 |   2 +-
 .orchestration/ISSUES.md                        |  12 ++++++------
 .orchestration/cleanup_progress.json            |   2 +-
 .orchestration/cleanup_reports/C05.md           |  61 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/external_research/MANIFEST.md    | 105 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/external_research/sha256sums.txt |   6 ++++++
 6 files changed, 180 insertions(+), 8 deletions(-)
```

### 7. `df61de2` — cleanup(C06): embed git_sha in tile_summary.json writes (R07-I3)
- SHA: `df61de28d79dd45682545424ec50ebebcb876838`
- Message: cleanup(C06): embed git_sha in tile_summary.json writes (R07-I3)
- Physical/numerical summary: Embedded the producing git SHA into tile_summary.json writes, adding scan provenance while preserving legacy readers.
- Files touched (`git show --stat`):
```text
 .orchestration/CLEANUP_QUEUE.md       |   2 +-
 .orchestration/ISSUES.md              |  12 ++++++------
 .orchestration/cleanup_progress.json  |   2 +-
 .orchestration/cleanup_reports/C06.md | 112 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 scripts/run_rs_anarchy.py             |  29 +++++++++++++++++++++++++++++
 tests/test_run_rs_anarchy.py          |  98 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 6 files changed, 247 insertions(+), 8 deletions(-)
```

### 8. `7eaa94a` — cleanup(C07): annotate master compile reports with consolidation reconciliation (R17-I2, R18-I3, R19-I4)
- SHA: `7eaa94a32658127892ead45464928cb8e0671c7e`
- Message: cleanup(C07): annotate master compile reports with consolidation reconciliation (R17-I2, R18-I3, R19-I4)
- Physical/numerical summary: Annotated v0.2/v0.3/v0.4 master-compile reports and added an aggregator that reconciles 102 total processes to 101 VERIFIED plus 1 PARTIAL.
- Files touched (`git show --stat`):
```text
 .orchestration/CLEANUP_QUEUE.md             |   2 +-
 .orchestration/ISSUES.md                    |  32 ++++++++++++--------------
 .orchestration/cleanup_progress.json        |   2 +-
 .orchestration/cleanup_reports/C07.md       | 177 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/master_compile_v02_report.md |  15 ++++++++++++
 flavor_catalog/master_compile_v03_report.md |  17 ++++++++++++++
 flavor_catalog/master_compile_v04_report.md |  18 +++++++++++++++
 tools/aggregate_factchecks.py               | 210 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 8 files changed, 453 insertions(+), 20 deletions(-)
```

### 9. `6e1fd9d` — cleanup(C08): close open_issues threads across wave-1 yamls (R10a-I2, R10a-I3, R10b-I3)
- SHA: `6e1fd9d736a6faa043281d95a7b9b5ddb8ad4cc2`
- Message: cleanup(C08): close open_issues threads across wave-1 yamls (R10a-I2, R10a-I3, R10b-I3)
- Physical/numerical summary: Closed Wave-1 open_issues strings in catalog YAMLs with explicit disposition markers, preserving the original notes.
- Files touched (`git show --stat`):
```text
 .orchestration/CLEANUP_QUEUE.md           |   2 +-
 .orchestration/ISSUES.md                  |  32 ++++++++++++++------------------
 .orchestration/cleanup_progress.json      |   2 +-
 .orchestration/cleanup_reports/C08.md     | 120 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/beauty/B002.yaml |   2 +-
 flavor_catalog/processes/beauty/B005.yaml |   2 +-
 flavor_catalog/processes/beauty/B009.yaml |   4 ++--
 flavor_catalog/processes/beauty/B011.yaml |   2 +-
 flavor_catalog/processes/beauty/B015.yaml |   4 ++--
 flavor_catalog/processes/kaon/K001.yaml   |   4 ++--
 flavor_catalog/processes/kaon/K003.yaml   |   2 +-
 11 files changed, 146 insertions(+), 30 deletions(-)
```

### 10. `c3ec1b7` — cleanup(C09): YAML typo + status_history + Belle II citation fixes (R10b-I1, R10b-I2, R11-I3)
- SHA: `c3ec1b7ec40a5f5eda2bcf633fc47150139073f0`
- Message: cleanup(C09): YAML typo + status_history + Belle II citation fixes (R10b-I1, R10b-I2, R11-I3)
- Physical/numerical summary: Fixed catalog metadata typos/status-history gaps, including B011 observable wording, B009 state, and B025 citation-policy documentation.
- Files touched (`git show --stat`):
```text
 .orchestration/CLEANUP_QUEUE.md                         |  2 +-
 .orchestration/ISSUES.md                                | 32 ++++++++++++++------------------
 .orchestration/cleanup_progress.json                    |  2 +-
 .orchestration/cleanup_reports/C09.md                   | 64 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/beauty/B009.yaml               |  3 ++-
 flavor_catalog/processes/beauty/B011.yaml               |  2 +-
 flavor_catalog/processes/beauty/B025.yaml               |  2 +-
 flavor_catalog/website/_data/citation_anchors/B025.yaml | 43 +++++++++++++++++++++++++++++++++++++++++++
 8 files changed, 127 insertions(+), 23 deletions(-)
```

### 11. `1e1428e` — cleanup(C10): T010/T011 cross-reference + sha256 backfill + worklog-alias notes (R10c-I1, R10c-I2, R10c-I3)
- SHA: `1e1428e14936d40ea40d7bb8ef9b181c86b93566`
- Message: cleanup(C10): T010/T011 cross-reference + sha256 backfill + worklog-alias notes (R10c-I1, R10c-I2, R10c-I3)
- Physical/numerical summary: Backfilled T010/T011 merge provenance, reference SHA256 ledgers, and wave/worklog aliases; no physics fields changed.
- Files touched (`git show --stat`):
```text
 .orchestration/CLEANUP_QUEUE.md                   |   2 +-
 .orchestration/ISSUES.md                          |  32 ++++++++++++++------------------
 .orchestration/cleanup_progress.json              |   2 +-
 .orchestration/cleanup_reports/C10.md             | 115 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/charged_lepton/L001.yaml |  25 +++++++++++++++++++++++++
 flavor_catalog/processes/charm/C001.yaml          |  19 +++++++++++++++++++
 flavor_catalog/processes/edm_neutrino/E001.yaml   |  19 +++++++++++++++++++
 flavor_catalog/processes/top_higgs_ew/T002.yaml   |  21 +++++++++++++++++++++
 flavor_catalog/processes/top_higgs_ew/T011.yaml   |  73 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/references/E001/sha256sums.txt     |   6 ++++++
 flavor_catalog/references/T002/sha256sums.txt     |   7 +++++++
 flavor_catalog/references/T010/sha256sums.txt     |   7 +++++++
 12 files changed, 308 insertions(+), 20 deletions(-)
```

### 12. `357828d` — cleanup(C11): DA-N worklog closure addenda (R11-I2, R12-I4, R15-I4)
- SHA: `357828dcd47299787625809f1779665b7e31ce25`
- Message: cleanup(C11): DA-N worklog closure addenda (R11-I2, R12-I4, R15-I4)
- Physical/numerical summary: Added closure addenda to DA discovery worklogs, resolving deferred process-discovery bookkeeping without changing catalog entries.
- Files touched (`git show --stat`):
```text
 .orchestration/CLEANUP_QUEUE.md                            |  2 +-
 .orchestration/ISSUES.md                                   | 32 ++++++++++++++------------------
 .orchestration/cleanup_progress.json                       |  2 +-
 .orchestration/cleanup_reports/C11.md                      | 92 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/worklogs/discovery/round_001_full_scope.md  | 22 ++++++++++++++++++++++
 flavor_catalog/worklogs/discovery/round_002_followup.md    | 18 ++++++++++++++++++
 flavor_catalog/worklogs/discovery/round_003_final_sweep.md | 19 +++++++++++++++++++
 flavor_catalog/worklogs/discovery/round_004_convergence.md | 26 ++++++++++++++++++++++++++
 8 files changed, 193 insertions(+), 20 deletions(-)
```

### 13. `c53ed43` — cleanup(C12): reconcile v0.X tag-annotation vs compile-report count drift (R18-I1, R19-I3, R21-I1, R22-I3)
- SHA: `c53ed43ca74a93e5dca3351b015854310e50b320`
- Message: cleanup(C12): reconcile v0.X tag-annotation vs compile-report count drift (R18-I1, R19-I3, R21-I1, R22-I3)
- Physical/numerical summary: Reconciled tag-annotation and compile-report count drift; canonical v0.4 remains 102 processes with 101 VERIFIED and 1 PARTIAL.
- Files touched (`git show --stat`):
```text
 .orchestration/CLEANUP_QUEUE.md                        |   2 +-
 .orchestration/ISSUES.md                               |  42 ++++++++++++++++++------------------------
 .orchestration/cleanup_progress.json                   |   2 +-
 .orchestration/cleanup_reports/C12.md                  | 138 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/audits/factcheck_status.md              |  46 ++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/website/WEBSITE_RUNBOOK.md              |   2 +-
 flavor_catalog/website/src/components/EntryTable.astro |   1 -
 7 files changed, 205 insertions(+), 28 deletions(-)
```

### 14. `0a545f4` — cleanup(C13): add spurion-seed provenance comment (R02-I1, R02-I2)
- SHA: `0a545f42197bd2587170e482de7c17c3ad2f25c1`
- Message: cleanup(C13): add spurion-seed provenance comment (R02-I1, R02-I2)
- Physical/numerical summary: Added in-source provenance comments for the PDG-derived spurion seed while leaving the pinned singular values unchanged.
- Files touched (`git show --stat`):
```text
 .orchestration/CLEANUP_QUEUE.md       |  2 +-
 .orchestration/ISSUES.md              | 12 ++++++------
 .orchestration/cleanup_progress.json  |  2 +-
 .orchestration/cleanup_reports/C13.md | 84 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 quarkConstraints/benchmarks.py        | 12 +++++++++++-
 5 files changed, 103 insertions(+), 9 deletions(-)
```

### 15. `9e9a9fe` — cleanup(C14): fix codex model version drift gpt-5.5 -> gpt-5.4 (R14-I3, R20-I2)
- SHA: `9e9a9fe0e5aabab1caa5ee62cc6ee9dadd08c2e8`
- Message: cleanup(C14): fix codex model version drift gpt-5.5 -> gpt-5.4 (R14-I3, R20-I2)
- Physical/numerical summary: Corrected prescriptive model-version text from gpt-5.5 to gpt-5.4 while leaving historical provenance stamps intact.
- Files touched (`git show --stat`):
```text
 .orchestration/CLEANUP_QUEUE.md                     |  2 +-
 .orchestration/ISSUES.md                            | 22 ++++++++++------------
 .orchestration/cleanup_progress.json                |  2 +-
 .orchestration/cleanup_reports/C14.md               | 77 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 docs/paper_execution_decisions.md                   |  2 +-
 docs/phase_logs/POST_COMPACTION_BRIEFING.md         |  6 +++---
 docs/phase_logs/flavor_catalog_codex_quota_pause.md |  2 +-
 flavor_catalog/AGENTIC_WORKFLOW.md                  | 10 +++++-----
 flavor_catalog/HANDOFF_PROMPT.md                    |  2 +-
 flavor_catalog/SESSION_NOTES.md                     |  4 ++--
 flavor_catalog/WEBSITE_BUILD_PROMPT.md              |  4 ++--
 flavor_catalog/audits/factcheck_status.md           |  2 +-
 flavor_catalog/website/README.md                    |  4 ++--
 13 files changed, 107 insertions(+), 32 deletions(-)
```

### 16. `77d587e` — cleanup(C15): backfill .orchestration/pytest_selection/<unit>.txt for 24 review units (R05-I2)
- SHA: `77d587eccf7e06296096c884891dd30119af5f4e`
- Message: cleanup(C15): backfill .orchestration/pytest_selection/<unit>.txt for 24 review units (R05-I2)
- Physical/numerical summary: Generated pytest-selection files for all 24 Phase-1 review units, distinguishing six test-touching units from eighteen no-test units.
- Files touched (`git show --stat`):
```text
 .orchestration/CLEANUP_QUEUE.md          |  2 +-
 .orchestration/ISSUES.md                 | 12 ++++++------
 .orchestration/cleanup_progress.json     |  2 +-
 .orchestration/cleanup_reports/C15.md    | 77 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 .orchestration/pytest_selection/R01.txt  |  9 +++++++++
 .orchestration/pytest_selection/R02.txt  |  3 +++
 .orchestration/pytest_selection/R03.txt  |  1 +
 .orchestration/pytest_selection/R04.txt  |  6 ++++++
 .orchestration/pytest_selection/R05.txt  |  1 +
 .orchestration/pytest_selection/R06.txt  |  1 +
 .orchestration/pytest_selection/R07.txt  |  1 +
 .orchestration/pytest_selection/R08.txt  |  1 +
 .orchestration/pytest_selection/R09.txt  |  1 +
 .orchestration/pytest_selection/R10a.txt |  1 +
 .orchestration/pytest_selection/R10b.txt |  1 +
 .orchestration/pytest_selection/R10c.txt |  1 +
 .orchestration/pytest_selection/R11.txt  |  1 +
 .orchestration/pytest_selection/R12.txt  |  1 +
 .orchestration/pytest_selection/R13.txt  |  1 +
 .orchestration/pytest_selection/R14.txt  |  1 +
 .orchestration/pytest_selection/R15.txt  |  1 +
 .orchestration/pytest_selection/R16.txt  |  1 +
 .orchestration/pytest_selection/R17.txt  |  1 +
 .orchestration/pytest_selection/R18.txt  |  1 +
 .orchestration/pytest_selection/R19.txt  |  1 +
 .orchestration/pytest_selection/R20.txt  |  1 +
 .orchestration/pytest_selection/R21.txt  |  1 +
 .orchestration/pytest_selection/R22.txt  |  1 +
 28 files changed, 124 insertions(+), 8 deletions(-)
```

### 17. `85d9ee2` — cleanup(C16): update SESSION_NOTES tally to post-v0.4 canonical numbers (R20-I1, R20-I3)
- SHA: `85d9ee29f783c3c976bff60b25427a2363dfc968`
- Message: cleanup(C16): update SESSION_NOTES tally to post-v0.4 canonical numbers (R20-I1, R20-I3)
- Physical/numerical summary: Updated SESSION_NOTES running tallies to the post-v0.4 canonical 102-process, 101 VERIFIED, 1 PARTIAL state.
- Files touched (`git show --stat`):
```text
 .orchestration/CLEANUP_QUEUE.md       |  2 +-
 .orchestration/ISSUES.md              | 28 ++++++++++++++++------------
 .orchestration/cleanup_progress.json  |  2 +-
 .orchestration/cleanup_reports/C16.md | 80 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/SESSION_NOTES.md       | 30 ++++++++++++++++++------------
 5 files changed, 116 insertions(+), 26 deletions(-)
```

### 18. `5c1a8f7` — cleanup(C17): remove redundant .gitkeep placeholders (R09-I2)
- SHA: `5c1a8f788591d63db05085f61a19d959abc80919`
- Message: cleanup(C17): remove redundant .gitkeep placeholders (R09-I2)
- Physical/numerical summary: Removed redundant .gitkeep placeholders; directory-tracking only, no numerical or physics effect.
- Files touched (`git show --stat`):
```text
 .orchestration/CLEANUP_QUEUE.md                  |  2 +-
 .orchestration/ISSUES.md                         | 12 ++++++------
 .orchestration/cleanup_progress.json             |  2 +-
 .orchestration/cleanup_reports/C17.md            | 47 +++++++++++++++++++++++++++++++++++++++++++++++
 flavor_catalog/processes/beauty/.gitkeep         |  0
 flavor_catalog/processes/charged_lepton/.gitkeep |  0
 flavor_catalog/processes/charm/.gitkeep          |  0
 flavor_catalog/processes/edm_neutrino/.gitkeep   |  0
 flavor_catalog/processes/kaon/.gitkeep           |  0
 flavor_catalog/processes/top_higgs_ew/.gitkeep   |  0
 flavor_catalog/references/.gitkeep               |  0
 flavor_catalog/signoff/by_process/.gitkeep       |  0
 flavor_catalog/worklogs/checker/.gitkeep         |  0
 flavor_catalog/worklogs/discovery/.gitkeep       |  0
 flavor_catalog/worklogs/pka/.gitkeep             |  0
 flavor_catalog/worklogs/writer/.gitkeep          |  0
 16 files changed, 55 insertions(+), 8 deletions(-)
```

### 19. `914814c` — cleanup(C18): MERGE_PLAN.md retroactive corrections (R12-I1..I3, R13-I1, R14-I1, R15-I1, R16-I3, R19-I1, R19-I2, R19-I5)
- SHA: `914814c7ee912f67a3dc0749486cf25fb167afb1`
- Message: cleanup(C18): MERGE_PLAN.md retroactive corrections (R12-I1..I3, R13-I1, R14-I1, R15-I1, R16-I3, R19-I1, R19-I2, R19-I5)
- Physical/numerical summary: Applied retroactive MERGE_PLAN SHA/schema/count corrections for review bookkeeping; no code or catalog numerics changed.
- Files touched (`git show --stat`):
```text
 .orchestration/CLEANUP_QUEUE.md       |   2 +-
 .orchestration/ISSUES.md              | 102 ++++++++++++++++++++++++++++++++++++++++++------------------------------------------------------------
 .orchestration/MERGE_PLAN.md          |  18 ++++++++++++------
 .orchestration/cleanup_progress.json  |   2 +-
 .orchestration/cleanup_reports/C18.md |  99 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 5 files changed, 155 insertions(+), 68 deletions(-)
```

### 20. `5a22b9f` — cleanup(C19): final docs sweep + methodology PDF rebuild (R08-I3, R04-I4)
- SHA: `5a22b9ff3d2844ae23c3cb5ad4b9ec9001bd8bec`
- Message: cleanup(C19): final docs sweep + methodology PDF rebuild (R08-I3, R04-I4)
- Physical/numerical summary: Closed final docs/methodology issues and rebuilt the methodology PDF with Wilson-RG cross-reference updates.
- Files touched (`git show --stat`):
```text
 .orchestration/CLEANUP_QUEUE.md       |   2 +-
 .orchestration/ISSUES.md              |  12 ++++++------
 .orchestration/cleanup_progress.json  |   2 +-
 .orchestration/cleanup_reports/C19.md |  37 +++++++++++++++++++++++++++++++++++++
 docs/audits/figure_prune_inventory.md |  29 +++++++++++++++++++++++++++++
 docs/quark_scan_methodology_note.pdf  | Bin 626853 -> 627505 bytes
 docs/quark_scan_methodology_note.tex  |   7 ++++++-
 7 files changed, 80 insertions(+), 9 deletions(-)
```

### 21. `759289a` — cleanup(C20): final ACCEPTED-RISK closures + verify pre-closed INFO items (R01-I2 + others)
- SHA: `759289afbdbcc9607ed1d0df632f1d1154ecc749`
- Message: cleanup(C20): final ACCEPTED-RISK closures + verify pre-closed INFO items (R01-I2 + others)
- Physical/numerical summary: Moved remaining INFO/LOW issues to accepted-risk or auto-resolved closure buckets, leaving only deferred cleanup pointers open.
- Files touched (`git show --stat`):
```text
 .orchestration/CLEANUP_QUEUE.md       |   2 +-
 .orchestration/ISSUES.md              | 343 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++--------------------------------------------------------------------------------------------
 .orchestration/cleanup_progress.json  |   2 +-
 .orchestration/cleanup_reports/C20.md | 126 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 4 files changed, 290 insertions(+), 183 deletions(-)
```

### 22. `6d305d8` — chore(orchestration): commit Phase 2 cleanup plan + 21 per-unit reports
- SHA: `6d305d8b1ef9b064ba74948674e3c91c2e4a6ee2`
- Message: chore(orchestration): commit Phase 2 cleanup plan + 21 per-unit reports
- Physical/numerical summary: Committed the cleanup plan, reviews, and per-unit reports as the Phase-2 ledger; bookkeeping only.
- Files touched (`git show --stat`):
```text
 .orchestration/CLEANUP_CHANGES_R2.md               |  189 +++++++++++++++++++++++++++
 .orchestration/CLEANUP_PLAN.md                     | 1135 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 .orchestration/CLEANUP_REVIEW_R1.md                |   53 ++++++++
 .orchestration/CLEANUP_REVIEW_R2.md                |   29 +++++
 .orchestration/cleanup_reports/C01.md              |  141 ++++++++++++++++++++
 .orchestration/cleanup_reports/C01_REVIEW.md       |  113 ++++++++++++++++
 .orchestration/cleanup_reports/C02a-code.md        |   99 ++++++++++++++
 .orchestration/cleanup_reports/C02a-code_REVIEW.md |  130 +++++++++++++++++++
 .orchestration/cleanup_reports/C03_REVIEW.md       |  107 +++++++++++++++
 .orchestration/cleanup_reports/C04_REVIEW.md       |   57 ++++++++
 .orchestration/cleanup_reports/C05_REVIEW.md       |   46 +++++++
 .orchestration/cleanup_reports/C06_REVIEW.md       |   77 +++++++++++
 .orchestration/cleanup_reports/C07_REVIEW.md       |   27 ++++
 .orchestration/cleanup_reports/C08_REVIEW.md       |   49 +++++++
 .orchestration/cleanup_reports/C09_REVIEW.md       |   41 ++++++
 .orchestration/cleanup_reports/C10_REVIEW.md       |   64 +++++++++
 .orchestration/cleanup_reports/C11_REVIEW.md       |   64 +++++++++
 .orchestration/cleanup_reports/C12_REVIEW.md       |   51 ++++++++
 .orchestration/cleanup_reports/C13_REVIEW.md       |   55 ++++++++
 .orchestration/cleanup_reports/C14_REVIEW.md       |   38 ++++++
 .orchestration/cleanup_reports/C15_REVIEW.md       |   52 ++++++++
 .orchestration/cleanup_reports/C16_REVIEW.md       |   34 +++++
 .orchestration/cleanup_reports/C18_REVIEW.md       |   52 ++++++++
 23 files changed, 2703 insertions(+)
```

### 23. `71a736b` — cleanup(C02c-scripts): create SLURM dispatcher + post-processing script for epsilon_K band reruns
- SHA: `71a736bd17e9bf9db7d2ff25aaa746f5f6644f1f`
- Message: cleanup(C02c-scripts): create SLURM dispatcher + post-processing script for epsilon_K band reruns
- Physical/numerical summary: Added the C02c SLURM dispatcher and post-processing script for epsilon_K budget-edge reruns; prepares measured band extraction but does not itself run scans.
- Files touched (`git show --stat`):
```text
 scripts/c02c_post_process.py       | 230 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 scripts/run_rs_anarchy_c02c.sbatch |  82 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 2 files changed, 312 insertions(+)
```

