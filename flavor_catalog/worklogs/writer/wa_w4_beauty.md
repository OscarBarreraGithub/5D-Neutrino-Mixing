# WA Worklog: wa_w4_beauty

Date: 2026-05-16.  Agent: WA.  Family: beauty.

Batch processes: B004, B006, B019, B026.

Scope followed: edited only `flavor_catalog/processes/beauty/` sidecars/drafts
and this writer worklog.  No process-local reference snapshots, PKA worklogs,
catalog indexes, templates, macros, or other families were modified by this WA
batch.  I did not prepare `references/catalog.bib` changes because the batch
hard rules forbid edits outside the allowed paths.

## Common writer pass

- Read plan v1 Section B and the WA deliverables/success criteria in Section D,
  plus the orchestrator decisions file.
- Read each assigned `.tex`, `.yaml`, PKA worklog, source manifest, and checked
  local reference snapshots for the quoted numerical values.
- Normalized Section B headings in all four TeX drafts and added source-key
  anchors near numerical claims.
- Appended `WRITER-DONE` status entries and updated `last_updated_at` in all
  four YAML sidecars without changing the `pdg_or_equivalent`,
  `code_coverage`, or `implementation_difficulty` blocks.

## Per-process changes

### B004: `phi_s` in `B_s -> J/psi phi`

- Tightened the process and post-2008 prose, and normalized the headings to
  `PDG or equivalent value`, `Relevance to RS / anarchic flavor`, and
  `Constraint validity / model dependence`.
- Added explicit source-key anchors for the HFLAV PDG 2025 averages, the LHCb
  2024 Run-2 input, and the SM `-2 beta_s` comparison.
- Kept the HFLAV all-combined `phi_s^ccs` average as the canonical headline and
  retained the strict `J/psi K+K-` combination as secondary provenance.

### B006: `B^0 -> mu+ mu-`

- Normalized headings and tightened the PDG/HFLAV/theory value paragraph.
- Added source-key anchors for the PDG live/API limit, CMS 2023 input, HFLAV
  ratio, ATLAS and LHCb limits, and Bobeth et al. SM prediction.
- Kept the HFLAV `B^0/B_s^0` dimuon ratio in the prose as an equivalent
  provenance item because it is documented in the sidecar.

### B019: `R_K*`

- Normalized headings, fixed `\mathcal{B}` usage, and removed extra prose around
  the neutral `K^*` notation.
- Added source-key anchors for the Dec. 2025 HFLAV values, the 2017 anomaly
  milestone, the 2023 superseding LHCb LFU result, and the BIP SM/QED context.
- Rewrote the key-reference list as process-local source keys rather than a raw
  prose list.

### B026: `R_D*`

- Normalized headings and tightened the HFLAV CKM 2025 value paragraph.
- Added source-key anchors for the HFLAV average/logfile, BaBar post-2008
  milestones, FLAG context, and the 2008 RS-flavor baseline.
- Rewrote the key-reference list using process-local source keys instead of
  filenames with `.txt` suffixes.

## Open issues for CA

- B004: Confirm that the all-combined HFLAV `phi_s^ccs` average should remain
  the canonical headline, with the strict `B_s -> J/psi K+K-` average secondary.
- B006: Confirm whether the HFLAV `B^0/B_s^0` ratio should stay in the final
  prose or move to sidecar-only provenance.
- B019: Confirm the citation convention for direct HFLAV Dec. 2025 EOS pages
  versus the broader rare-decay landing page.  The sidecar does not reconstruct
  HFLAV correlations or likelihoods.
- B026: The current HFLAV CKM 2025 page quotes a `3.8 sigma` combined
  discrepancy; this draft follows that source and does not mention older
  approximate `3.3 sigma` wording.  The Belle II CKM-2025 hadronic-tag input
  remains sourced through HFLAV rather than a tracked Indico PDF.
- Status lineage: B004, B006, and B019 sidecars entered this WA pass with
  `WRITER-INITIATED` as the final status rather than `PKA-DONE`; I appended
  `WRITER-DONE` from the observed current state and did not rewrite prior PKA
  history.

No CHECK markers were left in the TeX drafts.
