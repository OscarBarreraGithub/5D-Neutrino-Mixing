# `flavor_catalog/external_research/` — MANIFEST

Provenance log for external research artifacts and the in-house review memos
that compare them with the catalog. Closes R17-I1 (cleanup unit C05).

The directory holds two kinds of files:

1. **Imported artifacts** — external GPT Deep Research outputs (PDF + plain-text
   extract) on RS flavor constraints. These are read-only third-party reports;
   the catalog never depends on them for numerical claims (every Wave-7+
   promotion was independently re-verified against PDG / HFLAV / experiment
   sources in the per-family `factcheck_*.md` addenda).
2. **Review memos** — in-house markdown comparisons of each imported artifact
   against the catalog at the relevant snapshot (v0.1 for both memos). These
   feed Wave-7 deferred-scope reconciliation (`worklogs/discovery/
   round_004_addendum_deferred_scope.md`) and the R15/R16 review trail.

Cryptographic digests for all six files are recorded in `sha256sums.txt` and
can be re-validated with:

```
cd flavor_catalog/external_research && sha256sum -c sha256sums.txt
```

`MANIFEST.md` and `sha256sums.txt` are themselves intentionally excluded from
`sha256sums.txt` so the manifest can be edited without invalidating the
digest list.

---

## Imported artifacts (GPT Deep Research)

### `deepresearch_may15.pdf` / `deepresearch_may15.txt`

| field | value |
|-------|-------|
| Kind | imported artifact (PDF + plain-text extract) |
| Source | GPT Deep Research (external session, output exported to PDF) |
| Subject | RS flavor constraints — process scope review against catalog v0.1 |
| Date imported | 2026-05-16 (commit `022a20c`) |
| PDF size | 203806 bytes |
| Text size | 1109 lines |
| sha256 (PDF) | `008be1940ca03866715a2bd932263a9bda7503b62bd38fd0a410bd4ee949ec3a` |
| sha256 (TXT) | `d2f4dc68124e2e0b2fbf98391bc1464dd930ecd2fabb7735d19fa6df0cd04c11` |
| Description | Identifies catalog scope gaps in radiative-B chirality (`B -> K* gamma`, `B_s -> phi gamma`, photon polarization, `C7/C7'`), top FCNC photon modes (`t -> q gamma`), and related operator-diagnosis observables. Drives the Wave-7 promotion of T003/T004/T008/T012/B012 — see `signoff/round_003_index.md:10`. |
| Extraction provenance | `.txt` adjacent to the `.pdf`; both committed in `022a20c` with no explicit `pdftotext` invocation recorded in the commit message. Treat the `.txt` as a best-effort text mirror, not a canonical extract. |
| Session URL | not recorded at import time (R17-I1 documents the gap; no retroactive recovery) |

### `deepresearch_may16.pdf` / `deepresearch_may16.txt`

| field | value |
|-------|-------|
| Kind | imported artifact (PDF + plain-text extract) |
| Source | GPT Deep Research (external session, output exported to PDF) |
| Subject | RS flavor constraints — expanded process scope review against catalog v0.1 (followup to may15) |
| Date imported | 2026-05-16 (commit `022a20c`) |
| PDF size | 343307 bytes |
| Text size | 2092 lines |
| sha256 (PDF) | `d35d20c9eb43227acfb4c7d129d1b5d34ec28232cffa0f3bcb70da82e95e0ede` |
| sha256 (TXT) | `e99b950b46fbdde5fa6192f0a216342910eb5a60867d219c04278c4440d09781` |
| Description | Extends the may15 scope to LFV kaon tails (`K_L -> mu e`, etc.), LFV B and semileptonic LFV meson modes (`B -> K mu e`, `B^0 -> mu tau`), semileptonic asymmetries / `|q/p|_{d,s}`, rare-charm four-body angular modes, tau/W universality ratios, and pulls current PDG/HFLAV values for cross-check (some predate the catalog's HFLAV CKM 2025 sourcing). |
| Extraction provenance | `.txt` adjacent to the `.pdf`; both committed in `022a20c` with no explicit `pdftotext` invocation recorded. |
| Session URL | not recorded at import time |

---

## Review memos (in-house)

### `deepresearch_may15_review.md`

| field | value |
|-------|-------|
| Kind | review memo (markdown) |
| Author | in-house (Claude Code session) |
| Date | 2026-05-16 (commit `2c00d84`) |
| sha256 | `e27a190941638eb8a60b0bf96e11f2664dfc23bce3007a1ebfe719049e8c23ae` |
| Description | Compares `deepresearch_may15.{pdf,txt}` line-by-line against catalog v0.1 (75-row factcheck inventory + plan-v1 128 seed + DA-4 convergence). Sections: (1) processes missing from the catalog, (2) numeric disagreements, (3) judgment on deferral soundness. Cited by `signoff/round_003_index.md:10` as a Wave-7 input. |

### `deepresearch_may16_review.md`

| field | value |
|-------|-------|
| Kind | review memo (markdown) |
| Author | in-house (Claude Code session) |
| Date | 2026-05-16 (commit `8fb5f91`) |
| sha256 | `696b5ea97eda320dc2f03f8a0b45e1842e54e03a466f3e33bb75aecb18b2be07` |
| Description | Same structure as the may15 review but applied to `deepresearch_may16.{pdf,txt}`. Identifies eight missing-scope items (radiative-helicity B, top photon FCNCs, `t -> uH`, LFV kaon/B tails, semileptonic mixing asymmetries, rare-charm four-body, tau-universality ratios) and several numeric/source-vintage disagreements (`C004`, `B025/B026`, `B009`). Feeds the Wave-7 deferred-scope addendum. |

---

## Notes on remaining provenance gaps

R17-I1 noted three sub-items: (a) PDF binary digests, (b) GPT Deep Research
session URLs / dates / prompts, (c) text-extraction provenance.

- **(a) Closed** by `sha256sums.txt` plus the per-row digests above.
- **(b) Partially closed.** The import date is recorded from git history; the
  GPT session URLs and originating prompts were not captured at import time
  and are not retroactively recoverable. Future external-research imports
  should attach the session URL + prompt to the import commit message or to
  this MANIFEST in a `prompts/` sub-directory.
- **(c) Partially closed.** The `.txt` files were committed alongside the
  PDFs in `022a20c` without an explicit extraction recipe. They should be
  treated as best-effort text mirrors. If a future workflow re-extracts text
  from the PDFs, record the tool + flags (e.g. `pdftotext -layout`) here.
