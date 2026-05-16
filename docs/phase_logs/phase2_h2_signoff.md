# Phase 2 Hole #2 Sign-off

**Verdict**: PASS

This sign-off ratifies the APPROVE verdict in
`docs/phase_logs/phase2_h2_review_v2.md` after independent re-verification
of the six convergence checks below. The canonical run manifest for the
quark-sector paper RC1 (`quarkscan_paper_rc1`) is now sealed: tracked in
git, checksum-locked, bundled, cross-referenced from the methodology
note, and equipped with explicit CFW mapping rows.

## Verification

### 1. Clean commit chain since `ad2e3ee`

`git log --oneline ad2e3ee..HEAD` returns 16 commits, of which the
hole-#2 manifest work occupies a contiguous block of six commits, and
the hole-#9 closeout occupies four commits at the tip. The remaining
six commits cover hole-#8 (scope appendix and zero-pass wording) and
its prerequisite physics/stats helper. Hole-#2 commits, in order
applied:

```
80641cf chore(gitignore): track artifacts/ manifest and checksums
cb27631 docs(artifacts): seal canonical artifact manifest for paper
b86e66d docs(paper): cite artifact manifest in methodology note
2729659 chore(gitignore): track manifest bundle tarball + sidecar
b946abf docs(artifacts): add CFW mapping rows to canonical manifest
5e60bc0 docs(paper): add canonical-artifact-manifest anchor phrase in methodology note
```

Each commit has a single focused purpose, the bundle/sidecar tracking
follows the artifact creation rather than preceding it, and the
methodology-note edits sit on top of a sealed manifest. The chain is
clean: PASS.

### 2. Tracked artifact inventory: exactly 5 files

`git ls-files artifacts/` returns:

```
artifacts/README.md
artifacts/checksums.sha256
artifacts/quarkscan_paper_rc1_manifest.json
artifacts/quarkscan_paper_rc1_manifest_bundle.tar.gz
artifacts/quarkscan_paper_rc1_manifest_bundle.tar.gz.sha256
```

Both the release tarball and its `.sha256` sidecar are version-controlled
alongside the human-readable README, the per-input `checksums.sha256`,
and the JSON manifest. PASS.

### 3. Methodology note cites the canonical artifact manifest

`grep -nE "canonical artifact manifest|artifact_manifest" docs/quark_scan_methodology_note.tex`
returns:

```
980:All numerical claims in this note map to the canonical artifact manifest at
```

The anchor phrase is present in the methodology note itself (not just in
an auxiliary log), so the paper-side reader is given an explicit pointer
to `docs/artifact_manifest.md`. PASS.

### 4. CFW mapping rows for 23.4 TeV and 10.5 TeV

`grep -nE "23\.4|10\.5" docs/artifact_manifest.md` returns rows at lines
42 and 43:

- Line 42 maps `docs/quark_scan_methodology_note.tex:892` (CFW-matched p50
  projection quoted as `23.4 TeV`, rounded from `23.37 TeV`) to scan
  directory `rs_anarchy_runA_20260515T085316` and figure
  `results/figures/quark/rs_anarchy_cfw_comparison.pdf`
  (`b78b5d3c...da6f95`).
- Line 43 maps `docs/quark_scan_methodology_note.tex:894` (CFW no-UV-
  boundary-term `g_s^*=3` reference quoted as `10.5 TeV`) to the external
  literature value CFW 0804.1954 Section 2, with no local scan directory.

The locally generated comparison and the external literature value are
clearly separated. PASS.

### 5. Checksum spot-check (first 3 lines)

`sha256sum -c artifacts/checksums.sha256 | head -3` returns:

```
scan_outputs/rs_anarchy_runA_20260515T085316/tile_summary.json: OK
scan_outputs/rs_anarchy_runA_20260515T085316/draws.jsonl: OK
scan_outputs/rs_anarchy_run3_baseline_20260515T085324/tile_summary.json: OK
```

All three spot-checked entries verify OK against on-disk content. PASS.

### 6. JSON manifest has structure (> 5 keys)

`python -c "import json; d=json.load(open('artifacts/quarkscan_paper_rc1_manifest.json')); print('ok' if len(d) > 5 else 'short')"`
returns `ok`. The JSON manifest is non-trivially structured rather than
a placeholder stub. PASS.

## Recommendation

Proceed to Phase 3 hole #10 (figure prune + doc hygiene).

===PHASE_2_H2_SIGNOFF_END===
