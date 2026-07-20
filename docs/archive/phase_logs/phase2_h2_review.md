### Verdict
REJECT-WITH-REVISIONS

### Item-by-item (1-9)
1. PASS - `paper/quark-scan-2026q2` is at `b86e66d`, matches `origin/paper/quark-scan-2026q2`, `docs/canonical-manifest` points to the same SHA, and `760f23d..paper/quark-scan-2026q2` contains 3 commits.
2. PASS - `docs/artifact_manifest.md` lists all 9 canonical scan dirs with path, `tile_summary` SHA256, `draws` SHA256, size, and run command; it also has the paper-claim mapping table, audit dependency chain, and Python/key-package environment block.
3. PASS - full `sha256sum -c artifacts/checksums.sha256` completed with `OK` for all canonical scan files and the historical 800k cross-check entries.
4. PASS - JSON parse command completed without error and printed `15`.
5. PASS - `tar -tzf artifacts/quarkscan_paper_rc1_manifest_bundle.tar.gz | head` listed `docs/artifact_manifest.md`, manifest JSON, checksums, and README.
6. FAIL - `git ls-files artifacts/` shows only `README.md`, `checksums.sha256`, and `quarkscan_paper_rc1_manifest.json`; the manifest bundle tarball and `.sha256` sidecar exist locally but are ignored by `.gitignore` lines 66-67.
7. FAIL - the required grep `grep -nE "artifact_manifest|canonical artifact" docs/quark_scan_methodology_note.tex` exits 1; the note has `docs/artifact\_manifest.md` at lines 980-985, but no literal matching anchor.
8. FAIL - `git show --stat 80641cf cb27631 b86e66d` also touches `docs/phase_logs/phase2_h2_impl.md`, outside the stated allowlist. No source code drift was observed, but the literal check is not satisfied.
9. FAIL - mapping spot-check passes for `47.26`, `127.13`, `2.3e-6`, `45.73/120.39`, and RunB crossings, but the methodology note quotes `23.4` TeV and CFW `10.5` TeV while the manifest row lists `23.37 TeV` and omits `10.5`.

### Findings
1. High - Artifact tracking is incomplete for the requested small manifest bundle. Fix by tracking `artifacts/quarkscan_paper_rc1_manifest_bundle.tar.gz` and its `.sha256` sidecar, or revise the checklist if the intended policy is to keep even the small manifest bundle ignored. The heavy scan-output tarballs should remain untracked.
2. Medium - The methodology-note audit hook does not pass the exact command. Fix by adding a literal raw-source phrase such as `canonical artifact manifest` or `artifact_manifest` near the existing manifest pointer.
3. Medium - The paper-claim mapping table is not complete for the CFW comparison quote as written in the note. Fix by adding the rounded `23.4` TeV claim and the external CFW comparator value `10.5` TeV to the mapping row, with source scan dir and audit/plot pointer.
4. Low - The no-drift allowlist is exceeded by the implementation log. Fix by either excluding that file from the reviewed commit range or explicitly allowing phase-log implementation reports in the checklist.

### Final
Send back for the tracked bundle/sidecar, exact methodology-note anchor, and CFW mapping fixes before Opus sign-off.
===PHASE_2_H2_REVIEW_END===
