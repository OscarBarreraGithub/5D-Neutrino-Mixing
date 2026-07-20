# Quarkscan Paper RC1 Artifact Files

> **LEGACY (May 2026 ΔF=2-only, pre-audit).** These bless the `quarkscan-paper-rc1`
> ΔF=2 scan artifacts and predate the June 2026 audit; they do not reflect the
> corrected minimal-RS floors. See `docs/FLOOR_SUMMARY.md` and
> `reports/collaborator_2026-06/CONTENT.md`.

This directory tracks the small, reviewable files that bless the quark-scan paper artifacts:

- `quarkscan_paper_rc1_manifest.json`: machine-readable manifest matching `docs/archive/superseded/artifact_manifest.md`
- `checksums.sha256`: SHA256 checksums for canonical scan outputs, plus the historical 800k cross-check

The scan outputs are not committed. They remain under the cluster paths listed in `docs/archive/superseded/artifact_manifest.md`.

The local cluster bundle is:

```text
/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/artifacts/quarkscan_paper_rc1_manifest_bundle.tar.gz
/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/artifacts/quarkscan_paper_rc1_manifest_bundle.tar.gz.sha256
```

The tarball is intentionally ignored by git; it packages only the tracked manifest, README, and checksum files, not the heavy `scan_outputs/` tree.
