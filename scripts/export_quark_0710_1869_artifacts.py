#!/usr/bin/env python3
"""Write the canonical PR6 paper artifacts for the default kaon benchmark."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def _parse_args(*, default_output_dir: Path) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Export deterministic paper_0710_1869 artifacts for the default kaon benchmark."
        )
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=default_output_dir,
        help=(
            "Destination directory for wilsons.json, hadronic.json, "
            "observables.json, provenance.json."
        ),
    )
    return parser.parse_args()


def main() -> int:
    from quarkConstraints.paper_0710_1869.artifacts import (
        DEFAULT_KAON_ARTIFACT_DIR,
        write_default_paper_0710_1869_kaon_artifact_exports,
    )

    args = _parse_args(default_output_dir=DEFAULT_KAON_ARTIFACT_DIR)
    paths = write_default_paper_0710_1869_kaon_artifact_exports(root_dir=args.output_dir)
    print(paths.wilson_path)
    print(paths.hadronic_path)
    print(paths.observable_path)
    print(paths.provenance_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
