#!/usr/bin/env python3
"""Aggregate VERIFIED / PARTIAL / MISMATCH / FAILED counts across the
flavor-catalog per-family fact-check files.

Reads:
  flavor_catalog/audits/factcheck_status.md          (consolidated table)
  flavor_catalog/audits/factcheck_<family>.md        (per-family files)

Scans each file for process-row markdown like:
    | <ID> | <family> | VERIFIED | <date> | <agent> | <notes> |
or family-addendum table rows like:
    | <ID> | VERIFIED | 0 | <notes> |

and produces a unified per-family count plus an overall sum. This is a
read-only tool: it does not edit any catalog file. It exists to make
the per-master-compile reconciliation arithmetic trivially reproducible
(see C07 cleanup, ISSUES.md R17-I2 / R18-I3 / R19-I4).

Usage:
    python tools/aggregate_factchecks.py
    python tools/aggregate_factchecks.py --json   # machine-readable

Exit status: 0 on success; 1 if no factcheck files are found.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
AUDITS = REPO_ROOT / "flavor_catalog" / "audits"

# Status tokens recognised in fact-check tables.
STATUSES = ("VERIFIED", "PARTIAL", "MISMATCH", "FAILED")

# Match a row whose first column is a catalog process ID like K001, B012,
# CR001, EW001, etc. We deliberately allow 1-3 leading letters followed by
# digits to cover all current and any plausible future family prefixes.
_ID_RE = re.compile(r"^\|\s*([A-Z]{1,3}[0-9]{2,4})\s*\|")
_STATUS_RE = re.compile(r"\|\s*(VERIFIED|PARTIAL|MISMATCH|FAILED)\s*\|")
# Wave-7 / Wave-8 addenda sometimes record a single process verdict in a
# section header rather than a table row, e.g. "## B012 - VERIFIED" or
# "### T003 - VERIFIED". Match those too.
_HEADER_RE = re.compile(
    r"^#{2,4}\s+([A-Z]{1,3}[0-9]{2,4})\s*-\s*(VERIFIED|PARTIAL|MISMATCH|FAILED)\b"
)

# Map ID prefix -> family bucket used by the catalog. Anything not listed
# here will be reported under family "unknown".
_PREFIX_TO_FAMILY = {
    "K": "kaon",
    "C": "charm",
    "B": "beauty",
    "T": "top_higgs_ew",
    "EW": "top_higgs_ew",
    "L": "charged_lepton",
    "E": "edm_neutrino",
    "CR": "collider_rs",
}


def _family_for_id(pid: str) -> str:
    """Map a process ID like 'EW001' or 'CR014' to its family bucket."""
    for prefix_len in (2, 1):
        prefix = pid[:prefix_len]
        if prefix in _PREFIX_TO_FAMILY:
            return _PREFIX_TO_FAMILY[prefix]
    return "unknown"


def _scan_file(path: Path) -> dict[str, tuple[str, str]]:
    """Return {process_id: (status, source_path)} for rows whose first
    column is a known process-ID shape and which contain a status token.

    Later occurrences of the same process_id in the same file overwrite
    earlier ones — this preserves the convention that the consolidated
    summary block at the bottom of a file supersedes ad-hoc earlier
    references. Across files, the caller decides precedence.
    """
    rows: dict[str, tuple[str, str]] = {}
    try:
        text = path.read_text(encoding="utf-8")
    except OSError:
        return rows
    for line in text.splitlines():
        m_header = _HEADER_RE.match(line)
        if m_header:
            pid, status = m_header.group(1), m_header.group(2)
            # A later table-row entry should still override the header; but
            # if no table row exists for this ID in this file, the header
            # verdict stands.
            rows.setdefault(pid, (status, str(path.relative_to(REPO_ROOT))))
            continue
        m_id = _ID_RE.match(line)
        if not m_id:
            continue
        m_status = _STATUS_RE.search(line)
        if not m_status:
            continue
        pid = m_id.group(1)
        status = m_status.group(1)
        rows[pid] = (status, str(path.relative_to(REPO_ROOT)))
    return rows


def aggregate() -> dict[str, object]:
    """Scan factcheck_status.md (lowest precedence) then per-family files
    (higher precedence — addendum rows live there). Return a structured
    report dict.
    """
    if not AUDITS.is_dir():
        return {"error": f"audits dir not found: {AUDITS}"}

    consolidated = AUDITS / "factcheck_status.md"
    family_files = sorted(
        p for p in AUDITS.glob("factcheck_*.md") if p.name != "factcheck_status.md"
    )

    # Precedence: the consolidated `factcheck_status.md` table is the
    # authoritative final verdict for any row it lists (it folds in
    # follow-up corrections like the T020 MISMATCH -> VERIFIED fix). The
    # per-family files contribute *additional* rows from Wave-7+ addenda
    # that may not have been folded back into the consolidated table yet,
    # but should not override a row already present in consolidated.
    merged: dict[str, tuple[str, str]] = {}
    for f in family_files:
        for pid, row in _scan_file(f).items():
            merged.setdefault(pid, row)
    if consolidated.exists():
        # Consolidated rows take precedence: overwrite any family-file row.
        merged.update(_scan_file(consolidated))

    if not merged:
        return {"error": "no fact-check rows found"}

    by_family: dict[str, dict[str, int]] = defaultdict(
        lambda: {s: 0 for s in STATUSES} | {"total": 0}
    )
    totals = {s: 0 for s in STATUSES} | {"total": 0}
    for pid, (status, _src) in merged.items():
        fam = _family_for_id(pid)
        by_family[fam][status] += 1
        by_family[fam]["total"] += 1
        totals[status] += 1
        totals["total"] += 1

    return {
        "totals": totals,
        "by_family": dict(by_family),
        "process_count": len(merged),
        "consolidated_path": str(consolidated.relative_to(REPO_ROOT)),
        "family_files": [str(p.relative_to(REPO_ROOT)) for p in family_files],
    }


def _format_text(report: dict[str, object]) -> str:
    if "error" in report:
        return f"ERROR: {report['error']}\n"
    lines: list[str] = []
    totals = report["totals"]  # type: ignore[index]
    lines.append(
        f"Aggregated fact-check totals across {report['process_count']} processes:"
    )
    lines.append("")
    lines.append(f"  VERIFIED : {totals['VERIFIED']}")
    lines.append(f"  PARTIAL  : {totals['PARTIAL']}")
    lines.append(f"  MISMATCH : {totals['MISMATCH']}")
    lines.append(f"  FAILED   : {totals['FAILED']}")
    lines.append(f"  TOTAL    : {totals['total']}")
    lines.append("")
    lines.append("Per family:")
    by_family = report["by_family"]  # type: ignore[index]
    for fam in sorted(by_family.keys()):
        c = by_family[fam]
        lines.append(
            f"  {fam:<16s} V={c['VERIFIED']:>3d}  P={c['PARTIAL']:>2d}  "
            f"M={c['MISMATCH']:>2d}  F={c['FAILED']:>2d}  total={c['total']:>3d}"
        )
    lines.append("")
    lines.append(f"Consolidated source: {report['consolidated_path']}")
    lines.append("Per-family sources scanned (later wins on ID collision):")
    for f in report["family_files"]:  # type: ignore[index]
        lines.append(f"  - {f}")
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--json", action="store_true", help="emit JSON instead of text"
    )
    args = parser.parse_args(argv)
    report = aggregate()
    if "error" in report:
        print(report["error"], file=sys.stderr)
        return 1
    if args.json:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        print(_format_text(report), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
