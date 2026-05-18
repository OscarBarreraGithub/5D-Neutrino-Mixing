#!/usr/bin/env python3
"""Ingest flavor-catalog YAML sidecars + .tex prose into JSON for the Astro site.

Walks ``flavor_catalog/processes/**/*.yaml``, pairs each with its sibling
``.tex`` source, extracts the canonical prose sections (handling header-name
variants observed across the 102 entries), and writes:

  * ``src/content/entries/<process_id>.json`` -- full payload per entry
  * ``src/content/catalog_index.json``        -- summary used by the landing page
  * ``src/content/families.json``             -- per-family counts/labels

Re-runnable: every invocation regenerates everything from source.

Usage:
    python3 scripts/ingest_catalog.py

No third-party dependencies are required; YAML parsing uses a small purpose-built
parser tuned to the flavor-catalog sidecar dialect (block mappings/sequences,
scalars and quoted strings).  PyYAML is preferred if available.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

HERE = Path(__file__).resolve().parent
SITE_ROOT = HERE.parent
REPO_ROOT = SITE_ROOT.parents[1]
PROCESSES_DIR = REPO_ROOT / "flavor_catalog" / "processes"
OUT_DIR = SITE_ROOT / "src" / "content" / "entries"
INDEX_PATH = SITE_ROOT / "src" / "content" / "catalog_index.json"
FAMILIES_PATH = SITE_ROOT / "src" / "content" / "families.json"

# ---------------------------------------------------------------------------
# YAML loader: prefer PyYAML, fall back to a minimal hand-rolled parser.
# ---------------------------------------------------------------------------

try:  # pragma: no cover - environment-dependent
    import yaml  # type: ignore

    def load_yaml(text: str) -> Any:
        return yaml.safe_load(text)

except ImportError:  # pragma: no cover
    def load_yaml(text: str) -> Any:
        return _MinimalYAML(text).parse()


class _MinimalYAML:
    """Small YAML subset parser sufficient for flavor-catalog sidecars.

    Supports block mappings, block sequences, plain/quoted scalars, ints,
    floats (incl. scientific), booleans, null, and nested combinations.
    Comments (``# ...``) and blank lines are skipped.  Flow style and anchors
    are NOT supported -- the catalog does not use them.
    """

    _NUM_RE = re.compile(r"^-?\d+(\.\d+)?([eE][+-]?\d+)?$")
    _INT_RE = re.compile(r"^-?\d+$")

    def __init__(self, text: str) -> None:
        self.lines: list[tuple[int, str]] = []
        for raw in text.splitlines():
            if not raw.strip() or raw.lstrip().startswith("#"):
                continue
            stripped = raw.rstrip()
            indent = len(stripped) - len(stripped.lstrip(" "))
            self.lines.append((indent, stripped))
        self.pos = 0

    def parse(self) -> Any:
        if not self.lines:
            return {}
        return self._parse_block(0)

    # -- core -----------------------------------------------------------------

    def _peek(self) -> tuple[int, str] | None:
        if self.pos >= len(self.lines):
            return None
        return self.lines[self.pos]

    def _parse_block(self, indent: int) -> Any:
        """Parse a block node at the given indent.  Returns dict, list, or scalar."""
        node = self._peek()
        if node is None:
            return None
        if node[0] < indent:
            return None
        # Sequence?
        if node[1].lstrip().startswith("- "):
            return self._parse_sequence(node[0])
        return self._parse_mapping(node[0])

    def _parse_mapping(self, indent: int) -> dict:
        out: dict = {}
        while True:
            node = self._peek()
            if node is None or node[0] < indent:
                return out
            if node[0] > indent:
                # Defensive: shouldn't normally happen at mapping start.
                return out
            line = node[1][indent:]
            if line.startswith("- "):
                return out
            key, sep, rest = line.partition(":")
            if not sep:
                # Stray line; bail out.
                return out
            key = key.strip()
            self.pos += 1
            rest = rest.strip()
            if rest == "" or rest == "|" or rest == ">":
                # Nested block (mapping or sequence) starts on next line.
                nxt = self._peek()
                if nxt is None or nxt[0] <= indent:
                    out[key] = None
                else:
                    out[key] = self._parse_block(nxt[0])
            else:
                out[key] = self._scalar(rest)
        # unreachable

    def _parse_sequence(self, indent: int) -> list:
        out: list = []
        while True:
            node = self._peek()
            if node is None or node[0] < indent:
                return out
            line = node[1][indent:]
            if not line.startswith("- "):
                return out
            self.pos += 1
            item_body = line[2:]
            # Inline scalar item: "- foo" or "- key: value"
            if ":" in item_body and not item_body.startswith('"'):
                # Could be a mapping item with first key inline.
                first_key, sep, first_val = item_body.partition(":")
                first_val = first_val.strip()
                # Build a synthetic mapping using inline indent + any following lines.
                item: dict = {first_key.strip(): self._scalar(first_val) if first_val else None}
                # If the first value is empty (block), it may have nested children.
                if not first_val:
                    nxt = self._peek()
                    if nxt is not None and nxt[0] > indent:
                        item[first_key.strip()] = self._parse_block(nxt[0])
                # Then parse remaining sibling keys at indent+2 (relative to the dash).
                child_indent = indent + 2
                while True:
                    nxt = self._peek()
                    if nxt is None or nxt[0] < child_indent:
                        break
                    if nxt[0] == child_indent and not nxt[1][child_indent:].startswith("- "):
                        sub = self._parse_mapping(child_indent)
                        item.update(sub)
                    else:
                        break
                out.append(item)
            else:
                out.append(self._scalar(item_body))
        # unreachable

    def _scalar(self, raw: str) -> Any:
        s = raw.strip()
        if s == "" or s.lower() == "null" or s == "~":
            return None
        if s.lower() == "true":
            return True
        if s.lower() == "false":
            return False
        if len(s) >= 2 and s[0] == s[-1] and s[0] in ('"', "'"):
            inner = s[1:-1]
            if s[0] == '"':
                # decode simple escapes
                inner = inner.replace('\\"', '"').replace("\\\\", "\\")
            return inner
        if self._INT_RE.match(s):
            try:
                return int(s)
            except ValueError:
                return s
        if self._NUM_RE.match(s):
            try:
                return float(s)
            except ValueError:
                return s
        return s


# ---------------------------------------------------------------------------
# .tex section extraction
# ---------------------------------------------------------------------------

SECTION_RE = re.compile(r"^\\subsection\*\{(?P<title>[^}]+)\}\s*$", re.MULTILINE)

# Map of normalized prose-section keys to a list of regex patterns that match
# the various header strings observed in the catalog.
SECTION_PATTERNS: dict[str, list[re.Pattern[str]]] = {
    "relevance": [
        re.compile(r"^relevance to the rs / anarchic-flavor pipeline$", re.IGNORECASE),
        re.compile(r"^relevance to rs / anarchic flavor$", re.IGNORECASE),
        re.compile(r"^relevance to rs / anarchic-flavor pipeline$", re.IGNORECASE),
        re.compile(r"^relevance to rs$", re.IGNORECASE),
        re.compile(r"^relevance$", re.IGNORECASE),
    ],
    "post_2008": [
        re.compile(r"^post-2008 developments$", re.IGNORECASE),
        re.compile(r"^post-2010 developments$", re.IGNORECASE),
        re.compile(r"^post 2008 developments$", re.IGNORECASE),
        re.compile(r"^post 2010 developments$", re.IGNORECASE),
    ],
    "validity": [
        re.compile(r"^constraint validity and model dependence$", re.IGNORECASE),
        re.compile(r"^constraint validity / model dependence$", re.IGNORECASE),
        re.compile(r"^validity / model dependence$", re.IGNORECASE),
        re.compile(r"^validity$", re.IGNORECASE),
    ],
    "code_coverage_prose": [
        re.compile(r"^code coverage in this repo$", re.IGNORECASE),
        re.compile(r"^code coverage$", re.IGNORECASE),
    ],
    "implementation_difficulty_prose": [
        re.compile(r"^implementation difficulty$", re.IGNORECASE),
    ],
    "key_references": [
        re.compile(r"^key references$", re.IGNORECASE),
    ],
    "pdg_prose": [
        re.compile(r"^pdg / equivalent value\(s\)$", re.IGNORECASE),
        re.compile(r"^pdg or equivalent value$", re.IGNORECASE),
        re.compile(r"^pdg-or-equivalent value$", re.IGNORECASE),
        re.compile(r"^pdg-or-equivalent value\(s\)$", re.IGNORECASE),
        re.compile(r"^pdg value$", re.IGNORECASE),
        re.compile(r"^current best limit\(s\)$", re.IGNORECASE),
    ],
    "process_prose": [
        re.compile(r"^process$", re.IGNORECASE),
    ],
}


def classify_header(title: str) -> str | None:
    title_clean = title.strip()
    for key, patterns in SECTION_PATTERNS.items():
        for pat in patterns:
            if pat.match(title_clean):
                return key
    return None


def split_sections(tex: str) -> dict[str, str]:
    """Return a dict of normalized prose-section key -> raw LaTeX body."""
    matches = list(SECTION_RE.finditer(tex))
    sections: dict[str, str] = {}
    for i, m in enumerate(matches):
        key = classify_header(m.group("title"))
        if key is None:
            continue
        body_start = m.end()
        body_end = matches[i + 1].start() if i + 1 < len(matches) else len(tex)
        body = tex[body_start:body_end].strip()
        # If the same normalized key appears more than once, keep the first.
        sections.setdefault(key, body)
    return sections


# ---------------------------------------------------------------------------
# Per-entry transformation
# ---------------------------------------------------------------------------

# Display names + one-line descriptors for the 8 family cards on the landing page.
FAMILY_META: dict[str, dict[str, str]] = {
    "kaon": {
        "label": "Kaon",
        "descriptor": "Neutral kaon mixing, rare K decays, CP violation in the strange sector.",
    },
    "charm": {
        "label": "Charm",
        "descriptor": "D-meson mixing and rare charm decays (Delta C=2 and Delta C=1).",
    },
    "beauty": {
        "label": "Beauty",
        "descriptor": "B-meson mixing, b->s/d transitions, lepton-flavor-universality tests.",
    },
    "top_higgs_ew": {
        "label": "Top / Higgs / EW",
        "descriptor": "Electroweak precision, top/Higgs couplings, custodial sector observables.",
    },
    "charged_lepton": {
        "label": "Charged Lepton",
        "descriptor": "Lepton-flavor-violating dipole and conversion processes.",
    },
    "edm_neutrino": {
        "label": "EDM / Neutrino",
        "descriptor": "Electric dipole moments and neutrino oscillation constraints.",
    },
    "collider_rs": {
        "label": "Collider RS",
        "descriptor": "Direct LHC searches for KK gauge resonances and custodial top partners.",
    },
    "secondary": {
        "label": "Secondary (deferred)",
        "descriptor": "Deferred SECONDARY scope -- LFV-kaon trio and lepton-bulk extensions promoted in Wave-8.",
    },
}

FAMILY_ORDER = [
    "kaon",
    "charm",
    "beauty",
    "top_higgs_ew",
    "charged_lepton",
    "edm_neutrino",
    "collider_rs",
    "secondary",
]


def derive_family(yaml_path: Path, sidecar: dict) -> str:
    """Return the family bucket used for the landing page.

    SECONDARY entries live under processes/secondary/<original-family>/; they
    are grouped under the 'secondary' card regardless of their original family.
    """
    rel = yaml_path.relative_to(PROCESSES_DIR)
    parts = rel.parts
    if parts and parts[0] == "secondary":
        return "secondary"
    return sidecar.get("family") or (parts[0] if parts else "unknown")


def normalize_pdg_values(sidecar: dict) -> list[dict]:
    """Return a flat list of PDG-equivalent value rows for table rendering."""
    pdg = sidecar.get("pdg_or_equivalent") or {}
    rows: list[dict] = []

    # Many entries use a top-level pdg_or_equivalent block as a single value.
    # Newer entries also include a `values: [ ... ]` array; prefer that when
    # present, otherwise treat the top block as a single row.
    if isinstance(pdg, dict):
        if isinstance(pdg.get("values"), list):
            for v in pdg["values"]:
                if isinstance(v, dict):
                    rows.append(_row_from(v))
        # Sub-keyed block (canonical_experimental_value, standard_model_reference, ...)
        measurement_keys = {"value", "display", "display_value", "total_display_value", "observable"}
        if not rows and not measurement_keys.intersection(pdg.keys()):
            for sub_key, sub_val in pdg.items():
                if not isinstance(sub_val, dict):
                    continue
                # Skip sub-blocks that don't carry a measurable value -- those
                # are typically auxiliary lattice inputs or paper-era context.
                if not measurement_keys.intersection(sub_val.keys()):
                    continue
                row = _row_from(sub_val)
                row.setdefault("label", sub_key.replace("_", " "))
                rows.append(row)
        # Fallback: the top block itself is one row.
        if not rows and measurement_keys.intersection(pdg.keys()):
            rows.append(_row_from(pdg))
    return rows


def _row_from(v: dict) -> dict:
    return {
        "label": v.get("value_id") or v.get("source_key") or v.get("source") or v.get("observable"),
        "observable": v.get("observable"),
        "value": v.get("value"),
        "uncertainty": v.get("uncertainty"),
        "units": v.get("units"),
        "display": v.get("display") or v.get("display_value") or v.get("total_display_value"),
        "year": v.get("year"),
        "experiment": v.get("experiment") or v.get("source"),
        "source_url": v.get("source_url"),
        "snapshot_path": v.get("snapshot_path"),
        "sha256": v.get("sha256") or v.get("sha256_of_local_snapshot"),
        "access_date": v.get("access_date"),
        "cl": v.get("cl"),
        "limit_type": v.get("limit_type"),
    }


def derive_tier(sidecar: dict, yaml_path: Path) -> str:
    tier = sidecar.get("priority_tier")
    if tier:
        return tier
    rel = yaml_path.relative_to(PROCESSES_DIR)
    if rel.parts and rel.parts[0] == "secondary":
        return "SECONDARY"
    return "PRIMARY"


def latest_status_state(sidecar: dict) -> str | None:
    hist = sidecar.get("status_history")
    if isinstance(hist, list) and hist:
        last = hist[-1]
        if isinstance(last, dict):
            return last.get("state") or last.get("to")
    return None


def collect_access_dates(sidecar: dict) -> list[str]:
    """Walk the sidecar looking for access_date strings; dedupe + sort."""
    dates: set[str] = set()

    def walk(node: Any) -> None:
        if isinstance(node, dict):
            for k, v in node.items():
                if k == "access_date" and isinstance(v, str):
                    dates.add(v)
                walk(v)
        elif isinstance(node, list):
            for it in node:
                walk(it)

    walk(sidecar)
    return sorted(dates)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> int:
    if not PROCESSES_DIR.is_dir():
        print(f"ERROR: processes dir not found: {PROCESSES_DIR}", file=sys.stderr)
        return 2

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    # Wipe stale per-entry JSON so the ingest is authoritative.
    for stale in OUT_DIR.glob("*.json"):
        stale.unlink()

    yaml_paths = sorted(PROCESSES_DIR.rglob("*.yaml"))
    print(f"Found {len(yaml_paths)} entry sidecars under {PROCESSES_DIR.relative_to(REPO_ROOT)}")

    index_rows: list[dict] = []
    family_counts: dict[str, int] = {fam: 0 for fam in FAMILY_ORDER}
    verdict_counts: dict[str, int] = {}

    for yp in yaml_paths:
        try:
            sidecar = load_yaml(yp.read_text())
        except Exception as exc:  # noqa: BLE001
            print(f"  ! parse failed: {yp.relative_to(REPO_ROOT)} -- {exc}", file=sys.stderr)
            continue
        if not isinstance(sidecar, dict):
            print(f"  ! non-mapping YAML at {yp.relative_to(REPO_ROOT)}", file=sys.stderr)
            continue

        process_id = sidecar.get("process_id") or yp.stem
        family = derive_family(yp, sidecar)
        tex_path = yp.with_suffix(".tex")
        sections: dict[str, str] = {}
        if tex_path.is_file():
            sections = split_sections(tex_path.read_text())
        else:
            print(f"  ! missing .tex for {process_id}", file=sys.stderr)

        tier = derive_tier(sidecar, yp)
        verdict = sidecar.get("fact_check_verdict") or "UNKNOWN"
        difficulty = sidecar.get("implementation_difficulty") or "UNKNOWN"
        cov = sidecar.get("code_coverage") or {}
        cov_status = (cov.get("status") if isinstance(cov, dict) else None) or "UNKNOWN"
        pdg_values = normalize_pdg_values(sidecar)

        rel_yaml = yp.relative_to(REPO_ROOT).as_posix()
        rel_tex = tex_path.relative_to(REPO_ROOT).as_posix() if tex_path.is_file() else None
        worklog_path = f"flavor_catalog/worklogs/{process_id}/"

        payload = {
            "process_id": process_id,
            "family": family,
            "family_original": sidecar.get("family"),
            "tier": tier,
            "standard_notation": sidecar.get("standard_notation"),
            "process_name": sidecar.get("process_name"),
            "owner": sidecar.get("owner"),
            "fact_check_verdict": verdict,
            "implementation_difficulty": difficulty,
            "implementation_difficulty_reason": (
                sidecar.get("implementation_difficulty_reason")
                or sidecar.get("implementation_difficulty_rationale")
            ),
            "cycle_count": sidecar.get("cycle_count"),
            "priority_rationale": sidecar.get("priority_rationale"),
            "promoted_in_wave": sidecar.get("promoted_in_wave"),
            "code_coverage": cov if isinstance(cov, dict) else None,
            "code_coverage_status": cov_status,
            "pdg_values": pdg_values,
            "supporting_measurements": sidecar.get("supporting_measurements"),
            "paper_era_reference": sidecar.get("paper_era_reference"),
            "theory_context": sidecar.get("theory_context"),
            "auxiliary_theory_inputs": sidecar.get("auxiliary_theory_inputs"),
            "post_2010_developments": sidecar.get("post_2010_developments")
                or sidecar.get("post_2008_developments"),
            "source_shas": sidecar.get("source_shas"),
            "source_shas_count": (
                len(sidecar.get("source_shas")) if isinstance(sidecar.get("source_shas"), dict) else 0
            ),
            "access_dates": collect_access_dates(sidecar),
            "latest_status": latest_status_state(sidecar),
            "status_history": sidecar.get("status_history"),
            "sections": {
                "process": sections.get("process_prose"),
                "pdg_prose": sections.get("pdg_prose"),
                "relevance": sections.get("relevance"),
                "post_2008": sections.get("post_2008"),
                "validity": sections.get("validity"),
                "code_coverage_prose": sections.get("code_coverage_prose"),
                "implementation_difficulty_prose": sections.get("implementation_difficulty_prose"),
                "key_references": sections.get("key_references"),
            },
            "source_yaml": rel_yaml,
            "source_tex": rel_tex,
            "worklog_path": worklog_path,
        }

        out_path = OUT_DIR / f"{process_id}.json"
        out_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=False))

        family_counts[family] = family_counts.get(family, 0) + 1
        verdict_counts[verdict] = verdict_counts.get(verdict, 0) + 1
        index_rows.append({
            "process_id": process_id,
            "family": family,
            "tier": tier,
            "fact_check_verdict": verdict,
            "standard_notation": sidecar.get("standard_notation"),
            "process_name": sidecar.get("process_name"),
            "implementation_difficulty": difficulty,
            "code_coverage_status": cov_status,
        })

    index_rows.sort(key=lambda r: (FAMILY_ORDER.index(r["family"]) if r["family"] in FAMILY_ORDER else 99, r["process_id"]))

    INDEX_PATH.write_text(json.dumps({
        "total_entries": len(index_rows),
        "verdict_counts": verdict_counts,
        "family_counts": family_counts,
        "entries": index_rows,
    }, indent=2, ensure_ascii=False))

    families_payload = []
    for fam in FAMILY_ORDER:
        meta = FAMILY_META.get(fam, {"label": fam, "descriptor": ""})
        families_payload.append({
            "key": fam,
            "label": meta["label"],
            "descriptor": meta["descriptor"],
            "count": family_counts.get(fam, 0),
        })
    FAMILIES_PATH.write_text(json.dumps(families_payload, indent=2, ensure_ascii=False))

    print(f"  -> {len(index_rows)} entries written to {OUT_DIR.relative_to(SITE_ROOT)}")
    print(f"  -> catalog_index.json: {INDEX_PATH.relative_to(SITE_ROOT)}")
    print(f"  -> families.json:      {FAMILIES_PATH.relative_to(SITE_ROOT)}")
    print(f"  -> family counts:      {family_counts}")
    print(f"  -> verdict counts:     {verdict_counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
