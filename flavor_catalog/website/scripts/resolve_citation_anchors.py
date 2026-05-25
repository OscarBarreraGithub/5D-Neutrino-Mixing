#!/usr/bin/env python3
"""Resolve line-level citation anchors for flavor-catalog value claims.

This script is intentionally conservative: it records all matching line
contexts when a value is ambiguous, and marks a value unresolved when the local
snapshot does not contain a sufficiently clear anchor.

Usage:
    python3 flavor_catalog/website/scripts/resolve_citation_anchors.py \
        --family top_higgs_ew --include-secondary top_higgs_ew
    python3 flavor_catalog/website/scripts/resolve_citation_anchors.py --ids T001
"""
from __future__ import annotations

import argparse
import re
import sys
import unicodedata
from dataclasses import dataclass
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any

try:
    import yaml  # type: ignore
except ImportError as exc:  # pragma: no cover
    raise SystemExit("PyYAML is required for citation-anchor resolution") from exc


class LiteralString(str):
    pass


def literal_string_representer(dumper: yaml.SafeDumper, data: LiteralString) -> yaml.nodes.ScalarNode:
    return dumper.represent_scalar("tag:yaml.org,2002:str", data, style="|")


yaml.add_representer(LiteralString, literal_string_representer, Dumper=yaml.SafeDumper)

HERE = Path(__file__).resolve().parent
SITE_ROOT = HERE.parent
REPO_ROOT = SITE_ROOT.parents[1]
PROCESSES_DIR = REPO_ROOT / "flavor_catalog" / "processes"
REFERENCES_DIR = REPO_ROOT / "flavor_catalog" / "references"
OUT_DIR = SITE_ROOT / "_data" / "citation_anchors"

MEASUREMENT_KEYS = {
    "value",
    "display",
    "display_value",
    "total_display_value",
    "observable",
}

STOPWORDS = {
    "about",
    "access",
    "and",
    "atlas",
    "branching",
    "cms",
    "combined",
    "confidence",
    "context",
    "dataset",
    "direct",
    "equivalent",
    "experimental",
    "fit",
    "from",
    "global",
    "headline",
    "level",
    "limit",
    "observable",
    "observed",
    "pdg",
    "process",
    "projection",
    "source",
    "standard",
    "summary",
    "theory",
    "upper",
    "used",
    "value",
    "with",
    "year",
}


@dataclass(frozen=True)
class ValueBlock:
    process_id: str
    block_key: str
    raw: dict[str, Any]
    parent: dict[str, Any]
    index: int


@dataclass(frozen=True)
class MatchResult:
    anchor_field: str
    anchor_string: str
    status: str
    line_numbers: tuple[int, ...]


def load_yaml(path: Path) -> Any:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def clean_sha(value: Any) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def scalar_to_string(value: Any) -> str | None:
    if value is None:
        return None
    if isinstance(value, bool):
        return "true" if value else "false"
    text = str(value).strip()
    return text or None


def normalize_text(text: str) -> str:
    text = unicodedata.normalize("NFKC", text)
    replacements = {
        "\u00a0": " ",
        "\u2010": "-",
        "\u2011": "-",
        "\u2012": "-",
        "\u2013": "-",
        "\u2014": "-",
        "\u2212": "-",
        "\u2192": " -> ",
        "\u27f6": " -> ",
        "\u00b1": " +/- ",
        "\u00d7": " x ",
        "\u03c1": "rho",
        "\u03c3": "sigma",
        "\u0393": "Gamma",
        "\u03b3": "gamma",
        "\u039b": "Lambda",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    latex_replacements = {
        r"\pm": " +/- ",
        r"\times": " x ",
        r"\cdot": " x ",
        r"\to": " -> ",
        r"\rightarrow": " -> ",
        r"\sqrt": "sqrt",
        r"\mathrm": "",
        r"\text": "",
        r"\bar": "bar",
        r"\left": "",
        r"\right": "",
    }
    for old, new in latex_replacements.items():
        text = text.replace(old, new)
    text = text.replace("~", " ")
    text = text.replace("\\,", " ")
    text = text.replace("\\", "")
    text = re.sub(r"[{}]", "", text)
    text = re.sub(r"\s*/\s*-\s*", "/-", text)
    text = re.sub(r"\+\s*/\s*-", "+/-", text)
    text = re.sub(r"\bper\s+cent\b", "%", text, flags=re.IGNORECASE)
    text = re.sub(r"\binverse\s+femtobarns?\b", "fb^-1", text, flags=re.IGNORECASE)
    text = re.sub(r"\bfemtobarns?\s*\^\s*-?1\b", "fb^-1", text, flags=re.IGNORECASE)
    text = re.sub(r"\bbr\s*\(", "b(", text, flags=re.IGNORECASE)
    text = re.sub(r"\s+", " ", text)
    return text.strip().lower()


SCI_RE = re.compile(
    r"(?P<num>[+-]?\d+(?:\.\d+)?)\s*(?:x|times)\s*10\s*\^?\s*(?P<exp>[+-]?\d+)"
)
E_RE = re.compile(r"(?P<num>[+-]?\d+(?:\.\d+)?)\s*e\s*(?P<exp>[+-]?\d+)")
BARE_POWER_RE = re.compile(r"(?<![\w.])10\s*\^?\s*(?P<exp>[+-]?\d+)")


def decimal_scientific(match: re.Match[str]) -> str:
    num = match.group("num")
    exp = match.group("exp")
    try:
        value = Decimal(num) * (Decimal(10) ** int(exp))
    except (InvalidOperation, ValueError):
        return match.group(0)
    return format(value.normalize(), "f").rstrip("0").rstrip(".") or "0"


def canonical_for_match(text: str) -> str:
    text = normalize_text(text)
    text = SCI_RE.sub(decimal_scientific, text)
    text = E_RE.sub(decimal_scientific, text)
    text = BARE_POWER_RE.sub(lambda match: decimal_scientific(match_with_unit_mantissa(match)), text)
    text = re.sub(r"\s*([<>=(),:%^+\-/])\s*", r"\1", text)
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def match_with_unit_mantissa(match: re.Match[str]) -> re.Match[str]:
    class _Match:
        def group(self, key: str | int = 0) -> str:
            if key == "num":
                return "1"
            if key == "exp":
                return match.group("exp")
            return match.group(0)

    return _Match()  # type: ignore[return-value]


def compact(text: str) -> str:
    return re.sub(r"\s+", "", canonical_for_match(text))


def parenthetical_uncertainty(value: Any, uncertainty: Any) -> str | None:
    value_text = scalar_to_string(value)
    uncertainty_text = scalar_to_string(uncertainty)
    if not value_text or not uncertainty_text:
        return None
    if not re.fullmatch(r"[+-]?\d+(?:\.\d+)?", value_text):
        return None
    if not re.fullmatch(r"[+-]?\d+(?:\.\d+)?", uncertainty_text):
        return None
    if "." not in value_text:
        return None
    decimals = len(value_text.split(".", 1)[1])
    try:
        uncertainty_decimal = Decimal(uncertainty_text).copy_abs()
    except InvalidOperation:
        return None
    scaled = uncertainty_decimal * (Decimal(10) ** decimals)
    digits = str(int(scaled.to_integral_value()))
    if not digits or int(digits) == 0:
        return None
    return f"{value_text}({digits})"


def qualifier_variants(text: str) -> list[str]:
    variants = [text]
    stripped = re.sub(r"^(approximately|approx\.?|about|order of|order)\s+", "", text, flags=re.I).strip()
    if stripped and stripped != text:
        variants.append(stripped)
    variants.append(text.replace("approximately", "approx."))
    variants.append(text.replace("approximately", ""))
    variants.append(text.replace("about", ""))
    variants.append(text.replace("order ", ""))
    return unique_nonempty(variants)


def scientific_variants(text: str) -> list[str]:
    variants = [text]
    for match in re.finditer(r"(?P<num>[+-]?\d+(?:\.\d+)?)e(?P<exp>[+-]?\d+)", text, flags=re.I):
        num = match.group("num")
        exp = match.group("exp")
        replacement = f"{num} x 10^{exp}"
        variants.append(text[: match.start()] + replacement + text[match.end() :])
        try:
            dec = Decimal(num) * (Decimal(10) ** int(exp))
            variants.append(text[: match.start()] + format(dec.normalize(), "f") + text[match.end() :])
        except (InvalidOperation, ValueError):
            pass
    return unique_nonempty(variants)


def pdg_unit_variants(value: Any, pdg_units: Any) -> list[str]:
    value_text = scalar_to_string(value) or ""
    unit_text = scalar_to_string(pdg_units)
    if not unit_text:
        return []
    variants = [unit_text]
    if value_text.strip().startswith("<") and not unit_text.strip().startswith("<"):
        variants.append(f"< {unit_text}")
    first_number = re.search(r"[+-]?\d+(?:\.\d+)?", unit_text)
    if first_number:
        number = first_number.group(0)
        variants.append(number)
        if value_text.strip().startswith("<"):
            variants.append(f"< {number}")
    return unique_nonempty(variants)


def unique_nonempty(values: list[str]) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for value in values:
        text = value.strip()
        if not text:
            continue
        key = canonical_for_match(text)
        if key in seen:
            continue
        seen.add(key)
        out.append(text)
    return out


def candidate_anchors(block: dict[str, Any]) -> list[tuple[str, str]]:
    candidates: list[tuple[str, str]] = []

    for field in ("display_value", "display", "total_display_value"):
        text = scalar_to_string(block.get(field))
        if text:
            candidates.append((field, text))

    for field in sorted(k for k in block if k.startswith("related_")):
        value = block.get(field)
        if isinstance(value, list):
            for item in value:
                text = scalar_to_string(item)
                if text:
                    candidates.append((field, text))
        elif isinstance(value, dict):
            for item in value.values():
                text = scalar_to_string(item)
                if text:
                    candidates.append((field, text))
        else:
            text = scalar_to_string(value)
            if text:
                candidates.append((field, text))

    value = block.get("value")
    uncertainty = block.get("uncertainty")
    value_text = scalar_to_string(value)
    uncertainty_text = scalar_to_string(uncertainty)
    if value_text and uncertainty_text:
        candidates.append(("value_uncertainty", f"{value_text} +/- {uncertainty_text}"))
        paren = parenthetical_uncertainty(value, uncertainty)
        if paren:
            candidates.append(("value_uncertainty", paren))
    if value_text:
        candidates.append(("value", value_text))
    norm = scalar_to_string(block.get("normalized_value"))
    if norm:
        candidates.append(("normalized_value", norm))
    for variant in pdg_unit_variants(value, block.get("pdg_units")):
        candidates.append(("pdg_units", variant))

    deduped: list[tuple[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for field, text in candidates:
        key = (field, canonical_for_match(text))
        if key not in seen:
            seen.add(key)
            deduped.append((field, text))
    return deduped


def candidate_variants(anchor_string: str) -> list[str]:
    variants = [anchor_string]
    for text in list(variants):
        variants.extend(scientific_variants(text))
    for text in list(variants):
        variants.extend(qualifier_variants(text))
    normalized = normalize_text(anchor_string)
    luminosity_tokens = re.findall(r"\d+(?:\.\d+)?\s*fb\^-?1", normalized)
    variants.extend(luminosity_tokens)
    if not luminosity_tokens:
        variants.extend(re.findall(r"sqrt\s*\(?s\)?\s*=\s*\d+(?:\.\d+)?\s*(?:tev|gev)", normalized))
        variants.extend(re.findall(r"\d+(?:\.\d+)?\s*(?:tev|gev)", normalized))
    variants.extend([v.replace(" x 10^", " * 10^") for v in list(variants)])
    variants.extend([v.replace("sqrt(s)", "sqrt s") for v in list(variants)])
    variants.extend([v.replace("fb^-1", "inverse femtobarns") for v in list(variants)])
    variants.extend([v.replace("TeV", "tev").replace("GeV", "gev") for v in list(variants)])
    return unique_nonempty(variants)


def meaningful_hint_tokens(block: dict[str, Any]) -> set[str]:
    text = " ".join(
        str(block.get(field, ""))
        for field in ("observable", "value_id", "source", "limit_type", "role")
    )
    normalized = canonical_for_match(text)
    raw_tokens = set(re.findall(r"[a-zA-Z][a-zA-Z0-9_+\-]{1,}", normalized))
    tokens: set[str] = set()
    for token in raw_tokens:
        plain = token.strip("_-")
        if len(plain) < 2 or plain in STOPWORDS:
            continue
        tokens.add(plain)
        collapsed = plain.replace("_", "").replace("-", "")
        if len(collapsed) >= 2:
            tokens.add(collapsed)
    return tokens


def window_context(lines: list[str], line_number: int) -> str:
    start = max(1, line_number - 3)
    end = min(len(lines), line_number + 3)
    out: list[str] = []
    for i in range(start, end + 1):
        body = lines[i - 1].rstrip()
        out.append(f"L{i}: {body}\n" if body else f"L{i}:\n")
    return "".join(out)


def choose_line_for_window(
    line_numbers: list[int],
    canonical_lines: list[str],
    canonical_variant: str,
) -> int:
    tokens = re.findall(r"\d+(?:\.\d+)?|[a-zA-Z][a-zA-Z0-9_+\-]{1,}", canonical_variant)
    preferred = [token for token in tokens if token not in STOPWORDS]
    for token in preferred:
        for line_number in line_numbers:
            if token in canonical_lines[line_number - 1]:
                return line_number
    return line_numbers[0]


def search_snapshot(
    lines: list[str],
    anchor_string: str,
    block: dict[str, Any],
) -> tuple[int, ...]:
    canonical_lines = [canonical_for_match(line) for line in lines]
    compact_lines = [compact(line) for line in lines]
    line_matches: set[int] = set()

    for variant in candidate_variants(anchor_string):
        canonical_variant = canonical_for_match(variant)
        compact_variant = compact(variant)
        if not canonical_variant:
            continue
        for idx, canonical_line in enumerate(canonical_lines, start=1):
            if canonical_variant in canonical_line or (
                len(compact_variant) >= 3 and compact_variant in compact_lines[idx - 1]
            ):
                line_matches.add(idx)

        if line_matches:
            continue

        for width in range(2, 5):
            for start in range(0, len(lines) - width + 1):
                line_numbers = list(range(start + 1, start + width + 1))
                window = " ".join(canonical_lines[start : start + width])
                compact_window = "".join(compact_lines[start : start + width])
                if canonical_variant in window or (
                    len(compact_variant) >= 3 and compact_variant in compact_window
                ):
                    line_matches.add(
                        choose_line_for_window(line_numbers, canonical_lines, canonical_variant)
                    )

    if len(line_matches) <= 1:
        return tuple(sorted(line_matches))

    hints = meaningful_hint_tokens(block)
    if not hints:
        return tuple(sorted(line_matches))

    scored: list[tuple[int, int]] = []
    for line_number in sorted(line_matches):
        start = max(0, line_number - 2)
        end = min(len(lines), line_number + 1)
        context = canonical_for_match(" ".join(lines[start:end]))
        score = sum(1 for token in hints if token in context)
        scored.append((score, line_number))
    best_score = max(score for score, _ in scored)
    if best_score > 0:
        best_lines = [line_number for score, line_number in scored if score == best_score]
        if len(best_lines) == 1:
            return (best_lines[0],)
    return tuple(sorted(line_matches))


def resolve_anchor(lines: list[str], block: dict[str, Any]) -> MatchResult:
    first_candidate: tuple[str, str] | None = None
    first_ambiguous: MatchResult | None = None

    for field, anchor_string in candidate_anchors(block):
        if first_candidate is None:
            first_candidate = (field, anchor_string)
        matches = search_snapshot(lines, anchor_string, block)
        if len(matches) == 1:
            return MatchResult(field, anchor_string, "RESOLVED", matches)
        if len(matches) > 1 and first_ambiguous is None:
            first_ambiguous = MatchResult(field, anchor_string, "AMBIGUOUS", matches)

    if first_ambiguous:
        return first_ambiguous
    if first_candidate:
        return MatchResult(first_candidate[0], first_candidate[1], "UNRESOLVED", ())
    return MatchResult("value", "", "UNRESOLVED", ())


def source_manifest_maps(process_id: str) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, Any]]]:
    manifest_path = REFERENCES_DIR / process_id / "source_manifest.yaml"
    by_sha: dict[str, dict[str, Any]] = {}
    by_url: dict[str, dict[str, Any]] = {}
    if manifest_path.exists():
        manifest = load_yaml(manifest_path) or {}
        for source in manifest.get("sources", []) or []:
            if not isinstance(source, dict):
                continue
            sha = clean_sha(source.get("sha256_of_local_snapshot") or source.get("sha256"))
            url = scalar_to_string(source.get("url") or source.get("source_url"))
            if sha:
                by_sha[sha] = source
            if url:
                by_url[url] = source
    return by_sha, by_url


def source_shas_map(sidecar: dict[str, Any], process_id: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for filename, sha in (sidecar.get("source_shas") or {}).items():
        sha_text = clean_sha(sha)
        if sha_text:
            out[sha_text] = str(REFERENCES_DIR / process_id / filename)
    return out


def repo_relative(path: Path) -> str:
    try:
        return path.resolve().relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return path.as_posix()


def resolve_snapshot_metadata(
    sidecar: dict[str, Any],
    value: dict[str, Any],
    parent: dict[str, Any],
    process_id: str,
    manifest_by_sha: dict[str, dict[str, Any]],
    manifest_by_url: dict[str, dict[str, Any]],
    shas_to_paths: dict[str, str],
) -> dict[str, Any] | None:
    sha = clean_sha(value.get("sha256") or value.get("sha256_of_local_snapshot"))
    source_url = scalar_to_string(value.get("source_url"))
    snapshot_path = scalar_to_string(value.get("snapshot_path"))

    manifest_source = None
    if sha and sha in manifest_by_sha:
        manifest_source = manifest_by_sha[sha]
    elif source_url and source_url in manifest_by_url:
        manifest_source = manifest_by_url[source_url]

    if not snapshot_path and manifest_source:
        snapshot_path = scalar_to_string(manifest_source.get("snapshot_path"))
    if not snapshot_path and sha and sha in shas_to_paths:
        snapshot_path = repo_relative(Path(shas_to_paths[sha]))

    parent_sha = clean_sha(parent.get("sha256") or parent.get("sha256_of_local_snapshot"))
    parent_snapshot = scalar_to_string(parent.get("snapshot_path"))
    if not snapshot_path and parent_snapshot and (not sha or sha == parent_sha):
        snapshot_path = parent_snapshot
    if not snapshot_path:
        return None

    if not source_url and manifest_source:
        source_url = scalar_to_string(manifest_source.get("url") or manifest_source.get("source_url"))
    if not source_url:
        source_url = scalar_to_string(parent.get("source_url"))

    if not sha and manifest_source:
        sha = clean_sha(manifest_source.get("sha256_of_local_snapshot") or manifest_source.get("sha256"))
    if not sha:
        sha = parent_sha

    access_date = scalar_to_string(value.get("access_date"))
    if not access_date:
        access_date = scalar_to_string(parent.get("access_date"))

    return {
        "source_url": source_url,
        "snapshot_path": snapshot_path,
        "sha256": sha,
        "access_date": access_date,
    }


def merge_parent_value(parent: dict[str, Any], value: dict[str, Any]) -> dict[str, Any]:
    merged = dict(parent)
    merged.pop("values", None)
    merged.update(value)
    return merged


def block_key_for(value: dict[str, Any], fallback: str) -> str:
    for key in ("value_id", "source_key", "key", "observable"):
        text = scalar_to_string(value.get(key))
        if text:
            return text
    return fallback


def collect_value_blocks(process_id: str, sidecar: dict[str, Any]) -> list[ValueBlock]:
    pdg = sidecar.get("pdg_or_equivalent")
    blocks: list[ValueBlock] = []

    def add(value: dict[str, Any], parent: dict[str, Any], index: int, fallback: str) -> None:
        if not MEASUREMENT_KEYS.intersection(value.keys()):
            return
        blocks.append(
            ValueBlock(
                process_id=process_id,
                block_key=block_key_for(value, fallback),
                raw=value,
                parent=parent,
                index=index,
            )
        )

    if isinstance(pdg, list):
        for idx, item in enumerate(pdg):
            if not isinstance(item, dict):
                continue
            if isinstance(item.get("values"), list):
                for child_idx, child in enumerate(item["values"]):
                    if isinstance(child, dict):
                        add(child, item, child_idx, f"pdg_or_equivalent[{idx}].values[{child_idx}]")
            else:
                add(item, {}, idx, f"pdg_or_equivalent[{idx}]")
    elif isinstance(pdg, dict):
        if isinstance(pdg.get("values"), list):
            for idx, value in enumerate(pdg["values"]):
                if isinstance(value, dict):
                    add(value, pdg, idx, f"pdg_or_equivalent.values[{idx}]")
        elif MEASUREMENT_KEYS.intersection(pdg.keys()):
            add(pdg, {}, 0, "pdg_or_equivalent")
        else:
            idx = 0
            for key, value in pdg.items():
                if isinstance(value, dict):
                    add(value, pdg, idx, key)
                    idx += 1
    return blocks


def output_for_process(path: Path, generated_at: str) -> tuple[dict[str, Any], dict[str, int]]:
    sidecar = load_yaml(path)
    process_id = sidecar.get("process_id") or path.stem
    manifest_by_sha, manifest_by_url = source_manifest_maps(process_id)
    shas_to_paths = source_shas_map(sidecar, process_id)
    sources: list[dict[str, Any]] = []
    counts = {"anchors": 0, "RESOLVED": 0, "AMBIGUOUS": 0, "UNRESOLVED": 0}

    for value_block in collect_value_blocks(process_id, sidecar):
        value = merge_parent_value(value_block.parent, value_block.raw)
        metadata = resolve_snapshot_metadata(
            sidecar,
            value_block.raw,
            value_block.parent,
            process_id,
            manifest_by_sha,
            manifest_by_url,
            shas_to_paths,
        )
        if metadata is None:
            continue

        snapshot = REPO_ROOT / metadata["snapshot_path"]
        source_entry = {
            "block_key": value_block.block_key,
            "source_url": metadata.get("source_url"),
            "snapshot_path": metadata.get("snapshot_path"),
            "sha256": metadata.get("sha256"),
            "access_date": metadata.get("access_date"),
            "anchors": [],
        }

        if snapshot.exists():
            lines = snapshot.read_text(encoding="utf-8", errors="replace").splitlines()
            result = resolve_anchor(lines, value)
            matches = [
                {
                    "line_number": line_number,
                    "context": LiteralString(window_context(lines, line_number)),
                }
                for line_number in result.line_numbers
            ]
        else:
            result = MatchResult("value", scalar_to_string(value.get("value")) or "", "UNRESOLVED", ())
            matches = []

        source_entry["anchors"].append(
            {
                "anchor_field": result.anchor_field,
                "anchor_string": result.anchor_string,
                "status": result.status,
                "matches": matches,
            }
        )
        counts["anchors"] += 1
        counts[result.status] += 1
        sources.append(source_entry)

    return {
        "process_id": process_id,
        "generated_at": generated_at,
        "sources": sources,
    }, counts


def batch_paths(args: argparse.Namespace) -> list[Path]:
    paths: list[Path] = []
    if args.ids:
        requested = set(args.ids)
        for path in PROCESSES_DIR.glob("**/*.yaml"):
            try:
                sidecar = load_yaml(path)
            except Exception:
                continue
            process_id = sidecar.get("process_id") or path.stem
            if process_id in requested:
                paths.append(path)
    else:
        if args.family:
            paths.extend(sorted((PROCESSES_DIR / args.family).glob("*.yaml")))
        for family in args.include_secondary or []:
            paths.extend(sorted((PROCESSES_DIR / "secondary" / family).glob("*.yaml")))

    def sort_key(path: Path) -> tuple[str, str]:
        try:
            sidecar = load_yaml(path)
            process_id = sidecar.get("process_id") or path.stem
        except Exception:
            process_id = path.stem
        return (process_id, path.as_posix())

    return sorted(paths, key=sort_key)


def write_yaml(path: Path, payload: dict[str, Any]) -> None:
    text = "---\n" + yaml.safe_dump(
        payload,
        sort_keys=False,
        allow_unicode=False,
        width=120,
    )
    path.write_text(text, encoding="utf-8")


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--family", help="Primary family directory under flavor_catalog/processes")
    parser.add_argument(
        "--include-secondary",
        action="append",
        default=[],
        help="Secondary family directory under flavor_catalog/processes/secondary",
    )
    parser.add_argument("--ids", nargs="+", help="Resolve only these process IDs")
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = parser.parse_args(argv)

    paths = batch_paths(args)
    if not paths:
        raise SystemExit("No process YAML files matched the requested batch")

    generated_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    totals = {"entries": 0, "anchors": 0, "RESOLVED": 0, "AMBIGUOUS": 0, "UNRESOLVED": 0}
    per_entry: list[dict[str, Any]] = []
    for path in paths:
        payload, counts = output_for_process(path, generated_at)
        out_path = args.out_dir / f"{payload['process_id']}.yaml"
        write_yaml(out_path, payload)
        totals["entries"] += 1
        for key in ("anchors", "RESOLVED", "AMBIGUOUS", "UNRESOLVED"):
            totals[key] += counts[key]
        per_entry.append({"process_id": payload["process_id"], **counts})

    print(
        f"entries={totals['entries']} anchors={totals['anchors']} "
        f"resolved={totals['RESOLVED']} ambiguous={totals['AMBIGUOUS']} unresolved={totals['UNRESOLVED']}"
    )
    for row in per_entry:
        print(
            f"{row['process_id']}: anchors={row['anchors']} "
            f"resolved={row['RESOLVED']} ambiguous={row['AMBIGUOUS']} unresolved={row['UNRESOLVED']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
