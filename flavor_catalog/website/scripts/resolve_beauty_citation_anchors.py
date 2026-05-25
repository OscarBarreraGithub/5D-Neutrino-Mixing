#!/usr/bin/env python3
"""Resolve citation anchors for the beauty catalog batch.

The script walks pdg_or_equivalent blocks, builds conservative anchor strings,
searches local text snapshots, and writes website citation-anchor YAML files.
"""

from __future__ import annotations

import argparse
import html
import math
import re
import sys
import unicodedata
from dataclasses import dataclass
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any

import yaml


ROOT = Path(__file__).resolve().parents[3]
PRIMARY_DIR = ROOT / "flavor_catalog/processes/beauty"
SECONDARY_DIR = ROOT / "flavor_catalog/processes/secondary/beauty"
OUT_DIR = ROOT / "flavor_catalog/website/_data/citation_anchors"
SECONDARY_IDS = ("B007", "B008", "B013", "B014")

STATUS_RESOLVED = "RESOLVED"
STATUS_AMBIGUOUS = "AMBIGUOUS"
STATUS_UNRESOLVED = "UNRESOLVED"

SKIP_RECURSE_KEYS = {
    "snapshot_path",
    "supporting_snapshot_path",
    "logfile_snapshot_path",
}
META_KEYS = {
    "source_url",
    "access_date",
    "sha256",
    "sha256_of_text_snapshot",
    "sha256_of_local_snapshot",
}
DISPLAY_KEYS = ("display_value", "total_display_value", "display")
NESTED_VALUE_KEYS = (
    "observables",
    "values",
    "inputs",
    "average_inputs",
    "predictions",
    "dataset_numerical_claims",
    "reference_scales",
)
GENERIC_TOKENS = {
    "the",
    "and",
    "for",
    "from",
    "with",
    "into",
    "value",
    "values",
    "average",
    "averages",
    "input",
    "inputs",
    "source",
    "branching",
    "br",
    "fraction",
    "dimensionless",
    "units",
    "unit",
    "observable",
    "observables",
    "canonical",
    "experimental",
    "prediction",
    "standard",
    "model",
    "current",
    "combined",
    "equivalent",
    "display",
    "confidence",
    "level",
    "upper",
    "limit",
    "lower",
    "cl",
    "pdg",
    "hflav",
    "sm",
    "stat",
    "syst",
    "statistical",
    "systematic",
    "external",
    "total",
}


class LiteralStr(str):
    pass


def literal_representer(dumper: yaml.Dumper, data: LiteralStr) -> yaml.ScalarNode:
    return dumper.represent_scalar("tag:yaml.org,2002:str", data, style="|")


yaml.SafeDumper.add_representer(LiteralStr, literal_representer)


@dataclass
class AnchorCandidate:
    field: str
    text: str
    value: Any | None = None
    uncertainties: tuple[Any, ...] = ()
    tokens_source: str = ""


@dataclass
class SourceBlock:
    block_key: str
    block: dict[str, Any]
    inherited: dict[str, Any]
    anchors: list[AnchorCandidate]


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"{path} did not parse to a mapping")
    return data


def dump_yaml(path: Path, data: dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(
            data,
            handle,
            explicit_start=True,
            sort_keys=False,
            allow_unicode=True,
            width=120,
        )


def batch_paths(ids: set[str] | None = None) -> list[Path]:
    primary = sorted(PRIMARY_DIR.glob("*.yaml"))
    secondary = [SECONDARY_DIR / f"{process_id}.yaml" for process_id in SECONDARY_IDS]
    paths = primary + secondary
    if ids is not None:
        paths = [path for path in paths if path.stem in ids]
    return paths


def path_to_string(path: list[str | int]) -> str:
    parts: list[str] = []
    for item in path:
        if isinstance(item, int):
            if parts:
                parts[-1] = f"{parts[-1]}[{item}]"
            else:
                parts.append(f"[{item}]")
        else:
            parts.append(str(item))
    return ".".join(parts)


def block_key(path: list[str | int], block: dict[str, Any]) -> str:
    base = path_to_string(path)
    ident = block.get("id") or block.get("value_id")
    if ident:
        return f"{base}:{ident}"
    return base


def scalar(value: Any) -> bool:
    return isinstance(value, (str, int, float, bool)) and value is not None


def text_value(value: Any) -> str | None:
    if value is None:
        return None
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float, str)):
        text = str(value).strip()
        return text or None
    return None


def inherited_with(node: dict[str, Any], inherited: dict[str, Any]) -> dict[str, Any]:
    current = dict(inherited)
    for key in META_KEYS:
        if key in node and node[key] is not None:
            current[key] = node[key]
    return current


def walk_snapshot_blocks(
    node: Any,
    path: list[str | int] | None = None,
    inherited: dict[str, Any] | None = None,
) -> list[tuple[list[str | int], dict[str, Any], dict[str, Any]]]:
    if path is None:
        path = []
    if inherited is None:
        inherited = {}
    found: list[tuple[list[str | int], dict[str, Any], dict[str, Any]]] = []
    if isinstance(node, dict):
        current = inherited_with(node, inherited)
        if node.get("snapshot_path"):
            found.append((path, node, current))
        for key, value in node.items():
            if key in SKIP_RECURSE_KEYS:
                continue
            if isinstance(value, (dict, list)):
                found.extend(walk_snapshot_blocks(value, path + [key], current))
    elif isinstance(node, list):
        for index, value in enumerate(node):
            found.extend(walk_snapshot_blocks(value, path + [index], inherited))
    return found


def normalize_sha(value: Any) -> str | None:
    text = text_value(value)
    if text is None:
        return None
    if text.startswith("sha256:"):
        return text
    if re.fullmatch(r"[0-9a-fA-F]{64}", text):
        return f"sha256:{text.lower()}"
    return text


def source_sha(block: dict[str, Any], inherited: dict[str, Any]) -> str | None:
    for key in ("sha256_of_local_snapshot", "sha256_of_text_snapshot", "sha256"):
        if key in block:
            return normalize_sha(block[key])
    for key in ("sha256_of_local_snapshot", "sha256_of_text_snapshot", "sha256"):
        if key in inherited:
            return normalize_sha(inherited[key])
    return None


def source_url(block: dict[str, Any], inherited: dict[str, Any]) -> str | None:
    return text_value(block.get("source_url")) or text_value(inherited.get("source_url"))


def access_date(block: dict[str, Any], inherited: dict[str, Any]) -> str | None:
    return text_value(block.get("access_date")) or text_value(inherited.get("access_date"))


def decimal_or_none(value: Any) -> Decimal | None:
    if value is None:
        return None
    if isinstance(value, bool):
        return None
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "null"}:
        return None
    text = text.replace("−", "-").replace("+/-", "")
    match = re.search(r"[-+]?\d+(?:\.\d+)?(?:e[-+]?\d+)?", text, flags=re.I)
    if not match:
        return None
    try:
        return Decimal(match.group(0))
    except InvalidOperation:
        return None


def decimal_exponent(value: Any) -> int | None:
    dec = decimal_or_none(value)
    if dec is None or dec == 0:
        return None
    return int(math.floor(math.log10(abs(float(dec)))))


def scaled_text(value: Any, exponent: int) -> str | None:
    dec = decimal_or_none(value)
    if dec is None:
        return None
    scaled = dec / (Decimal(10) ** exponent)
    return pretty_decimal(scaled)


def pretty_decimal(value: Decimal) -> str:
    text = format(value.normalize(), "f")
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    if text == "-0":
        text = "0"
    return text


def is_small_or_large(value: Any) -> bool:
    dec = decimal_or_none(value)
    if dec is None or dec == 0:
        return False
    abs_dec = abs(dec)
    return abs_dec < Decimal("0.001") or abs_dec >= Decimal("10000")


def field_label(block: dict[str, Any]) -> str:
    for key in ("observable", "name", "label", "type", "source_key", "id", "value_id"):
        label = text_value(block.get(key))
        if label:
            return label
    return ""


def uncertainty_values(block: dict[str, Any]) -> list[Any]:
    keys = (
        "uncertainty",
        "uncertainties",
        "stat_uncertainty",
        "syst_uncertainty",
        "statistical_uncertainty",
        "systematic_uncertainty",
        "uncertainty_stat",
        "uncertainty_syst",
        "uncertainty_positive",
        "uncertainty_negative",
        "uncertainty_plus",
        "uncertainty_minus",
        "upper_uncertainty",
        "lower_uncertainty",
        "upper_uncertainty_degrees",
        "lower_uncertainty_degrees",
        "uncertainty_stat_plus",
        "uncertainty_stat_minus",
        "uncertainty_syst_plus",
        "uncertainty_syst_minus",
        "stat_only_uncertainty",
    )
    values: list[Any] = []
    for key in keys:
        if key in block and block[key] not in (None, "not quoted"):
            values.append(block[key])
    for key, value in block.items():
        if key.endswith("_uncertainty") and value not in (None, "not quoted") and value not in values:
            values.append(value)
    return values


def value_key(block: dict[str, Any]) -> str | None:
    for key in (
        "value",
        "upper_limit",
        "upper_limit_90cl",
        "upper_limit_95cl",
        "upper_limit_95cl_run2",
        "upper_limit_95cl_run1_run2_combined",
        "lower_limit",
        "value_degrees",
        "bound_abs_delta_phi_d_degrees",
        "phi_s_value",
        "lambda_value",
        "a_u_value",
        "a_v_value",
        "s_value",
        "a_value",
        "c_cp_value",
        "sin2beta_eff_value",
        "signal_events",
        "n_bb_millions",
        "integrated_luminosity_fb_inverse",
        "approximate_signal_decays",
        "b0_to_phiphi_limit_value",
    ):
        if key in block and block[key] is not None:
            return key
    for key, value in block.items():
        if key.endswith("_value") and value is not None:
            return key
    return None


def format_value_uncertainty(block: dict[str, Any], key: str) -> tuple[str, tuple[Any, ...]]:
    value = block[key]
    uncertainties = uncertainty_values(block)
    units = text_value(block.get("units") or block.get(f"{key.rsplit('_', 1)[0]}_units"))
    cl = text_value(block.get("cl") or block.get("confidence_level"))
    relation = "<" if "limit" in key or block.get("limit_type") == "upper" else "="
    exponent = decimal_exponent(value) if is_small_or_large(value) else None
    value_text = str(value)

    numeric_uncertainties = [unc for unc in uncertainties if decimal_or_none(unc) is not None]
    if exponent is not None:
        scaled_value = scaled_text(value, exponent)
        if scaled_value is not None:
            value_text = f"{scaled_value}e{exponent}"
            if numeric_uncertainties:
                scaled_unc = scaled_text(numeric_uncertainties[0], exponent)
                if scaled_unc is not None and relation == "=":
                    value_text = f"({scaled_value} +- {scaled_unc})e{exponent}"
    elif numeric_uncertainties and relation == "=":
        value_text = f"{value} +- {numeric_uncertainties[0]}"

    if uncertainties and not numeric_uncertainties and relation == "=":
        value_text = f"{value} {uncertainties[0]}"

    parts = [field_label(block), relation, value_text]
    if units and units not in {"dimensionless", "branching fraction"}:
        parts.append(units)
    if cl:
        if cl in {"0.9", "0.90"}:
            parts.append("90% CL")
        elif cl in {"0.95"}:
            parts.append("95% CL")
        elif "cl" in cl.lower() or "%" in cl:
            parts.append(cl)
    return " ".join(part for part in parts if part), tuple(uncertainties)


def collect_related_anchors(block: dict[str, Any], prefix: str = "") -> list[AnchorCandidate]:
    anchors: list[AnchorCandidate] = []
    for key, value in block.items():
        if not key.startswith("related_"):
            continue
        field = f"{prefix}{key}" if prefix else key
        if scalar(value):
            anchors.append(AnchorCandidate(field, str(value), tokens_source=str(value)))
        elif isinstance(value, list):
            for index, item in enumerate(value):
                if scalar(item):
                    anchors.append(
                        AnchorCandidate(f"{field}[{index}]", str(item), tokens_source=str(item))
                    )
        elif isinstance(value, dict):
            for subkey, item in value.items():
                if scalar(item):
                    anchors.append(
                        AnchorCandidate(f"{field}.{subkey}", str(item), tokens_source=str(item))
                    )
    return anchors


def anchors_for_single_value(block: dict[str, Any], prefix: str = "") -> list[AnchorCandidate]:
    for key in DISPLAY_KEYS:
        value = text_value(block.get(key))
        if value:
            field = f"{prefix}{key}" if prefix else key
            return [AnchorCandidate(field, value, tokens_source=value)]

    related = collect_related_anchors(block, prefix)
    if related:
        return related

    key = value_key(block)
    if key is None:
        return []
    anchor_text, uncertainties = format_value_uncertainty(block, key)
    if not anchor_text:
        return []
    field = f"{prefix}fallback_value_uncertainty" if prefix else "fallback_value_uncertainty"
    return [
        AnchorCandidate(
            field,
            anchor_text,
            value=block.get(key),
            uncertainties=uncertainties,
            tokens_source=" ".join(
                str(part)
                for part in (
                    block.get("observable"),
                    block.get("name"),
                    block.get("label"),
                    block.get("q2_region"),
                )
                if part is not None
            ),
        )
    ]


def anchors_for_nested_item(value: Any, prefix: str) -> list[AnchorCandidate]:
    anchors: list[AnchorCandidate] = []
    if isinstance(value, dict):
        if value.get("snapshot_path"):
            return []
        anchors.extend(anchors_for_single_value(value, prefix))
        if anchors:
            return anchors
        for key in NESTED_VALUE_KEYS:
            child = value.get(key)
            if isinstance(child, list):
                for index, item in enumerate(child):
                    anchors.extend(anchors_for_nested_item(item, f"{prefix}{key}[{index}]."))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            anchors.extend(anchors_for_nested_item(item, f"{prefix}[{index}]."))
    return anchors


def anchors_for_block(block: dict[str, Any]) -> list[AnchorCandidate]:
    direct = anchors_for_single_value(block)
    if direct:
        return direct
    anchors: list[AnchorCandidate] = []
    for key in NESTED_VALUE_KEYS:
        value = block.get(key)
        if isinstance(value, list):
            for index, item in enumerate(value):
                anchors.extend(anchors_for_nested_item(item, f"{key}[{index}]."))
        elif isinstance(value, dict):
            anchors.extend(anchors_for_nested_item(value, f"{key}."))
    return anchors


def normalize_text(text: str) -> str:
    text = html.unescape(text)
    text = re.sub(r"<[^>]+>", " ", text)
    text = unicodedata.normalize("NFKC", text)
    replacements = {
        "±": "+-",
        "∓": "-+",
        "−": "-",
        "–": "-",
        "—": "-",
        "×": "x",
        "⋅": "x",
        "·": "x",
        "→": "->",
        "⟶": "->",
        "⇒": "=>",
        "ν": "nu",
        "μ": "mu",
        "τ": "tau",
        "π": "pi",
        "ρ": "rho",
        "ω": "omega",
        "φ": "phi",
        "ϕ": "phi",
        "Φ": "phi",
        "γ": "gamma",
        "Γ": "gamma",
        "Δ": "delta",
        "δ": "delta",
        "β": "beta",
        "χ": "chi",
        "λ": "lambda",
        "Λ": "lambda",
        "η": "eta",
        "ψ": "psi",
        "Υ": "upsilon",
        "√": "sqrt",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    text = text.lower()
    text = text.replace("\\pm", "+-")
    text = text.replace("+/-", "+-")
    text = text.replace("\\times", " x ")
    text = re.sub(r"\\(?:mathrm|text|rm|mathcal|bar|overline)\s*\{([^{}]*)\}", r"\1", text)
    text = re.sub(r"\\(?:to|rightarrow)", " -> ", text)
    text = re.sub(r"\\[a-zA-Z]+", " ", text)
    text = text.replace("{", "").replace("}", "")
    text = text.replace("^", "")
    text = text.replace("_", " ")
    text = re.sub(r"\)\s*(?:x|\*)\s*10\s*([+-]?\d+)", r")e\1", text)
    text = re.sub(
        r"(?P<num>[-+]?\d+(?:\.\d+)?)\s*(?:x|\*)\s*10\s*(?P<exp>[+-]?\d+)",
        r"\g<num>e\g<exp>",
        text,
    )
    text = re.sub(r"\be([+-])0+(\d+)", r"e\1\2", text)
    text = re.sub(r"(?<=\d)e\+?(-?\d+)", lambda m: "e" + m.group(1), text)
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def token_set(text: str) -> set[str]:
    normalized = normalize_text(text)
    tokens = set(re.findall(r"[a-z][a-z0-9]+|[a-z]", normalized))
    return {token for token in tokens if token not in GENERIC_TOKENS and len(token) > 1}


def decimal_close(left: Decimal, right: Decimal) -> bool:
    if left == right:
        return True
    scale = max(abs(left), abs(right), Decimal("1"))
    return abs(left - right) <= scale * Decimal("0.0000001")


def normalized_numbers(text: str) -> list[Decimal]:
    normalized = normalize_text(text)
    values: list[Decimal] = []
    for match in re.finditer(r"[-+]?\d+(?:\.\d+)?(?:e[-+]?\d+)?", normalized, flags=re.I):
        token = match.group(0)
        try:
            values.append(Decimal(token))
        except InvalidOperation:
            continue
    return values


def number_variants(value: Any, scale_exponent: int | None = None) -> set[str]:
    dec = decimal_or_none(value)
    if dec is None:
        return set()
    variants: set[str] = set()
    plain = pretty_decimal(dec)
    variants.add(normalize_text(plain))
    variants.add(normalize_text(str(value)))
    if dec != 0:
        exponent = decimal_exponent(dec)
        if exponent is not None:
            mantissa = scaled_text(dec, exponent)
            if mantissa is not None:
                variants.add(normalize_text(f"{mantissa}e{exponent}"))
                variants.add(normalize_text(f"{mantissa} x 10^{exponent}"))
        if scale_exponent is not None:
            scaled = scaled_text(dec, scale_exponent)
            if scaled is not None:
                variants.add(normalize_text(scaled))
                variants.add(normalize_text(f"{scaled}e{scale_exponent}"))
                variants.add(normalize_text(f"{scaled} x 10^{scale_exponent}"))
    return {variant for variant in variants if variant}


def number_present(value: Any, text: str, scale_exponent: int | None = None) -> bool:
    dec = decimal_or_none(value)
    if dec is None:
        return False
    normalized = normalize_text(text)
    for variant in number_variants(value, scale_exponent):
        if variant and variant in normalized:
            return True
    for candidate in normalized_numbers(text):
        if decimal_close(dec, candidate):
            return True
    if scale_exponent is not None:
        scaled = scaled_text(dec, scale_exponent)
        if scaled and normalize_text(scaled) in normalized:
            exp_token = f"e{scale_exponent}"
            if exp_token in normalized or f"10{scale_exponent}" in normalized:
                return True
    return False


def exact_anchor_present(anchor: str, text: str) -> bool:
    normalized_anchor = normalize_text(anchor)
    normalized_text = normalize_text(text)
    if not normalized_anchor:
        return False
    if normalized_anchor in normalized_text:
        return True
    anchor_no_space = re.sub(r"[^a-z0-9.+-]+", "", normalized_anchor)
    text_no_space = re.sub(r"[^a-z0-9.+-]+", "", normalized_text)
    return bool(anchor_no_space and anchor_no_space in text_no_space)


def anchor_numeric_values(anchor: AnchorCandidate) -> list[Any]:
    values: list[Any] = []
    if anchor.value is not None:
        values.append(anchor.value)
    for uncertainty in anchor.uncertainties:
        if isinstance(uncertainty, str):
            matches = re.findall(r"[-+]?\d+(?:\.\d+)?(?:e[-+]?\d+)?", normalize_text(uncertainty), flags=re.I)
            values.extend(match.lstrip("+-") for match in matches)
            if not matches:
                values.append(uncertainty)
        else:
            values.append(uncertainty)
    if values:
        return values
    found: list[str] = []
    for match in re.finditer(r"[-+]?\d+(?:\.\d+)?(?:e[-+]?\d+)?", normalize_text(anchor.text), flags=re.I):
        found.append(match.group(0))
    return found[:4]


def candidate_matches(anchor: AnchorCandidate, line: str, window: str) -> bool:
    if exact_anchor_present(anchor.text, line):
        return True

    numeric_values = [value for value in anchor_numeric_values(anchor) if decimal_or_none(value) is not None]
    if not numeric_values:
        return False

    scale_exponent = decimal_exponent(numeric_values[0]) if is_small_or_large(numeric_values[0]) else None
    primary_present, matched_scale_exponent = number_present_near(
        numeric_values[0], line, window, scale_exponent
    )
    if not primary_present:
        return False

    secondary_values = numeric_values[1:3]
    if secondary_values and not any(
        number_present(value, window, matched_scale_exponent)
        or compact_uncertainty_present(numeric_values[0], value, window, matched_scale_exponent)
        for value in secondary_values
    ):
        return False
    if secondary_values:
        return True
    if "<" in anchor.text and "<" in line:
        return True

    tokens = token_set(anchor.tokens_source or anchor.text)
    if not tokens:
        return True
    window_tokens = token_set(window)
    overlap = tokens & window_tokens
    if len(tokens) <= 2:
        return bool(overlap)
    return len(overlap) >= min(2, len(tokens))


def candidate_scale_exponents(value: Any, preferred: int | None) -> list[int | None]:
    candidates: list[int | None] = [preferred]
    if preferred is not None:
        candidates.extend(range(preferred - 2, preferred + 3))
    dec = decimal_or_none(value)
    if dec is not None and abs(dec) < Decimal("1"):
        candidates.extend(range(-2, -13, -1))
    deduped: list[int | None] = []
    for candidate in candidates:
        if candidate not in deduped:
            deduped.append(candidate)
    return deduped


def number_present_near(
    value: Any,
    line: str,
    window: str,
    preferred_scale_exponent: int | None,
) -> tuple[bool, int | None]:
    if number_present(value, line, preferred_scale_exponent):
        return True, preferred_scale_exponent
    normalized_line = normalize_text(line)
    normalized_window = normalize_text(window)
    for exponent in candidate_scale_exponents(value, preferred_scale_exponent):
        if exponent is None:
            continue
        scaled = scaled_text(value, exponent)
        if not scaled:
            continue
        normalized_scaled = normalize_text(scaled)
        exp_tokens = {f"e{exponent}", f"10{exponent}"}
        if normalized_scaled in normalized_line and (
            exp_tokens & set(normalized_window.split()) or any(token in normalized_window for token in exp_tokens)
        ):
            return True, exponent
        dec = decimal_or_none(value)
        if dec is not None and abs(dec) < Decimal("0.001") and normalized_scaled in normalized_line:
            return True, exponent
    return False, preferred_scale_exponent


def compact_uncertainty_present(
    primary_value: Any,
    uncertainty: Any,
    text: str,
    scale_exponent: int | None,
) -> bool:
    if scale_exponent is None:
        return False
    primary = scaled_text(primary_value, scale_exponent)
    unc = scaled_text(uncertainty, scale_exponent)
    if not primary or not unc:
        return False
    normalized = normalize_text(text)
    unc_digits = normalize_text(unc).replace("0.", "").replace(".", "")
    pattern = re.escape(normalize_text(primary)) + r"\s*\(\s*" + re.escape(unc_digits) + r"\s*\)"
    return bool(unc_digits and re.search(pattern, normalized))


def line_window(lines: list[str], index: int, radius: int = 6) -> str:
    start = max(0, index - radius)
    end = min(len(lines), index + radius + 1)
    return "\n".join(lines[start:end])


def context_block(lines: list[str], index: int) -> LiteralStr:
    start = max(0, index - 3)
    end = min(len(lines), index + 4)
    rendered: list[str] = []
    for line_no in range(start + 1, end + 1):
        line = lines[line_no - 1].rstrip()
        rendered.append(f"L{line_no}: {line}" if line else f"L{line_no}:")
    context = "\n".join(rendered)
    return LiteralStr(context + "\n")


def resolve_anchor(anchor: AnchorCandidate, snapshot_path: Path) -> tuple[str, list[dict[str, Any]]]:
    if not snapshot_path.exists():
        return STATUS_UNRESOLVED, []
    lines = snapshot_path.read_text(encoding="utf-8", errors="replace").splitlines()
    matching_indices: list[int] = []
    for index, line in enumerate(lines):
        window = line_window(lines, index)
        if candidate_matches(anchor, line, window):
            matching_indices.append(index)
    deduped: list[int] = []
    for index in matching_indices:
        if index not in deduped:
            deduped.append(index)
    matches = [
        {"line_number": index + 1, "context": context_block(lines, index)}
        for index in deduped
    ]
    if len(matches) == 1:
        return STATUS_RESOLVED, matches
    if len(matches) > 1:
        return STATUS_AMBIGUOUS, matches
    return STATUS_UNRESOLVED, []


def build_source_blocks(pdg: Any) -> list[SourceBlock]:
    sources: list[SourceBlock] = []
    for path, block, inherited in walk_snapshot_blocks(pdg):
        anchors = anchors_for_block(block)
        if not anchors:
            continue
        sources.append(SourceBlock(block_key(path, block), block, inherited, anchors))
    return sources


def output_for_process(process_path: Path, generated_at: str) -> tuple[dict[str, Any], dict[str, int]]:
    data = load_yaml(process_path)
    process_id = process_path.stem
    pdg = data.get("pdg_or_equivalent")
    source_blocks = build_source_blocks(pdg)
    output: dict[str, Any] = {
        "process_id": process_id,
        "generated_at": generated_at,
        "sources": [],
    }
    counts = {STATUS_RESOLVED: 0, STATUS_AMBIGUOUS: 0, STATUS_UNRESOLVED: 0}
    for source in source_blocks:
        snapshot_rel = text_value(source.block.get("snapshot_path"))
        if not snapshot_rel:
            continue
        snapshot_path = ROOT / snapshot_rel
        source_entry: dict[str, Any] = {
            "block_key": source.block_key,
            "source_url": source_url(source.block, source.inherited),
            "snapshot_path": snapshot_rel,
            "sha256": source_sha(source.block, source.inherited),
            "access_date": access_date(source.block, source.inherited),
            "anchors": [],
        }
        for anchor in source.anchors:
            status, matches = resolve_anchor(anchor, snapshot_path)
            counts[status] += 1
            source_entry["anchors"].append(
                {
                    "anchor_field": anchor.field,
                    "anchor_string": anchor.text,
                    "status": status,
                    "matches": matches,
                }
            )
        output["sources"].append(source_entry)
    return output, counts


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ids", nargs="*", help="Optional process IDs to process")
    args = parser.parse_args(argv)

    ids = set(args.ids) if args.ids else None
    paths = batch_paths(ids)
    if ids is None and len(paths) != 26:
        raise SystemExit(f"Expected 26 beauty entries, found {len(paths)}")
    missing = [path for path in paths if not path.exists()]
    if missing:
        raise SystemExit("Missing process files: " + ", ".join(str(path) for path in missing))

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    generated_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")
    total_counts = {STATUS_RESOLVED: 0, STATUS_AMBIGUOUS: 0, STATUS_UNRESOLVED: 0}
    per_entry: dict[str, dict[str, int]] = {}

    for process_path in paths:
        output, counts = output_for_process(process_path, generated_at)
        process_id = process_path.stem
        dump_yaml(OUT_DIR / f"{process_id}.yaml", output)
        per_entry[process_id] = counts
        for status, value in counts.items():
            total_counts[status] += value

    total = sum(total_counts.values())
    resolved = total_counts[STATUS_RESOLVED]
    percent = 100.0 * resolved / total if total else 0.0
    print(f"entries {len(paths)}/26")
    print(
        "total anchors "
        f"{total}; RESOLVED {resolved} ({percent:.1f}%); "
        f"AMBIGUOUS {total_counts[STATUS_AMBIGUOUS]}; "
        f"UNRESOLVED {total_counts[STATUS_UNRESOLVED]}"
    )
    for process_id in sorted(per_entry):
        counts = per_entry[process_id]
        total_entry = sum(counts.values())
        unresolved = counts[STATUS_UNRESOLVED]
        flag = " FLAG >50% UNRESOLVED" if total_entry and unresolved / total_entry > 0.5 else ""
        print(
            f"{process_id}: {counts[STATUS_RESOLVED]} RES, "
            f"{counts[STATUS_AMBIGUOUS]} AMB, {counts[STATUS_UNRESOLVED]} UNR{flag}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
