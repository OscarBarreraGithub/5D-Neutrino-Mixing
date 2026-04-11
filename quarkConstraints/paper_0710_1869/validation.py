"""Validation helpers for the dedicated paper-mode 0710.1869 package."""

from __future__ import annotations

import ast
import math
from pathlib import Path
from typing import Iterable


def require_known_schema_id(name: str, value: str, *, expected: str) -> str:
    """Require an exact, versioned schema identifier."""
    if value != expected:
        raise ValueError(f"{name} must be exactly {expected!r}")
    return value


def require_nonempty_identifier(name: str, value: str) -> str:
    """Require a non-empty string identifier."""
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty string")
    return value.strip()


def require_member(name: str, value: str, allowed: Iterable[str]) -> str:
    """Require membership in a closed set of identifiers."""
    allowed_values = tuple(allowed)
    if value not in allowed_values:
        raise ValueError(f"{name} must be one of {allowed_values!r}")
    return value


def require_positive_finite(name: str, value: float) -> float:
    """Normalize and validate a positive finite float."""
    numeric = float(value)
    if not math.isfinite(numeric) or numeric <= 0.0:
        raise ValueError(f"{name} must be a positive finite float")
    return numeric


def normalize_optional_positive_finite(name: str, value: float | None) -> float | None:
    """Normalize an optional positive finite float."""
    if value is None:
        return None
    return require_positive_finite(name, value)


def module_has_forbidden_import(
    module_path: str | Path, forbidden_modules: Iterable[str]
) -> bool:
    """Return ``True`` if a Python module imports any forbidden module path."""
    path = Path(module_path)
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    forbidden = {item for item in forbidden_modules}

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                if _matches_forbidden_module(alias.name, forbidden):
                    return True
        if isinstance(node, ast.ImportFrom):
            imported_from = node.module or ""
            if imported_from and _matches_forbidden_module(imported_from, forbidden):
                return True
            for alias in node.names:
                qualified = f"{imported_from}.{alias.name}" if imported_from else alias.name
                if _matches_forbidden_module(qualified, forbidden):
                    return True

    return False


def _matches_forbidden_module(candidate: str, forbidden_modules: set[str]) -> bool:
    return candidate in forbidden_modules or any(
        candidate.startswith(f"{name}.") for name in forbidden_modules
    )
