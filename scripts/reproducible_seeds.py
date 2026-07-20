"""Stable seed derivation for scan and research scripts.

Python's built-in ``hash()`` is intentionally salted per process.  Seeds that
must reproduce across workers, login sessions, and Python installations should
derive label offsets through this module instead.
"""

from __future__ import annotations

import hashlib


def stable_seed_offset(
    label: str,
    *,
    modulus: int,
    namespace: str,
) -> int:
    """Return a deterministic integer offset in ``range(modulus)``.

    ``namespace`` separates independent random streams that happen to use the
    same label.  The frozen SHA-256/first-eight-bytes convention is deliberately
    simple and language-independent.
    """

    if not isinstance(label, str) or not label:
        raise ValueError("label must be a non-empty string")
    if not isinstance(namespace, str) or not namespace:
        raise ValueError("namespace must be a non-empty string")
    if isinstance(modulus, bool) or not isinstance(modulus, int) or modulus <= 0:
        raise ValueError("modulus must be a positive integer")
    payload = f"{namespace}\0{label}".encode("utf-8")
    prefix = hashlib.sha256(payload).digest()[:8]
    return int.from_bytes(prefix, byteorder="big", signed=False) % modulus


__all__ = ["stable_seed_offset"]
