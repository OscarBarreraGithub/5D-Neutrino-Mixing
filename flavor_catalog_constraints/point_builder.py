"""Construction paths for :class:`ParameterPoint`.

Three entry points, in increasing order of "how much the caller already
has computed":

``make_point(raw=None, **extras)``
    Hand-build a point from explicit pieces. This is the primary,
    always-usable path: tests and callers pass whatever extras they
    have as keyword arguments. Unknown extra keys raise immediately, so
    a typo is caught at construction rather than silently surfacing as a
    ``None`` inside a constraint.

``empty_point()``
    A point with no extras — for smoke tests of constraints that need no
    pre-computed input.

``build_from_quark_couplings(couplings, *, raw=None)``
    Real, working builder for the quark-sector constraints: it takes a
    :class:`quarkConstraints.couplings.QuarkMassBasisCouplings` (the
    object the fit machinery already produces) and packages it into the
    ``"quark_mass_basis_couplings"`` extra. This is deliberately NOT a
    ``NotImplementedError`` stub — the previous scaffold left its
    primary path unusable, which this fixes. Full scan-driver wiring
    (mapping a raw scan row -> couplings) is a separate later task, but
    the construction path itself is usable today.

The set of recognized extra keys lives in :data:`KNOWN_EXTRA_KEYS`; the
contract test asserts every key a constraint reads is declared here, and
:func:`make_point` enforces it at runtime.
"""

from __future__ import annotations

from typing import Any

from .base import ParameterPoint

__all__ = [
    "KNOWN_EXTRA_KEYS",
    "make_point",
    "empty_point",
    "build_from_quark_couplings",
]

# The single registry of valid ``ParameterPoint.extras`` keys. Adding a
# new physics input means: (1) add the key here with a one-line comment,
# (2) populate it in a builder, (3) read it via point.get_extra(key).
KNOWN_EXTRA_KEYS: frozenset[str] = frozenset(
    {
        "quark_mass_basis_couplings",  # quarkConstraints.couplings.QuarkMassBasisCouplings
        "kk_gluon_mass_gev",           # float, KK gluon mass (GeV)
        "kk_ew_mass_gev",              # float, KK electroweak gauge boson mass (GeV)
        "lepton_mass_basis_couplings", # future lepton-sector couplings object
    }
)


def make_point(raw: Any = None, **extras: Any) -> ParameterPoint:
    """Hand-build a :class:`ParameterPoint` from explicit extras.

    Every keyword must be a key in :data:`KNOWN_EXTRA_KEYS`; an unknown
    key raises :class:`KeyError` so typos fail loudly at construction.
    """
    unknown = set(extras) - KNOWN_EXTRA_KEYS
    if unknown:
        raise KeyError(
            f"unknown ParameterPoint extra(s): {sorted(unknown)}; "
            f"declared keys are {sorted(KNOWN_EXTRA_KEYS)}"
        )
    return ParameterPoint(raw=raw, extras=dict(extras))


def empty_point() -> ParameterPoint:
    """Return a :class:`ParameterPoint` with no extras."""
    return ParameterPoint(raw=None, extras={})


def build_from_quark_couplings(couplings: Any, *, raw: Any = None) -> ParameterPoint:
    """Build a point carrying quark mass-basis couplings.

    ``couplings`` is expected to be a
    :class:`quarkConstraints.couplings.QuarkMassBasisCouplings` (the
    object the existing fit code produces). It is stored under the
    ``"quark_mass_basis_couplings"`` extra and exposes a usable
    ``kk_gluon_mass_gev`` when the couplings object carries an ``M_KK``.
    """
    extras: dict[str, Any] = {"quark_mass_basis_couplings": couplings}
    m_kk = getattr(couplings, "M_KK", None)
    if m_kk is not None:
        extras["kk_gluon_mass_gev"] = float(m_kk)
    return make_point(raw=raw, **extras)
