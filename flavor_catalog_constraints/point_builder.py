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

``build_from_rs_ew_inputs(quark_fit_result, *, Lambda_IR, ...)``
    Electroweak builder. It constructs the RS-EW spectrum, mass-basis
    light-Z shifts, neutral semileptonic contacts, and Wilson bundles as
    extras. It deliberately does not rewire any constraint consumer.

The set of recognized extra keys lives in :data:`KNOWN_EXTRA_KEYS`; the
contract test asserts every key a constraint reads is declared here, and
:func:`make_point` enforces it at runtime.
"""

from __future__ import annotations

from typing import Any, Callable

from .base import ParameterPoint

__all__ = [
    "KNOWN_EXTRA_KEYS",
    "make_point",
    "empty_point",
    "build_from_quark_couplings",
    "build_from_rs_ew_inputs",
]

# The single registry of valid ``ParameterPoint.extras`` keys. Adding a
# new physics input means: (1) add the key here with a one-line comment,
# (2) populate it in a builder, (3) read it via point.get_extra(key).
KNOWN_EXTRA_KEYS: frozenset[str] = frozenset(
    {
        "quark_mass_basis_couplings",  # quarkConstraints.couplings.QuarkMassBasisCouplings
        "kk_gluon_mass_gev",           # float, KK gluon mass (GeV)
        "kk_ew_mass_gev",              # float, KK electroweak gauge boson mass (GeV)
        "lepton_mass_basis_couplings", # quarkConstraints.rs_ew_couplings.RSLeptonMassBasisCouplings
        "rs_ew_spectrum",              # quarkConstraints.rs_ew_spectrum.RSEWSpectrum
        "rs_ew_couplings",             # quarkConstraints.rs_ew_couplings.RSEWMassBasisCouplings
        "rs_charged_current",          # quarkConstraints.rs_charged_current.RSChargedCurrentCouplings
        "rs_higgs_yukawas",            # quarkConstraints.rs_higgs_yukawas.RSHiggsYukawaCouplings
        "rs_semileptonic_wilsons",     # quarkConstraints.rs_semileptonic_wilsons.RSSemileptonicWilsonBundle
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


def build_from_rs_ew_inputs(
    quark_fit_result: Any,
    *,
    Lambda_IR: float,
    k: float | None = None,
    n_gauge_modes: int | None = None,
    ew_model: str = "minimal_rs",
    raw: Any = None,
    spectrum: Any | None = None,
    a_of_c: Callable[[float], float] | None = None,
    omega_of_c: Callable[[float], Any] | None = None,
    rs_ew_cache: Any | None = None,
    **kwargs: Any,
) -> ParameterPoint:
    """Build a point carrying RS electroweak extras.

    The fit result must expose ``bulk_state.c_Q``, ``bulk_state.c_u``,
    ``bulk_state.c_d``, and the four mass-basis rotations
    ``U_L_u``, ``U_L_d``, ``U_R_u``, ``U_R_d``.
    """

    from quarkConstraints.rs_ew_spectrum import DEFAULT_N_GAUGE_MODES
    from warpConfig.baseParams import MPL

    from .rs_ew_builder import build_rs_ew_extras

    extras = build_rs_ew_extras(
        quark_fit_result,
        Lambda_IR=float(Lambda_IR),
        k=float(MPL if k is None else k),
        n_gauge_modes=int(DEFAULT_N_GAUGE_MODES if n_gauge_modes is None else n_gauge_modes),
        ew_model=ew_model,
        spectrum=spectrum,
        a_of_c=a_of_c,
        omega_of_c=omega_of_c,
        rs_ew_cache=rs_ew_cache,
        **kwargs,
    )
    return make_point(raw=raw, **extras)
