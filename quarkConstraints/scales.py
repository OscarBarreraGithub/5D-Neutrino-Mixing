"""Scale conventions for the quark-sector MFV module.

Two scales coexist and must stay orthogonal:

* ``DEFAULT_QUARK_FIT_SCALE_GEV`` — the renormalization scale at which
  PDG 2024 MS-bar quark mass targets are reported and the fitter scores
  residuals. Chosen as ``qcd.constants.M_TOP_MS = 162.5 GeV`` so all
  quarks use the PDG-consistent top MS-bar reference scale.

* ``DEFAULT_QUARK_TARGET_SCALE_GEV`` — the *Wilson-coefficient*
  matching/reference scale at which alpha_s is anchored for Delta-F=2
  running. Stays at 3 TeV. **Do not** repurpose this for mass targets.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from qcd.constants import M_TOP_MS

DEFAULT_QUARK_XI_KK = 1.0
GAUGE_KK_ROOT_NN = 2.448687135269161
KK_GLUON_PHYSICAL_MASS_CONVENTION_ID = (
    "kk_gluon_m1_physical.xi_kk_times_lambda_ir.v1"
)
KK_GLUON_LEGACY_GEOMETRIC_ALIAS_MASS_CONVENTION_ID = (
    "legacy_quark_geometric_alias.M_KK_equals_Lambda_IR.v1"
)
SPIN2_GRAVITON_KK_ROOT = 3.8317059702075125
DEFAULT_QUARK_BENCHMARK_XI_KK = DEFAULT_QUARK_XI_KK
# Wilson-coefficient reference scale (UNCHANGED at 3 TeV).
DEFAULT_QUARK_TARGET_SCALE_GEV = 3000.0
# Mass-target / fitter scoring scale (NEW: PDG 2024 -> mu = m_t(m_t)).
DEFAULT_QUARK_FIT_SCALE_GEV = M_TOP_MS
DEFAULT_QUARK_PAPER_H_RS_MAX = 0.3
DEFAULT_QUARK_BENCHMARK_H_RS_MAX = 1.0
DEFAULT_QUARK_MAX_MISALIGNMENT_STRESS = 6.0


class KKGluonMassConventionError(ValueError):
    """Raised when a row mixes physical and geometric KK-gluon masses."""


@dataclass(frozen=True)
class KKGluonMassConvention:
    """Explicit KK-gluon mass convention payload.

    ``M_KK`` in scan outputs is the physical first KK-gluon mass.  The
    geometric IR scale remains ``lambda_ir_gev`` and must satisfy
    ``m_kk_physical_gev == xi_kk * lambda_ir_gev`` at boundaries where the
    KK-gluon scale flows into constraints.
    """

    lambda_ir_gev: float
    m_kk_physical_gev: float
    xi_kk: float
    mass_convention_id: str = KK_GLUON_PHYSICAL_MASS_CONVENTION_ID

    def __post_init__(self) -> None:
        _validate_positive_finite("lambda_ir_gev", self.lambda_ir_gev)
        _validate_positive_finite("m_kk_physical_gev", self.m_kk_physical_gev)
        _validate_positive_finite("xi_kk", self.xi_kk)
        if not isinstance(self.mass_convention_id, str) or not self.mass_convention_id:
            raise KKGluonMassConventionError("mass_convention_id must be a non-empty string")
        if (
            self.mass_convention_id == KK_GLUON_PHYSICAL_MASS_CONVENTION_ID
            and math.isclose(
                float(self.xi_kk),
                DEFAULT_QUARK_XI_KK,
                rel_tol=0.0,
                abs_tol=1e-15,
            )
        ):
            raise KKGluonMassConventionError(
                "physical KK-gluon mass_convention_id cannot be paired with "
                "xi_kk=1.0; use the legacy geometric-alias convention id or pass "
                "the physical first-gauge-KK xi_KK"
            )
        if (
            self.mass_convention_id
            == KK_GLUON_LEGACY_GEOMETRIC_ALIAS_MASS_CONVENTION_ID
            and not math.isclose(
                float(self.xi_kk),
                DEFAULT_QUARK_XI_KK,
                rel_tol=0.0,
                abs_tol=1e-15,
            )
        ):
            raise KKGluonMassConventionError(
                "legacy geometric-alias mass_convention_id requires xi_kk=1.0"
            )
        expected = float(self.xi_kk * self.lambda_ir_gev)
        if not math.isclose(
            float(self.m_kk_physical_gev),
            expected,
            rel_tol=1e-12,
            abs_tol=1e-9,
        ):
            raise KKGluonMassConventionError(
                "KK-gluon mass convention mismatch: "
                f"m_kk_physical_gev={float(self.m_kk_physical_gev):.17g} "
                f"but xi_kk*lambda_ir_gev={expected:.17g} "
                f"(xi_kk={float(self.xi_kk):.17g}, "
                f"lambda_ir_gev={float(self.lambda_ir_gev):.17g})"
            )

    @property
    def M_KK(self) -> float:
        """Backward-compatible alias for the physical first KK-gluon mass."""

        return float(self.m_kk_physical_gev)

    def as_params(self) -> dict[str, float | str]:
        """Return the canonical scan-parameter payload."""

        return {
            "Lambda_IR": float(self.lambda_ir_gev),
            "lambda_ir_gev": float(self.lambda_ir_gev),
            "M_KK": float(self.m_kk_physical_gev),
            "m_kk_physical_gev": float(self.m_kk_physical_gev),
            "xi_KK": float(self.xi_kk),
            "mass_convention_id": self.mass_convention_id,
        }


def _validate_positive_finite(name: str, value: float) -> float:
    numeric = float(value)
    if not math.isfinite(numeric) or numeric <= 0.0:
        raise KKGluonMassConventionError(f"{name} must be a positive finite float")
    return numeric


def mass_convention_id_for_xi_kk(xi_kk: float) -> str:
    """Return a loud convention ID for a chosen ``xi_KK`` value."""

    if math.isclose(float(xi_kk), DEFAULT_QUARK_XI_KK, rel_tol=0.0, abs_tol=1e-15):
        return KK_GLUON_LEGACY_GEOMETRIC_ALIAS_MASS_CONVENTION_ID
    return KK_GLUON_PHYSICAL_MASS_CONVENTION_ID


def assert_kk_gluon_mass_convention(
    *,
    lambda_ir_gev: float,
    m_kk_physical_gev: float,
    xi_kk: float,
    mass_convention_id: str | None = None,
    context: str = "KK-gluon mass boundary",
) -> KKGluonMassConvention:
    """Validate and return the explicit KK-gluon mass convention payload."""

    try:
        return KKGluonMassConvention(
            lambda_ir_gev=float(lambda_ir_gev),
            m_kk_physical_gev=float(m_kk_physical_gev),
            xi_kk=float(xi_kk),
            mass_convention_id=(
                mass_convention_id
                if mass_convention_id is not None
                else mass_convention_id_for_xi_kk(float(xi_kk))
            ),
        )
    except KKGluonMassConventionError as exc:
        raise KKGluonMassConventionError(f"{context}: {exc}") from exc


def default_quark_m_kk_from_lambda_ir(
    Lambda_IR: float,
    xi_KK: float = DEFAULT_QUARK_XI_KK,
) -> float:
    """Map the geometric IR scale to the quark-sector KK convention.

    Parameters
    ----------
    Lambda_IR : float
        Geometric IR scale, ``Lambda_IR = 1 / z_v``.
    xi_KK : float, optional
        Explicit KK-mass convention. The low-level quark helpers keep the repo's
        bookkeeping default ``xi_KK = 1.0`` (``M_KK = Lambda_IR``), while the
        benchmark/validation helpers can opt into a physical first-KK
        convention such as ``DEFAULT_QUARK_BENCHMARK_XI_KK``.
    """
    if Lambda_IR <= 0:
        raise ValueError("Lambda_IR must be positive")
    if xi_KK <= 0:
        raise ValueError("xi_KK must be positive")
    return float(xi_KK * Lambda_IR)


def spin2_graviton_mass_from_lambda_ir(
    Lambda_IR: float,
    root: float = SPIN2_GRAVITON_KK_ROOT,
) -> float:
    """Return the first spin-2 KK-graviton mass from the geometric IR scale."""
    if Lambda_IR <= 0:
        raise ValueError("Lambda_IR must be positive")
    if root <= 0:
        raise ValueError("root must be positive")
    return float(root * Lambda_IR)
