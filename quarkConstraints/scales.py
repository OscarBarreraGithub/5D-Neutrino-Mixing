"""Scale conventions for the quark-sector MFV module.

Two scales coexist and must stay orthogonal:

* ``DEFAULT_QUARK_FIT_SCALE_GEV`` — the renormalization scale at which
  PDG 2024 MS-bar quark mass targets are reported and the fitter scores
  residuals. Chosen as ``qcd.constants.M_TOP_MS = 163.5 GeV`` so all
  quarks (light, charm, bottom, top) live in MS-bar n_f=5 between m_b
  and m_t (top runs from 162.5 to 163.5 in n_f=6, per-mille effect).

* ``DEFAULT_QUARK_TARGET_SCALE_GEV`` — the *Wilson-coefficient*
  matching/reference scale at which alpha_s is anchored for Delta-F=2
  running. Stays at 3 TeV. **Do not** repurpose this for mass targets.
"""

from __future__ import annotations

from qcd.constants import M_TOP_MS

DEFAULT_QUARK_XI_KK = 1.0
GAUGE_KK_ROOT_NN = 2.448687135269161
SPIN2_GRAVITON_KK_ROOT = 3.8317059702075125
DEFAULT_QUARK_BENCHMARK_XI_KK = DEFAULT_QUARK_XI_KK
# Wilson-coefficient reference scale (UNCHANGED at 3 TeV).
DEFAULT_QUARK_TARGET_SCALE_GEV = 3000.0
# Mass-target / fitter scoring scale (NEW: PDG 2024 -> mu = m_t(m_t)).
DEFAULT_QUARK_FIT_SCALE_GEV = M_TOP_MS
DEFAULT_QUARK_PAPER_H_RS_MAX = 0.3
DEFAULT_QUARK_BENCHMARK_H_RS_MAX = 1.0
DEFAULT_QUARK_MAX_MISALIGNMENT_STRESS = 6.0


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
