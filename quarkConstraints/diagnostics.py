"""Diagnostics for the quark-sector accepted-points CSVs.

Helpers in this module recompute MS-bar quark masses from the explicit
Yukawa matrices and bulk-mass-parameter columns stored in
``data/accepted_points_with_yukawas*.csv`` and run the resulting singular
values to a single common scale via :func:`qcd.mass_running.run_msbar_mass`.

The intended use is post-hoc auditing: after a scan completes, we need a
way to confirm that *if* the Yukawa columns are taken at face value,
the resulting fitted masses really do reproduce the PDG MS-bar targets at
``mu_common = m_t(m_t)``. This is an independent re-derivation; the
fitter's own ``masses_up`` / ``masses_down`` columns are the production
truth, but the diagnostics let us cross-check that the saved Yukawas
encode that truth without information loss.

Per plan v3 §6 (NEW: ``quarkConstraints/diagnostics.py``).
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

from qcd.mass_running import run_msbar_mass
from warpConfig.baseParams import V_EWSB, get_warp_params
from warpConfig.wavefuncs import f_IR

_UP_FLAVORS: tuple[str, str, str] = ("u", "c", "t")
_DOWN_FLAVORS: tuple[str, str, str] = ("d", "s", "b")
_NF_AT_TARGET = 5  # masses are reported in MS-bar n_f=5 between m_b and m_t.


def _row_lookup(row: Mapping[str, object], key: str) -> float:
    if key not in row:
        raise KeyError(
            f"diagnostics row is missing required column {key!r}; "
            f"got columns: {sorted(row.keys()) if hasattr(row, 'keys') else type(row)}"
        )
    return float(row[key])


def _read_yukawa(row: Mapping[str, object], block: str) -> np.ndarray:
    """Read a 3x3 complex Yukawa from columns ``Y_<block>_<i><j>_re/im``."""
    out = np.zeros((3, 3), dtype=np.complex128)
    for i in range(3):
        for j in range(3):
            re = _row_lookup(row, f"Y_{block}_{i+1}{j+1}_re")
            im = _row_lookup(row, f"Y_{block}_{i+1}{j+1}_im")
            out[i, j] = complex(re, im)
    return out


def _read_c_triplet(row: Mapping[str, object], prefix: str) -> np.ndarray:
    return np.array(
        [_row_lookup(row, f"{prefix}{i+1}") for i in range(3)],
        dtype=float,
    )


def extract_msbar_masses_from_yukawa_row(
    row: Mapping[str, object],
    scale_GeV: float,
    *,
    k: float = 1.2209e19,
    Lambda_IR: float = 3000.0,
    n_loops: int = 4,
) -> dict[str, float]:
    """Recompute MS-bar quark masses at ``scale_GeV`` from a CSV row.

    The row must expose the columns laid down by
    ``scripts/export_accepted_quark_scan_with_yukawas.py``:

    * ``Y_u_<ij>_re``, ``Y_u_<ij>_im`` for i,j in 1..3
    * ``Y_d_<ij>_re``, ``Y_d_<ij>_im`` for i,j in 1..3
    * ``c_Q1..c_Q3``, ``c_u1..c_u3``, ``c_d1..c_d3``

    Parameters
    ----------
    row
        Mapping from column names to scalars (e.g. a ``pd.Series`` or dict).
    scale_GeV
        Common scale at which to report masses.
    k, Lambda_IR
        Geometry parameters used to build the warp factor; defaults match
        the repo benchmark (Planck-mass curvature, 3 TeV IR brane).
    n_loops
        Loop order forwarded to :func:`qcd.mass_running.run_msbar_mass`.

    Returns
    -------
    dict[str, float]
        Mapping ``flavor -> m_q(scale_GeV)`` in GeV, ordered
        ``(u, c, t, d, s, b)``.

    Notes
    -----
    The Yukawa matrices stored in the CSV are taken to live at the
    fitter's natural matching scale (``m_t(m_t)`` in the PDG-2024 layout,
    or ``mu = 3 TeV`` in the legacy layout). Singular values are run from
    that natural scale to ``scale_GeV`` in MS-bar n_f=5; threshold matching
    handles crossing m_t and m_b if ``scale_GeV`` is below them.
    """
    Y_u = _read_yukawa(row, "u")
    Y_d = _read_yukawa(row, "d")
    c_Q = _read_c_triplet(row, "c_Q")
    c_u = _read_c_triplet(row, "c_u")
    c_d = _read_c_triplet(row, "c_d")

    params = get_warp_params(k=k, Lambda_IR=Lambda_IR)
    epsilon = float(params["epsilon"])

    F_Q = np.asarray(f_IR(c_Q, epsilon), dtype=float)
    F_u = np.asarray(f_IR(c_u, epsilon), dtype=float)
    F_d = np.asarray(f_IR(c_d, epsilon), dtype=float)

    prefactor = 2.0 * V_EWSB
    M_u = prefactor * np.diag(F_Q) @ Y_u @ np.diag(F_u)
    M_d = prefactor * np.diag(F_Q) @ Y_d @ np.diag(F_d)

    sv_u = np.sort(np.linalg.svd(M_u, compute_uv=False))
    sv_d = np.sort(np.linalg.svd(M_d, compute_uv=False))

    # The CSV-stored Yukawas are taken at the fitter's working scale; for
    # the PDG-2024 layout that scale is m_t(m_t). Run each SV to the user
    # requested ``scale_GeV`` in MS-bar n_f=5.
    from qcd.constants import M_TOP_MS

    fitter_scale = M_TOP_MS

    out: dict[str, float] = {}
    for flavor, m_at_fit in zip(_UP_FLAVORS, sv_u):
        out[flavor] = run_msbar_mass(
            m_ref=float(m_at_fit),
            mu_ref=fitter_scale,
            mu_target=scale_GeV,
            n_f_ref=_NF_AT_TARGET,
            n_loops=n_loops,
        )
    for flavor, m_at_fit in zip(_DOWN_FLAVORS, sv_d):
        out[flavor] = run_msbar_mass(
            m_ref=float(m_at_fit),
            mu_ref=fitter_scale,
            mu_target=scale_GeV,
            n_f_ref=_NF_AT_TARGET,
            n_loops=n_loops,
        )
    return out
