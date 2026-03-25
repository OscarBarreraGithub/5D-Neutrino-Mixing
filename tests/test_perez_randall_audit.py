import importlib.util
from pathlib import Path

import numpy as np

_SCRIPT_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "audit_perez_randall_consistency.py"
)
_SPEC = importlib.util.spec_from_file_location("audit_perez_randall_consistency", _SCRIPT_PATH)
_MODULE = importlib.util.module_from_spec(_SPEC)
assert _SPEC is not None and _SPEC.loader is not None
_SPEC.loader.exec_module(_MODULE)

masses_from_kyn = _MODULE.masses_from_kyn
paper_like_geometry = _MODULE.paper_like_geometry
required_kyn = _MODULE.required_kyn


def test_paper_like_geometry_matches_quoted_overlap_scales():
    """Paper-like geometry should reproduce Eq. (11)/Table I overlap scales."""
    state = paper_like_geometry()

    assert np.isclose(state["f_L"], 0.016, rtol=0.06)
    assert np.isclose(state["f_N"], 0.48, rtol=0.02)
    assert np.isclose(state["f_N_UV"], 1.6e-4, rtol=0.15)


def test_displayed_eq10_neutrino_yukawas_do_not_reproduce_eq7_masses():
    """Eq. (10) kY_N values are not consistent with Eq. (7) under Eq. (6)."""
    state = paper_like_geometry()

    eq7_masses = np.array([0.002, 0.009, 0.05])
    eq10_kyn = np.array([0.02, 0.03, 0.07])

    implied_masses = masses_from_kyn(
        eq10_kyn,
        f_L=state["f_L"],
        f_N=state["f_N"],
        f_N_UV=state["f_N_UV"],
        M_N=state["M_N"],
    )
    needed_kyn = required_kyn(
        eq7_masses,
        f_L=state["f_L"],
        f_N=state["f_N"],
        f_N_UV=state["f_N_UV"],
        M_N=state["M_N"],
    )

    # The displayed Eq. (10) values undershoot the required neutrino masses badly.
    assert np.all(implied_masses < 0.2 * eq7_masses)
    assert implied_masses[2] < 0.1 * eq7_masses[2]
    # The kY_N values required by Eq. (6) are much larger than Eq. (10).
    assert np.all(needed_kyn > 2.0 * eq10_kyn)
    assert np.isclose(needed_kyn[2], 0.24626311, rtol=0.03)
