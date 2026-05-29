"""Adapter boundary for CKM extraction arithmetic.

Constraint modules import this adapter rather than reaching into
``quarkConstraints.ckm_extraction`` directly.
"""

from __future__ import annotations

from quarkConstraints.ckm_extraction import (
    KL3VusExtraction,
    VusConsistencyResult,
    extract_vus_from_kl3,
    vus_consistency_pull,
)

K018_RS_MATCHING_GAP = (
    "NEEDS-HUMAN-PHYSICS: rigorous RS matching for K_l3 |V_us| extraction "
    "requires charged-current W/W'/KK electroweak quark couplings, possible "
    "right-handed-current and lepton-sector/G_F shifts, and radiative/isospin "
    "convention bookkeeping that are not available on ParameterPoint; v1 uses "
    "the YAML data/lattice extraction with no point-dependent RS shift."
)

__all__ = [
    "K018_RS_MATCHING_GAP",
    "KL3VusExtraction",
    "VusConsistencyResult",
    "extract_vus_from_kl3",
    "vus_consistency_pull",
]
