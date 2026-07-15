"""Physical constants for QCD coupling computations.

This module contains PDG values for the strong coupling at M_Z,
quark masses, and flavor threshold definitions.
"""

# Strong coupling constant at M_Z, PDG 2024 world average (MS-bar, n_f=5)
ALPHA_S_MZ = 0.1180
ALPHA_S_MZ_UNCERTAINTY = 0.0009

# Z boson mass (GeV), PDG 2024
M_Z = 91.1876

# Quark masses (GeV)
# charm and bottom: MS-bar mass evaluated at its own scale
# top: MS-bar mass at m_t from PDG 2024 cross-section/MS-bar extraction
M_CHARM = 1.27
M_BOTTOM = 4.18
M_TOP_MS = 162.5

# Top direct-reconstruction mass (GeV), PDG 2024, kept for reference.
M_TOP_POLE = 172.57

# Default alias — points to MS-bar mass (changed from pole mass in v0.2)
M_TOP = M_TOP_MS

# Ordered flavor thresholds: (mass_GeV, n_f_below, n_f_above)
THRESHOLD_LIST = [
    (M_CHARM, 3, 4),
    (M_BOTTOM, 4, 5),
    (M_TOP_MS, 5, 6),
]
