"""Physical constants for Yukawa computations.

This module contains PDG values for charged lepton masses and unit conversion factors.
"""

# Charged lepton masses (GeV) - PDG 2024
M_ELECTRON = 0.51099895e-3   # 0.511 MeV
M_MUON = 105.6583755e-3      # 105.7 MeV
M_TAU = 1.77686              # 1.777 GeV

# Tuple for convenience (e, μ, τ)
LEPTON_MASSES = (M_ELECTRON, M_MUON, M_TAU)

# Unit conversion
EV_TO_GEV = 1e-9
