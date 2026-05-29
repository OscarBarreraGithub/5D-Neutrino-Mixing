# Implement E002 — muon EDM |d_μ|. family=edm_neutrino. REUSE edm module.
REUSE quarkConstraints/edm.py (built for E001) with lepton='mu'. SM negligible → pure-NP bound vs the experimental limit (from E002.yaml), HARD. RS one-loop CP-odd dipole = documented proxy NEEDS-HUMAN-PHYSICS. Append-only (don't modify E001's electron path). Anchor guard: require the canonical limit block.
