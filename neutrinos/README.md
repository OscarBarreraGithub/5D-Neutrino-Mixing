This directory contains tools for computing the PMNS matrix and allowed neutrino masses from experimental constraints.

The functions allow either normal or inverted mass ordering.

### `neutrinoValues.py`

Here we set the default neutrino parameters derived from experiments (mass-squared splittings, mixing angles, etc..) and helper functions to compute neutrino masses and the PMNS matrix.

### `massConstraints.py`

Performs parameter-space sweeps over neutrino masses and collects allowed parameter sets.

### `allowedMass.ipynb`

Notebook demonstrating how to import `neutrinoValues` and `massConstraints`, compute allowed mass ranges, and produce plots of the allowed neutrino masses.

### `PMNS.ipynb`

Notebook that builds the PMNS matrix from model parameters. Note: there may be a subtlety in the sign/convention of the Majorana phase.
