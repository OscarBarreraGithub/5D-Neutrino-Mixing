This directory contains a small set of tools for computing neutrino masses and finding the allowed parameter range under a cosmological sum constraint.

We fit for for both normal and inverted ordering.

### massConfig.py:

Experimental constraints and a function: `compute_masses`. The function returns $m_1, m_2, m_3, M_{tot}$. 


### massConstraints.py:
Sweeps the parameter space and returns arrays of the allowed values, subject to the configured sum-of-masses constraint.


### allowedMass.ipynb:
Example notebook that imports the above modules, computes allowed ranges, and makes plots.