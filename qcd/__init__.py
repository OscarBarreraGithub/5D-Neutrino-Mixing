"""QCD running coupling module.

Provides numerical computation of the strong coupling constant alpha_s(mu)
at arbitrary energy scales via MS-bar RG evolution with flavor threshold
crossings.  Primarily used for neutron EDM calculations at TeV scales.

Example
-------
>>> from qcd import alpha_s
>>>
>>> # alpha_s at 1 TeV (4-loop running, 3-loop matching by default)
>>> a_s = alpha_s(1000.0)
>>> print(f"alpha_s(1 TeV) = {a_s:.4f}")
>>> 
>>> # alpha_s at 3 TeV (2-loop)
>>> a_s_2loop = alpha_s(3000.0, n_loops=2)
"""

from .running import alpha_s, alpha_s_array
from .decoupling import match_alpha_s
from .beta_function import (
    beta_coefficients,
    beta_rhs,
    beta_0,
    beta_1,
    beta_2,
    beta_3,
)
from .constants import (
    ALPHA_S_MZ,
    ALPHA_S_MZ_UNCERTAINTY,
    M_Z,
    M_CHARM,
    M_BOTTOM,
    M_TOP,
    M_TOP_MS,
    M_TOP_POLE,
    THRESHOLD_LIST,
)

__all__ = [
    # Main API
    'alpha_s',
    'alpha_s_array',
    'match_alpha_s',
    # Beta function
    'beta_coefficients',
    'beta_rhs',
    'beta_0',
    'beta_1',
    'beta_2',
    'beta_3',
    # Constants
    'ALPHA_S_MZ',
    'ALPHA_S_MZ_UNCERTAINTY',
    'M_Z',
    'M_CHARM',
    'M_BOTTOM',
    'M_TOP',
    'M_TOP_MS',
    'M_TOP_POLE',
    'THRESHOLD_LIST',
]
