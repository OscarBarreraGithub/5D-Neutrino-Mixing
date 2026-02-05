"""QCD running coupling module.

Provides numerical computation of the strong coupling constant alpha_s(mu)
at arbitrary energy scales via MS-bar RG evolution with flavor threshold
crossings.

Example
-------
>>> from qcd import alpha_s
>>>
>>> # Default: 4-loop running, 3-loop decoupling
>>> a_s = alpha_s(1000.0)
>>>
>>> # Precision presets
>>> alpha_s(3000.0, precision='low')   # 3-loop, continuous matching
>>> alpha_s(3000.0, precision='high')  # 4-loop, 3-loop decoupling
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
