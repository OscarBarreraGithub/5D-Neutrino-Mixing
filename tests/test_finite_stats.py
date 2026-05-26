from __future__ import annotations

import math

import pytest

from quarkConstraints.finite_stats import wilson_upper_limit


def test_zero_pass_limit_matches_large_n_rule_of_thumb():
    assert wilson_upper_limit(0, 1_000) == pytest.approx(0.003673, rel=1e-4)


def test_zero_pass_limit_for_run3_sample_size():
    assert wilson_upper_limit(0, 1_600_000) == pytest.approx(2.304e-6, rel=1e-4)


def test_wilson_upper_limit_validates_counts():
    with pytest.raises(ValueError):
        wilson_upper_limit(1, 0)
    with pytest.raises(ValueError):
        wilson_upper_limit(11, 10)


def _wilson_upper_limit_closed_form(k: int, n: int, z: float) -> float:
    """Reference closed-form Wilson upper limit, independent of scipy.

    Used as a second cross-check for k>0 cases.  Matches the formula in
    Brown, Cai, DasGupta (Stat. Sci. 2001) §3:

        UL = (p_hat + z^2/(2n) + z * sqrt(p_hat*(1-p_hat)/n + z^2/(4n^2)))
             / (1 + z^2/n)

    where ``p_hat = k/n``.  This is the upper endpoint of the Wilson-score
    confidence interval at confidence level ``1 - 2 * (1 - Phi(z))`` (so
    z=1.96 -> 95%, z=1.92 -> ~94.5% which is the convention used by
    ``quarkConstraints.finite_stats.wilson_upper_limit``; see R07-I1 for
    the z=1.92 vs z=1.96 footnote).
    """
    p_hat = k / n
    z2 = z * z
    denom = 1.0 + z2 / n
    center = (p_hat + z2 / (2.0 * n)) / denom
    margin = z * math.sqrt(p_hat * (1.0 - p_hat) / n + z2 / (4.0 * n * n)) / denom
    return center + margin


@pytest.mark.parametrize("k,n", [(1, 100), (1, 1000), (5, 1000), (10, 1000), (50, 1000)])
def test_wilson_upper_limit_k_gt_zero_matches_closed_form(k, n):
    """R07-I2: extend Wilson-UL coverage from k=0 to k > 0.

    Cross-check the repo helper (which uses z=1.92, the local convention
    documented in R07-I1) against the closed-form formula at the same z,
    for k in {1, 5, 10, 50} and n in {100, 1000}.
    """
    z_local = 1.92
    expected = _wilson_upper_limit_closed_form(k, n, z=z_local)
    actual = wilson_upper_limit(k, n, z=z_local)
    assert actual == pytest.approx(expected, rel=1e-12, abs=1e-15)


@pytest.mark.parametrize("k,n", [(1, 100), (5, 1000), (50, 1000)])
def test_wilson_upper_limit_k_gt_zero_matches_scipy_at_z_eq_1p96(k, n):
    """R07-I2: cross-validate k>0 against scipy's standard 95% Wilson UL.

    ``scipy.stats.binomtest(...).proportion_ci(method='wilson', ...)``
    returns the standard Wilson score interval (z = 1.959964..., 95% CL).
    We feed that exact z into the repo helper and assert agreement.

    Requires ``scipy >= 1.10`` for ``BinomTestResult.proportion_ci``.  If
    that floor is ever lowered the closed-form parametrized test above is
    sufficient on its own; this case is the "scipy as oracle" cross-check
    requested in CLEANUP_PLAN.md §C C04.
    """
    pytest.importorskip("scipy", minversion="1.10")
    from scipy.stats import binomtest

    result = binomtest(k=k, n=n)
    ci = result.proportion_ci(method="wilson", confidence_level=0.95)
    scipy_upper = float(ci.high)

    # Reproduce the same Wilson interval at z corresponding to the 95% CL.
    # statsmodels-style two-sided 95% uses z = Phi^{-1}(0.975) = 1.959964...
    from scipy.stats import norm

    z_95 = float(norm.ppf(0.975))
    helper_upper = wilson_upper_limit(k, n, z=z_95)
    assert helper_upper == pytest.approx(scipy_upper, rel=1e-12, abs=1e-15)


def test_wilson_upper_limit_k_gt_zero_spot_check_k5_n1000():
    """Spot-check tying down a concrete numerical k>0 value (R07-I2).

    For k=5, n=1000 at the repo's z=1.92 convention the Wilson upper
    limit is approximately 0.01146.  (The C04 dispatch-prompt hand-
    computed value of ~0.0089 was off by the ``+z*sqrt(...)`` margin
    term; the correct closed-form evaluates to 0.011463... — see the
    parametrized closed-form cross-check above.)
    """
    value = wilson_upper_limit(5, 1000, z=1.92)
    assert value == pytest.approx(0.01146, abs=5e-5)
