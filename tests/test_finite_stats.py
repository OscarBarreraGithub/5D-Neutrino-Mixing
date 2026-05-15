from __future__ import annotations

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
