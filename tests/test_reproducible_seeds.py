"""Regression tests for process-independent research-script seeds."""

from __future__ import annotations

import pytest

from scripts.reproducible_seeds import stable_seed_offset


def test_stable_seed_offset_has_frozen_cross_process_values():
    assert stable_seed_offset(
        "flat_typical",
        modulus=1000,
        namespace="yukawa_perturbation.base_point.v1",
    ) == 564
    assert stable_seed_offset(
        "nelson_barr",
        modulus=97,
        namespace="yukawa_perturbation.gradient_field.v1",
    ) == 95


def test_stable_seed_offset_separates_namespaces_and_validates_inputs():
    first = stable_seed_offset("u2", modulus=1000, namespace="stream-a")
    second = stable_seed_offset("u2", modulus=1000, namespace="stream-b")

    assert first != second
    assert stable_seed_offset("u2", modulus=1000, namespace="stream-a") == first
    with pytest.raises(ValueError, match="label"):
        stable_seed_offset("", modulus=10, namespace="stream-a")
    with pytest.raises(ValueError, match="namespace"):
        stable_seed_offset("u2", modulus=10, namespace="")
    with pytest.raises(ValueError, match="modulus"):
        stable_seed_offset("u2", modulus=0, namespace="stream-a")
