"""Shared fixtures for the flavor_catalog_constraints test suite."""

from __future__ import annotations

import pytest

import flavor_catalog_constraints as fcc


@pytest.fixture(scope="session")
def registry():
    """The discovered ``{process_id: instance}`` registry (session-scoped)."""
    return fcc.all_constraints()


@pytest.fixture(scope="session")
def import_failures():
    """Per-module import failures captured during discovery."""
    return fcc.import_failures()
