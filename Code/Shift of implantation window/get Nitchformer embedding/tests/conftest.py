"""Test configuration for loading pipeline scripts as modules."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))


@pytest.fixture
def load_script():
    """Return a helper that imports a numbered script by filename."""
    def _load(filename: str):
        name = "testloaded_" + filename.replace(".py", "").replace("-", "_")
        spec = importlib.util.spec_from_file_location(name, SCRIPTS / filename)
        if spec is None or spec.loader is None:
            raise ImportError(filename)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    return _load
