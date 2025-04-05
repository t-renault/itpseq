# module_dir/tests/conftest.py
import os
import shutil
import tempfile
from pathlib import Path

import pytest


@pytest.fixture(scope='module')
def data_dir():
    """Fixture that returns the path to the test data directory."""
    return Path(__file__).parent / 'test_data'


@pytest.fixture(scope='class')
def tmp_outdir():
    """
    Fixture that creates a temporary directory for test outputs.
    Automatically cleans up after tests are done.
    """
    temp_dir = tempfile.mkdtemp()
    yield Path(temp_dir)
    shutil.rmtree(temp_dir)
