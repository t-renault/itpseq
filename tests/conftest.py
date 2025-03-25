# module_dir/tests/conftest.py
import os
import pytest
import tempfile
import shutil

from pathlib import Path



@pytest.fixture(scope='module')
def data_dir():
    """Fixture that returns the path to the test data directory."""
    return Path(__file__).parent / 'test_data/'

@pytest.fixture(scope='class')
def tmp_outdir():
    """
    Fixture that creates a temporary directory for test outputs.
    Automatically cleans up after tests are done.
    """
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    shutil.rmtree(temp_dir)

@pytest.fixture
def parsed_output_dir(test_data_dir, temp_output_dir):
    """
    Fixture that runs the parsing.parse_all function on test files
    and returns the output directory.
    This is useful for tests that depend on the output of the parsing.
    """
    from itpseq.parsing import parse_all
    
    fastq_files = [
        os.path.join(test_data_dir, f) 
        for f in os.listdir(test_data_dir) 
        if f.endswith('.fastq')
    ]
    
    output_dir = os.path.join(temp_output_dir, 'parsed_output')
    os.makedirs(output_dir, exist_ok=True)
    
    parse_all(fastq_files, output_dir)
    return output_dir
