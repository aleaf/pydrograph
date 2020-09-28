import os
from pathlib import Path
import shutil
import pytest
from .test_examples import nwis_instance, extent_poly, field_sites


@pytest.fixture(scope="session")
def project_root_path():
    filepath = os.path.split(os.path.abspath(__file__))[0]
    return os.path.normpath(os.path.join(filepath, '../../'))


@pytest.fixture(scope="session")
def test_data_path(project_root_path):
    """Root folder for the project (with setup.py),
    two levels up from the location of this file.
    """
    return Path(project_root_path, 'pydrograph/tests/data')


@pytest.fixture(scope="session", autouse=True)
def tmpdir(project_root_path):
    folder = project_root_path + '/pydrograph/tests/tmp'
    if os.path.isdir(folder):
        shutil.rmtree(folder)
    os.makedirs(folder)
    return folder