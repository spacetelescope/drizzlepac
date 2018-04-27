"""Project default for pytest"""
import os
import pytest

from astropy.tests.plugins.display import PYTEST_HEADER_MODULES
from astropy.tests.helper import enable_deprecations_as_exceptions


# Uncomment the following line to treat all DeprecationWarnings as exceptions
enable_deprecations_as_exceptions()

"""
try:
    PYTEST_HEADER_MODULES['Astropy'] = 'astropy'
    PYTEST_HEADER_MODULES['asdf'] = 'asdf'
    del PYTEST_HEADER_MODULES['h5py']
except (NameError, KeyError):
    pass

pytest_plugins = [
    'asdf.tests.schema_tester'
]
"""

def pytest_addoption(parser):
    # Add option to run slow tests
    parser.addoption(
        "--runslow",
        action="store_true",
        help="run slow tests"
    )

    # Add option to use big data sets
    parser.addoption(
        "--bigdata",
        action="store_true",
        help="use big data sets (intranet)"
    )


@pytest.fixture(scope='function', autouse=True)
def _jail(tmpdir):
    """ Perform test in a pristine temporary working directory
    """
    os.chdir(tmpdir.strpath)
    yield
