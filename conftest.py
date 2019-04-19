"""Project default for pytest"""

from astropy.tests.helper import enable_deprecations_as_exceptions


# Uncomment the following line to treat all DeprecationWarnings as exceptions
enable_deprecations_as_exceptions()

# Utilities to support command line passing of arguments to test_randomlist.py
def pytest_addoption(parser):
    """Support options for the pytest test_randomlist.py."""
    parser.addoption("--start_row", action="store", default=0)
    parser.addoption("--num_rows", action="store", default=50)

def pytest_generate_tests(metafunc):
    """Get the command line option."""
    option_value1 = metafunc.config.option.start_row
    if 'start_row' in metafunc.fixturenames and option_value1 is not None:
        metafunc.parametrize("start_row", [option_value1])

    option_value2 = metafunc.config.option.num_rows
    if 'num_rows' in metafunc.fixturenames and option_value2 is not None:
        metafunc.parametrize("num_rows", [option_value2])
