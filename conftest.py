"""Project default for pytest"""
from astropy.tests.helper import enable_deprecations_as_exceptions

# Uncomment the following line to treat all DeprecationWarnings as exceptions
enable_deprecations_as_exceptions()


# Utilities to support command line passing of arguments to test_randomlist.py
def pytest_addoption(parser):
    """Support options for the pytest test_randomlist.py."""
    parser.addoption("--start_row", action="store", default=0)
    parser.addoption("--num_rows", action="store", default=5)
    parser.addoption("--master_list", action="store", default="ACSWFC3ListDefault50.csv")
