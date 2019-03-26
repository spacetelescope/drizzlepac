"""Project default for pytest"""
import os
import pytest

from astropy.tests.plugins.display import PYTEST_HEADER_MODULES
from astropy.tests.helper import enable_deprecations_as_exceptions


# Uncomment the following line to treat all DeprecationWarnings as exceptions
enable_deprecations_as_exceptions()
