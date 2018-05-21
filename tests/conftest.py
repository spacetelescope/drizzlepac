"""Custom ``pytest`` configurations."""
import pytest


def pytest_addoption(parser):
    parser.addoption("--slow", action="store_true",
                     default=False, help="Run slow tests.")


def pytest_configure(config):
    config.getini('markers').append(
        'slow: Run tests that are resource intensive')


def pytest_runtest_setup(item):
    if 'slow' in item.keywords and not item.config.getvalue("slow"):
        pytest.skip("need --slow option to run")
