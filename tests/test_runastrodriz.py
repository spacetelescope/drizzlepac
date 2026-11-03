"""Tests for the SKIP_APOSTERIORI_IF_VISIT_LEVEL_ALIGNMENT environment
variable used by runastrodriz to keep visit-level (FIT_SVM) a posteriori
WCS solutions and skip further alignment.
"""
import os

import pytest

from drizzlepac import runastrodriz, util

ENVVAR_NAME = 'SKIP_APOSTERIORI_IF_VISIT_LEVEL_ALIGNMENT'


@pytest.fixture(autouse=True)
def clean_envvar():
    """Keep the process environment free of the switch around every test."""
    os.environ.pop(ENVVAR_NAME, None)
    yield
    os.environ.pop(ENVVAR_NAME, None)


def read_switch():
    """Read the switch exactly as runastrodriz.process() does."""
    return util.get_envvar_switch(
        runastrodriz.envvar_skip_visit_level_aposteriori_alignment,
        description="'Skip visit-level a posteriori alignment'",
        default=False,
    )


def test_envvar_constant_matches_documented_name():
    """The module constant must match the documented variable name."""
    assert (runastrodriz.envvar_skip_visit_level_aposteriori_alignment
            == ENVVAR_NAME)
    assert ENVVAR_NAME in runastrodriz.__doc__


def test_defaults_to_off():
    """When unset, the switch must default to False (no behavior change)."""
    assert read_switch() is False


@pytest.mark.parametrize("env_value,expected", [
    ("on", True), ("yes", True), ("true", True), ("ON", True),
    ("off", False), ("no", False), ("false", False),
])
def test_envvar_values(env_value, expected):
    """The switch must honor all documented on/off values."""
    os.environ[ENVVAR_NAME] = env_value
    assert read_switch() is expected


def test_invalid_value_raises():
    """An unrecognized value must raise a ValueError."""
    os.environ[ENVVAR_NAME] = 'maybe'
    with pytest.raises(ValueError, match=ENVVAR_NAME):
        read_switch()


@pytest.mark.parametrize("skip_switch,wcsnames,expected_skip", [
    # switch on + all visit-level solutions -> skip
    (True, ['IDC_2731450pi-FIT_SVM_GAIADR2', 'IDC_2731450pi-FIT_SVM_GAIAeDR3'], True),
    # switch on but one exposure lacks a visit-level solution -> no skip
    (True, ['IDC_2731450pi-FIT_SVM_GAIADR2', 'IDC_2731450pi-HSC30'], False),
    # other a posteriori or a priori solutions do not count as visit-level
    (True, ['IDC_2731450pi-FIT_REL_GAIAeDR3'], False),
    (True, ['IDC_2731450pi-GSC240'], False),
    # switch off -> never skip, even when all are visit-level
    (False, ['IDC_2731450pi-FIT_SVM_GAIADR2'], False),
])
def test_skip_decision_logic(skip_switch, wcsnames, expected_skip):
    """Mirror of the runastrodriz.process() gating expression: alignment is
    only skipped when the switch is on AND every input exposure's primary
    WCSNAME contains 'FIT_SVM'."""
    _wnames_calfiles = [(f'file{i}_flt.fits', w) for i, w in enumerate(wcsnames)]
    svm_wcsname = all('FIT_SVM' in _w for _c, _w in _wnames_calfiles)
    assert (skip_switch and svm_wcsname) is expected_skip
