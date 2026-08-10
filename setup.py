#!/usr/bin/env python
import sys
import sysconfig
from glob import glob
from pathlib import Path

import numpy
from astropy import wcs
from setuptools import Extension, setup

FREE_THREADED_PYTHON = sysconfig.get_config_var("Py_GIL_DISABLED") == 1

# Setup C module include directories
include_dirs = []
numpy_includes = [numpy.get_include()]
wcs_include_path = Path(wcs.get_include())
wcs_includes = [
    str(wcs_include_path / "astropy_wcs"),
    str(wcs_include_path / "wcslib"),
]

include_dirs.extend(numpy_includes)
include_dirs.extend(wcs_includes)

# Setup C module macros
define_macros = [("NUMPY", "1")]
if not FREE_THREADED_PYTHON:
    define_macros.append(("Py_LIMITED_API", 0x030B0000))  # PY_VERSION_HEX for 3.11

# Handle MSVC `wcsset` redefinition
if sys.platform == "win32":
    define_macros += [("_CRT_SECURE_NO_WARNING", None), ("__STDC__", 1)]

# importing these extension modules is tested in `.github/workflows/build.yml`;
# when adding new modules here, make sure to add them to the `test_command` entry there
ext_modules = [
    Extension(
        "drizzlepac.cdriz",
        glob("src/*.c"),
        include_dirs=include_dirs,
        define_macros=define_macros,
        py_limited_api=not FREE_THREADED_PYTHON,
    )
]

SETUPTOOLS_OPTIONS = {}
if not FREE_THREADED_PYTHON:
    SETUPTOOLS_OPTIONS["bdist_wheel"] = {"py_limited_api": "cp311"}

setup(
    ext_modules=ext_modules,
    options=SETUPTOOLS_OPTIONS,
)
