#!/usr/bin/env python
import sys
from glob import glob
from pathlib import Path

import numpy
from astropy import wcs
from setuptools import setup, Extension

# Setup C module include directories
include_dirs = []
numpy_includes = [numpy.get_include()]
wcs_include_path = Path(wcs.get_include())
wcs_includes = [
    str(wcs_include_path / 'astropy_wcs'),
    str(wcs_include_path / 'wcslib'),
]

include_dirs.extend(numpy_includes)
include_dirs.extend(wcs_includes)

# Setup C module macros
define_macros = []

# Handle MSVC `wcsset` redefinition
if sys.platform == 'win32':
    define_macros += [
        ('_CRT_SECURE_NO_WARNING', None),
        ('__STDC__', 1)
    ]

# importing these extension modules is tested in `.github/workflows/build.yml`; 
# when adding new modules here, make sure to add them to the `test_command` entry there
ext_modules = [
    Extension(
        'drizzlepac.cdriz',
        glob('src/*.c'),
        include_dirs=include_dirs,
        define_macros=define_macros,
    )
]

setup(
    ext_modules=ext_modules,
)
