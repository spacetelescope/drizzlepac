#!/usr/bin/env python
import os.path
import sys

import numpy
from astropy import wcs

from glob import glob
from setuptools import setup, find_packages, Extension


# Setup C module include directories
include_dirs = []
numpy_includes = [numpy.get_include()]
wcs_includes = [os.path.join(wcs.get_include(), 'astropy_wcs'),
                os.path.join(wcs.get_include(), 'wcslib')]

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

setup(
    use_scm_version=True,
    setup_requires=["setuptools_scm"],
    packages=find_packages(),
    package_data={
        '': ['README.md', 'LICENSE.txt'],
        'drizzlepac': [
            'pars/*',
            'pars/hap_pars/*',
            'pars/hap_pars/mvm_parameters/*',
            'pars/hap_pars/mvm_parameters/acs/hrc/*',
            'pars/hap_pars/mvm_parameters/acs/sbc/*',
            'pars/hap_pars/mvm_parameters/acs/wfc/*',
            'pars/hap_pars/mvm_parameters/any/*',
            'pars/hap_pars/mvm_parameters/wfc3/ir/*',
            'pars/hap_pars/mvm_parameters/wfc3/uvis/*',
            'pars/hap_pars/svm_parameters/*',
            'pars/hap_pars/svm_parameters/acs/hrc/*',
            'pars/hap_pars/svm_parameters/acs/sbc/*',
            'pars/hap_pars/svm_parameters/acs/wfc/*',
            'pars/hap_pars/svm_parameters/any/*',
            'pars/hap_pars/svm_parameters/wfc3/ir/*',
            'pars/hap_pars/svm_parameters/wfc3/uvis/*',
            'pars/psfs/*',
            'pars/psfs/acs/hrc/*',
            'pars/psfs/acs/sbc/*',
            'pars/psfs/acs/wfc/*',
            'pars/psfs/wfc3/ir/*',
            'pars/psfs/wfc3/uvis/*',
            '*.help',
            'htmlhelp/*',
            'htmlhelp/_*/*',
            'htmlhelp/_*/*/*',
        ]
    },
    ext_modules=[
        Extension('drizzlepac.cdriz',
                  glob('src/*.c'),
                  include_dirs=include_dirs,
                  define_macros=define_macros),
    ],
)
