#!/usr/bin/env python
import os
import pkgutil
import sys

import numpy
from astropy import wcs

from glob import glob
from setuptools import setup, find_packages, Extension
from subprocess import check_call, CalledProcessError


if not pkgutil.find_loader('relic'):
    relic_local = os.path.exists('relic')
    relic_submodule = (relic_local and
                       os.path.exists('.gitmodules') and
                       not os.listdir('relic'))
    try:
        if relic_submodule:
            check_call(['git', 'submodule', 'update', '--init', '--recursive'])
        elif not relic_local:
            check_call(['git', 'clone', 'https://github.com/spacetelescope/relic.git'])

        sys.path.insert(1, 'relic')
    except CalledProcessError as e:
        print(e)
        exit(1)

import relic.release

PACKAGENAME = 'drizzlepac'
version = relic.release.get_info()
relic.release.write_template(version, PACKAGENAME)

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

TESTS_REQUIRE = [
    'ci_watson',
    'crds',
    'pytest',
    'pytest-remotedata'
]

setup(
    name=PACKAGENAME,
    version=version.pep386,
    author='Megan Sosey, Warren Hack, Christopher Hanley, '
           'Chris Sontag, Mihai Cara',
    author_email='help@stsci.edu',
    description='drizzle tools: combines astronomical images, including '
                'modeling distortion, removing cosmic rays, and generally '
                'improving fidelity of data in the final image',
    url='https://github.com/spacetelescope/drizzlepac',
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    python_requires='>=3.5',
    install_requires=[
        'astropy>=3.1.2',
        'fitsblender',
        'nictools',
        'nose',
        'numpy>=1.16',
        'scipy',
        'matplotlib',
        'acstools',
        'spherical-geometry>=1.2.17',
        'stsci.tools>=3.6',
        'stsci.image>=2.3.0',
        'stsci.imagestats',
        'stsci.skypac',
        'stsci.stimage',
        'stwcs',
        'tweakwcs>=0.5.0',
        'stregion',
        'requests',
        # HAP-pipeline specific:
        'astroquery',
        'photutils>=0.7',
        'pysynphot',
        'lxml',
        'scikit-image>=0.14.2',
    ],
    extras_require={
        'test': TESTS_REQUIRE
    },
    packages=find_packages(),
    package_data={
        '': ['README.md', 'LICENSE.txt'],
        'drizzlepac': [
            'pars/*',
            'pars/hap_pars/*',
            'pars/hap_pars/default_parameters/acs/hrc/*',
            'pars/hap_pars/default_parameters/acs/sbc/*',
            'pars/hap_pars/default_parameters/acs/wfc/*',
            'pars/hap_pars/default_parameters/any/*',
            'pars/hap_pars/default_parameters/wfc3/ir/*',
            'pars/hap_pars/default_parameters/wfc3/uvis/*',
            'pars/hap_pars/user_parameters/acs/hrc/*',
            'pars/hap_pars/user_parameters/acs/sbc/*',
            'pars/hap_pars/user_parameters/acs/wfc/*',
            'pars/hap_pars/user_parameters/any/*',
            'pars/hap_pars/user_parameters/wfc3/ir/*',
            'pars/hap_pars/user_parameters/wfc3/uvis/*',
            '*.help',
            'htmlhelp/*',
            'htmlhelp/_*/*',
            'htmlhelp/_*/*/*',
        ]
    },
    entry_points={
        'console_scripts': [
            'mdriz=drizzlepac.mdriz:main',
            'resetbits=drizzlepac.resetbits:main',
            'updatenpol=drizzlepac.updatenpol:main',
            'runastrodriz=drizzlepac.runastrodriz:main',
            'runsinglehap=drizzlepac.runsinglehap:main'
        ],
    },
    ext_modules=[
        Extension('drizzlepac.cdriz',
                  glob('src/*.c'),
                  include_dirs=include_dirs,
                  define_macros=define_macros),
    ],
    project_urls={
        'Bug Reports': 'https://github.com/spacetelescope/drizzlepac/issues/',
        'Source': 'https://github.com/spacetelescope/drizzlepac/',
        'Help': 'https://hsthelp.stsci.edu/',
    },
)
