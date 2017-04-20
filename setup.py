#!/usr/bin/env python
import os
import numpy
import pandokia
import subprocess
import sys
from glob import glob
from astropy import wcs
from setuptools import setup, find_packages, Extension


if os.path.exists('relic'):
    sys.path.insert(1, 'relic')
    import relic.release
else:
    try:
        import relic.release
    except ImportError:
        try:
            subprocess.check_call(['git', 'clone',
                'https://github.com/jhunkeler/relic.git'])
            sys.path.insert(1, 'relic')
            import relic.release
        except subprocess.CalledProcessError as e:
            print(e)
            exit(1)


version = relic.release.get_info()
relic.release.write_template(version, 'drizzlepac')

numpy_includes = numpy.get_include()
wcs_includes = [
    os.path.join(wcs.get_include(), 'astropy_wcs/'),
    os.path.join(wcs.get_include(), 'wcslib/')
]
fctx_includes = [
    os.path.join(os.path.dirname(pandokia.__file__),
                'runners', 'maker')
]

include_dirs = []
include_dirs.append(numpy_includes)
include_dirs.append(fctx_includes)
[ include_dirs.append(i) for i in wcs_includes ]


setup(
    name = 'drizzlepac',
    version = version.pep386,
    author = 'Megan Sosey, Warren Hack, Christopher Hanley, Chris Sontag, Mihai Cara',
    author_email = 'help@stsci.edu',
    description = 'drizzle tools: combines astronomical images, including modeling distortion, removing cosmic rays, and generally improving fidelity of data in the final image',
    url = 'https://github.com/spacetelescope/drizzlepac',
    classifiers = [
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    install_requires = [
        'astropy',
        'fitsblender',
        'nictools',
        'nose',
        'numpy',
        'pandokia',
        'scipy',
        'sphinx',
        'stsci.tools',
        'stsci.convolve',
        'stsci.image',
        'stsci.imagemanip',
        'stsci.imagestats',
        'stsci.ndimage',
        'stsci.stimage',
        'stwcs',
    ],

    packages = find_packages(),
    package_data = {
        '': [
            'pars/*',
            '*.help',
        ]
    },
    entry_points = {
        'console_scripts': [
            'mdriz=drizzlepac.mdriz:main',
            'resetbits=drizzlepac.resetbits:main',
            'updatenpol=drizzlepac.updatenpol:main',
            'runastrodriz=drizzlepac.runastrodriz:main'
        ],
    },
    ext_modules=[
        Extension('drizzlepac.cdriz',
            glob('src/*.c'),
            include_dirs=include_dirs,
            define_macros=[('PYDRIZZLE','1')]),
    ],
)
