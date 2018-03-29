#!/usr/bin/env python
import os
import numpy
import shutil
import subprocess
import sys

try:
    _mnfe = ModuleNotFoundError
except NameError:
    ModuleNotFoundError = ImportError

try:
    import pandokia
except (ImportError, NameError, ModuleNotFoundError):
    pandokia = False

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


NAME = 'drizzlepac'
version = relic.release.get_info()
relic.release.write_template(version, NAME)

# Setup C module include directories
include_dirs = []
numpy_includes = [numpy.get_include()]
wcs_includes = [os.path.join(wcs.get_include(), 'astropy_wcs'),
                os.path.join(wcs.get_include(), 'wcslib')]

include_dirs.extend(numpy_includes)
include_dirs.extend(wcs_includes)

# Deprecation warning:
#    Pandokia integration will be removed in a later release.
if pandokia:
    fctx_includes = [os.path.join(os.path.dirname(pandokia.__file__),
                                  'runners', 'maker')]
    include_dirs.extend(fctx_includes)

# Distribute compiled documentation alongside the installed package
docs_compiled_src = os.path.normpath('build/sphinx/html')
docs_compiled_dest = os.path.normpath('{0}/htmlhelp'.format(NAME))

if os.path.exists(docs_compiled_src):
    if os.path.exists(docs_compiled_dest):
        shutil.rmtree(docs_compiled_dest)

    shutil.copytree(docs_compiled_src, docs_compiled_dest)
else:
    if len(sys.argv) > 1 and 'build_sphinx' not in sys.argv[1]:
        print('\nwarning: SPHINX DOCUMENTATION WILL NOT BE INSTALLED!\n'
              '         Please run: python {0} build_sphinx\n'
              ''.format(sys.argv[0]), file=sys.stderr)


setup(
    name=NAME,
    version=version.pep386,
    author='Megan Sosey, Warren Hack, Christopher Hanley, Chris Sontag, Mihai Cara',
    author_email='help@stsci.edu',
    description='drizzle tools: combines astronomical images, including modeling distortion, removing cosmic rays, and generally improving fidelity of data in the final image',
    url='https://github.com/spacetelescope/drizzlepac',
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    setup_requires=[
        'graphviz',
        'sphinx_rtd_theme',
        'stsci_rtd_theme',
    ],
    install_requires=[
        'astropy',
        'fitsblender',
        'nictools',
        'nose',
        'numpy',
        'scipy',
        'spherical_geometry',
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
    packages=find_packages(),
    package_data={
        'drizzlepac': [
            'pars/*',
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
            'runastrodriz=drizzlepac.runastrodriz:main'
        ],
    },
    ext_modules=[
        Extension('drizzlepac.cdriz',
                  glob('src/*.c'),
                  include_dirs=include_dirs),
    ],
)
