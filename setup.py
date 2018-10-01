#!/usr/bin/env python
from __future__ import print_function


import inspect
import os
import pkgutil
import shutil
import sys
import importlib


try:
    _mnfe = ModuleNotFoundError
except NameError:
    ModuleNotFoundError = ImportError

try:
    import pandokia
except (ImportError, NameError, ModuleNotFoundError):
    pandokia = False

from glob import glob
from setuptools import setup, find_packages, Extension, _install_setup_requires
from setuptools.command.install import install
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
SETUP_REQUIRES = [
    'numpy',
    'astropy',
    'sphinx',
]

# Due to overriding `install` and `build_sphinx` we need to download
# setup_requires dependencies before reaching `setup()`. This allows
# `sphinx` to exist before the `BuildSphinx` class is injected.
_install_setup_requires(dict(setup_requires=SETUP_REQUIRES))

for dep_pkg in SETUP_REQUIRES:
    try:
        importlib.import_module(dep_pkg)
    except ImportError:
        print("{0} is required in order to install '{1}'.\n"
              "Please install {0} first.".format(dep_pkg, PACKAGENAME),
              file=sys.stderr)
        exit(1)

import numpy
from astropy import wcs
from sphinx.cmd.build import build_main
from sphinx.setup_command import BuildDoc

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

# Deprecation warning:
#    Pandokia integration will be removed in a later release.
if pandokia:
    fctx_includes = [os.path.join(os.path.dirname(pandokia.__file__),
                                  'runners', 'maker')]
    include_dirs.extend(fctx_includes)

# Distribute compiled documentation alongside the installed package
docs_compiled_src = os.path.normpath('build/sphinx/html')
docs_compiled_dest = os.path.normpath('{0}/htmlhelp'.format(PACKAGENAME))

class InstallCommand(install):
    """Ensure drizzlepac's C extensions are available when imported relative
    to the documentation, instead of relying on `site-packages`. What comes
    from `site-packages` may not be the same drizzlepac that was *just*
    compiled.
    """
    def run(self):
        build_cmd = self.reinitialize_command('build_ext')
        build_cmd.inplace = 1
        self.run_command('build_ext')

        # Explicit request for old-style install?  Just do it
        if self.old_and_unmanageable or self.single_version_externally_managed:
            install.run(self)
        elif not self._called_from_setup(inspect.currentframe()):
            # Run in backward-compatibility mode to support bdist_* commands.
            install.run(self)
        else:
            self.do_egg_install()

        if not os.path.exists(docs_compiled_dest):
            print('\nwarning: Sphinx "htmlhelp" documentation was NOT bundled!\n'
                  '         Execute the following then reinstall:\n\n'
                  '         $ python setup.py build_sphinx\n\n',
                  file=sys.stderr)



class BuildSphinx(BuildDoc):
    """Build Sphinx documentation after compiling C extensions"""

    description = 'Build Sphinx documentation'

    def initialize_options(self):
        BuildDoc.initialize_options(self)

    def finalize_options(self):
        BuildDoc.finalize_options(self)

    def run(self):
        build_cmd = self.reinitialize_command('build_ext')
        build_cmd.inplace = 1
        self.run_command('build_ext')
        build_main(['-b', 'html', 'doc/source', 'build/sphinx/html'])

        # Bundle documentation inside of drizzlepac
        if os.path.exists(docs_compiled_src):
            if os.path.exists(docs_compiled_dest):
                shutil.rmtree(docs_compiled_dest)

            shutil.copytree(docs_compiled_src, docs_compiled_dest)


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
    setup_requires=SETUP_REQUIRES,
    install_requires=[
        'astropy',
        'fitsblender',
        'nictools',
        'nose',
        'numpy',
        'scipy',
        'spherical-geometry',
        'stsci.tools',
        'stsci.convolve',
        'stsci.image>=2.3.0',
        'stsci.imagemanip',
        'stsci.imagestats',
        'stsci.ndimage',
        'stsci.skypac',
        'stsci.stimage',
        'stwcs',
        'stregion',
    ],
    packages=find_packages(),
    package_data={
        '': ['README.md', 'LICENSE.txt'],
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
                  include_dirs=include_dirs,
                  define_macros=define_macros),
    ],
    cmdclass={
        'install': InstallCommand,
        'build_sphinx': BuildSphinx,
    },
    project_urls={
        'Bug Reports': 'https://github.com/spacetelescope/drizzlepac/issues/',
        'Source': 'https://github.com/spacetelescope/drizzlepac/',
        'Help': 'https://hsthelp.stsci.edu/',
    },
)
