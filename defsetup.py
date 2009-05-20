from distutils.core import Extension
import sys, os.path, os
from distutils import sysconfig

try:
    import numpy
    import numpy.numarray as nn
except ImportError:
    "Numpy was not found. It may not be installed or it may not be on your PYTHONPATH. Pydrizzle requires numpy v 1.0.2 or later.\n"

try:
    #import pywcs
    pywcs_path = ['/user/hack/dev/release/lib/python/pywcs']
    pywcslib = pywcs_path[0]
except ImportError:
    "PyWCS was not found. It may not be installed or it may not be on your PYTHONPATH. \nPydrizzle requires numpy v 1.0.2 or later.\n"

if numpy.__version__ < "1.0.2":
    raise SystemExit, "Numpy 1.0.2 or later required to build pydrizzle."

print "Building C extensions using NUMPY."

numpyinc = numpy.get_include()
numpynumarrayinc = nn.get_numarray_include_dirs()

pythonlib = sysconfig.get_python_lib(plat_specific=1)
pythoninc = sysconfig.get_python_inc()
ver = sysconfig.get_python_version()
pythonver = 'python' + ver

if 'CFITSIO_LIB' in os.environ :
    cfitsio = os.environ['CFITSIO_LIB']
else :
    cfitsio = '/usr/stsci/cfitsio'

if os.path.exists(cfitsio+"/include") :
    # pointing at the installed library
    cfitsio_inc = cfitsio + "/include"
    cfitsio_lib = cfitsio + "/lib"
else :
    # pointing at a source distribution
    # (still needs to be compiled with "make" but not installed with "make install")
    cfitsio_inc = cfitsio 
    cfitsio_lib = cfitsio 

if sys.platform != 'win32':
    pydrizzle_libraries = ['m']
    cfitsioinc = [ cfitsio_inc ]
    EXTRA_LINK_ARGS = ['-L'+cfitsio_lib, pywcslib+'/_pywcs.so']
else:
    raise Exception("Nobody ever wrote Windows support for linking with CFITSIO")
    pydrizzle_libraries = []
    EXTRA_LINK_ARGS = ['/NODEFAULTLIB:MSVCRT']

cfitsioinc += [os.path.join(pywcslib, 'include'), os.path.join(pywcslib, 'include', 'wcslib')]


def getNumpyExtensions():
    ext = [Extension("betadrizzle.cdriz",['src/arrdrizmodule.c',
                                          'src/cdrizzleblot.c',
                                          'src/cdrizzlebox.c',
                                          'src/cdrizzleio.c',
                                          'src/cdrizzlemap.c',
                                          'src/cdrizzleutil.c',
                                          'src/cdrizzlewcs.c'],
                     define_macros=[('NUMPY', '1')],
                     # undef_macros=['NDEBUG'],
                     include_dirs=[pythoninc] + [numpyinc]+ cfitsioinc + numpynumarrayinc,
                     # library_dirs=[],
                     extra_link_args=EXTRA_LINK_ARGS,
                     libraries=['m', 'cfitsio'],
                     extra_compile_args=['-funroll-loops', '-DPYDRIZZLE','-g'] # , '-fno-inline', '-O0']
                     )]

    return ext


pkg = "betadrizzle"

setupargs = {

    'version' :         '0.1',
    'description' :     "C-based MultiDrizzle",
    'author' :          "Megan Sosey, Warren Hack, Christopher Hanley",
    'author_email' :    "help@stsci.edu",
    'license' :         "http://www.stsci.edu/resources/software_hardware/pyraf/LICENSE",
    'platforms' :       ["Linux","Solaris","Mac OS X","Win"],
    'data_files' :        [( pkg+"/pars", ['lib/pars/*']),( pkg, ['lib/*.help'])],
    'scripts' :         ["lib/mdriz.py"] ,
    'ext_modules' :     getNumpyExtensions()
    }

