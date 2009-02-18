from distutils.core import Extension
import sys, os.path
from distutils import sysconfig

try:
    import numpy
    import numpy.numarray as nn
except ImportError:
    "Numpy was not found. It may not be installed or it may not be on your PYTHONPATH. Pydrizzle requires numpy v 1.0.2 or later.\n"

if numpy.__version__ < "1.0.2":
    raise SystemExit, "Numpy 1.0.2 or later required to build pydrizzle."

print "Building C extensions using NUMPY."

numpyinc = numpy.get_include()
numpynumarrayinc = nn.get_numarray_include_dirs()

pythonlib = sysconfig.get_python_lib(plat_specific=1)
pythoninc = sysconfig.get_python_inc()
ver = sysconfig.get_python_version()
pythonver = 'python' + ver

cfitsioinc = []

if sys.platform != 'win32':
    pydrizzle_libraries = ['m']
    cfitsioinc = ['/user/dencheva/mdriz-dev/soft-build/cfitsio']
    EXTRA_LINK_ARGS = ['-L/user/dencheva/mdriz-dev/soft-build/cfitsio']
else:
    pydrizzle_libraries = []
    EXTRA_LINK_ARGS = ['/NODEFAULTLIB:MSVCRT']


 
def getNumpyExtensions():
    ext = [Extension("pydrizzle.arrdriz",['src/arrdrizmodule.c',
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
                     extra_compile_args=['-funroll-loops', '-DPYDRIZZLE'] # , '-fno-inline', '-O0']
                     )]

    return ext


pkg = "BigBlackBox"

setupargs = {

    'version' :         '0.1',
    'description' :     "C-based MultiDrizzle",
    'author' :          "Megan Sosey, Warren Hack, Christopher Hanley",
    'author_email' :    "help@stsci.edu",
    'license' :         "http://www.stsci.edu/resources/software_hardware/pyraf/LICENSE",
    'platforms' :       ["Linux","Solaris","Mac OS X","Win"],
    'scripts' :         [] ,
    'ext_modules' :     getNumpyExtensions()
    }

