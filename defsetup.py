from distutils.core import Extension
import sys, os.path, os
from distutils import sysconfig

# BUILD should be 'debug', 'profile' or 'release'
BUILD = 'release'

try:
    import numpy
except ImportError:
    print "Numpy was not found. It may not be installed or it may not be on your PYTHONPATH. Multidrizzle requires numpy v 1.0.2 or later.\n"
    raise

# This is the case for building as part of stsci_python
if os.path.exists('pywcs'):
    pywcsincludes = [os.path.join('pywcs', 'src')]
    candidates = []
    for path in os.listdir('pywcs'):
        if path.startswith('wcslib'):
            candidates.append(path)
    if len(candidates) == 1:
        pywcsincludes.append(os.path.join('pywcs', candidates[0], 'C'))
    else:
        raise SystemExit, "No suitable version of wcslib found in the current distribution of pywcs"
else:
    try:
        import pywcs
        pywcslib = pywcs.__path__[0]
        pywcsincludes = [os.path.join(pywcslib, 'include'),
                         os.path.join(pywcslib, 'include', 'wcslib')]
    except ImportError:
        raise ImportError("PyWCS was not found. It may not be installed or it may not be on your PYTHONPATH. \nPydrizzle requires pywcs 1.4 or later.\n")

if numpy.__version__ < "1.0.2":
    raise SystemExit, "Numpy 1.0.2 or later required to build Multidrizzle."

print "Building C extensions using NUMPY."

numpyinc = numpy.get_include()

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
    cfitsioinc = [ cfitsio_inc ]
    EXTRA_LINK_ARGS = []
else:
    raise Exception("Nobody ever wrote Windows support for linking with CFITSIO")
    EXTRA_LINK_ARGS = ['/NODEFAULTLIB:MSVCRT', pywcslib+'/_pywcs.dll']

def getNumpyExtensions():
    define_macros = [('PYDRIZZLE', None)]
    undef_macros = []
    EXTRA_COMPILE_ARGS = []
    if BUILD.lower() == 'debug':
        define_macros.append(('DEBUG', None))
        undef_macros.append('NDEBUG')
        if not sys.platform.startswith('sun') and \
           not sys.platform == 'win32':
            EXTRA_COMPILE_ARGS.extend(["-fno-inline", "-O0", "-g"])
    elif BUILD.lower() == 'profile':
        define_macros.append(('NDEBUG', None))
        undef_macros.append('DEBUG')
        if not sys.platform.startswith('sun') and \
           not sys.platform == 'win32':
            EXTRA_COMPILE_ARGS.extend(["-O3", "-g"])
    elif BUILD.lower() == 'release':
        # Define ECHO as nothing to prevent spurious newlines from
        # printing within the libwcs parser
        define_macros.append(('NDEBUG', None))
        undef_macros.append('DEBUG')
    else:
        raise ValueError("BUILD should be one of 'debug', 'profile', or 'release'")


    ext = [Extension("betadrizzle.cdriz",['src/arrdrizmodule.c',
                                          'src/cdrizzleblot.c',
                                          'src/cdrizzlebox.c',
                                          'src/cdrizzleio.c',
                                          'src/cdrizzlemap.c',
                                          'src/cdrizzleutil.c',
                                          'src/cdrizzlewcs.c'],
                     define_macros=define_macros,
                     undef_macros=undef_macros,
                     include_dirs=[pythoninc] + [numpyinc] + cfitsioinc + \
                         pywcsincludes,
                     extra_link_args=EXTRA_LINK_ARGS,
                     extra_compile_args=EXTRA_COMPILE_ARGS,
                     library_dirs=[cfitsio_lib],
                     libraries=['m', 'cfitsio']
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
    'scripts' :         ["lib/mdriz.py","lib/resetbits","lib/updatenpol"] ,
    'ext_modules' :     getNumpyExtensions()
    }

