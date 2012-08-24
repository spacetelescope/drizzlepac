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
    # use the not-installed-yet source tree for pywcs
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

if sys.platform != 'win32':
    EXTRA_LINK_ARGS = []
else:
    EXTRA_LINK_ARGS = ['/NODEFAULTLIB:MSVCRT' ] # , pywcslib+'/_pywcs.dll']
    EXTRA_LINK_ARGS = []
    

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


    ext = [Extension("drizzlepac.cdriz",['src/arrdrizmodule.c',
                                          'src/cdrizzleblot.c',
                                          'src/cdrizzlebox.c',
                                          'src/cdrizzleio.c',
                                          'src/cdrizzlemap.c',
                                          'src/cdrizzleutil.c',
                                          'src/cdrizzlewcs.c'],
                     define_macros=define_macros,
                     undef_macros=undef_macros,
                     include_dirs=[pythoninc] + [numpyinc] + \
                         pywcsincludes,
                     extra_link_args=EXTRA_LINK_ARGS,
                     extra_compile_args=EXTRA_COMPILE_ARGS,
                     # not needed on windows; not needed on mac
                     # libraries=['m']
                     )]

    return ext


pkg = "drizzlepac"

setupargs = {

    'version' :         '1.1.0dev',
    'description' :     "C-based MultiDrizzle",
    'author' :          "Megan Sosey, Warren Hack, Christopher Hanley, Chris Sontag",
    'author_email' :    "help@stsci.edu",
    'license' :         "http://www.stsci.edu/resources/software_hardware/pyraf/LICENSE",
    'platforms' :       ["Linux","Solaris","Mac OS X","Win"],
    'data_files' :        [( pkg+"/pars", ['lib/drizzlepac/pars/*']),
                            ( pkg+"/htmlhelp/_images", ['lib/drizzlepac/htmlhelp/_images/*']),
                            ( pkg+"/htmlhelp/_sources", ['lib/drizzlepac/htmlhelp/_sources/*']),
                            ( pkg+"/htmlhelp/_static", ['lib/drizzlepac/htmlhelp/_static/*']),
                            ( pkg+"/htmlhelp", ['lib/drizzlepac/htmlhelp/*.html']),
                            ( pkg, ['lib/drizzlepac/*.help'])],
    'scripts' :         ["scripts/mdriz","scripts/resetbits","scripts/updatenpol"] ,
    'ext_modules' :     getNumpyExtensions(),
    'package_dir' :     { 'drizzlepac' : 'lib/drizzlepac', },

    }

