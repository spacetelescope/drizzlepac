import glob
import os
import sys


# BUILD should be 'debug', 'profile' or 'release'
# TODO: How often is this actually mucked with? Would it be worth adding a
# custom command that adds a command-line option for this?
BUILD = 'release'


def setup_hook(config):
    # First get the pywcs includes path
    # This is the case for building as part of stsci_python
    pywcs_dir = os.path.abspath(os.path.join(os.path.pardir, 'pywcs'))
    if os.path.exists(pywcs_dir):
        pywcsincludes = [os.path.join(pywcs_dir, 'src')]
        wcslibs = glob.glob(os.path.join(pywcs_dir, 'wcslib*'))
        if len(wcslibs) == 1:
            pywcsincludes.append(os.path.join(pywcs_dir, wcslibs[0], 'C'))
        else:
            raise SystemExit('No suitable version of wcslib found in the '
                             'pywcs distribution at %s' % pywcs_dir)
    else:
        # If pywcs is otherwise already installed...
        # TODO: Maybe we can eventually make pywcs a setup requirement for
        # drizzlepac, so long as pywcs itself installs easily enough...
        try:
            import pywcs
            # TODO: It would be nice if pywcs had a get_includes() function a
            # la numpy
            pywcslib = pywcs.__path__[0]
            pywcsincludes = [os.path.join(pywcslib, 'include'),
                             os.path.join(pywcslib, 'include', 'wcslib')]
            # TODO: Institute version check?
        except ImportError:
            raise ImportError('PyWCS was not found.  It may not be installed '
                              'or it may not be on your PYTHONPATH.\n'
                              'drizzlepac requires pywcs 1.4 or later.')

    # Add/remove macros and compile args based on the build type
    define_macros = []
    undef_macros = []
    extra_compile_args = []

    if BUILD.lower() == 'debug':
        define_macros.append(('DEBUG', None))
        undef_macros.append('NDEBUG')
        if not (sys.platform.startswith('sun') or sys.platform == 'win32'):
            extra_compile_args.extend(['-fno-inline', '-O0', '-g'])
    elif BUILD.lower() == 'profile':
        define_macros.append(('NDEBUG', None))
        undef_macros.append('DEBUG')
        if not (sys.platform.startswith('sun') or sys.platform == 'win32'):
            extra_compile_args.extend(['-O3', '-g'])
    elif BUILD.lower() == 'release':
        define_macros.append(('NDEBUG', None))
        undef_macros.append('DEBUG')
    else:
        raise ValueError("BUILD should be one of 'debug', 'profile', or "
                         "'release'; got %s" % BUILD)

    for idx, m in enumerate(define_macros):
        if m[1] is not None:
            define_macros[idx] = '%s = %s' % m
        else:
            define_macros[idx] = m[0]

    ext_opts = [('include_dirs', pywcsincludes),
                ('define_macros', define_macros),
                ('undef_macros', undef_macros),
                ('extra_compile_args', extra_compile_args)]

    ext = config['extension=drizzlepac.cdriz']
    for opt, value in ext_opts:
        if opt in ext:
            ext[opt] += '\n' + '\n'.join(value)
        else:
            ext[opt] = '\n'.join(value)
