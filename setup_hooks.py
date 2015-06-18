
'''
To compile another package that wants to link with the C code in astropy.wcs:

[build_ext]
pre-hook.wcs = drizzlepac.hooks.setup

[extension=drizzlepac.cdriz]
include_dirs = numpy wcs
libraries = wcs m

'''


import os
import os.path
from distutils import log

def wcs_include(command_obj) :

    # warning if we somehow get called when not building a c extension
    command_name = command_obj.get_command_name()
    if command_name != 'build_ext':
        log.warn('%s is meant to be used with the build_ext command only; '
                'it is not for use with the %s command.' %
                (__name__, command_name))

    # get information from astropy.wcs
    from astropy import wcs
    include_path = wcs.get_include()
    includes = [os.path.join(include_path, 'astropy_wcs'),
                os.path.join(include_path, 'wcslib')]

    for extension in command_obj.extensions:
            extension.include_dirs.extend(includes)


#####
#
# add includes for compiling tests that use fctx.  Use
#   [build_ext]
#   pre-hook.pdk_fctx = setup_hooks.pdk_fctx
# and list the directory
#   pdk_fctx
# in your includes
#

def pdk_fctx(command_obj) :

    command_name = command_obj.get_command_name()
    if command_name != 'build_ext':
        log.warn('%s is meant to be used with the build_ext command only; '
                 'it is not for use with the %s command.' %
                 (__name__, command_name))
    try :
        import pandokia.runners as x
    except ImportError :
        log.warn('pandokia import failed - no fctx include files')
        return

    # knowledge about internals of pandokia - bummer, but pandokia does
    # not report this any other way except running an external program.
    includes = [ os.path.join( os.path.dirname(x.__file__), 'maker' ) ]
    #includes = [numpy.get_numarray_include(), numpy.get_include()]
    for extension in command_obj.extensions:
        if 'pdk_fctx' not in extension.include_dirs:
            continue
        idx = extension.include_dirs.index('pdk_fctx')
        for inc in includes:
            extension.include_dirs.insert(idx, inc)
        extension.include_dirs.remove('pdk_fctx')

