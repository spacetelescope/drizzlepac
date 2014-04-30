
'''
To compile another package that wants to link with the C code in astropy.wcs:

[build_ext]
pre-hook.wcs = drizzlepac.hooks.setup

[extension=drizzlepac.cdriz]
include_dirs = numpy wcs
libraries = wcs m

'''


import os

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
