
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
	includes = [wcs.get_include()]

	# each place where the include directory is named exactly "wcs",
	# replace it with the include directories computed above
	for extension in command_obj.extensions:
		if 'wcs' not in extension.include_dirs:
			continue
		idx = extension.include_dirs.index('wcs')
		for inc in includes:
			extension.include_dirs.insert(idx, inc)
		extension.include_dirs.remove('wcs')





