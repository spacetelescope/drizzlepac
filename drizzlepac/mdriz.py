#!/usr/bin/env python
"""
Main program for running MultiDrizzle from the command line.

:Authors: Warren Hack

:License: :doc:`/LICENSE`

"""
import getopt
import sys
from drizzlepac.astrodrizzle import AstroDrizzle
from drizzlepac import __version__


ruler = '-' * 80
print('{0:s}\nDEPRECATION WARNING:\n'
      '    This module (mdriz) will be removed in a future release.\n'
      '    Please execute "runastrodriz" instead.\n{0:s}\n'.format(ruler), file=sys.stderr)

#-------------------------------------------------------------------------------
# a main program for running MultiDrizzle from the command line


def main() :

    try:
        optlist, args = getopt.getopt(sys.argv[1:], '?hg')
    except getopt.error as e:
        print(str(e))
        print(__doc__)
        print("\t", __version__)

    # read options
    help = False
    editpars = False # deprecated parameter left over from TEAL
    long_help = False
    for opt, value in optlist:
        if opt == '-h':
            help = True
        elif opt == '-g':
            editpars = True
            raise ValueError('The -g option and TEAL GUI are no longer supported.')
        elif opt == '-?':
            long_help = True

    if long_help:
        print(drizzlepac.getHelpAsString())

    if help:
        print('Syntax: mdriz.py -[h|g|?] [name=value,...]')
        print('  Options: -h: help')
        print('           -g: deprecated parameter left over from TEAL')
        print('           -?: extended help message')
        print('  Parameters should be given as "name=value" pairs for all parameters')
        print('      understood by MultiDrizzle. These values will ALWAYS override')
        print('      any parameter value set in the TEAL if TEAL is used at all.')

    if len(args) < 1 and not editpars:
        print('No input specified... Use -h for help.')
    else:
        input_dict = {}
        for a in args:
            avar,aval = a.split('=')
            input_dict[avar] = aval
        AstroDrizzle(editpars=editpars,**input_dict)

if __name__ == '__main__':
    main()
