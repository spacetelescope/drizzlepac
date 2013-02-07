#!/usr/bin/env python
from __future__ import division  # confidence high

import sys
from drizzlepac import AstroDrizzle
from drizzlepac.version import __version__


#-------------------------------------------------------------------------------
# a main program for running MultiDrizzle from the command line


def main() :

    try:
        optlist, args = getopt.getopt(sys.argv[1:], '?hg')
    except getopt.error, e:
        print str(e)
        print __doc__
        print "\t", __version__

    # read options
    help = False
    editpars = False
    long_help = False
    for opt, value in optlist:
        if opt == '-h':
            help = True
        elif opt == '-g':
            editpars = True
        elif opt == '-?':
            long_help = True

    if long_help:
        print drizzlepac.getHelpAsString()

    if help:
        print 'Syntax: mdriz.py -[h|g|?] [name=value,...]'
        print '  Options: -h: help'
        print '           -g: edit parameters with TEAL'
        print '           -?: extended help message'
        print '  Parameters should be given as "name=value" pairs for all parameters'
        print '      understood by MultiDrizzle. These values will ALWAYS override'
        print '      any parameter value set in the TEAL if TEAL is used at all.'

    if len(args) < 1 and not editpars:
        print 'No input specified... Use -h for help.'
    else:
        input_dict = {}
        for a in args:
            avar,aval = a.split('=')
            input_dict[avar] = aval
        AstroDrizzle(editpars=editpars,**input_dict)

if __name__ == '__main__':
    main()
