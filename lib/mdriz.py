#!/usr/bin/env python

import sys, getopt
from BigBlackBox import MultiDrizzle,util,wcs_functions
from BigBlackBox import __version__
    

#-------------------------------------------------------------------------------
# a main program for running MultiDrizzle from the command line


def main() :

    import getopt

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
        if opt == 'h':
            help = True
        elif opt == 'g':
            editpars = True
        elif opt == '?':
            long_help = True

    if long_help:
        print BigBlackBox.getHelpAsString()

    if help:
        print 'Syntax: mdriz.py [-h|-g] input ...'
        print '  Options: -h: help'
        print '           -g: edit parameters with GUI'
        
    if len(args) < 1 and not editpars:
        print 'No input specified... Use -h for help.'
    else:
        MultiDrizzle(editpars=editpars,input_dict=args)
    
if __name__=='__main__':
    main()
    