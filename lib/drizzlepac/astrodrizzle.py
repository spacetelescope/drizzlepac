""" AstroDrizzle - Python implementation of MultiDrizzle

AstroDrizzle automates the process of aligning images in an output frame,
identifying cosmic-rays, removing distortion, and then combining the images
after removing the identified cosmic-rays.

This process involves a number of steps, namely:
  1.  Processing the input images and input parameters
  2.  Creating a static mask
  3.  Performing sky subtraction
  4.  Drizzling onto separate output images
  5.  Creating the median image
  6.  Blotting the median image
  7.  Identifying and flagging cosmic-rays
  8.  Final combination

A full description of this process can be found in the DrizzlePac Handbook
available online at:

http://drizzlepac.stsci.edu

**Output**: The primary output from this task is the distortion-corrected,
cosmic-ray cleaned, and combined image as a FITS file.

This task requires numerous user-settable parameters to control the primary
aspects of each of the processing steps.

"""

from __future__ import division  # confidence high


import os
import sys

from stsci.tools import teal, logutil, textutil

from . import adrizzle
from . import ablot
from . import createMedian
from . import drizCR
from . import processInput
from . import sky
from . import staticMask
from . import util
from . import wcs_functions
from .version import *


__taskname__ = "astrodrizzle"

# Definitions for flags on when to raise an EXCEPTION
RAISE = 1
DO_NOT_RAISE = 0

# Pointer to the included Python class for WCS-based coordinate transformations
PYTHON_WCSMAP = wcs_functions.WCSMap

log = logutil.create_logger(__name__)


#
#### Interactive user interface (functional form)
#
def AstroDrizzle(input=None, mdriztab=False, editpars=False, configobj=None,
                 wcsmap=None, **input_dict):
    """ AstroDrizzle command-line interface
    """
    # support input of filenames from command-line without a parameter name
    # then copy this into input_dict for merging with TEAL ConfigObj parameters
    if input_dict is None:
        input_dict = {}

    if input is None and configobj is None:
        raise TypeError('AstroDrizzle() needs either "input" or "configobj" arg')

    if input and not util.is_blank(input):
        input_dict['input'] = input

    # input_dict['mdriztab'] = mdriztab

    # Load any user-specified configobj
    if isinstance(configobj, str) and configobj != 'defaults':
        if not os.path.exists(configobj):
            raise RuntimeError('Cannot find .cfg file: '+configobj)
        configobj = teal.load(configobj, strict=False)

    # If called from interactive user-interface, configObj will not be
    # defined yet, so get defaults using EPAR/TEAL.
    #
    # Also insure that the input_dict (user-specified values) are folded in
    # with a fully populated configObj instance.
    try:
        configObj = util.getDefaultConfigObj(__taskname__, configobj,
                                             input_dict,
                                             loadOnly=(not editpars))
    except ValueError:
        print >> sys.stderr, "Problem with input parameters. Quitting..."
        return

    if configObj is None:
        return

    configObj['mdriztab'] = mdriztab
    # If 'editpars' was set to True, util.getDefaultConfigObj() will have
    # already called 'run()'.
    if editpars == False:
        run(configObj, wcsmap=wcsmap)

#
#### Interfaces used by TEAL
#


def getHelpAsString(docstring=False):
    """
    return useful help from a file in the script directory called
    __taskname__.help
    """

    install_dir = os.path.dirname(__file__)
    htmlfile = os.path.join(install_dir, 'htmlhelp', __taskname__ + '.html')
    helpfile = os.path.join(install_dir, __taskname__ + '.help')
    if docstring or (not docstring and not os.path.exists(htmlfile)):
        helpString = ' '.join([__taskname__, 'Version', __version__,
                               ' updated on ', __vdate__]) + '\n\n'
        if os.path.exists(helpfile):
            helpString += teal.getHelpFileAsString(__taskname__, __file__)
    else:
        helpString = 'file://' + htmlfile

    return helpString


AstroDrizzle.__doc__ = getHelpAsString(docstring=True)


@util.with_logging
def run(configobj, wcsmap=None):
    """
    Initial example by Nadia ran MD with configobj EPAR using:
    It can be run in one of two ways:

        from stsci.tools import teal

        1. Passing a config object to teal

        teal.teal('drizzlepac/pars/astrodrizzle.cfg')


        2. Passing a task  name:

        teal.teal('astrodrizzle')

        The example config files are in drizzlepac/pars

    """
    raise_status = RAISE
    #
    # turn on logging, redirecting stdout/stderr messages to a log file
    # while also printing them out to stdout as well
    # also, initialize timing of processing steps
    #
    # We need to define a default logfile name from the user's parameters
    input_list, output, ivmlist, odict = \
            processInput.processFilenames(configobj['input'])

    if output is not None:
        def_logname = output
    elif len(input_list) > 0:
        def_logname = input_list[0]
    else:
        print >> sys.stderr, textutil.textbox(
            'ERROR:\nNo valid input files found!   Please restart the task '
            'and check the value for the "input" parameter.')
        def_logname = None
        return

    stateObj = configobj['STATE OF INPUT FILES']
    procSteps = util.ProcSteps()

    print ('AstroDrizzle Version %s(%s) started at: %s\n' %
           (__version__, __vdate__, util._ptime()[0]))
    util.print_pkg_versions(log=log)

    #try:
    try:
        # Define list of imageObject instances and output WCSObject instance
        # based on input paramters
        procSteps.addStep('Initialization')
        imgObjList = None
        imgObjList, outwcs = processInput.setCommonInput(configobj)
        procSteps.endStep('Initialization')

        if not imgObjList:
            errmsg = "No valid images found for processing!\n"
            errmsg += "Check log file for full details.\n"
            errmsg += "Exiting AstroDrizzle now..."
            print textutil.textbox(errmsg,width=65)
            raise_status = DO_NOT_RAISE
            raise ValueError

        log.info("USER INPUT PARAMETERS common to all Processing Steps:")
        util.printParams(configobj, log=log)

        # Call rest of MD steps...
        #create static masks for each image
        staticMask.createStaticMask(imgObjList, configobj,
                                    procSteps=procSteps)

        #subtract the sky
        sky.subtractSky(imgObjList, configobj, procSteps=procSteps)

        #drizzle to separate images
        adrizzle.drizSeparate(imgObjList, outwcs, configobj, wcsmap=wcsmap,
                              procSteps=procSteps)

        #create the median images from the driz sep images
        createMedian.createMedian(imgObjList, configobj,
                                  procSteps=procSteps)

        #blot the images back to the original reference frame
        ablot.runBlot(imgObjList, outwcs, configobj, wcsmap=wcsmap,
                      procSteps=procSteps)

        #look for cosmic rays
        drizCR.rundrizCR(imgObjList, configobj,
                         procSteps=procSteps)

        #Make your final drizzled image
        adrizzle.drizFinal(imgObjList, outwcs, configobj, wcsmap=wcsmap,
                           procSteps=procSteps)

        print
        print ' '.join(['AstroDrizzle Version', __version__,
                        'is finished processing at ',
                        util._ptime()[0]]) + '!\n'
    except:
        print >> sys.stderr, textutil.textbox(
            'ERROR:\nAstroDrizzle Version %s encountered a problem!  '
            'Processing terminated at %s.' %
            (__version__, util._ptime()[0]))
        procSteps.reportTimes()
        if imgObjList:
            for image in imgObjList:
                image.close()
            del imgObjList
            del outwcs

        # Raise an exception ONLY if requested...
        if raise_status == RAISE:
            raise
        else:
            return

    procSteps.reportTimes()

    if imgObjList:
        for image in imgObjList:
            if stateObj['clean']:
                image.clean()
            image.close()

        del imgObjList
        del outwcs


def help(file=None):
    """
    Print out syntax help for running astrodrizzle

    Parameters
    ----------
    file : str (Default = None)
        If given, write out help to the filename specified by this parameter
        Any previously existing file with this name will be deleted before
        writing out the help.
    """
    helpstr = getHelpAsString(docstring=True)
    if file is None:
        print helpstr
    else:
        if os.path.exists(file): os.remove(file)
        f = open(file,mode='w')
        f.write(helpstr)
        f.close()
