"""
`AstroDrizzle` - Python implementation of `MultiDrizzle`

`AstroDrizzle` automates the process of aligning images in an output frame,
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

A full description of this process can be found
in the `DrizzlePac Handbook <http://drizzlepac.stsci.edu>`_\ .

:Output: The primary output from this task is the distortion-corrected, \
cosmic-ray cleaned, and combined image as a FITS file.

This task requires numerous user-settable parameters to control the primary
aspects of each of the processing steps.

:Authors: Warren Hack

:License: `<http://www.stsci.edu/resources/software_hardware/pyraf/LICENSE>`_

"""
from __future__ import absolute_import, division, print_function  # confidence high


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

    if configobj is None or configobj == 'defaults':
        configobj = teal.load(__taskname__)

    if 'updatewcs' in input_dict: # user trying to explicitly turn on updatewcs
        configobj['updatewcs'] = input_dict['updatewcs']
        del input_dict['updatewcs']

    # If called from interactive user-interface, configObj will not be
    # defined yet, so get defaults using EPAR/TEAL.
    #
    # Also insure that the input_dict (user-specified values) are folded in
    # with a fully populated configObj instance.
    try:
        configObj = util.getDefaultConfigObj(__taskname__, configobj,
                                             input_dict,
                                             loadOnly=(not editpars))
        print("INPUT_DICT: {}".format(input_dict))
        # If user specifies optional parameter for final_wcs specification in input_dict,
        #    insure that the final_wcs step gets turned on
        util.applyUserPars_steps(configObj, input_dict, step='3a')
        util.applyUserPars_steps(configObj, input_dict, step='7a')

    except ValueError:
        print("Problem with input parameters. Quitting...", file=sys.stderr)
        return

    if not configObj:
        return

    configObj['mdriztab'] = mdriztab
    # If 'editpars' was set to True, util.getDefaultConfigObj() will have
    # already called 'run()'.
    if not editpars:
        run(configObj, wcsmap=wcsmap)

#
#### Interfaces used by TEAL
#


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
        print(textutil.textbox(
            'ERROR:\nNo valid input files found!   Please restart the task '
            'and check the value for the "input" parameter.'), file=sys.stderr)
        def_logname = None
        return

    stateObj = configobj['STATE OF INPUT FILES']
    procSteps = util.ProcSteps()

    print('AstroDrizzle Version %s(%s) started at: %s\n' %
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
            print(textutil.textbox(errmsg,width=65))
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

#       _dbg_dump_virtual_outputs(imgObjList)

        #drizzle to separate images
        adrizzle.drizSeparate(imgObjList, outwcs, configobj, wcsmap=wcsmap,
                              procSteps=procSteps)

#       _dbg_dump_virtual_outputs(imgObjList)

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

        print()
        print(' '.join(['AstroDrizzle Version', __version__,
                        'is finished processing at ',
                        util._ptime()[0]]) + '!\n')
    except:
        print(textutil.textbox(
            'ERROR:\nAstroDrizzle Version %s encountered a problem!  '
            'Processing terminated at %s.' %
            (__version__, util._ptime()[0])), file=sys.stderr)
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
    helpstr = getHelpAsString(docstring=True, show_ver = True)
    if file is None:
        print(helpstr)
    else:
        if os.path.exists(file): os.remove(file)
        f = open(file, mode = 'w')
        f.write(helpstr)
        f.close()


def getHelpAsString(docstring = False, show_ver = True):
    """
    return useful help from a file in the script directory called
    __taskname__.help

    """
    install_dir = os.path.dirname(__file__)
    taskname = util.base_taskname(__taskname__, '')
    htmlfile = os.path.join(install_dir, 'htmlhelp', taskname + '.html')
    helpfile = os.path.join(install_dir, taskname + '.help')

    if docstring or (not docstring and not os.path.exists(htmlfile)):
        if show_ver:
            helpString = os.linesep + \
                ' '.join([__taskname__, 'Version', __version__,
                ' updated on ', __vdate__]) + 2*os.linesep
        else:
            helpString = ''
        if os.path.exists(helpfile):
            helpString += teal.getHelpFileAsString(taskname, __file__)
        else:
            if __doc__ is not None:
                helpString += __doc__ + os.linesep
    else:
        helpString = 'file://' + htmlfile

    return helpString


_fidx = 0
def _dbg_dump_virtual_outputs(imgObjList):
    """ dump some helpful information.  strictly for debugging """
    global _fidx
    tag = 'virtual'
    log.info((tag+'  ')*7)
    for iii in imgObjList:
        log.info('-'*80)
        log.info(tag+'  orig nm: '+iii._original_file_name)
        log.info(tag+'  names.data: '+str(iii.outputNames["data"]))
        log.info(tag+'  names.orig: '+str(iii.outputNames["origFilename"]))
        log.info(tag+'  id: '+str(id(iii)))
        log.info(tag+'  in.mem: '+str(iii.inmemory))
        log.info(tag+'  vo items...')
        for vok in sorted(iii.virtualOutputs.keys()):
            FITSOBJ = iii.virtualOutputs[vok]
            log.info(tag+': '+str(vok)+' = '+str(FITSOBJ))
            if vok.endswith('.fits'):
              if not hasattr(FITSOBJ, 'data'):  FITSOBJ = FITSOBJ[0] # list of PrimaryHDU ?
              if not hasattr(FITSOBJ, 'data'):  FITSOBJ = FITSOBJ[0] # was list of HDUList ?
              dbgname = 'DEBUG_%02d_'%(_fidx,); dbgname+=os.path.basename(vok); _fidx+=1
              FITSOBJ.writeto(dbgname)
              log.info(tag+'  wrote: '+dbgname)
              log.info('\n'+vok)
              if hasattr(FITSOBJ, 'data'):
                log.info(str(FITSOBJ._summary()))
                log.info('min and max are: '+str( (FITSOBJ.data.min(),
                                                   FITSOBJ.data.max()) ))
                log.info('avg and sum are: '+str( (FITSOBJ.data.mean(),
                                                   FITSOBJ.data.sum()) ))
#               log.info(str(FITSOBJ.data)[:75])
              else:
                log.info(vok+' has no .data attr')
                log.info(str(type(FITSOBJ)))
              log.info(vok+'\n')
    log.info('-'*80)

AstroDrizzle.__doc__ = getHelpAsString(docstring = True, show_ver = False)
