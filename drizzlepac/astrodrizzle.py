"""
``AstroDrizzle`` - Python implementation of ``MultiDrizzle``

``AstroDrizzle`` automates the process of aligning images in an output frame,
identifying cosmic-rays, removing distortion, and then combining the images
after removing the identified cosmic-rays.

This process involves a number of steps, such as:
  *  Processing the input images and input parameters
  *  Creating a static mask
  *  Performing sky subtraction
  *  Drizzling onto separate output images
  *  Creating the median image
  *  Blotting the median image
  *  Identifying and flagging cosmic-rays
  *  Final combination
  *  Cleaning-up of temporary files (when applicable)

A full description of this process can be found
in the `DrizzlePac Handbook <http://drizzlepac.stsci.edu>`_\ .

The primary output from this task is the distortion-corrected,
cosmic-ray cleaned, and combined image as a FITS file.

This task requires numerous user-settable parameters to control the primary
aspects of each of the processing steps.

:Authors: Warren Hack

:License: :doc:`/LICENSE`

"""
import os
import sys
import logging

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
from . import __version__


__taskname__ = "astrodrizzle"


# Pointer to the included Python class for WCS-based coordinate transformations
PYTHON_WCSMAP = wcs_functions.WCSMap

log = logutil.create_logger(__name__, level=logutil.logging.NOTSET)


def AstroDrizzle(input=None, mdriztab=False, editpars=False, configobj=None,
                 wcsmap=None, **input_dict):
    """ AstroDrizzle command-line interface """
    # Support input of filenames from command-line without a parameter name
    # then copy this into input_dict for merging with TEAL ConfigObj
    # parameters.

    # Load any user-specified configobj
    if isinstance(configobj, (str, bytes)):
        if configobj == 'defaults':
            # load "TEAL"-defaults (from ~/.teal/):
            configobj = teal.load(__taskname__)
        else:
            if not os.path.exists(configobj):
                raise RuntimeError('Cannot find .cfg file: ' + configobj)
            configobj = teal.load(configobj, strict=False)
    elif configobj is None:
        # load 'astrodrizzle' parameter defaults as described in the docs:
        configobj = teal.load(__taskname__, defaults=True)

    if input and not util.is_blank(input):
        input_dict['input'] = input
    elif configobj is None:
        raise TypeError("AstroDrizzle() needs either 'input' or "
                        "'configobj' arguments")

    if 'updatewcs' in input_dict:  # user trying to explicitly turn on updatewcs
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
        log.debug('')
        log.debug("INPUT_DICT:")
        util.print_cfg(input_dict, log.debug)
        log.debug('')
        # If user specifies optional parameter for final_wcs specification in input_dict,
        #    insure that the final_wcs step gets turned on
        util.applyUserPars_steps(configObj, input_dict, step='3a')
        util.applyUserPars_steps(configObj, input_dict, step='7a')

    except ValueError:
        log.error("Problem with input parameters. Quitting...", file=sys.stderr)
        return

    # add flag to configObj to indicate whether or not to use mdriztab
    configObj['mdriztab'] = mdriztab

    run(configObj, wcsmap=wcsmap, input_dict=input_dict)

##############################
#   Interfaces used by TEAL  #
##############################
@util.with_logging
def run(configobj, wcsmap=None, input_dict=None):
    """
    Initial example by Nadia ran MD with configobj EPAR using:
    It can be run in one of two ways:

        from stsci.tools import teal

        1. Passing a config object to teal

        teal.teal('drizzlepac/pars/astrodrizzle.cfg')


        2. Passing a task  name:

        teal.teal('astrodrizzle')

        The example config files are in drizzlepac/pars

    input_dict is a dictionary of user-specified parameters

    """
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
        log.error(textutil.textbox(
            "No valid input files found!   Please restart the task "
            "and check the value for the 'input' parameter."), file=sys.stderr)
        def_logname = None
        return

    # Build name of output trailer file
    logging_handlers = logging.getLogger().handlers
    log_name = [lh.name for lh in logging_handlers if lh.level > 0]
    logfile = log_name[0] if log_name else "{}.tra".format(def_logname)
    log.info("AstroDrizzle log file: {}".format(logfile))

    clean = configobj['STATE OF INPUT FILES']['clean']
    procSteps = util.ProcSteps()

    log.info("AstroDrizzle Version {:s} started at: {:s}\n"
          .format(__version__, util._ptime()[0]))
    util.print_pkg_versions(log=log)

    log.debug('')
    log.debug(
        "==== AstroDrizzle was invoked with the following parameters: ===="
    )
    log.debug('')
    util.print_cfg(configobj, log.debug)

    try:
        # Define list of imageObject instances and output WCSObject instance
        # based on input paramters
        imgObjList = None
        procSteps.addStep('Initialization')
        imgObjList, outwcs = processInput.setCommonInput(
            configobj, overwrite_dict=input_dict
        )
        procSteps.endStep("Initialization")

        if imgObjList is None or not imgObjList:
            errmsg = "No valid images found for processing!\n"
            errmsg += "Check log file for full details.\n"
            errmsg += "Exiting AstroDrizzle now..."
            log.error(textutil.textbox(errmsg, width=65))
            log.error(textutil.textbox(
                '\nAstroDrizzle Version {:s} encountered a problem!  '
                'Processing terminated at {:s}.'
                .format(__version__, util._ptime()[0])), file=sys.stderr)
            return

        log.info("USER INPUT PARAMETERS common to all Processing Steps:")
        util.printParams(configobj, log=log)

        step_name_single = util.getSectionName(configobj, adrizzle.STEP_NUM_SINGLE)
        do_single = configobj[step_name_single]["driz_separate"]

        step_name_crrej = util.getSectionName(configobj, drizCR.STEP_NUM)
        do_crrej = configobj[step_name_crrej]["driz_cr"]
        skip_crrej = False

        step_name_median = util.getSectionName(configobj, createMedian.STEP_NUM)
        do_median = configobj[step_name_median]["median"]
        skip_median = False

        step_name_blot = util.getSectionName(configobj, ablot.STEP_NUM)
        do_blot = configobj[step_name_blot]["blot"]
        skip_blot = False

        if len(imgObjList) > 1:
            if do_crrej and not do_blot:
                log.warning(
                    "Turning blot step on as it is required by 'driz_cr'."
                )
                configobj[step_name_blot]["blot"] = True
                do_blot = True

            if do_blot and not do_median:
                log.warning(
                    "Turning median step on as it is required by 'blot'."
                )
                configobj[step_name_median]["median"] = True
                do_median = True

            if do_median and not do_single:
                log.warning(
                    "Turning single drizzle step on as it is required by "
                    "'median'."
                )
                configobj[step_name_single]["driz_separate"] = True
                do_single = True

        else:
            if do_crrej:
                log.warning(
                    "Turning CR rejection step off as it requires two or more "
                    "input images."
                )
                configobj[step_name_crrej]["driz_cr"] = False
                do_crrej = False
                skip_crrej = True

            if do_blot:
                log.warning(
                    "Turning blot step off as it is requires two or more "
                    "input images."
                )
                configobj[step_name_blot]["blot"] = False
                do_blot = False
                skip_blot = True

            if do_median:
                log.warning(
                    "Turning median step off as it requires two or more "
                    "input images."
                )
                configobj[step_name_median]["median"] = False
                do_median = False
                skip_median = True

        # Call rest of MD steps...
        # create static masks for each image
        staticMask.createStaticMask(imgObjList, configobj,
                                    procSteps=procSteps)

        # subtract the sky
        sky.subtractSky(imgObjList, configobj, procSteps=procSteps)

        #       _dbg_dump_virtual_outputs(imgObjList)

        # drizzle to separate images
        adrizzle.drizSeparate(imgObjList, outwcs, configobj, wcsmap=wcsmap,
                              logfile=logfile,
                              procSteps=procSteps)

        #       _dbg_dump_virtual_outputs(imgObjList)

        # create the median images from the driz sep images
        try:
            createMedian.createMedian(
                imgObjList,
                configobj,
                procSteps=procSteps
            )

            if skip_median:
                procSteps.endStep(createMedian.PROCSTEPS_NAME, reason="skipped")
            elif not do_median:
                procSteps.endStep(createMedian.PROCSTEPS_NAME, reason="off")

        except util.StepAbortedError as e:
            if str(e).startswith("Rejecting all pixels"):
                log.warning(
                    "Create median step was aborted due the following error:"
                )
                log.warning(
                    f"ERROR: {str(e)}"
                )

                if do_blot:
                    log.warning(
                        "Turning blot step off due to aborted median step."
                    )
                    configobj[step_name_blot]['median'] = False
                    skip_blot = True
                    do_blot = False

                if do_crrej:
                    log.warning(
                        "Turning CR rejection step off due to aborted "
                        "median step."
                    )
                    configobj[step_name_crrej]['driz_cr'] = False
                    skip_crrej = True
                    do_crrej = False
            else:
                raise e

        # blot the images back to the original reference frame
        ablot.runBlot(imgObjList, outwcs, configobj, wcsmap=wcsmap,
                      procSteps=procSteps)
        if skip_blot:
            procSteps.endStep(ablot.PROCSTEPS_NAME, reason="skipped")
        elif not do_blot:
            procSteps.endStep(ablot.PROCSTEPS_NAME, reason="off")

        # look for cosmic rays
        drizCR.rundrizCR(imgObjList, configobj, procSteps=procSteps)
        if skip_crrej:
            procSteps.endStep(drizCR.PROCSTEPS_NAME, reason="skipped")
        elif not do_crrej:
            procSteps.endStep(drizCR.PROCSTEPS_NAME, reason="off")

        # Make your final drizzled image
        adrizzle.drizFinal(imgObjList, outwcs, configobj, wcsmap=wcsmap,
                           logfile=logfile,
                           procSteps=procSteps)

        log.info("AstroDrizzle Version {:s} is finished processing at {:s}.\n"
              .format(__version__, util._ptime()[0]))

    except Exception:
        clean = False
        log.error(textutil.textbox(
            "AstroDrizzle Version {:s} encountered a problem!  "
            "Processing terminated at {:s}."
            .format(__version__, util._ptime()[0])), file=sys.stderr)
        procSteps.endStep(None, reason="aborted")
        raise

    finally:
        procSteps.reportTimes()
        if imgObjList:
            for image in imgObjList:
                if clean:
                    image.clean()
                image.close()
            del imgObjList
            del outwcs


AstroDrizzle.__doc__ = util._def_help_functions(
    locals(), module_file=__file__, task_name=__taskname__, module_doc=__doc__
)

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
                if not hasattr(FITSOBJ, 'data'):
                    FITSOBJ = FITSOBJ[0] # list of PrimaryHDU ?
                if not hasattr(FITSOBJ, 'data'):
                    FITSOBJ = FITSOBJ[0] # was list of HDUList ?
                dbgname = 'DEBUG_%02d_'%(_fidx,)
                dbgname+=os.path.basename(vok)
                _fidx+=1
                FITSOBJ.writeto(dbgname)
                log.info(tag+'  wrote: '+dbgname)
                log.info('\n'+vok)
                if hasattr(FITSOBJ, 'data'):
                    log.info(str(FITSOBJ._summary()))
                    log.info('min and max are: '+str( (FITSOBJ.data.min(),
                                                       FITSOBJ.data.max()) ))
                    log.info('avg and sum are: '+str( (FITSOBJ.data.mean(),
                                                       FITSOBJ.data.sum()) ))
#                    log.info(str(FITSOBJ.data)[:75])
                else:
                    log.info(vok+' has no .data attr')
                    log.info(str(type(FITSOBJ)))
                log.info(vok+'\n')
    log.info('-'*80)
