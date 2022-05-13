#!/usr/bin/env python

""" runastrodriz.py - Module to control operation of astrodrizzle to
        remove distortion and combine HST images in the pipeline.

:License: :doc:`LICENSE`

USAGE: runastrodriz.py [-fhdaibng] inputFilename [newpath]

Alternative USAGE:
    python
    from acstools import runastrodriz
    runastrodriz.process(inputFilename,force=False,newpath=None)

If a value has been provided for the newpath parameter, all processing will be
performed in that directory/ramdisk.  The steps involved are:
   - create a temporary directory under that directory named after the input file
   - copy all files related to the input to that new directory
   - change to that new directory and run astrodrizzle
   - change back to original directory
   - move (not copy) ALL files from temp directory to original directory
   - delete temp sub-directory

The '-b' option will run this task in BASIC mode without creating headerlets
for each input image.

The '-n' option allows the user to specify the number of cores to be used in
running AstroDrizzle.

The '-g' option allows the user to TURN OFF alignment of the images to an external
astrometric catalog, such as GAIA, as accessible through the MAST interface.

Additional control over whether or not to attempt to align to an external
astrometric catalog, such as GAIA, is provided through the use of the
environment variables:

    - ASTROMETRY_COMPUTE_APOSTERIORI : Turn on/off alignment step.
      This environment variable will ALWAYS override any setting of the '-g' switch.
      Values (case-insensitive) can be 'on', 'off', 'yes', 'no'.

    - ASTROMETRY_APPLY_APRIORI : Replaces/resets ASTROMETRY_STEP_CONTROL
      variable used by `stwcs.updatewcs` to control whether or not a priori WCS
      solutions from the astrometry database should be applied to the data.
      If this is set, it will override any value set in the old variable.
      Values (case-insensitive) can be 'on','off','yes','no'.

Additionally, the output products can be evaluated to determine the quality of
the alignment and output data through the use of the environment variable:

    - PIPELINE_QUALITY_TESTING : Turn on quality assessment processing.
      This environment variable, if found with any value, will turn on
      processing to generate a JSON file which contains the results of
      evaluating the quality of the generated products.

This processing can also insure that the IDCTAB reference file in the FLT/FLC
files are as up-to-date as the IDCTAB specified in the RAW file using:

    - PIPELINE_RESET_IDCTAB : Turn on automatic reset of IDCTAB in FLT/FLC
      files so that they are identical to those in the RAW files.  This
      does nothing if they are already in sync.


*** INITIAL VERSION
W.J. Hack  12 Aug 2011: Initial version based on Version 1.2.0 of
                        STSDAS$pkg/hst_calib/wfc3/runwf3driz.py
W.J. Hack  27 Jun 2012: Implement support to process in different directory

W.J. Hack  24 Aug 2012: Provided interface for in-memory option

W.J. Hack  26 Nov 2012: Option to write out headerlets added and debugged

W.J. Hack  18 Oct 2019: Impelemented multi-stage alignment with verification

"""
# Import standard Python modules
import glob
import os
import shutil
import sys
import time
import logging
import json
import traceback
import stat
import errno
from collections import OrderedDict
import datetime
import fnmatch
try:
    from psutil import Process
except ImportError:
    Process = None

# THIRD-PARTY
import numpy as np
from astropy.io import fits

import stwcs
from stwcs import wcsutil
from stwcs.wcsutil import HSTWCS
from stwcs import updatewcs
from stwcs.wcsutil import headerlet, altwcs

from stsci.tools import fileutil, asnutil
import tweakwcs

import drizzlepac
from drizzlepac import processInput  # used for creating new ASNs for _flc inputs

from drizzlepac import align
from drizzlepac import resetbits
from drizzlepac.haputils import astrometric_utils as amutils
from drizzlepac.haputils import cell_utils
from drizzlepac.haputils import processing_utils
from drizzlepac import util
from drizzlepac import mdzhandler
from drizzlepac import updatehdr
from drizzlepac.haputils import quality_analysis as qa
from drizzlepac import wcs_functions
from . import __version__


__taskname__ = "runastrodriz"

# Local variables

# Implement WIN specific check
RM_LOGFILES = False if sys.platform.startswith('win') else True

# Define parameters which need to be set specifically for
#    pipeline use of astrodrizzle
# The parameter for resetbits resets DQ values of :
#  - 4096: pixels previously flagged as cosmic-rays.
#
PIPELINE_PARS = {'mdriztab': True,
                 'in_memory': True,
                 'stepsize': 10,
                 'output': '',
                 'preserve': False,
                 'clean': False,
                 'resetbits': 4096}

# Values of good_bits are set to treat these DQ bit values as 'good':
#  - 1024: sink pixel (ACS), charge trap (WFC3/UVIS)
#  -  256: saturated pixel (ACS), full-well saturation (WFC3)
#  -   64: warm pixel (ACS, WFC3)
#  -   16: hot pixel (ACS, WFC3)
#  -  512: bad reference file pixel (ACS), bad flat pixel (WFC3)

focus_pars = {"WFC3/IR": {'sigma': 2.0, 'good_bits': 512},
              "WFC3/UVIS": {'sigma': 1.5, 'good_bits': 1360},
              "ACS/WFC": {'sigma': 1.5, 'good_bits': 1360},
              "ACS/SBC": {'sigma': 2.0, 'good_bits': 0},
              "ACS/HRC": {'sigma': 1.5, 'good_bits': 1360},
              "WFPC2/PC": {'sigma': 1.5, 'good_bits': 1360}}
sub_dirs = ['OrIg_files', 'pipeline-default']
valid_alignment_modes = ['apriori', 'aposteriori', 'default-pipeline']
gsc240_date = '2017-10-01'
apriori_priority = ['HSC', 'GSC', '']


# default marker for trailer files
__trlmarker__ = '*** astrodrizzle Processing Version ' + __version__ + '***\n'

envvar_bool_dict = {'off': False, 'on': True, 'no': False, 'yes': True, 'false': False, 'true': True}
envvar_dict = {'off': 'off', 'on': 'on', 'yes': 'on', 'no': 'off', 'true': 'on', 'false': 'off'}

envvar_compute_name = 'ASTROMETRY_COMPUTE_APOSTERIORI'
# Replace ASTROMETRY_STEP_CONTROL with this new related name
envvar_new_apriori_name = "ASTROMETRY_APPLY_APRIORI"
envvar_old_apriori_name = "ASTROMETRY_STEP_CONTROL"
envvar_qa_stats_name = "PIPELINE_QUALITY_TESTING"
envvar_reset_idctab_name = "PIPELINE_RESET_IDCTAB"

# Order of preference for common WCS solutions
wcs_preference = ['IDC_?????????-FIT_REL_GAIA*', 'IDC_?????????-FIT_IMG_GAIA*', 'IDC_?????????-GSC240', 'IDC_?????????']

# History:
# Version 1.0.0 - Derived from v1.2.0 of wfc3.runwf3driz to run astrodrizzle


# Primary user interface
def process(inFile, force=False, newpath=None, num_cores=None, inmemory=True,
            headerlets=True, align_to_gaia=True, force_alignment=False, debug=False):
    """ Run astrodrizzle on input file/ASN table
        using default values for astrodrizzle parameters.
    """
    init_time = time.time()
    trlmsg = "{}: Calibration pipeline processing of {} started.\n".format(init_time, inFile)
    trlmsg += __trlmarker__
    trlmsg += "    drizzlepac version {}".format(drizzlepac.__version__)
    trlmsg += "    tweakwcs version {}".format(tweakwcs.__version__)
    trlmsg += "    stwcs version {}".format(stwcs.__version__)
    trlmsg += "    numpy version {}".format(np.__version__)
    pipeline_pars = PIPELINE_PARS.copy()
    _verify = True  # Switch to control whether to verify alignment or not

    # interpret envvar variable, if specified
    if envvar_compute_name in os.environ:
        val = os.environ[envvar_compute_name].lower()
        if val not in envvar_bool_dict:
            msg = "ERROR: invalid value for {}.".format(envvar_compute_name)
            msg += "  \n    Valid Values: on, off, yes, no, true, false"
            raise ValueError(msg)
        align_to_gaia = envvar_bool_dict[val]

    if envvar_new_apriori_name in os.environ:
        # Reset ASTROMETRY_STEP_CONTROL based on this variable
        # This provides backward-compatibility until ASTROMETRY_STEP_CONTROL
        # gets removed entirely.
        val = os.environ[envvar_new_apriori_name].lower()
        if val not in envvar_dict:
            msg = "ERROR: invalid value for {}.".format(envvar_new_apriori_name)
            msg += "  \n    Valid Values: on, off, yes, no, true, false"
            raise ValueError(msg)

        os.environ[envvar_old_apriori_name] = envvar_dict[val]
    else:
        # Insure os.environ ALWAYS contains an entry for envvar_new_apriori_name
        # and it will default to being 'on'
        if envvar_old_apriori_name in os.environ:
            val = os.environ[envvar_old_apriori_name].lower()
        else:
            val = 'on'
        os.environ[envvar_new_apriori_name] = val

    align_with_apriori = True
    if envvar_new_apriori_name in os.environ:
        val = os.environ[envvar_new_apriori_name].lower()
        align_with_apriori = envvar_bool_dict[val]

    # Add support for environment variable switch to automatically
    # reset IDCTAB in FLT/FLC files if different from IDCTAB in RAW files.
    reset_idctab_switch = False
    if envvar_reset_idctab_name in os.environ:
        val = os.environ[envvar_reset_idctab_name].lower()
        reset_idctab_switch = envvar_bool_dict[val]

    if headerlets or align_to_gaia:
        from stwcs.wcsutil import headerlet

    # Open the input file
    try:
        # Make sure given filename is complete and exists...
        inFilename = fileutil.buildRootname(inFile, ext=['.fits'])
        if not os.path.exists(inFilename):
            print("ERROR: Input file - %s - does not exist." % inFilename)
            return
    except TypeError:
        print("ERROR: Inappropriate input file.")
        return

    # If newpath was specified, move all files to that directory for processing
    if newpath:
        orig_processing_dir = os.getcwd()
        new_processing_dir = _createWorkingDir(newpath, inFilename)
        _copyToNewWorkingDir(new_processing_dir, inFilename)
        os.chdir(new_processing_dir)

    # Initialize for later use...
    _mname = None
    _new_asn = None
    _calfiles = []

    # Identify WFPC2 inputs to account for differences in WFPC2 inputs
    infile_inst = fits.getval(inFilename, 'instrume')
    infile_det = fits.getval(inFilename, 'detector')
    wfpc2_input = infile_inst == 'WFPC2'
    cal_ext = None

    # Check input file to see if [DRIZ/DITH]CORR is set to PERFORM
    if '_asn' in inFilename:
        # We are working with an ASN table.
        # Use asnutil code to extract filename
        inFilename = _lowerAsn(inFilename)
        _new_asn = [inFilename]
        _asndict = asnutil.readASNTable(inFilename, None, prodonly=False)
        _cal_prodname = _asndict['output'].lower()
        # _fname = fileutil.buildRootname(_cal_prodname,ext=['_drz.fits'])

        # Retrieve the first member's rootname for possible use later
        _fimg = fits.open(inFilename, memmap=False)
        for name in _fimg[1].data.field('MEMNAME'):
            if name[-1] != '*':
                _mname = name.split('\0', 1)[0].lower()
                break
        _fimg.close()
        del _fimg

    else:
        # Check to see if input is a _RAW file
        # If it is, strip off the _raw.fits extension...
        _indx = inFilename.find('_raw')
        if _indx < 0: _indx = len(inFilename)
        # ... and build the CALXXX product rootname.
        if wfpc2_input:
            # force code to define _c0m file as calibrated product to be used
            cal_ext = ['_c0m.fits']
        _mname = fileutil.buildRootname(inFilename[:_indx], ext=cal_ext)

        _cal_prodname = inFilename[:_indx]
        # Reset inFilename to correspond to appropriate input for
        # drizzle: calibrated product name.
        inFilename = _mname

        if _mname is None:
            errorMsg = 'Could not find calibrated product!'
            raise Exception(errorMsg)

    # Create trailer filenames based on ASN output filename or
    # on input name for single exposures
    if '_raw' in inFile:
        # Output trailer file to RAW file's trailer
        _trlroot = inFile[:inFile.find('_raw')]
    elif '_asn' in inFile:
        # Output trailer file to ASN file's trailer, not product's trailer
        _trlroot = inFile[:inFile.find('_asn')]
    else:
        # Default: trim off last suffix of input filename
        # and replacing with .tra
        _indx = inFile.rfind('_')
        if _indx > 0:
            _trlroot = inFile[:_indx]
        else:
            _trlroot = inFile

    _trlfile = _trlroot + '.tra'
    _alignlog = _trlroot + '_align.log'
    _calfiles_flc = []

    # Write message out to temp file and append it to full trailer file
    _updateTrlFile(_trlfile, trlmsg)

    # Open product and read keyword value
    # Check to see if product already exists...
    dkey = 'DRIZCORR'
    # ...if product does NOT exist, interrogate input file
    # to find out whether 'dcorr' has been set to PERFORM
    # Check if user wants to process again regardless of DRIZCORR keyword value
    if force:
        dcorr = 'PERFORM'
    else:
        if _mname:
            _fimg = fits.open(fileutil.buildRootname(_mname, ext=['_raw.fits']), memmap=False)
            _phdr = _fimg['PRIMARY'].header
            if dkey in _phdr:
                dcorr = _phdr[dkey]
            else:
                dcorr = None
            _fimg.close()
            del _fimg
        else:
            dcorr = None

    if '_asn.fits' not in inFilename:
        # Working with a singleton
        # However, we always want to make sure we always use
        # a calibrated product as input, if available.
        _infile = fileutil.buildRootname(_cal_prodname, ext=cal_ext)
        _infile_flc = fileutil.buildRootname(_cal_prodname, ext=['_flc.fits'])

        _cal_prodname = _infile
        _calfiles = [_infile]

        _inlist = [_infile]

        print("_calfiles initialized as: {}".format(_calfiles))
        if len(_calfiles) == 1 and "_raw" in _calfiles[0]:
            _verify = False

        # Add CTE corrected filename as additional input if present
        if os.path.exists(_infile_flc) and _infile_flc != _infile:
            _calfiles_flc = [_infile_flc]
            _inlist = [_infile, _infile_flc]

    else:
        # Working with an ASN table...
        _infile = inFilename
        flist, duplist = processInput.checkForDuplicateInputs(_asndict['order'])
        _calfiles = flist
        if len(duplist) > 0:
            origasn = processInput.changeSuffixinASN(inFilename, 'flt')
            dupasn = processInput.changeSuffixinASN(inFilename, 'flc')
            _inlist = [origasn, dupasn]
        else:
            _inlist = [_infile]
        # We want to keep the original specification of the calibration
        # product name, though, not a lower-case version...
        _cal_prodname = inFilename
        _new_asn.extend(_inlist)  # kept so we can delete it when finished

        # check to see whether FLC files are also present, and need to be updated
        # generate list of FLC files
        _calfiles_flc = [f.replace('_flt.fits', '_flc.fits')
                         for f in _calfiles
                         if os.path.exists(f.replace('_flt.fits', '_flc.fits'))]

    # If specified, insure that IDCTAB in FLT/FLC files are the same
    # as the IDCTAB found in the RAW files
    if reset_idctab_switch:
        reset_idctab_kw(_calfiles, _calfiles_flc, logfile=_trlfile)

    # Add S_REGION keyword to input files regardless of whether DRIZCORR is turned on
    for f in _calfiles+_calfiles_flc:
        processing_utils.compute_sregion(f)

    if dcorr == 'PERFORM':

        """
        Start updating the data and verifying that the new alignment is valid.
            1. Run updatewcs without astrometry database update on all input exposures (FLCs? and FLTs)
            2. Generate initial default products and perform verification
                a. perform cosmic-ray identification and generate
                    drizzle products using astrodrizzle for all sets of inputs
                b. verify relative alignment with focus index after masking out CRs
                c. copy all drizzle products to parent directory
                d. if alignment fails, update trailer file with failure information
                e. if alignment verified, copy updated input exposures to parent directory
            3. If alignment is verified,
                0. copy inputs to separate sub-directory for processing
                a. run updatewcs to get a priori updates
                a.1.  apply 'best' apriori (not aposteriori) solution
                b. generate drizzle products for all sets of inputs (FLC and/or FLT) without CR identification
                c. verify alignment using focus index on FLC or, if no FLC, FLT products
                d. if alignment fails, update trailer file with info on failure
                e. if product alignment verified,
                    - copy all drizzle products to parent directory
                    - copy updated input exposures to parent directory
            4. If a posteriori correction enabled,
                0. copy all inputs to separate sub-directory for processing
                a. run align to align the images
                b. generate drizzle products for all sets of inputs (FLC and/or FLT) without CR identification
                c. verify alignment using focus index on FLC or, if no FLC, FLT products
                d. determine similarity index relative to pipeline default product
                e. if either focus or similarity indicates a problem, update trailer file with info on failure
                f. if product alignment verified,
                    - copy all drizzle products to parent directory
                    - copy updated input exposures to parent directory
            5. Remove all processing sub-directories
        """
        inst_mode = "{}/{}".format(infile_inst, infile_det)
        _good_images = [f for f in _calfiles if fits.getval(f, 'exptime') > 0.]
        _good_images = [f for f in _good_images if fits.getval(f, 'ngoodpix', ext=("SCI", 1)) > 0.]
        if len(_good_images) == 0:
            _good_images = _calfiles
        adriz_pars = mdzhandler.getMdriztabParameters(_good_images)
        adriz_pars.update(pipeline_pars)
        adriz_pars['mdriztab'] = False
        adriz_pars['final_fillval'] = "INDEF"
        adriz_pars['driz_sep_kernel'] = 'turbo'
        adriz_pars['driz_sep_fillval'] = 0.0
        adriz_pars['num_cores'] = num_cores
        adriz_pars['resetbits'] = 0

        exptimes = np.array([fits.getval(flt, 'exptime') for flt in _calfiles])
        if exptimes.max() / exptimes.min() > 2:
            adriz_pars['combine_type'] = 'median'
            adriz_pars['combine_nhigh'] = 0

        # Run updatewcs on each list of images to define pipeline default WCS
        # based on latest distortion models
        # Always apply latest distortion to replace pipeline-default OPUS WCS
        # for successful creation of updated headerlets for the cases where
        # all inputs having EXPTIME==0 (for example)
        updatewcs.updatewcs(_calfiles, use_db=False, checkfiles=False)
        if _calfiles_flc:
            updatewcs.updatewcs(_calfiles_flc, use_db=False, checkfiles=False)

        # Integrate user-specified drizzle parameters into pipeline_pars
        _trlmsg = _timestamp('Starting alignment with bad-pixel identification')
        _trlmsg += __trlmarker__
        _updateTrlFile(_trlfile, _trlmsg)

        if align_with_apriori or force_alignment or align_to_gaia:
            # Generate initial default products and perform verification
            align_dicts, align_table = verify_alignment(_inlist,
                                             _calfiles, _calfiles_flc,
                                             _trlfile,
                                             tmpdir=None, debug=debug,
                                             force_alignment=force_alignment,
                                             find_crs=True, **adriz_pars)

        if align_with_apriori:
            _trlmsg = _timestamp('Starting alignment with a priori solutions')
            _trlmsg += __trlmarker__
            if align_dicts is not None:
                find_crs = not align_dicts[0]['alignment_verified']
            else:
                find_crs = False

            # run updatewcs with use_db=True to insure all products have
            # have a priori solutions as extensions
            # FIX: This should probably only be done in the apriori sub-directory!
            updatewcs.updatewcs(_calfiles)
            _trlmsg += "Adding apriori WCS solutions to {}\n".format(_calfiles)
            _trlmsg += verify_gaia_wcsnames(_calfiles) + '\n'
            _wnames_calfiles = [(c, fits.getval(c, 'wcsname', ext=1)) for c in _calfiles]
            _trlmsg += "Verifying apriori WCSNAMEs:\n"
            for (_cname, _wname) in _wnames_calfiles:
                _trlmsg += "   {}: {}\n".format(_cname, _wname)
            if _calfiles_flc:
                _trlmsg += "Adding apriori WCS solutions to {}\n".format(_calfiles_flc)
                updatewcs.updatewcs(_calfiles_flc)
                _trlmsg += verify_gaia_wcsnames(_calfiles_flc) + '\n'

            try:
                tmpname = "_".join([_trlroot, 'apriori'])
                sub_dirs.append(tmpname)
                # Generate initial default products and perform verification
                align_apriori, apriori_table = verify_alignment(_inlist,
                                                 _calfiles, _calfiles_flc,
                                                 _trlfile,
                                                 tmpdir=tmpname, debug=debug,
                                                 good_bits=focus_pars[inst_mode]['good_bits'],
                                                 alignment_mode='apriori',
                                                 force_alignment=force_alignment,
                                                 find_crs=find_crs,
                                                 **adriz_pars)
            except Exception:
                # Reset to state prior to applying a priori solutions
                traceback.print_exc()
                align_apriori = None
                _trlmsg += "ERROR in applying a priori solution.\n"

            if align_apriori is None or (not align_apriori[0]['alignment_verified']):
                _trlmsg += "Resetting WCS to pipeline-default solutions..."
                # This operation replaces the PRIMARY WCS with one from the attached
                # headerlet extensions that corresponds to the distortion-model
                # solution created in the first place with 'updatewcs(use_db=False)'
                # Doing so, retains all solutions added from the astrometry database
                # while resetting to use whatever solution was defined by the instrument
                # calibration, since 'updatewcs' does not by default replace solutions
                # it finds in the files.
                restore_pipeline_default(_calfiles)
                if _calfiles_flc:
                    restore_pipeline_default(_calfiles_flc)

            else:
                align_dicts = align_apriori
                if align_dicts[0]['alignment_quality'] == 0:
                    _trlmsg += 'A priori alignment SUCCESSFUL.\n'
                if align_dicts[0]['alignment_quality'] == 1:
                    _trlmsg += 'A priori alignment potentially compromised.  Please review final product!\n'
                if align_dicts[0]['alignment_quality'] > 1:
                    _trlmsg += 'A priori alignment FAILED! No a priori astrometry correction applied.\n'
            _updateTrlFile(_trlfile, _trlmsg)

        aposteriori_table=None
        if align_to_gaia:
            _trlmsg = _timestamp('Starting a posteriori alignment')
            _trlmsg += __trlmarker__

            #
            # Start by creating the 'default' product using a priori/pipeline WCS
            # This product will be used as the final output if alignment fails
            # and will be used as the reference to compare to the aligned
            # product to determine whether alignment was ultimately successful or not.
            #
            # Call astrodrizzle to create the drizzle products
            if align_dicts is not None:
                find_crs = not align_dicts[0]['alignment_verified']
            else:
                find_crs = False
            tmpname = "_".join([_trlroot, 'aposteriori'])
            sub_dirs.append(tmpname)
            align_aposteriori, aposteriori_table = verify_alignment(_inlist,
                                             _calfiles, _calfiles_flc,
                                             _trlfile,
                                             tmpdir=tmpname, debug=debug,
                                             good_bits=focus_pars[inst_mode]['good_bits'],
                                             alignment_mode='aposteriori',
                                             force_alignment=force_alignment,
                                             find_crs=find_crs,
                                             **adriz_pars)
            if align_aposteriori:
                align_dicts = align_aposteriori
                align_qual = align_dicts[0]['alignment_quality']
                if align_qual == 0:
                    _trlmsg += 'A posteriori alignment SUCCESSFUL.\n'
                elif align_qual == 1:
                    _trlmsg += 'A posteriori alignment potentially COMPROMISED with bad focus.\n'
                    _trlmsg += '  Please review final product for alignment!\n'
                elif align_dicts[0]['alignment_quality'] == 2:
                    _trlmsg += 'A posteriori alignment potentially COMPROMISED.\n'
                    _trlmsg += 'Please review final product!\n'
                else:
                    _trlmsg += 'A posteriori alignment FAILED! No a posteriori astrometry correction applied.\n'
            _updateTrlFile(_trlfile, _trlmsg)

        _trlmsg = _timestamp('Creating final combined,corrected product based on best alignment')
        _trlmsg += __trlmarker__
        _updateTrlFile(_trlfile, _trlmsg)

        # Generate final pipeline products based on 'best' alignment
        pipeline_pars['in_memory'] = inmemory
        pipeline_pars['clean'] = True
        pipeline_pars['num_cores'] = num_cores

        drz_products, asn_dicts, diff_dicts = run_driz(_inlist, _trlfile, _calfiles,
                                             verify_alignment=False,
                                             good_bits=focus_pars[inst_mode]['good_bits'],
                                             **pipeline_pars)

        # Save this for when astropy.io.fits can modify a file 'in-place'
        # Update calibration switch
        _fimg = fits.open(_cal_prodname, mode='update', memmap=False)
        _fimg['PRIMARY'].header[dkey] = 'COMPLETE'
        _fimg.close()
        del _fimg

        # Enforce pipeline convention of all lower-case product
        # names
        for drz_product in drz_products:
            _plower = drz_product.lower()
            if drz_product != _plower: os.rename(drz_product, _plower)

    else:
        # Create default trailer file messages when astrodrizzle is not
        # run on a file.  This will typically apply only to BIAS,DARK
        # and other reference images.
        # Start by building up the message...
        _trlmsg = _timestamp('astrodrizzle skipped ')
        _trlmsg += __trlmarker__
        _trlmsg += '%s: astrodrizzle processing not requested for %s.\n' % (_getTime(), inFile)
        _trlmsg += '       astrodrizzle will not be run at this time.\n'

        # Write message out to temp file and append it to full trailer file
        _updateTrlFile(_trlfile, _trlmsg)

    # If we created a new ASN table, we need to remove it
    if _new_asn is not None:
        for _name in _new_asn: fileutil.removeFile(_name)

    # Insure all input FLC/FLT files have updated WCSTYPE* keywords
    for fname in _calfiles + _calfiles_flc:
        with fits.open(fname, mode='update') as fhdu:
            numsci = fileutil.countExtn(fhdu)
            for extn in range(1, numsci + 1):
                scihdr = fhdu[('sci', extn)].header
                keys = altwcs.wcskeys(fhdu, ('sci', extn))
                for key in keys:
                    wname = 'wcsname' + key.strip()
                    wtype = 'wcstype' + key.strip()
                    if wname in scihdr:
                        wval = scihdr[wname]
                    else:
                        # IF WCSNAME is missing, update with default value of OPUS
                        wval = "OPUS"
                        scihdr[wname] = wval

                    scihdr[wtype] = updatehdr.interpret_wcsname_type(wval)

    # If headerlets have already been written out by alignment code,
    # do NOT write out this version of the headerlets
    if headerlets:
        # Generate headerlets for each updated FLT image
        hlet_msg = _timestamp("Writing Headerlets started")
        for fname in _calfiles:
            hlet_msg += "Creating new headerlet from {}".format(fname)
            frootname = fileutil.buildNewRootname(fname)
            hname = "%s_flt_hlet.fits" % frootname
            # Write out headerlet file used by astrodrizzle, however,
            # do not overwrite any that was already written out by align
            if not os.path.exists(hname):
                hlet_msg += "Created Headerlet file %s \n" % hname
                try:
                    wcsname = fits.getval(fname, 'wcsname', ext=1)
                    wcstype = updatehdr.interpret_wcsname_type(wcsname)
                    hdrname = "{}_{}-hlet.fits".format(fname.replace('.fits', ''), wcsname)
                    headerlet.write_headerlet(fname, hdrname, output='flt',
                                              wcskey='PRIMARY',
                                              author="OPUS",
                                              descrip=wcstype,
                                              attach=False,
                                              clobber=True,
                                              logging=False)
                except ValueError:
                    hlet_msg += _timestamp("SKIPPED: Headerlet not created for %s \n" % fname)
                    # update trailer file to log creation of headerlet files
        hlet_msg += _timestamp("Writing Headerlets completed")
        ftrl = open(_trlfile, 'a')
        ftrl.write(hlet_msg)
        ftrl.close()

    # Remove secondary log files for good...
    logging.shutdown()

    for _olog in [_alignlog]:
        if os.path.exists(_olog):
            os.remove(_olog)

    # If processing was done in a temp working dir, restore results to original
    # processing directory, return to original working dir and remove temp dir
    if newpath:
        _restoreResults(new_processing_dir, orig_processing_dir)
        os.chdir(orig_processing_dir)
        _removeWorkingDir(new_processing_dir)

    if debug and Process is not None:
        print("Files still open for this process include: ")
        print([ofile.path for ofile in Process().open_files()])

    if not debug:
        try:
            # Remove all temp sub-directories now that we are done
            for sd in sub_dirs:
                if os.path.exists(sd): rmtree2(sd)
        except Exception:
            # If we are unable to remove any of these sub-directories,
            # leave them for the user or calling routine/pipeline to clean up.
            print("WARNING: Unable to remove any or all of these sub-directories: \n{}\n".format(sub_dirs))
            if Process is not None:
                print("Files still open at this time include: ")
                print([ofile.path for ofile in Process().open_files()])
            pass

    # Append final timestamp to trailer file...
    end_time = _getTime()
    _delta_time = time.time() - init_time
    _final_msg = '%s: Finished processing %s in %.2f seconds \n' % (end_time, inFilename, _delta_time)
    _final_msg += _timestamp('astrodrizzle completed ')

    _updateTrlFile(_trlfile, _final_msg)

    # Provide feedback to user
    print(_final_msg)

    # Look to see whether we have products which can be evaluated
    # wcsname = fits.getval(drz_products[0], 'wcsname', ext=1)

    # interpret envvar variable, if specified
    qa_switch = _get_envvar_switch(envvar_qa_stats_name)

    if qa_switch and dcorr == 'PERFORM':

        # Generate quality statistics for astrometry if specified
        calfiles = _calfiles_flc if _calfiles_flc else _calfiles
        qa.run_all(inFile, calfiles, catalogs=aposteriori_table)


def run_driz(inlist, trlfile, calfiles, mode='default-pipeline', verify_alignment=True,
            debug=False, good_bits=512, **pipeline_pars):

    import drizzlepac
    pyver = drizzlepac.astrodrizzle.__version__
    drz_products = []
    focus_dicts = []
    diff_dicts = OrderedDict()

    pipeline_pars['runfile'] = trlfile.replace('.tra', '_pydriz')
    drizlog = pipeline_pars['runfile'] + ".log"  # the '.log' gets added automatically by astrodrizzle
    for infile in inlist:  # Run astrodrizzle for all inputs
        asndict, ivmlist, drz_product = processInput.process_input(infile, updatewcs=False,
                                                        preserve=False,
                                                        overwrite=False)
        del ivmlist
        calfiles = asndict['original_file_names'] if asndict is not None else calfiles

        drz_products.append(drz_product)

        # Create trailer marker message for start of astrodrizzle processing
        _trlmsg = _timestamp('astrodrizzle started ')
        _trlmsg += __trlmarker__
        _trlmsg += '%s: Processing %s with astrodrizzle Version %s\n' % (_getTime(), infile, pyver)
        print(_trlmsg)
        _updateTrlFile(trlfile, _trlmsg)

        _pyd_err = trlfile.replace('.tra', '_pydriz.stderr')

        try:
            drizzlepac.astrodrizzle.AstroDrizzle(input=infile, configobj=None,
                                                 **pipeline_pars)

            # Edit trailer file name since 'runastrodriz' copies what astrodrizzle used
            # to another file...
            if os.path.exists(drz_product):
                fits.setval(drz_product, 'DRIZPARS', value=trlfile)

            util.end_logging(drizlog)

        except Exception as errorobj:
            _appendTrlFile(trlfile, drizlog)
            _appendTrlFile(trlfile, _pyd_err)
            _ftrl = open(trlfile, 'a')
            _ftrl.write('ERROR: Could not complete astrodrizzle processing of %s.\n' % infile)
            _ftrl.write(str(sys.exc_info()[0]) + ': ')
            _ftrl.writelines(str(errorobj))
            _ftrl.write('\n')
            _ftrl.close()
            print('ERROR: Could not complete astrodrizzle processing of %s.' % infile)
            raise Exception(str(errorobj))

        # For singletons, there is no need to perform focus check since there is only 1 input exposure
        if len(calfiles) == 1:
            verify_alignment = False

        if verify_alignment and asndict is not None:
            # Evaluate generated products: single_sci vs drz/drc
            # FLT files are always first, and we want FLC when present
            cal_suffix = '_flt' if calfiles[0].endswith('_flt.fits') else '_flc'
            single_files = [calfile.replace(cal_suffix, '_single_sci') for calfile in calfiles]
            sfile = single_files[0]
            if not os.path.exists(sfile):
                # Working with data where CR is turned off by default (ACS/SBC, for example)
                # Reset astrodrizzle parameters to generate single_sci images
                reset_mdriztab_nocr(pipeline_pars, good_bits, pipeline_pars['skysub'])

                drizzlepac.astrodrizzle.AstroDrizzle(input=infile, configobj=None,
                                                    **pipeline_pars)

            instr_det = "{}/{}".format(fits.getval(sfile, 'instrume'), fits.getval(sfile, 'detector'))
            focus_sigma = focus_pars[instr_det]['sigma']
            print("Measuring similarity and focus for: \n{} \n    {}".format(single_files, drz_product))
            focus_dicts.append(amutils.build_focus_dict(single_files, drz_product, sigma=focus_sigma))
            if debug:
                json_name = drz_product.replace('.fits', '_{}_focus.json'.format(mode))
                with open(json_name, mode='w') as json_file:
                    json.dump(focus_dicts, json_file)

            # Compute additional verification based on Hamming distances of gradients in overlap region
            drz_wcs = HSTWCS(drz_product, ext=("SCI", 1))
            drz_footprint = cell_utils.SkyFootprint(drz_wcs)
            drz_footprint.build(single_files)
            diff_dicts[drz_product] = amutils.max_overlap_diff(drz_footprint.total_mask,
                                                               single_files, drz_product,
                                                               scale=1)

        else:
            focus_dicts = None

    # Now, append comments created by PyDrizzle to CALXXX trailer file
    print('Updating trailer file %s with astrodrizzle comments.' % trlfile)
    drizlog_copy = drizlog.replace('.log', '_copy.log')
    if os.path.exists(drizlog):
        shutil.copy(drizlog, drizlog_copy)
    _appendTrlFile(trlfile, drizlog_copy)
    # clean up log files
    if RM_LOGFILES and os.path.exists(drizlog):
        os.remove(drizlog)
    # Clean up intermediate files generated by astrodrizzle
    if not debug:
        for ftype in ['*mask*.fits']:
            [os.remove(file) for file in glob.glob(ftype)]

    return drz_products, focus_dicts, diff_dicts

def reset_mdriztab_nocr(pipeline_pars, good_bits, skysub):
    # Need to turn off MDRIZTAB if any other parameters are to be set
    pipeline_pars['mdriztab'] = False
    pipeline_pars['build'] = True
    pipeline_pars['resetbits'] = 0
    pipeline_pars['static'] = False
    pipeline_pars['skysub'] = skysub
    pipeline_pars['driz_separate'] = True
    pipeline_pars['driz_sep_bits'] = good_bits
    pipeline_pars['driz_sep_fillval'] = 0.0
    pipeline_pars['median'] = False
    pipeline_pars['blot'] = False
    pipeline_pars['driz_cr'] = False
    pipeline_pars['final_fillval'] = "INDEF"


def verify_alignment(inlist, calfiles, calfiles_flc, trlfile,
                     find_crs=True, tmpdir=None, debug=False, good_bits=512,
                     alignment_mode=None, force_alignment=False,
                     **pipeline_pars):

    headerlet_files = []
    for infile in inlist:
        asndict, ivmlist, drz_product = processInput.process_input(infile, updatewcs=False,
                                                    preserve=False,
                                                    overwrite=False)
        del ivmlist
        # If there are no products to be generated, there is nothing to align...
        if asndict is None:
            return None, None

    # if tmpdir is turned off (== None), tmpname set to 'default-pipeline'
    tmpname = tmpdir if tmpdir else 'default-pipeline'
    tmpmode = None
    for mode in valid_alignment_modes:
        if mode in tmpname:
            tmpmode = mode
            break
    if tmpmode is None:
        print("Invalid alignment mode {} requested.".format(tmpdir))
        raise ValueError

    full_table = None
    fraction_matched = 1.0
    num_sources = -1
    try:
        if not find_crs:
            # Need to turn off MDRIZTAB if any other parameters are to be set
            reset_mdriztab_nocr(pipeline_pars, good_bits, pipeline_pars['skysub'])

        if tmpdir:
            # Create tmp directory for processing
            if not os.path.exists(tmpdir):
                os.makedirs(tmpdir)

            # Now, copy all necessary files to tmpdir
            _ = [shutil.copy(f, tmpdir) for f in inlist + calfiles]
            if calfiles_flc:
                _ = [shutil.copy(f, tmpdir) for f in calfiles_flc]

            parent_dir = os.getcwd()
            os.chdir(tmpdir)

        # insure these files exist, if not, blank them out
        # Also pick out what files will be used for additional alignment to GAIA
        if not calfiles_flc or not os.path.exists(calfiles_flc[0]):
            calfiles_flc = None

        alignfiles = calfiles_flc if calfiles_flc else calfiles
        align_update_files = calfiles if calfiles_flc else None
        inst = fits.getval(alignfiles[0], 'instrume').lower()
        det = fits.getval(alignfiles[0], 'detector').lower()

        if find_crs:
            trlmsg = _timestamp("Resetting CRs ")
            # reset all DQ flags associated with CRs assuming previous attempts were inaccurate
            for f in alignfiles:
                trlmsg += "Resetting CR DQ bits for {}\n".format(f)
                resetbits.reset_dq_bits(f, "4096,8192")
                sat_flags = 256 + 2048
        else:
            sat_flags = 256 + 2048 + 4096 + 8192

        # Perform any requested alignment here...
        if alignment_mode == 'aposteriori':
            # Create trailer marker message for start of align_to_GAIA processing
            trlmsg = _timestamp("Align_to_GAIA started ")
            _updateTrlFile(trlfile, trlmsg)
            # Evaluate all input exposures and, if necessary, reset WCSs to a common WCS
            # This is necessary in order to avoid imprinting zero point differences between
            # coordinate systems into the final fit which would result in mis-alignment of the
            # sources in the final combined output images.
            if len(alignfiles) > 1:
                update_wcs_in_list(alignfiles, logfile=trlfile)

            instdet_pars = align.get_default_pars(inst, det)

            alignlog = trlfile.replace('.tra', '_align.log')
            alignlog_copy = alignlog.replace('_align', '_align_copy')
            try:

                full_table = align.perform_align(alignfiles,
                                                 catalog_list=instdet_pars['run_align']['catalog_list'],
                                                 num_sources=instdet_pars['general']['MAX_SOURCES_PER_CHIP'],
                                                 update_hdr_wcs=True, runfile=alignlog,
                                                 clobber=False, output=debug,
                                                 debug=debug, sat_flags=sat_flags)
                if full_table is None:
                    raise Exception("No successful aposteriori fit determined.")

                align_table = full_table.filtered_table

                num_sources = align_table['matchSources'][0]
                fraction_matched = num_sources / align_table['catalogSources'][0]

                for row in align_table:
                    if row['status'] == 0:
                        if row['compromised'] == 0:
                            trlstr = "Successfully aligned {} to {} astrometric frame\n"
                            trlmsg += trlstr.format(row['imageName'], row['catalog'])
                        else:
                            trlstr = "Alignment only partially successful for {}\n"
                            trlmsg += trlstr.format(row['imageName'])
                    else:
                        trlstr = "Could not align {} to absolute astrometric frame\n"
                        trlmsg += trlstr.format(row['imageName'])
                        print(trlmsg)
                        _updateTrlFile(trlfile, trlmsg)
                        return None, None
            except Exception as err:
                # Something went wrong with alignment to GAIA, so report this in
                # trailer file
                _trlmsg = "EXCEPTION encountered in align...\n"
                _trlmsg += "   No correction to absolute astrometric frame applied!\n"
                print(_trlmsg)
                _updateTrlFile(trlfile, _trlmsg)
                if 'aposteriori' not in repr(err):
                    traceback.print_exc()
                else:
                    print("WARNING: {}".format(err))
                return None, None

            _updateTrlFile(trlfile, trlmsg)
            # Write the perform_align log to the trailer file...(this will delete the _alignlog)
            if os.path.exists(alignlog):
                shutil.copy(alignlog, alignlog_copy)
                _appendTrlFile(trlfile, alignlog_copy)

            _trlmsg = ""
            # Check to see whether there are any additional input files that need to
            # be aligned (namely, FLT images)
            if align_update_files and align_table:
                # Apply headerlets from alignment to FLT version of the files
                for fltfile, flcfile in zip(align_update_files, alignfiles):
                    # Update non-headerlet-based keywords in fltfile
                    _update_wcs_fit_keywords(fltfile, flcfile)
                    row = align_table[align_table['imageName'] == flcfile]
                    headerlet_file = row['headerletFile'][0]
                    if headerlet_file not in ["None", '']:
                        apply_headerlet(fltfile, headerlet_file, flcfile=flcfile)
                        # headerlet.apply_headerlet_as_primary(fltfile, headerlet_file,
                        #                                     attach=True, archive=True)
                        headerlet_files.append(headerlet_file)
                        # append log file contents to _trlmsg for inclusion in trailer file
                        _trlstr = "Applying headerlet {} as Primary WCS to {}\n"
                        _trlmsg += _trlstr.format(headerlet_file, fltfile)
                    else:
                        _trlmsg += "No absolute astrometric headerlet applied to {}\n".format(fltfile)

            # Finally, append any further messages associated with alignement from this calling routine
            _trlmsg += _timestamp('Align_to_GAIA completed ')
            _updateTrlFile(trlfile, _trlmsg)

        if find_crs:
            drz, fdicts, ddicts = run_driz(inlist, trlfile, calfiles, mode=tmpmode,
                                        verify_alignment=False, debug=debug,
                                        good_bits=good_bits, **pipeline_pars)

        # Run astrodrizzle in desired mode
        drz_products, focus_dicts, diff_dicts = run_driz(inlist, trlfile, calfiles,
                                                         mode=tmpmode, verify_alignment=True,
                                                         debug=debug, good_bits=good_bits,
                                                         **pipeline_pars)

        # Start verification of alignment using focus and similarity indices
        _trlmsg = _timestamp('Verification of {} alignment started '.format(tmpmode))

        if focus_dicts is not None:
            # Only check focus on CTE corrected, when available
            align_focus = focus_dicts[-1] if 'drc' in focus_dicts[-1]['prodname'] else focus_dicts[0]

            pscale = HSTWCS(alignfiles[0], ext=1).pscale

            det_pars = align.get_default_pars(inst, det)['generate_source_catalogs']
            default_fwhm = det_pars['fwhmpsf'] / pscale
            align_fwhm = amutils.get_align_fwhm(align_focus, default_fwhm)

            if align_fwhm:
                _trlmsg += "align_fwhm: {}[{},{}]={:0.4f}pix\n".format(align_focus['prodname'],
                                                            align_focus['prod_pos'][1],
                                                            align_focus['prod_pos'][0],
                                                            align_fwhm)

            # Interpret the overlap differences computed for this alignment
            dkeys = [k for k in diff_dicts.keys()]
            diff_verification, max_diff = amutils.evaluate_overlap_diffs(diff_dicts[dkeys[-1]])
            _trlmsg += "Fraction of sources matched: {}  out of {} sources\n".format(fraction_matched, num_sources)

            # For any borderline situation with alignment, perform an extra check on alignment
            if fraction_matched < 0.1 or -1 < num_sources < 10:
                focus_verification = amutils.evaluate_focus(align_focus)
                alignment_verified = True if (diff_verification and focus_verification) else False
                alignment_quality = 0 if alignment_verified else 3
            else:
                alignment_verified = diff_verification
                alignment_quality = 0 if diff_verification else 3

            if alignment_verified:
                _trlmsg += "Focus verification indicated that {} alignment SUCCEEDED.\n".format(tmpmode)
            else:
                _trlmsg += "Focus verification indicated that {} alignment FAILED.\n".format(tmpmode)
                _trlmsg += "  Reverting to previously determined WCS alignment.\n"

            prodname = align_focus['prodname']
        else:
            fd = amutils.FOCUS_DICT.copy()
            fd['expnames'] = calfiles
            fd['prodname'] = drz_products[0]
            alignment_verified = True
            alignment_quality = 0
            focus_dicts = [fd]

            prodname = drz_products[0]

        # For default pipeline alignment, we have nothing else to compare
        # similarity to, so skip this step...
        if alignment_mode:
            alignprod = fits.getdata(prodname, ext=1)

            # compute similarity_index as well and fold into alignment_verified state
            refname = os.path.abspath(os.path.join('..', prodname))
            align_ref = fits.getdata(refname, ext=1)

            print("Computing sim_indx for: {} ".format(os.path.join(tmpdir, prodname)))
            sim_indx = amutils.compute_similarity(alignprod, align_ref)
            align_sim_fail = sim_indx > 1

            if not align_sim_fail and alignment_verified:
                _trlmsg += "Alignment appeared to SUCCEED based on similarity index of {:0.4f} \n".format(sim_indx)
            else:
                _trlmsg += "Alignment appeared to FAIL based on similarity index of {:0.4f} \n".format(sim_indx)
                # _trlmsg += "  Reverting to previously determined WCS alignment.\n"
                # alignment_verified = False
                alignment_quality += 3

        for fd in focus_dicts:
            fd['alignment_verified'] = alignment_verified
            fd['alignment_quality'] = alignment_quality

        # If CRs were identified, copy updated input files to main directory
        if tmpdir and alignment_verified:
            _trlmsg += "Saving products with new alignment.\n"
            _ = [shutil.copy(f, parent_dir) for f in calfiles]
            if calfiles_flc:
                _ = [shutil.copy(f, parent_dir) for f in calfiles_flc]
            # Copy drizzle products to parent directory to replace 'less aligned' versions
            _ = [shutil.copy(f, parent_dir) for f in headerlet_files]

        _trlmsg += _timestamp('Verification of alignment completed ')
        _updateTrlFile(trlfile, _trlmsg)

    finally:
        if tmpdir:
            _appendTrlFile(os.path.join(parent_dir, trlfile), trlfile)
            # Return to main processing dir
            os.chdir(parent_dir)

    return focus_dicts, full_table

def apply_headerlet(filename, headerlet_file, flcfile=None):

    # Use headerlet module to apply headerlet as PRIMARY WCS
    headerlet.apply_headerlet_as_primary(filename, headerlet_file,
                                        attach=True, archive=True)
    # Verify that all keywords from headerlet got applied
    if flcfile is not None:
        with fits.open(filename, mode='update') as fhdu:
            num_sci = fileutil.countExtn(fhdu)
            for sciext in range(1, num_sci + 1):
                extn = ('sci', sciext)
                fhdu[extn].header['wcstype'] = fits.getval(flcfile, 'wcstype', extn)


def verify_gaia_wcsnames(filenames, catalog_name='GSC240', catalog_date=gsc240_date):
    """Insure that data taken with GAIA has WCSNAME reflecting that"""
    gsc240 = catalog_date.split('-')
    gdate = datetime.date(int(gsc240[0]), int(gsc240[1]), int(gsc240[2]))
    msg = ''

    for f in filenames:
        wcsnames = None
        with fits.open(f, mode='update') as fhdu:
            # Check to see whether a RAW/uncalibrated file has been provided
            # If so, skip it since updatewcs has not been run on it yet.
            if '_raw' in f or 'wcsname' not in fhdu[('sci', 1)].header:
                continue
            num_sci = fileutil.countExtn(fhdu)
            dateobs = fhdu[0].header['date-obs'].split('-')
            # convert to datetime object
            fdate = datetime.date(int(dateobs[0]), int(dateobs[1]), int(dateobs[2]))
            for sciext in range(num_sci):
                wcsname = fhdu['sci', sciext + 1].header['wcsname']
                if fdate > gdate and '-' not in wcsname:
                    wcsname = "{}-{}".format(wcsname, catalog_name)
                    fhdu['sci', sciext + 1].header['wcsname'] = wcsname
                    msg += "Updating WCSNAME of {}[sci,{}] for use of {} catalog \n".format(f,
                            sciext + 1, catalog_name)
                    continue
                # Check to see whether it is an aposteriori solution
                # If so, replace it with an apriori solution instead
                if '-FIT' in wcsname:
                    if wcsnames is None:
                        wcsnames = headerlet.get_headerlet_kw_names(fhdu, kw='WCSNAME')
                        hdrnames = headerlet.get_headerlet_kw_names(fhdu)
                        extvers = headerlet.get_headerlet_kw_names(fhdu, kw='extver')

                        # Remove duplicate hdrlet extensions
                        c = [hdrnames.count(h) for h in hdrnames]
                        extdict = {}
                        for num, hname, ev in zip(c, hdrnames, extvers):
                            if hname not in extdict:
                                extdict[hname] = []
                            if num > 1:
                                extdict[hname].append(ev)
                        for extvals in extdict.values():
                            if extvals:
                                extvals.sort()
                                remove_e = []
                                for e in extvals[:-1]:
                                    del fhdu[('hdrlet'), e]
                                    remove_e.append(e)
                                for e in remove_e:
                                    del extvers[extvers.index(e)]
                        # Remove OPUS based solutions
                        opus_indx = []
                        for i, w in enumerate(wcsnames):
                            if 'OPUS' in w: opus_indx.append(i)
                        opus_indx.reverse()
                        for i in opus_indx:
                            del wcsnames[i]
                            del hdrnames[i]
                        gscwcs = any(['GSC' in w for w in wcsnames])
                        if not gscwcs:
                            priwcs = fhdu['sci', sciext + 1].header['wcsname']
                            if priwcs not in wcsnames:
                                wcsnames.append(priwcs)
                                hdrnames.append(priwcs)
                                # Build IDC_* only WCSNAME
                                defwcs = priwcs.split("-")[0]
                                # delete this from list of WCSNAMEs since it was
                                # already replaced by GSC240 WCSNAME
                                indx = wcsnames.index(defwcs)
                                del wcsnames[indx]
                                del hdrnames[indx]

                    # This actually returns the date for the IDCTAB itself,
                    # not just when the specific WCS was created using it.
                    hlets = [fhdu[('hdrlet',e)].headerlet for e in extvers]
                    idctabs = [h[0].header['idctab'] for h in hlets]
                    wcsdates = [fits.getval(fileutil.osfn(idc), 'date') for idc in idctabs]

                    # Look for priority apriori WCS
                    restored = False
                    most_recent_wcs = None
                    for apriori_type in apriori_priority:
                        # For each WCSNAME/HDRNAME in the file...
                        for wname, hname, wdate in zip(wcsnames, hdrnames, wcsdates):
                            # Look for apriori_type (HSC, GSC,...) based on newest IDCTAB
                            if apriori_type in wname:
                                if (most_recent_wcs and wdate > most_recent_wcs[0]) \
                                    or most_recent_wcs is None:
                                    most_recent_wcs = (wdate, hname)

                        if most_recent_wcs:
                            break

                    if most_recent_wcs:
                        # restore this WCS
                        msg += 'Restoring apriori WCS {} as primary WCS in {}\n'.format(wname, f)
                        headerlet.restore_from_headerlet(fhdu,
                                                         force=True,
                                                         hdrname=most_recent_wcs[1],
                                                         archive=False)
                        # insure IDCSCALE is still present
                        if 'idcscale' not in fhdu[('sci', 1)].header:
                            msg += "Headerlet {} was missing IDCSCALE keyword".format(most_recent_wcs[1])
                            # get IDCTAB name
                            itabroot = fhdu[0].header['idctab'].split('$')[1].split('_')[0]
                            fhdu_idscale = None
                            # pull a value from one of the other headerlet extensions
                            for extn in extvers:
                                if itabroot in fhdu[('HDRLET', extn)].header['wcsname']:
                                    _hdrlet = fhdu[('HDRLET', extn)].headerlet
                                    if 'idcscale' in _hdrlet[('SIPWCS', 1)].header:
                                        fhdu_idscale = _hdrlet[('SIPWCS', 1)].header['idcscale']
                                        break
                            if fhdu_idscale is None:
                                cd11 = fhdu[('sci', sciext + 1)].header['CD1_1']
                                cd21 = fhdu[('sci', sciext + 1)].header['CD2_1']
                                fhdu_idscale = round(np.sqrt(np.power(cd11, 2) + np.power(cd21, 2)) * 3600., 3)
                            # Set the value of the IDCSCALE keyword
                            for extn in range(num_sci):
                                msg +=  'Adding IDCSCALE {} to {}[sci,{}]'.format(fhdu_idscale, fhdu.filename(), extn + 1)
                                fhdu[('sci', extn + 1)].header['idcscale'] = fhdu_idscale
    return msg

def restore_pipeline_default(files):
    """Restore pipeline-default IDC_* WCS as PRIMARY WCS in all input files"""
    print("Restoring pipeline-default WCS as PRIMARY WCS... using updatewcs.")
    updatewcs.updatewcs(files, use_db=False)
    # Remove HDRNAME, if added by some other code.
    #  This keyword only needs to be included in the headerlet file itself.
    for f in files:
        with fits.open(f, mode='update') as fhdu:
            num_sci = fileutil.countExtn(fhdu)
            for sciext in range(num_sci):
                if 'hdrname' in fhdu[('sci', sciext + 1)].header:
                    del fhdu[('sci', sciext + 1)].header['hdrname']


def reset_idctab_kw(files, files_flc, logfile=None):
    """Insure IDCTAB in files are the same as those in RAW files"""
    trlmsg = "Insuring IDCTAB keywords are up-to-date\n"

    raw_files = [f.replace('_flt', '_raw') for f in files]
    # determine what IDCTAB should be in FLT files
    raw_idctab = fits.getval(raw_files[0], 'idctab')
    raw_wcsname = 'IDC_{}'.format(raw_idctab.split('_')[0].split('$')[1])
    trlmsg += "IDCTAB from RAW file: {}\n".format(raw_idctab)

    if logfile:
        # Write message out to temp file and append it to full trailer file
        _updateTrlFile(logfile, trlmsg)
    else:
        print(trlmsg)

    # Now, restore headerlet with this WCSNAME as primary WCS
    for flt in files+files_flc:
        flt_idctab = fits.getval(flt, 'idctab')
        if flt_idctab == raw_idctab:
            # We don't need to update anything
            continue
        newmsg = "Updating IDCTAB {} in {}\n".format(flt_idctab, flt)
        if logfile:
            # Write message out to temp file and append it to full trailer file
            _updateTrlFile(logfile, newmsg)
        else:
            print(newmsg)

        # Get info from all hdrlet extensions in file
        wnames = headerlet.get_headerlet_kw_names(flt, 'wcsname')
        hnames = headerlet.get_headerlet_kw_names(flt, 'hdrname')
        # find the HDRNAME kw value for the hlet extension that has the desired WCSNAME
        idc_hdrname = hnames[wnames.index(raw_wcsname)]
        # Find the actual extension number for the headerlet with this HDRNAME
        hdrname_ext = headerlet.find_headerlet_HDUs(flt, hdrname=idc_hdrname)[0]
        # Restore the WCS from this HDRLET extension
        headerlet.restore_from_headerlet(flt, hdrext=hdrname_ext, force=True)


# ------------------------------------------------------------------------------


def update_wcs_in_list(exp_list, logfile=None):
    """Examine entries in the exposure list for inconsistent WCS solutions

    Note: The direct exposure image active WCS solutions will be modified in-place.

    Parameters
    ----------
    exp_list : list
        List of filenames for exposures to be evaluated.

    logfile : str, optional
        Name of file to write out processing messages generated by this function.

    Returns
    -------
    Nothing.

    """

    msg = "\n***** Processing List for Consistent WCS's *****"
    print(msg)
    update_msg = msg
    primary_wcsnames = set([fits.getval(fname, 'wcsname', ext=('SCI',1)) for fname in exp_list])

    if len(primary_wcsnames) == 1:
        msg = "\nAll Primary WCS's confirmed as consistent as {}.".format(primary_wcsnames)
        print(msg)
        update_msg += msg
        # Update trailer file with log messages
        if logfile:
            _updateTrlFile(logfile, update_msg)
        return

    # We have detected inconsistent WCSs...
    # Loop over all the direct images for this detector in the visit to update the WCS
    primary_idctabs = []
    for filename in exp_list:
        hdu = fits.open(filename)
        idctab = hdu[0].header['idctab'].split('$')[1]
        idcroot = idctab.split('_')[0]
        primary_idctabs.append('IDC_{}*'.format(idcroot))  # Add all WCS that start with IDC_<idctab> as desired WCSs

        # Make sure the WCS is up-to-date - doing the update for the
        # direct images here so there is no need to modify a
        # pre-existing class which could cause an issue
        drizcorr = hdu[0].header['DRIZCORR']
        if drizcorr == "OMIT":
            updatewcs.updatewcs(filename, use_db=True)

        hdu.close()
        del hdu
        # Insure HDRNAME keywords are properly populated in SCI extensions.
        wcs_functions.verify_sci_hdrname(filename)

    final_wcs_set, skip_direct_list, d_keyword_wcs_names_dict, direct_dict = collect_wcs_names(exp_list, 'DIRECT', logfile=logfile)
    msg = "WCS solutions common to all viable direct images: {}\n".format(final_wcs_set)
    print(msg)
    update_msg += msg

    # At this point, the IDCTAB points to the most current reference file
    # for these exposures.
    primary_idctabs = set(primary_idctabs)
    # Remove all WCSs from final_wcs_set that do not use the most recent IDCTAB
    primary_wcs_set = []
    for pidc in primary_idctabs:
        match_list = fnmatch.filter(final_wcs_set, pidc)
        if match_list:
            primary_wcs_set.extend(match_list)
    # Redefine using only those WCSs based on most recent distortion model
    # Should never need this logic, but just to be safe...
    if len(primary_wcs_set) == 0:
        # Simply rely on WCS solutions already in the headers
        msg = "NO Common WCS solutions for images: {}\n".format(final_wcs_set)
        msg += "  for the latest calibration IDCTAB: {}\n".format(primary_idctabs)
        print(msg)
        update_msg += msg
        if logfile:
            _updateTrlFile(logfile, update_msg)
        return
    final_wcs_set = set(primary_wcs_set)

    # There is a preference for the active WCS for the viable images in the visit
    # Check the final_wcs_set for the preferential solutions
    match_list = []
    final_wcsname = ''
    for wcs_item in wcs_preference:
        match_list = fnmatch.filter(final_wcs_set, wcs_item)
        if match_list:
            final_wcsname = match_list[0]
            msg = "Final WCS solution to use for all images: {}\n".format(final_wcsname)
            print(msg)
            update_msg += msg
            break

    if final_wcsname:
        # Finally, if the image is not in a skip list, reset the primary WCS in all the images
        for filename in exp_list:
            if filename not in skip_direct_list:
                msg = "\nSetting the primary WCS for direct image {} to {}.".format(filename, final_wcsname)
                print(msg)
                update_msg += msg
                update_active_wcs(filename, final_wcsname, logfile=logfile)
    else:
        # Do nothing
        pass
    # Update trailer file with log messages
    if logfile:
        _updateTrlFile(logfile, update_msg)

# ------------------------------------------------------------------------------


def collect_wcs_names(exp_list, image_type, logfile=None):
    """
    Utility to collect all the WCS solution names common to the input image list

    Parameters
    ----------
    exp_list: str list
        List containing the SVM FLT/FLC filenames

    image_type: string
        String containing either 'GRISM' or 'DIRECT' to use as a switch for
        output information

    Returns
    -------
    image_wcs_set: set of WCS solutions
        The set contains the WCS solution names common to all of the
        input exposures

    skip_image_list: list
        This is a list of exposures in the input list which should be
        skipped/ignored when updating the active WCS solution

    keyword_wcs_names_dict: dictionary {filename: list}
        The dictionary is used to associate an individual image/filename with
        a list of WCS solution names in the file stored as keywords (not headerlets)

    image_dict: dictionary {filename: list}
        The dictionary is used to associate an individual image/filename with
        a list of *all* WCS solution names in the file

    """
    update_msg = ""
    image_wcs_set = set()
    skip_image_list = []
    exist_image_set = False
    image_dict = {}
    keyword_wcs_names_dict = {}
    # Loop over all the Grism/Prism images for this detector in the visit
    for filename in exp_list:

        # Get all the WCS names which are common to all of the images
        # Note that WCS solutions may be represented as FITS keyword values in the
        # SCI extension and/or as headerlets in the HDRLET extensions.
        # Get the keyword WCS solution names.
        keyword_wcs_names = list(wcsutil.altwcs.wcsnames(filename, ext=1).values())

        # Get the headerlet WCS solution names
        headerlet_wcs_names = wcsutil.headerlet.get_headerlet_kw_names(filename, kw="WCSNAME")
        all_wcs_names = keyword_wcs_names + headerlet_wcs_names
        keyword_wcs_names_dict[filename] = keyword_wcs_names
        image_dict[filename] = all_wcs_names
        if all_wcs_names:
            msg = "WCS solutions for file {} are {}.".format(filename, all_wcs_names)
            print(msg)
            update_msg += msg

            # Initialize a set with wcsnames
            if not exist_image_set:
                image_wcs_set = set(all_wcs_names)
                exist_image_set = True
            # Generate the intersection with the set of existing wcsnames and the wcsnames
            # from the current image
            else:
                image_wcs_set &= set(all_wcs_names)

            # Oops...no common wcsnames
            if not image_wcs_set:
                msg = "ERROR: There are no common WCS solutions with this image {} and previously processed images\n".format(filename)
                msg += "       There is a problem with this image/visit.\n"
                msg += "       Make sure the input data are not *_raw.fits files.\n"
                print(msg)
                sys.exit(1)
        # If there are no WCS solutions, the image could be bad (e.g., EXPTIME=0 or EXPFLAG="TDF-DOWN...")
        else:
            msg = "WARNING: There are no WCS solutions in the image {} in this visit.".format(filename)
            print(msg)
            update_msg += msg
            skip_image_list.append(filename)
            if image_type == 'GRISM':
                msg = "WARNING:    Skip and delete this image."
                print(msg)
                update_msg += msg
                # Delete the SVM FLT/FlC image as it has no updated WCS
                try:
                    os.remove(filename)
                    msg = "WARNING: Deleted image {}.".format(filename)
                    print(msg)
                    update_msg += msg
                except OSError:
                    pass
            else:
                msg = "WARNING:    Skip this image."
                print(msg)
                update_msg += msg

    if logfile:
        _updateTrlFile(logfile, update_msg)

    return image_wcs_set, skip_image_list, keyword_wcs_names_dict, image_dict

# ------------------------------------------------------------------------------


def update_active_wcs(filename, wcsname, logfile=None):
    """
    Utility to update the active/primary WCS solution

    This small utility updates the active/primary WCS solution for the input
    file with the WCS solution indicted by the input parameter "wcsname"

    Parameters
    ----------
    filename : str
        Input/Output SVM FLT/FLC filename - the file is updated in-place

    wcsname : str
        Name of the desired WCS active/primary solution to be set for the filename

    Returns
    -------
    None

    """
    update_msg = ""
    # For exposures with multiple science extensions (multiple chips),
    # generate a combined WCS
    num_sci_ext, extname = util.count_sci_extensions(filename)
    extname_list = [(extname, x + 1) for x in range(num_sci_ext)]

    hdu = fits.open(filename)

    # Check if the desired WCS solution is already the active solution
    # whereupon there is nothing to do
    key = wcsutil.altwcs.getKeyFromName(hdu['SCI', 1].header, wcsname)
    keyword_wcs_list = [key]  # Initialize to a default value

    # No need to keep this file handle open anymore
    hdu.close()
    del hdu

    # Case where the desired active solution is not the current active solution
    if key != ' ':
        # Get the distortion model identification of the desired active WCS solution
        tmp_wcsname = wcsname.split('-')[0]
        index = tmp_wcsname.upper().find('IDC_')
        idc_new_string = ''
        if index > -1:
            idc_new_string = tmp_wcsname[index:]

        # Get the headerlet HDRNAMES for comparison to the alternate WCS solutions
        headerlet_hdrnames = wcsutil.headerlet.get_headerlet_kw_names(filename, kw="HDRNAME")

        if headerlet_hdrnames:
            # Examine the alternate WCS solutions to determine if they will be auto-archived due
            # to a distortion model change.  The auto-archiving will happen when the desired WCS
            # solution is installed as the active solution - just deleting duplicates here pro-actively.
            wcs_key_dict = wcsutil.altwcs.wcsnames(filename, ext=1)
            for wkey, wname in wcs_key_dict.items():
                if wkey == ' ':
                    continue

                index = wname.upper().find(idc_new_string.upper())

                # No match so solution will be copied to a headerlet automatically when the new primary is set
                if index == -1 and wkey.upper() != 'O':
                    msg = "Archiving alternate WCS solution as a headerlet as necessary: {}\n".format(wname)
                    print(msg)
                    update_msg += msg

                    # Now check if the HDRNAME between this solution and a headerlet already exists
                    hdr_keyword = fits.getval(filename, 'HDRNAME{}'.format(wkey.upper()), ext=1)

                    # Solution already exists as a headerlet extension, so just delete it
                    if hdr_keyword in headerlet_hdrnames:
                        wcsutil.altwcs.deleteWCS(filename, extname_list, wcskey=wkey)

            # Get all the WCS solution names
            headerlet_wcsnames = wcsutil.headerlet.get_headerlet_kw_names(filename, kw="WCSNAME")
            keyword_wcs_list = list(wcs_key_dict.values())

            # Prepare to install a new active WCS, but need to do some checking first
            #
            # This returns the first matching instance
            hdrname = headerlet_hdrnames[headerlet_wcsnames.index(wcsname)]
            extensions = []
            extensions = wcsutil.headerlet.find_headerlet_HDUs(filename, hdrname=hdrname)

            # It is possible the hdrname is not unique, so need to delete the dups
            msg = ''
            for ext in reversed(extensions[1:]):
                wcsutil.headerlet.delete_headerlet(filename, hdrext=ext)
                msg += "Delete duplicate headerlet extension {} in filename {}.\n".format(ext, filename)

            msg += "Desired active WCS solution {} has an HDRNAME of {}.\n".format(wcsname, hdrname)
            print(msg)
            update_msg += msg

            # Finally, install the desired WCS as the active WCS solution
            # Is the source of the wcsname for this image from a headerlet extension
            # or from the alternate solutions in the header as the source dictates how
            # the WCS will be made the active WCS. If available, restore a WCS solution
            # from the headerlet extension.
            try:
                wcsutil.headerlet.restore_from_headerlet(filename, hdrname=hdrname, force=True)
                # Update value of nmatches based on headerlet
                fhdu = fits.open(filename, mode='update')
                for sciext in range(1, num_sci_ext+1):
                    nm = fhdu[extensions[0]].header['nmatch'] if 'nmatch' in fhdu[extensions[0]].header else 0
                    fhdu[(extname, sciext)].header['nmatches'] = nm
                fhdu.close()
                del fhdu
            except ValueError as err:
                msg = "WARNING: Trapped ValueError - attempting recovery: {}\n".format(str(err))
                print(msg)
                update_msg += msg

                found_string = [i for i in keyword_wcs_list if wcsname == i]
                if found_string:
                    wcsutil.altwcs.restoreWCS(filename, ext=extname_list, wcsname=found_string[0])
                else:
                    msg = "WARNING: Could not restore the common WCS, {}, as the active WCS in this file {}.\n".format(wcsname, filename)
                    print(msg)
                    update_msg += msg
            except AssertionError:
                _, _, tb = sys.exc_info()
                tb_info = traceback.extract_tb(tb)
                _, _, _, text = tb_info[-1]
                msg = "WARNING: Trapped AssertionError: {}.\n".format(text)
                msg += "         Could not restore the common WCS, {}, as the active WCS in this file {}.\n".format(wcsname, filename)
                print(msg)
                update_msg += msg
        else:
            found_string = [i for i in keyword_wcs_list if wcsname == i]
            if found_string:
                wcsutil.altwcs.restoreWCS(filename, ext=extname_list, wcsname=found_string[0])
            else:
                msg = "WARNING: Could not restore the common WCS from alternate WCS solutions, {},\n".format(wcsname)
                msg += "as the active WCS in this file {}.\n".format(filename)
                print(msg)
                update_msg += msg
    else:
        msg = "No need to update active WCS solution of {} for {} as it is already the active solution.\n".format(wcsname, filename)
        print(msg)
        update_msg += msg

    if logfile:
        _updateTrlFile(logfile, update_msg)

def _update_wcs_fit_keywords(fltfile, flcfile):
    """Update the header of the FLT file with the a posteriori fit results"""
    fit_kws_sci = [('RMS_RA', -1.0), ('RMS_DEC', -1.0),
                   ('CRDER1', -1.0), ('CRDER2', -1.0),
                   ('NMATCHES', 0), ('FITGEOM', 'N/A'),
                   ('HDRNAME', '')]

    hdulist = fits.open(fltfile, mode='update')
    hdulist_flc = fits.open(flcfile)  # source header

    if 'HISTORY' in hdulist[0].header:
        after_kw = None
        before_kw = 'HISTORY'
    elif 'ASN_MTYP' in hdulist[0].header:
        after_kw = 'ASN_MTYP'
        before_kw = None
    else:
        after_kw = hdulist[0].header.cards[-1][0]
        before_kw = None

    hdulist[0].header.set('UPWCSVER', value=hdulist_flc[0].header['UPWCSVER'],
                          comment="Version of STWCS used to update the WCS",
                          after=after_kw, before=before_kw)
    hdulist[0].header.set('PYWCSVER', value=hdulist_flc[0].header['PYWCSVER'],
                          comment="Version of Astropy used to update the WCS",
                          after='UPWCSVER')

    num_sci_ext = amutils.countExtn(hdulist)
    for extnum in range(1, num_sci_ext+1):
        sci_extn = ('SCI', extnum)
        for kw in fit_kws_sci:
            src_hdr = hdulist_flc[sci_extn].header
            # Account for situations where FLC was not updated due to
            # EXPTIME=0 or other condition
            src_val = src_hdr[kw[0]] if kw[0] in src_hdr else kw[1]
            hdulist[sci_extn].header.set(kw[0], value=src_val, after='WCSNAME')

    hdulist.flush()
    hdulist.close()
    hdulist_flc.close()
    del hdulist
    del hdulist_flc

def _lowerAsn(asnfile):
    """ Create a copy of the original asn file and change
        the case of all members to lower-case.
    """
    # Start by creating a new name for the ASN table
    _indx = asnfile.find('_asn.fits')
    _new_asn = asnfile[:_indx] + '_pipeline' + asnfile[_indx:]
    if os.path.exists(_new_asn):
        os.remove(_new_asn)
    # copy original ASN table to new table
    shutil.copy(asnfile, _new_asn)

    # Open up the new copy and convert all MEMNAME's to lower-case
    fasn = fits.open(_new_asn, mode='update', memmap=False)
    for i in range(len(fasn[1].data)):
        fasn[1].data[i].setfield('MEMNAME', fasn[1].data[i].field('MEMNAME').lower())
    fasn.close()

    return _new_asn

def _updateTrlFile(trlfile, trl_lines):
    tmptrl = trlfile.replace('.tra', '_tmp.tra')

    print(trl_lines)

    # Write message out to temp file and append it to full trailer file
    ftmp = open(tmptrl, 'w')
    ftmp.writelines(trl_lines)
    ftmp.close()
    _appendTrlFile(trlfile, tmptrl)


def _appendTrlFile(trlfile, drizfile):
    """ Append drizfile to already existing trlfile from CALXXX.
    """
    if not os.path.exists(drizfile):
        return
    # Open already existing CALWF3 trailer file for appending
    ftrl = open(trlfile, 'a')
    # Open astrodrizzle trailer file
    fdriz = open(drizfile)

    # Read in drizzle comments
    _dlines = fdriz.readlines()

    # Append them to CALWF3 trailer file
    ftrl.writelines(_dlines)

    # Close all files
    ftrl.close()
    fdriz.close()

    try:
        # Now, clean up astrodrizzle trailer file
        os.remove(drizfile)
    except Exception:
        pass

def _timestamp(_process_name):
    """Create formatted time string recognizable by OPUS."""
    _prefix = time.strftime("%Y%j%H%M%S-I-----", time.localtime())
    _lenstr = 60 - len(_process_name)
    return _prefix + _process_name + (_lenstr * '-') + '\n'

def _getTime():
    # Format time values for keywords IRAF-TLM, and DATE
    _ltime = time.localtime(time.time())
    time_str = time.strftime('%H:%M:%S (%d-%b-%Y)', _ltime)

    return time_str


# Functions used to manage processing in a separate directory/ramdisk
def _createWorkingDir(rootdir, input):
    """
    Create a working directory based on input name under the parent directory specified as rootdir
    """
    # extract rootname from input
    rootname = input[:input.find('_')]
    newdir = os.path.join(rootdir, rootname)
    if not os.path.exists(newdir):
        os.mkdir(newdir)
    return newdir

def _copyToNewWorkingDir(newdir, input):
    """ Copy input file and all related files necessary for processing to the new working directory.

        This function works in a greedy manner, in that all files associated
        with all inputs(have the same rootname) will be copied to the new
        working directory.
    """
    flist = []
    if '_asn.fits' in input:
        asndict = asnutil.readASNTable(input, None)
        flist.append(input[:input.find('_')])
        flist.extend(asndict['order'])
        flist.append(asndict['output'])
    else:
        flist.append(input[:input.find('_')])
    # copy all files related to these rootnames into new dir
    for rootname in flist:
        for fname in glob.glob(rootname + '*'):
            shutil.copy(fname, os.path.join(newdir, fname))

def _get_envvar_switch(envvar_name):
    # interpret envvar variable, if specified
    if envvar_name in os.environ:
        val = os.environ[envvar_name].lower()
        if val not in envvar_bool_dict:
            msg = "ERROR: invalid value for {}.".format(envvar_name)
            msg += "  \n    Valid Values: on, off, yes, no, true, false"
            raise ValueError(msg)
        switch_val = envvar_bool_dict[val]
    else:
        switch_val = None

    return switch_val

def _restoreResults(newdir, origdir):
    """ Move (not copy) all files from newdir back to the original directory
    """
    for fname in glob.glob(os.path.join(newdir, '*')):
        shutil.move(fname, os.path.join(origdir, os.path.basename(fname)))

def _removeWorkingDir(newdir):
    """ Delete working directory
    """
    os.rmdir(newdir)

def rmtree2(path, n=3):
    """Wrapper around shutil.rmtree to make it more robust when used on NFS mounted file systems."""
    ok = False
    for i in range(0, n):
        try:
            shutil.rmtree(path, ignore_errors=False, onerror=handle_remove_readonly)
            ok = True
            break
        except OSError as err:
            print("Failed to remove path %s with shutil.rmtree at attempt %d: %s" % (path, n, err))
        time.sleep(3)

    if not ok:
        print("Failed to remove path %s with shutil.rmtree, even after %d attempts.".format(path, n))
        raise OSError
    else:
        print("Path %s successfully removed." % path)

def handle_remove_readonly(func, path, exc):
    excvalue = exc[1]
    if func in (os.rmdir, os.remove) and excvalue.errno == errno.EACCES:
        os.chmod(path, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)  # 0777
        func(path)
    else:
        raise

# Functions to support execution from the shell.
def main():

    import getopt

    try:
        optlist, args = getopt.getopt(sys.argv[1:], 'bdahfgin:')
    except getopt.error as e:
        print(str(e))
        print(__doc__)
        print("\t", __version__)

    # initialize default values
    help = 0
    force = False
    newdir = None
    inmemory = False
    num_cores = None
    headerlets = True
    align_to_gaia = True
    debug = False
    force_alignment = False

    # read options
    for opt, value in optlist:
        if opt == "-g":
            align_to_gaia = False
        if opt == "-d":
            debug = True
        if opt == "-a":
            force_alignment = False
        if opt == "-h":
            help = 1
        if opt == "-f":
            force = True
        if opt == "-i":
            inmemory = True
        if opt == '-n':
            if not value.isdigit():
                print('ERROR: num_cores value must be an integer!')
                raise ValueError
            num_cores = int(value)
        if opt == '-b':
            # turn off writing headerlets
            headerlets = False
    if len(args) < 1:
        print("syntax: runastrodriz.py [-fhibng] inputFilename [newpath]")
        sys.exit()
    if len(args) > 1:
        newdir = args[-1]
    if (help):
        print(__doc__)
        print("\t", __version__)
    else:
        try:
            process(args[0], force=force, newpath=newdir, num_cores=num_cores,
                    inmemory=inmemory, headerlets=headerlets,
                    align_to_gaia=align_to_gaia, force_alignment=force_alignment, debug=debug)

        except Exception as errorobj:
            print(str(errorobj))
            print("ERROR: Cannot run astrodrizzle on %s." % sys.argv[1])
            raise Exception(str(errorobj))

    sys.exit()


if __name__ == "__main__":
    main()
