#!/usr/bin/env python

""" runastrodriz.py - Module to control operation of astrodrizzle to
        remove distortion and combine HST images in the pipeline.

:License: :doc:`/LICENSE`

USAGE: runastrodriz.py [-bdahfginv] inputFilename [newpath]

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

See the main() at the end of this module for all available options and a
terse description.

Additional control over whether or not to attempt to align to an external
astrometric catalog, such as GAIA, is provided through the use of the
environment variables:

    - ASTROMETRY_COMPUTE_APOSTERIORI : Turn on/off alignment step.
      This environment variable will ALWAYS override any setting of the '-g' switch.
      Values (case-insensitive) can be 'on', 'off', 'yes', 'no'.

    - ASTROMETRY_APPLY_APRIORI : Replaces the obsolete ASTROMETRY_STEP_CONTROL
      variable used by ``stwcs.updatewcs`` to control whether or not a priori WCS
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

W.J. Hack  18 Oct 2019: Implemented multi-stage alignment with verification

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
import astropy.units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
import photutils

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
from drizzlepac.haputils import analyze
from drizzlepac.haputils import astrometric_utils as amutils
from drizzlepac.haputils import cell_utils
from drizzlepac.haputils import processing_utils
from drizzlepac import util
from drizzlepac import mdzhandler
from drizzlepac import updatehdr
from drizzlepac.haputils import quality_analysis as qa
from drizzlepac import wcs_functions
# for WFPC2 support
from drizzlepac.haputils import config_utils
from drizzlepac import wfpc2Data
from drizzlepac import photeq

from drizzlepac import __version__


__taskname__ = "runastrodriz"
package_level_logger = logging.getLogger('drizzlepac')
log = logging.getLogger(f'drizzlepac.{__taskname__}')

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

FILTER_NAMES = {'WFPC2': ['FILTNAM1', 'FILTNAM2'],
                'ACS': ['FILTER1', 'FILTER2'],
                'WFC3': ['FILTER']}

__trlmarker__ = '*** Drizzlepac Processing Version ' + drizzlepac.__version__ + '***'

envvar_compute_name = 'ASTROMETRY_COMPUTE_APOSTERIORI'
# ASTROMETRY_APPLY_APRIORI supersedes a previously existing environment variable
envvar_new_apriori_name = "ASTROMETRY_APPLY_APRIORI"
envvar_qa_stats_name = "PIPELINE_QUALITY_TESTING"
envvar_reset_idctab_name = "PIPELINE_RESET_IDCTAB"

# Order of preference for common WCS solutions - brute force choice of the GAIA catalog preference.
# Enforce the chosen order for the WCS preferential solution for all exposures in a dataset.
# REL fit method is better than IMG, and GAIA*3 catalog is better than GAIADR2 or GAIADR1.
wcs_preference = ['IDC_?????????-FIT_REL_GAIA*3', 'IDC_?????????-FIT_IMG_GAIA*3', 'IDC_?????????-FIT_REL_GAIADR2',
                  'IDC_?????????-FIT_IMG_GAIADR2', 'IDC_?????????-FIT_REL_GAIADR1', 'IDC_?????????-FIT_IMG_GAIADR1',
                  'IDC_?????????-GSC240', 'IDC_?????????']

# Primary user interface
def process(inFile, force=False, newpath=None, num_cores=None, inmemory=True,
            headerlets=True, align_to_gaia=True, force_alignment=False,
            do_verify_guiding=False, debug=False, make_manifest=False):
    """
    Run astrodrizzle on input file/ASN table 
    using default values for astrodrizzle parameters.

    Parameters
    ----------
    inFile : str
        Input.
    force : bool
        If True make a product even if verification fails.
    newpath : str
        If provided processing will be done in ``newpath`` directory.
    num_cores : int
        Number of cores to use.
    inmemory : bool
        Process the inputs in memory (default=True).
    headerlets : bool
        Whether to generate headerlets or not.
    align_to_gaia : bool
        Whether to align to an external catalog.
    force_alignment : bool
        When ``False`` (default) turn off alignment to GAIA if there were problems with guiding.
    do_verify_guiding : bool
        Whether to check for guiding problems. If True and force=True,
        a product will be generated even if issues with guiding are found.
    debug : bool
        Debug logging on/off. Additional debug logging is saved in a json file.
    make_manifest : bool
        Whether to generate a MANIFEST file.

    """

    # determine log name
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
    _manifest_filename = _trlroot + '_manifest.txt'
    _calfiles_flc = []

    # remove previous file and stream handlers
    log.handlers.clear()
    log.parent.handlers.clear()
    
    if debug:
        default_log_level = logging.DEBUG
        formatter = logging.Formatter('[%(levelname)s:%(name)s] %(message)s')
    else:
        default_log_level = logging.INFO
        formatter = logging.Formatter('[%(levelname)-8s] %(message)s')
    
    file_handler = logging.FileHandler(f'{_trlfile}')
    stream_handler = logging.StreamHandler(sys.stdout)
    
    file_handler.setLevel(default_log_level)
    stream_handler.setLevel(default_log_level)
    file_handler.setFormatter(formatter)
    stream_handler.setFormatter(formatter)
    package_level_logger.addHandler(file_handler)
    package_level_logger.addHandler(stream_handler)

    msg = (f"""Calibration pipeline processing of {inFile} started.
                 {__trlmarker__} 
                 drizzlepac version {drizzlepac.__version__}
                 tweakwcs version {tweakwcs.__version__}
                 stwcs version {stwcs.__version__}
                 numpy version {np.__version__}
                 photutils version {photutils.__version__}""")
    log.debug(msg)

    init_time = time.time()
    pipeline_pars = PIPELINE_PARS.copy()
    _verify = True  # Switch to control whether to verify alignment or not
    manifest_list = []
    
    # interpret envvar variable, if specified
    align_to_gaia = util.get_envvar_switch(
        envvar_compute_name, default=align_to_gaia, description="'align to gaia'"
    )

    # Insure os.environ ALWAYS contains an entry for envvar_new_apriori_name
    # and it will default to being 'on'
    align_with_apriori = util.get_envvar_switch(
        envvar_new_apriori_name, default=True, description="'align with apriori'"
    )

    # Add support for environment variable switch to automatically
    # reset IDCTAB in FLT/FLC files if different from IDCTAB in RAW files.
    reset_idctab_switch = util.get_envvar_switch(
        envvar_reset_idctab_name,
        default=False,
        description="'reset idctab in flt if different from raw'",
    )

    if headerlets or align_to_gaia:
        from stwcs.wcsutil import headerlet

    # Open the input file
    try:
        # Make sure given filename is complete and exists...
        inFilename = fileutil.buildRootname(inFile, ext=['.fits'])
        if not os.path.exists(inFilename):
            log.error(f"ERROR: Input file - {inFilename} - does not exist.")
            return
    except TypeError:
        log.error("ERROR: Appropriate input file could not be found.")
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
    wfpc2_input = infile_inst == 'WFPC2'
    if not wfpc2_input:
        raw_suffix = '_raw.fits'
        goodpix_name = 'NGOODPIX'
    else:
        # start by turning off 'verify_guiding' since loss of lock
        # was exceedingly rare and noted in the quality keywords
        # which get checked in '_analyze_exposure'.
        do_verify_guiding = False
        # Convert input c0m file into compatible flt file
        if 'd0m' in inFilename:
            inFilename = inFilename.replace('d0m', 'c0m')
        # This returns the name of the _flt.fits file that was created
        inFilename = wfpc2Data.wfpc2_to_flt(inFilename)
        if _analyze_exposure(inFilename):
            # Update header of WFPC2 data to use latest reference
            # files from CRDS
            log.debug(f"Updating distortion reference files for: {inFilename}")
            wfpc2Data.apply_bestrefs(inFilename)
            try:
                photeq.photeq(files=inFilename, ref_phot_ext=3, readonly=False)
            except Exception as err:
                log.error(err)
                log.warning(f"WARNING: PHOTEQ was unable to run on {inFilename}")

            raw_suffix = '_d0m.fits'
            goodpix_name = 'GPIXELS'
        else:
            # Remove FLT file created here, since the calibration file can NOT be aligned or drizzled
            os.remove(inFilename)
            log.error("ERROR: Inappropriate input file.  Deleting converted WFPC2 FLT file.")

            # write out manifest file, if requested
            mani_filename = inFile[:inFile.find('_d0m')] + '_manifest.txt'
            if make_manifest:
                if os.path.exists(mani_filename):
                    os.remove(mani_filename)
                with open(mani_filename, 'w') as fout:
                    _ = [fout.write(f"{fname}\n") for fname in manifest_list]
                log.debug(f"Created manifest file: {mani_filename}")
            sys.exit(analyze.Ret_code.NO_VIABLE_DATA.value)

    infile_det = fits.getval(inFilename, 'detector')
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
        _indx = inFilename.find(raw_suffix)
        if _indx < 0: _indx = len(inFilename)
        # ... and build the CALXXX product rootname.
        _mname = fileutil.buildRootname(inFilename[:_indx], ext=cal_ext)

        _cal_prodname = inFilename[:_indx]
        # Reset inFilename to correspond to appropriate input for
        # drizzle: calibrated product name.
        inFilename = _mname

        if _mname is None:
            errorMsg = 'Could not find calibrated product!'
            raise Exception(errorMsg)

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
            _fimg = fits.open(fileutil.buildRootname(_mname, ext=[raw_suffix]), memmap=False)
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
        reset_idctab_kw(_calfiles, _calfiles_flc)

    # Add S_REGION keyword to input files regardless of whether DRIZCORR is turned on
    for f in _calfiles+_calfiles_flc:
        processing_utils.compute_sregion(f)

    # If we no longer have any valid images to process due to guiding problems,
    # set drizcorr to OMIT and finish processing gracefully.
    if len(_calfiles) == 0:
        dcorr = 'OMIT'    
    if dcorr == 'PERFORM':
        # Run updatewcs on each list of images to define pipeline default WCS
        # based on latest distortion models
        # Always apply latest distortion to replace pipeline-default OPUS WCS
        # for successful creation of updated headerlets for the cases where
        # all inputs having EXPTIME==0 (for example) or guiding is bad.
        updatewcs.updatewcs(_calfiles, use_db=False, checkfiles=False)
        if _calfiles_flc:
            updatewcs.updatewcs(_calfiles_flc, use_db=False, checkfiles=False)

        # adds skycell keyword to science header of all flt(c) and drz(c,w) files.
        # the SKYCELL value for IPPPSSOOT and SVM products may be different as the
        # current computation is based upon the WCSNAME of the input images.
        for filename in _calfiles+_calfiles_flc:
            processing_utils.add_skycell_to_header(filename)

        # Check to see whether or not guide star failure affected these observations
        # They would show up as images all sources streaked as if taken in SCAN mode or with a GRISM
        #
        # Note: This functionality is intentionally turned off for pipeline processing at this time.
        # However, a user may wish to invoke this functionality which is controlled by the
        # parameter "do_verify_guiding".
        if do_verify_guiding:
            for fltimg in _calfiles:
                flcimg = fltimg.replace('_flt.fits', '_flc.fits')
                guide_img = flcimg if os.path.exists(flcimg) else fltimg
                # We want to use the FLC image, if possible, to avoid any
                # possible detection of CTE tails as false guide-star trailing lines
                bad_guiding = analyze.verify_guiding(guide_img)
                if bad_guiding:
                    # If the user did not specify they wanted a drizzle product no matter what...
                    if not force:
                        # Remove the affected image(s) from further processing
                        # except when user specifies they want to make a product no matter what (force=True)
                        _calfiles.remove(fltimg)
                        if os.path.exists(flcimg):
                            _calfiles_flc.remove(flcimg)
                        # If ANY input exposure has bad guiding, none of the data can be
                        # trusted.  However, only allow alignment if the user forces the alignment.
                    if not force_alignment:
                        # Turn off any alignment to GAIA
                        # After all, if one exposure has bad guiding,
                        # none of the WCS information can be trusted.
                        align_to_gaia = False
                        align_with_apriori = False

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
        inst_mode = f"{infile_inst}/{infile_det}"
        _good_images = [f for f in _calfiles if fits.getval(f, 'exptime') > 0.]
        _good_images = [f for f in _good_images if fits.getval(f, goodpix_name, ext=("SCI", 1)) > 0.]
        if len(_good_images) == 0:
            if len(_calfiles)==0:
                log.error("ERROR: No valid data found in input exposures.")
                sys.exit(analyze.Ret_code.NO_VIABLE_DATA.value)
            else:
                _good_images = _calfiles
                # Turn off alignment to GAIA since there is no valid data in the exposure.
                align_to_gaia = False

        if not wfpc2_input:
            adriz_pars = mdzhandler.getMdriztabParameters(_good_images)
        else:
            default_pars = config_utils.get_wfpc2_pars(_good_images)
            adriz_pars = default_pars['astrodrizzle'].outpars

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

        # Integrate user-specified drizzle parameters into pipeline_pars
        log.info('Starting alignment with bad-pixel identification')
        log.info(__trlmarker__)

        if align_with_apriori or force_alignment or align_to_gaia:
            # Generate initial default products and perform verification
            align_dicts, align_table = verify_alignment(_inlist,
                                                        _calfiles, _calfiles_flc,
                                                        _trlfile,
                                                        tmpdir=None, debug=debug,
                                                        force_alignment=force_alignment,
                                                        find_crs=True, **adriz_pars)
        if align_with_apriori:
            log.info('Starting alignment with a priori solutions')
            log.info(__trlmarker__)
            if align_dicts is not None:
                find_crs = not align_dicts[0]['alignment_verified']
            else:
                find_crs = False

            # run updatewcs with use_db=True to insure all products have
            # have a priori solutions as extensions
            # FIX: This should probably only be done in the apriori sub-directory!
            updatewcs.updatewcs(_calfiles, verbose=True)
            for _file in _calfiles:
                confirm_aposteriori_hdrlets(_file)

                # verify WCS solution (CRVALs) near target coordinates
                warning_separation_threshold = 0.05*u.deg # value determine by Rick White from experience
                header_ex0 = fits.getheader(_file, ext=0)
                header_ex1 = fits.getheader(_file, ext=1)
                targ_pos = SkyCoord(header_ex0['RA_TARG']*u.deg, header_ex0['DEC_TARG']*u.deg)
                wcs_pos = SkyCoord(header_ex1['CRVAL1']*u.deg, header_ex1['CRVAL2']*u.deg)
                sep = wcs_pos.separation(targ_pos)
                if sep > warning_separation_threshold:
                    log.warning(f'WARNING: WCS reference pixel is {sep.value:.2f} degrees '+
                                   'from target position. The astrometry database solution is suspect!')

            log.debug(f"Adding apriori WCS solutions to {_calfiles}")
            log.debug(verify_gaia_wcsnames(_calfiles))
            _wnames_calfiles = [(c, fits.getval(c, 'wcsname', ext=1)) for c in _calfiles]
            log.debug("Verifying apriori WCSNAMEs:")
            for (_cname, _wname) in _wnames_calfiles:
                log.debug(f"   {_cname}: {_wname}")
            if _calfiles_flc:
                log.debug(f"Adding apriori WCS solutions to {_calfiles_flc}")
                updatewcs.updatewcs(_calfiles_flc)
                for _file in _calfiles_flc:
                    confirm_aposteriori_hdrlets(_file)

                log.debug(verify_gaia_wcsnames(_calfiles_flc))

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
                log.debug("ERROR in applying a priori solution.")
            if align_apriori is None or (not align_apriori[0]['alignment_verified']):
                log.debug("Resetting WCS to pipeline-default solutions...")
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
                if align_dicts[0]["alignment_quality"] == 0:
                    log.debug("A priori alignment SUCCESSFUL.")
                if align_dicts[0]["alignment_quality"] == 1:
                    log.debug("A priori alignment potentially compromised.")
                    log.debug("Please review final product!")
                if align_dicts[0]["alignment_quality"] > 1:
                    log.debug(
                        "A priori alignment FAILED! No a priori astrometry correction applied."
                    )

        aposteriori_table=None
        if align_to_gaia:
            log.info("Starting a posteriori alignment")
            log.info(__trlmarker__)

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
                align_qual = align_dicts[0]["alignment_quality"]
                if align_qual == 0:
                    log.debug("A posteriori alignment SUCCESSFUL.")
                elif align_qual == 1:
                    log.debug(
                        "A posteriori alignment potentially COMPROMISED with bad focus."
                    )
                    log.debug("  Please review final product for alignment!")
                elif align_dicts[0]["alignment_quality"] == 2:
                    log.debug(
                        "A posteriori alignment potentially COMPROMISED."
                    )
                    log.debug("Please review final product!")
                else:
                    log.debug(
                        "A posteriori alignment FAILED! No a posteriori astrometry correction applied."
                    )

        log.info(
            "Creating final combined,corrected product based on best alignment"
        )
        log.info(__trlmarker__)

        # Generate final pipeline products based on 'best' alignment
        pipeline_pars['in_memory'] = inmemory
        pipeline_pars['clean'] = True
        pipeline_pars['num_cores'] = num_cores
        if wfpc2_input:
            pipeline_pars.update(adriz_pars)
            pipeline_pars['mdriztab'] = False

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
            if drz_product != _plower:
                os.rename(drz_product, _plower)
                drz_output = _plower
            else:
                drz_output = drz_product
            # keep track of output files
            manifest_list.append(drz_output)
            # update s_region header keyword only if the drz_output
            # contains actual data, not just header information.
            if asn_dicts is not None:
                processing_utils.compute_sregion(drz_output)

    else:
        # Create default trailer file messages when astrodrizzle is not
        # run on a file.  This will typically apply only to BIAS,DARK
        # and other reference images.
        # Start by building up the message...
        log.info("astrodrizzle skipped ")
        log.info(__trlmarker__)
        log.info(f"{_getTime()}: astrodrizzle processing not requested for {inFile}.")
        log.info("       astrodrizzle will not be run at this time.")

    # If we created a new ASN table, we need to remove it
    if _new_asn is not None:
        for _name in _new_asn:
            fileutil.removeFile(_name)

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
    # Keep track of input files that have been updated
    manifest_list.extend(_calfiles + _calfiles_flc)

    # If headerlets have already been written out by alignment code,
    # do NOT write out this version of the headerlets
    if headerlets:
        # Generate headerlets for each updated FLT image
        log.debug("Writing Headerlets started")
        for fname in _calfiles:
            log.debug(f"Creating new headerlet from {fname}")
            frootname = fileutil.buildNewRootname(fname)
            hname = f"{frootname}_flt_hlet.fits"
            # Write out headerlet file used by astrodrizzle, however,
            # do not overwrite any that was already written out by align
            if not os.path.exists(hname):
                log.debug(f"Created Headerlet file {hname} ")
                try:
                    wcsname = fits.getval(fname, 'wcsname', ext=1)
                    wcstype = updatehdr.interpret_wcsname_type(wcsname)
                    hdrname = f"{fname.replace('.fits', '')}_hlet.fits"
                    headerlet.write_headerlet(fname, hdrname, output='flt',
                                              wcskey='PRIMARY',
                                              author="OPUS",
                                              descrip=wcstype,
                                              attach=False,
                                              clobber=True,
                                              logging=False)
                    # Keep track of headerlet files written out to disk
                    manifest_list.append(hdrname)
                except ValueError:
                    log.debug(f"SKIPPED: Headerlet not created for {fname}")
                    # update trailer file to log creation of headerlet files

        log.debug("Writing Headerlets completed")

    # Keep track of headerlet files written out to disk.
    # Those headerlets would have been written out by 'updatewcs.updatewcs()'
    # using a SHA256 hash for the WCSNAME portion of the headerlet filename.
    # The headerlets written out with that naming convention are the only ones
    # which need to be added to the manifest, since they represent the only WCS
    # solutions which are new for the exposure and, thus, the ones which need to
    # be saved in the astrometry database during pipeline processing.
    for fname in _calfiles:
        # Look for HLET files new to astrometry database written out by STWCS
        updatewcs_hdrlets = sorted(glob.glob(f'{fname.replace(".fits", "")}_??????_hlet.fits'))
        manifest_list.extend(updatewcs_hdrlets)
        # Look for HLET representing the WCS used in the updated FLT/FLC file.
        final_hdrlet = glob.glob(f'{fname.replace(".fits", "")}_hlet.fits')
        manifest_list.extend(final_hdrlet)

    # add trailer file to list of output products
    manifest_list.append(_trlfile)

    # Remove secondary log files for good...
    logging.shutdown()

    # write out manifest file, if requested
    if make_manifest:
        if os.path.exists(_manifest_filename):
            os.remove(_manifest_filename)
        with open(_manifest_filename, 'w') as fout:
            _ = [fout.write(f"{fname}\n") for fname in manifest_list]
        log.debug(f"Created manifest file: {_manifest_filename}")

    # delete log files generated by alignment code
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
        log.debug("Files still open for this process include: ")
        log.debug([ofile.path for ofile in Process().open_files()])

    if not debug:
        # clean up any temporary ref files created for WFPC2 data
        if wfpc2_input:
            d2im_file = f"{inFilename.split('_')[0]}_flt_d2im.fits"
            if os.path.exists(d2im_file):
                os.remove(d2im_file)
        try:
            # Remove all temp sub-directories now that we are done
            for sd in sub_dirs:
                if os.path.exists(sd):
                    rmtree2(sd)
        except Exception:
            # If we are unable to remove any of these sub-directories,
            # leave them for the user or calling routine/pipeline to clean up.
            log.warning(f"WARNING: Unable to remove any or all of these sub-directories: \n{sub_dirs}")
            if Process is not None:
                log.debug("Files still open at this time include: ")
                log.debug([ofile.path for ofile in Process().open_files()])
            pass

    end_time = _getTime()
    _delta_time = time.time() - init_time
    log.debug(f"{end_time}: Finished processing {inFilename} in {_delta_time:.2f} seconds \n")
    log.debug("astrodrizzle completed")

    # Look to see whether we have products which can be evaluated
    # wcsname = fits.getval(drz_products[0], 'wcsname', ext=1)

    # interpret envvar variable, if specified
    qa_switch = util.get_envvar_switch(envvar_qa_stats_name, description="'QA statistics'", default=False)

    if qa_switch and dcorr == 'PERFORM':

        # Generate quality statistics for astrometry if specified
        calfiles = _calfiles_flc if _calfiles_flc else _calfiles
        qa.run_all(inFile, calfiles, catalogs=aposteriori_table)


def run_driz(inlist, trlfile, calfiles, mode='default-pipeline', verify_alignment=True,
            debug=False, good_bits=512, **pipeline_pars):
    """
    Parameters
    ----------

    inlist : list
        List of input files to process
    trlfile : str
        Trailer file.
    calfiles : list of str
        Original file names of input
    mode : str
        One of ['apriori', 'aposteriori', 'default-pipeline']
    verify_alignment : bool
    debug : bool
        Debug logging on/off
    good_bits : int
        Good bits to use in drizzle
    pipeline_pars : dict
        Drizzle parameters
    """

    import drizzlepac
    pyver = drizzlepac.astrodrizzle.__version__
    drz_products = []
    focus_dicts = []
    diff_dicts = OrderedDict()

    pipeline_pars['runfile'] = trlfile.replace('.tra', '_pydriz')
    drizlog = pipeline_pars['runfile'] + ".log"  # the '.log' gets added automatically by astrodrizzle
    _trlmsg = "In run_driz"
    _trlmsg += "inlist is, %s\n" % (inlist,)
    print(_trlmsg)
    for infile in inlist:  # Run astrodrizzle for all inputs

        asndict, ivmlist, drz_product = processInput.process_input(infile, updatewcs=False,
                                                        preserve=False,
                                                        overwrite=False)
        del ivmlist
        calfiles = asndict['original_file_names'] if asndict is not None else calfiles

        drz_products.append(drz_product)

        # Create trailer marker message for start of astrodrizzle processing
        log.info('astrodrizzle started ')
        log.info(f'{_getTime()}: Processing {infile} with astrodrizzle Version {pyver}')
        log.info(__trlmarker__)

        try:
            drizzlepac.astrodrizzle.AstroDrizzle(input=infile, configobj=None,
                                                 **pipeline_pars)

            # Edit trailer file name since 'runastrodriz' copies what astrodrizzle used
            # to another file...
            if os.path.exists(drz_product):
                fits.setval(drz_product, 'DRIZPARS', value=trlfile)

            # util.end_logging(drizlog)

        except Exception as errorobj:
            log.error(f"ERROR: Could not complete astrodrizzle processing of {infile}.")
            log.error(str(errorobj))
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

            instr_det = f"{fits.getval(sfile, 'instrume')}/{fits.getval(sfile, 'detector')}"
            focus_sigma = focus_pars[instr_det]['sigma']
            log.debug(f"Measuring similarity and focus for: \n{single_files} \n    {drz_product}")
            focus_dicts.append(amutils.build_focus_dict(single_files, drz_product, sigma=focus_sigma))
            if debug:
                json_name = drz_product.replace('.fits', f'_{mode}_focus.json')
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
    log.debug(f"Updating trailer file {trlfile} with astrodrizzle comments.")
    drizlog_copy = drizlog.replace('.log', '_copy.log')
    if os.path.exists(drizlog):
        shutil.copy(drizlog, drizlog_copy)

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
    """

    Parameters
    ----------
    inlist : list
        List of files to process
    calfiles : list
    calfiles_flc : list
    trlfile : str
        Trailer file name
    find_crs : bool
        Whether to run code detecting cosmic rays. Default is True.
    tmpdir : str
        Temporary directory for processing mode.
    debug : bool
        Debug logging on/off.
    good_bits : int
        Good bits to use with drizzle.
    alignment_mode : str
        Alignment mode, one of 'apriori', 'aposteriori', 'default-pipeline'.
    force_alignment : bool
        If False (default) turn off alignment to GAIA if there were problems with guiding.
    pipeline_pars : dict
        Drizzle parameters.
    """

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
        log.error(f"Invalid alignment mode {tmpdir} requested.")
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
        hdr0 = fits.getheader(alignfiles[0])
        inst = hdr0.get('instrume').lower()
        if inst == 'wfpc2' and 'detector' not in hdr0:
            det = 'pc'
        else:
            det = hdr0.get('detector').lower()

        if find_crs:
            log.debug("Resetting CRs ")
            # reset all DQ flags associated with CRs assuming previous attempts were inaccurate
            for f in alignfiles:
                log.debug(f"Resetting CR DQ bits for {f}")
                resetbits.reset_dq_bits(f, "4096,8192")
                sat_flags = 256 + 2048
        else:
            sat_flags = 256 + 2048 + 4096 + 8192

        # Perform any requested alignment here...
        if alignment_mode == 'aposteriori':
            # Create trailer marker message for start of align_to_GAIA processing
            log.debug("Align_to_GAIA started ")
            # Evaluate all input exposures and, if necessary, reset WCSs to a common WCS
            # This is necessary in order to avoid imprinting zero point differences between
            # coordinate systems into the final fit which would result in mis-alignment of the
            # sources in the final combined output images.

            if len(alignfiles) > 1:
                update_wcs_in_list(alignfiles)

            instdet_pars = align.get_default_pars(inst, det)

            # alignlog = trlfile.replace('.tra', '_align.log')
            try:

                full_table = align.perform_align(alignfiles,
                                                 catalog_list=instdet_pars['run_align']['catalog_list'],
                                                 num_sources=instdet_pars['general']['MAX_SOURCES_PER_CHIP'],
                                                 update_hdr_wcs=True, runfile=trlfile,
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
                            log.debug(f"""Successfully aligned {row['imageName']} 
                                               to {row['catalog']} astrometric frame""")
                        else:
                            log.debug(f"""Alignment only partially 
                                               successful for {row['imageName']}""")
                    else:
                        log.debug(f"""Could not align {row['imageName']} 
                                           to absolute astrometric frame""")
                        return None, None
            except Exception as err:
                # Something went wrong with alignment to GAIA, so report this in
                # trailer file
                log.warning("EXCEPTION encountered in align...")
                log.warning("No correction to absolute astrometric frame applied!")
                if 'aposteriori' not in repr(err):
                    traceback.print_exc()
                else:
                    log.warning(f"WARNING: {err}")
                return None, None

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
                        log.debug(f"""Applying headerlet {headerlet_file} 
                                           as Primary WCS to {fltfile}""")
                    else:
                        log.debug(f"""No absolute astrometric headerlet 
                                           applied to {fltfile}""")

            # Finally, append any further messages associated with alignement from this calling routine
            log.debug('Align_to_GAIA completed ')

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
        log.debug(f"Verification of {tmpmode} alignment started ")

        if focus_dicts is not None:
            # Only check focus on CTE corrected, when available
            align_focus = focus_dicts[-1] if 'drc' in focus_dicts[-1]['prodname'] else focus_dicts[0]

            pscale = HSTWCS(alignfiles[0], ext=1).pscale

            det_pars = align.get_default_pars(inst, det)['generate_source_catalogs']
            default_fwhm = det_pars['fwhmpsf'] / pscale
            align_fwhm = amutils.get_align_fwhm(align_focus, default_fwhm)

            if align_fwhm:
                log.debug(f"""align_fwhm: {align_focus['prodname']}
                                   [{align_focus['prod_pos'][1]},
                                   {align_focus['prod_pos'][0]}]
                                   ={align_fwhm:0.4f}pix""")

            # Interpret the overlap differences computed for this alignment
            dkeys = [k for k in diff_dicts.keys()]
            diff_verification, max_diff = amutils.evaluate_overlap_diffs(diff_dicts[dkeys[-1]])
            log.debug(f"""Fraction of sources matched: {fraction_matched} 
                               out of {num_sources} sources""")

            # For any borderline situation with alignment, perform an extra check on alignment
            if fraction_matched < 0.1 or -1 < num_sources < 10:
                focus_verification = amutils.evaluate_focus(align_focus)
                alignment_verified = True if (diff_verification and focus_verification) else False
                alignment_quality = 0 if alignment_verified else 3
            else:
                alignment_verified = diff_verification
                alignment_quality = 0 if diff_verification else 3

            if alignment_verified:
                log.debug(f"""Focus verification indicated that {tmpmode} 
                                   alignment SUCCEEDED.""")
            else:
                log.debug(f"""Focus verification indicated that {tmpmode} 
                                   alignment FAILED.""")
                log.debug(f"  Reverting to previously determined WCS alignment.")

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

            log.debug(f"Computing sim_indx for: {os.path.join(tmpdir, prodname)}")
            sim_indx = amutils.compute_similarity(alignprod, align_ref)
            align_sim_fail = sim_indx > 1

            if not align_sim_fail and alignment_verified:
                log.debug(f"Alignment appeared to SUCCEED based on similarity index of {sim_indx:0.4f}")
            else:
                log.debug(f"Alignment appeared to FAIL based on similarity index of {sim_indx:0.4f}")
                # alignment_verified = False
                alignment_quality += 3

        for fd in focus_dicts:
            fd['alignment_verified'] = alignment_verified
            fd['alignment_quality'] = alignment_quality

        # If CRs were identified, copy updated input files to main directory
        if tmpdir and alignment_verified:
            log.debug("Saving products with new alignment.")
            _ = [shutil.copy(f, parent_dir) for f in calfiles]
            if calfiles_flc:
                _ = [shutil.copy(f, parent_dir) for f in calfiles_flc]
            # Copy drizzle products to parent directory to replace 'less aligned' versions
            _ = [shutil.copy(f, parent_dir) for f in headerlet_files]

        log.debug('Verification of alignment completed ')
    finally:
        if tmpdir:
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
                    log.debug(f"""Updating WCSNAME of {f}[sci,{sciext + 1}] 
                                       for use of {catalog_name} catalog.""")
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
                                if (
                                    most_recent_wcs and wdate > most_recent_wcs[0]
                                ) or most_recent_wcs is None:
                                    most_recent_wcs = (wdate, hname)

                        if most_recent_wcs:
                            break

                    if most_recent_wcs:
                        # restore this WCS
                        log.debug(f"Restoring apriori WCS {wname} as primary WCS in {f}")
                        headerlet.restore_from_headerlet(fhdu,
                                                         force=True,
                                                         hdrname=most_recent_wcs[1],
                                                         archive=False)
                        # Ensure IDCSCALE is still present
                        if 'idcscale' not in fhdu[('sci', 1)].header:
                            msg += "Headerlet {} was missing IDCSCALE keyword".format(most_recent_wcs[1])
                            _update_idcscale(fhdu)
    return msg


def restore_pipeline_default(files):
    """Restore pipeline-default IDC_* WCS as PRIMARY WCS in all input files"""
    log.debug("Restoring pipeline-default WCS as PRIMARY WCS... using updatewcs.")
    updatewcs.updatewcs(files, use_db=False)
    # Remove HDRNAME, if added by some other code.
    #  This keyword only needs to be included in the headerlet file itself.
    for f in files:
        with fits.open(f, mode='update') as fhdu:
            num_sci = fileutil.countExtn(fhdu)
            for sciext in range(num_sci):
                if 'hdrname' in fhdu[('sci', sciext + 1)].header:
                    del fhdu[('sci', sciext + 1)].header['hdrname']


def reset_idctab_kw(files, files_flc):
    """Insure IDCTAB in files are the same as those in RAW files"""

    raw_files = [f.replace('_flt', '_raw') for f in files]
    # determine what IDCTAB should be in FLT files
    raw_idctab = fits.getval(raw_files[0], 'idctab')
    raw_wcsname = 'IDC_{}'.format(raw_idctab.split('_')[0].split('$')[1])
    log.debug(f"IDCTAB from RAW file: {raw_idctab}")

    log.debug("Insuring IDCTAB keywords are up-to-date")

    # Now, restore headerlet with this WCSNAME as primary WCS
    for flt in files+files_flc:
        flt_idctab = fits.getval(flt, 'idctab')
        if flt_idctab == raw_idctab:
            # We don't need to update anything
            continue
        log.debug(f"Updating IDCTAB {flt_idctab} in {flt}")

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


def update_wcs_in_list(exp_list):
    """Examine entries in the exposure list for inconsistent WCS solutions

    Note: The direct exposure image active WCS solutions will be modified in-place.

    Parameters
    ----------
    exp_list : list
        List of filenames for exposures to be evaluated.

    Returns
    -------
    Nothing.

    """

    log.debug("\n***** Processing List for Consistent WCS's *****")
    primary_wcsnames = set([fits.getval(fname, 'wcsname', ext=('SCI',1)) for fname in exp_list])
    if len(primary_wcsnames) == 1:
        log.debug(f"\nAll Primary WCS's confirmed as consistent as {primary_wcsnames}.")

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

    final_wcs_set, skip_direct_list, d_keyword_wcs_names_dict, direct_dict = collect_wcs_names(exp_list, 'DIRECT')
    log.debug(f"WCS solutions common to all viable direct images: {final_wcs_set}")


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
        log.debug(f"NO Common WCS solutions for images: {final_wcs_set}")
        log.debug(f" for the latest calibration IDCTAB: {primary_idctabs}")
    final_wcs_set = set(primary_wcs_set)

    # There is a preference for the active WCS for the viable images in the visit
    # Check the final_wcs_set for the preferential solutions
    match_list = []
    final_wcsname = ''
    for wcs_item in wcs_preference:
        match_list = fnmatch.filter(final_wcs_set, wcs_item)
        if match_list:
            final_wcsname = match_list[0]
            log.debug(f"Final WCS solution to use for all images: {final_wcsname}")
            break

    if final_wcsname:
        # Finally, if the image is not in a skip list, reset the primary WCS in all the images
        for filename in exp_list:
            if filename not in skip_direct_list:
                log.debug(f"Setting the primary WCS for direct image {filename} to {final_wcsname}.")
                update_active_wcs(filename, final_wcsname, logfile=logfile)
    else:
        # Do nothing
        pass

# ------------------------------------------------------------------------------


def collect_wcs_names(exp_list, image_type):
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
            log.debug(f"WCS solutions for file {filename} are {all_wcs_names}.")

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
                msg = """ERROR: There are no common WCS solutions with this image 
                {filename} and previously processed images. There is a problem with 
                this image/visit. Make sure the input data are not *_raw.fits files."""
                log.error(msg)
                sys.exit(1)
        # If there are no WCS solutions, the image could be bad (e.g., EXPTIME=0 or EXPFLAG="TDF-DOWN...")
        else:
            log.warning(f"WARNING: There are no WCS solutions in the image {filename} in this visit.")
            skip_image_list.append(filename)
            if image_type == 'GRISM':
                log.warning("WARNING:    Skip and delete this image.")
                # Delete the SVM FLT/FlC image as it has no updated WCS
                try:
                    os.remove(filename)
                    log.warning(f"WARNING: Deleted image {filename}.")
                except OSError:
                    pass
            else:
                log.warning(f"WARNING:    Skip this image.")

    return image_wcs_set, skip_image_list, keyword_wcs_names_dict, image_dict

# ------------------------------------------------------------------------------


def update_active_wcs(filename, wcsname):
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
                    log.debug(f"Archiving alternate WCS solution as a headerlet as necessary: {wname}")

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
            extensions = wcsutil.headerlet.find_headerlet_HDUs(filename, hdrname=hdrname)

            # It is possible the hdrname is not unique, so need to delete the dups
            msg = ''
            for ext in reversed(extensions[1:]):
                wcsutil.headerlet.delete_headerlet(filename, hdrext=ext)
                log.debug(f"Delete duplicate headerlet extension {ext} in filename {filename}.")

            log.debug(f"Desired active WCS solution {wcsname} has an HDRNAME of {hdrname}.")
            
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
                log.warning(f"WARNING: Trapped ValueError - attempting recovery: {str(err)}")

                found_string = [i for i in keyword_wcs_list if wcsname == i]
                if found_string:
                    wcsutil.altwcs.restoreWCS(filename, ext=extname_list, wcsname=found_string[0])
                else:
                    log.debug(f"""WARNING: Could not restore the common WCS, 
                                 {wcsname}, as the active WCS in this file {filename}.""")

            except AssertionError:
                _, _, tb = sys.exc_info()
                tb_info = traceback.extract_tb(tb)
                _, _, _, text = tb_info[-1]
                msg = f"""WARNING: Trapped AssertionError: {text}. Could not restore 
                the common WCS, {wcsname}, as the active WCS in this file {filename}."""
                log.warning(msg)
        else:
            found_string = [i for i in keyword_wcs_list if wcsname == i]
            if found_string:
                wcsutil.altwcs.restoreWCS(filename, ext=extname_list, wcsname=found_string[0])
            else:
                msg = f"""WARNING: Could not restore the common WCS from alternate 
                WCS solutions, {wcsname}, as the active WCS in this file {filename}."""
                log.warning(msg)
    else:
        msg = f"""No need to update active WCS solution of {wcsname} for {filename} 
        as it is already the active solution."""
        log.debug(msg)
        update_msg += msg


def confirm_aposteriori_hdrlets(filename):
    """Confirm that all the a posteriori headerlets are valid, and remove any that are invalid."""
    update_msg = ""
    num_sci_ext, extname = util.count_sci_extensions(filename)
    extname_list = [(extname, x + 1) for x in range(num_sci_ext)]

    hdu = fits.open(filename, mode='update')
    # Get a list of all HDRLET extensions
    hdrlet_extensions = wcsutil.headerlet.get_extname_extver_list(hdu, sciext='HDRLET')
    hdrlet_wcsnames = wcsutil.headerlet.get_headerlet_kw_names(hdu, 'wcsname')

    invalid_extns = []
    for wname, extname in zip(hdrlet_wcsnames, hdrlet_extensions):
        # Only evaluate a posteriori headerlets -- those with '-FIT' in WCSNAME
        if '-FIT' in wname:
            hdrlet = hdu[extname].headerlet
            valid_ctype = '-SIP' in hdrlet[1].header['ctype1']
            valid_dist_kws = 'A_ORDER' in hdrlet[1].header
            if not valid_ctype or not valid_dist_kws:
                key = wcsutil.altwcs.getKeyFromName(hdu['SCI', 1].header, wname)
                # Guard against the case where the headerlet WCS is not an alternate WCS.
                keyname = key.rstrip() if key is not None else key
                invalid_extns.append({'extname': extname, 'wcsname': wname, 'key': keyname})

    hdu.close()
    # If any invalid headerlets are found, remove them from
    if invalid_extns:
        for extn in reversed(invalid_extns):
            wcsutil.headerlet.delete_headerlet(filename, hdrext=extn['extname'])
            # also remove this solution from SCI headers
            if extn['key']:
                wcsutil.altwcs.deleteWCS(filename, extname_list, wcskey=extn['key'])
            log.debug(f"Delete duplicate headerlet extension {extn} in filename {filename}.")


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
            log.warning(f"""Failed to remove path {path} with 
                                 shutil.rmtree at attempt {n}: {err}""")
        time.sleep(3)

    if not ok:
        log.error(f"""Failed to remove path {path} with shutil.rmtree, 
                             even after {n} attempts.""")
        raise OSError
    else:
        log.debug(f"Path {path} successfully removed.")


def handle_remove_readonly(func, path, exc):
    excvalue = exc[1]
    if func in (os.rmdir, os.remove) and excvalue.errno == errno.EACCES:
        os.chmod(path, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)  # 0777
        func(path)
    else:
        raise


def _analyze_exposure(filename):
    """Evaluate whether or not this exposure should be processed at all."""

    process_exposure = True

    fhdu = fits.getheader(filename)
    targname = fhdu['targname']
    filts = fhdu['filt*']
    instrument = fhdu['instrume']

    filters = [filts[filtname].strip() for filtname in FILTER_NAMES[instrument]]

    # Let's start by checking whether the header indicates any problems with
    # the guide stars or tracking.
    # Using .get() insures that this check gets done even if keyword is missing.
    gs_quality = fhdu.get('quality', default="")
    if 'gsfail' in gs_quality.lower() or 'tdf-down' in gs_quality.lower():
        log.error(f"ERROR: Image {filename}'s QUALITY keywords: '{gs_quality}'")
        log.error("        GUIDING == BAD.  Skipping processing ")
        process_exposure = False  # Yes, there was bad guiding...

    badtab, _ = analyze.analyze_data([filename])
    if badtab['doProcess'][0] == 0:
        process_exposure = False

    # Also check to see whether this observation was taken with a blank filter name.
    if all(filter == '' for filter in filters):
        log.warning(f"ERROR: Inappropriate filter for exposure of {filters}")
        process_exposure = False

    return process_exposure


def _update_idcscale(filename):
    """
    Update IDCSCALE if a headerlet does not include it.

    Parameters
    ----------
    filename : str or fits.HDUList
        The science file to be updated.
        If HDUList it must be opened in "update" mode.

    """

    if isinstance(filename, str):
        hdul = fits.open(filename, mode='update')
    else:
        # Assume it's an HDUList object in "update" mode.
        hdul = filename
    num_sci = fileutil.countExtn(hdul)
    # get IDCTAB name
    itabroot = hdul[0].header['idctab'].split('$')[1].split('_')[0]
    fhdu_idscale = None
    # pull a value from one of the other headerlet extensions
    extvers = headerlet.get_headerlet_kw_names(hdul, kw='extver')
    for sciext in range(num_sci):
        for extn in extvers:
            if itabroot in hdul[('HDRLET', extn)].header['wcsname']:
                _hdrlet = hdul[('HDRLET', extn)].headerlet
                if 'idcscale' in _hdrlet[('SIPWCS', sciext + 1)].header:
                    fhdu_idscale = _hdrlet[('SIPWCS', sciext + 1)].header['idcscale']
                    break
        if fhdu_idscale is None:
            cd11 = hdul[('sci', sciext + 1)].header['CD1_1']
            cd21 = hdul[('sci', sciext + 1)].header['CD2_1']
            fhdu_idscale = np.sqrt(np.power(cd11, 2) + np.power(cd21, 2)) * 3600
        # Set the value of the IDCSCALE keyword
        for extn in range(num_sci):
            msg =  'Adding IDCSCALE {} to {}[sci,{}]'.format(fhdu_idscale, hdul.filename(), extn + 1)
            hdul[('sci', extn + 1)].header['idcscale'] = fhdu_idscale
            print(msg)
            update_msg = msg
    # No need to keep this file handle open anymore
    if isinstance(filename, str):
        hdul.close()
    del hdul


# Functions to support execution from the shell.
def main():

    import getopt

    try:
        optlist, args = getopt.getopt(sys.argv[1:], 'bdahfgimn:v:')
    except getopt.error as e:
        log.error(str(e))
        log.error(__doc__)
        log.error("\t", drizzlepac.__version__)

    # initialize default values
    help = 0
    force = False
    newdir = None
    inmemory = False
    num_cores = None
    headerlets = True
    align_to_gaia = True
    debug = True
    force_alignment = False
    do_verify_guiding = False
    make_manifest = False

    # read options
    for opt, value in optlist:
        if opt == "-g":
            align_to_gaia = False
        if opt == "-d":
            debug = True
        if opt == "-a":
            force_alignment = True
        if opt == "-h":
            help = 1
        if opt == "-f":
            force = True
        if opt == "-i":
            inmemory = True
        # The "-m" option is specific for WFPC2 data
        if opt == "-m":
            make_manifest = True
        if opt == "-v":
            do_verify_guiding = True
        if opt == '-n':
            if not value.isdigit():
                log.error('ERROR: num_cores value must be an integer!')
                raise ValueError
            num_cores = int(value)
        if opt == '-b':
            # turn off writing headerlets
            headerlets = False
    if len(args) < 1:
        log.error("syntax: runastrodriz.py [-bdahfginv] inputFilename [newpath]")
        sys.exit()
    if len(args) > 1:
        newdir = args[-1]
    if (help):
        log.debug(__doc__)
        log.debug("\t", drizzlepac.__version__)
    else:
        try:
            process(args[0], force=force, newpath=newdir, num_cores=num_cores,
                    inmemory=inmemory, headerlets=headerlets,
                    align_to_gaia=align_to_gaia, force_alignment=force_alignment,
                    do_verify_guiding=do_verify_guiding, debug=debug,
                    make_manifest=make_manifest)

        except Exception as errorobj:
            log.error(str(errorobj))
            log.error(f"ERROR: Cannot run astrodrizzle on {' '.join(sys.argv)}.")
            raise Exception(str(errorobj))

        # This except handles sys.exit() which raises the SystemExit exception which inherits from BaseException.
        except BaseException:
            exc_type, exc_value, exc_tb = sys.exc_info()
            log.debug(f"Return Value: {exc_value}")

if __name__ == "__main__":
    main()
