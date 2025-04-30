""" Utilities to analyze an input and determine whether the input is viable for a given process

The fundamental function analyze_data opens an input list containing FLT and/or FLC FITS
filenames in order to access the primary header data.  Based upon the values of specific
FITS keywords, the function determines whether or not each file within this dataset
can or should be reconciled against an astrometric catalog and, for multiple images, used
to create a mosaic.

The functions, mvm_analyze_wrapper and analyze_wrapper, are thin wrappers around analyze_data
to accommodate special case uses.  The mvm_analyze_wrapper takes a filename as input and
returns a boolean indicator.  The analyze_wrapper has the same function signature as
analyze_data, but it returns a list instead of an astropy table.
"""
import math
import sys

from enum import Enum
from astropy.io import fits
from astropy.io.fits import getheader, getdata
from astropy.table import Table
from astropy.stats import sigma_clipped_stats
import numpy as np

from photutils.segmentation import SourceCatalog, detect_sources

from scipy import ndimage
from skimage.transform import probabilistic_hough_line
from skimage.feature import canny

from stsci.tools import logutil
from stsci.tools.bitmask import bitfield_to_boolean_mask

from .astrometric_utils import classify_sources
from ..util import count_sci_extensions

__taskname__ = 'analyze'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)

__all__ = ['analyze_data', 'analyze_wrapper', 'mvm_analyze_wrapper']

# Define global default keyword names for these fields
"""
OBSKEY = 'OBSTYPE'
MTKEY = 'MTFLAG'
SCNKEY = 'SCAN_TYP'
FILKEY = 'FILTER'
FILKEY1 = 'FILTER1'
FILKEY2 = 'FILTER2'
APKEY = 'APERTURE'
TARKEY = 'TARGNAME'
EXPKEY = 'EXPTIME'
FGSKEY = 'FGSLOCK'
CHINKEY = 'CHINJECT'
DRIZKEY = 'DRIZCORR'
TYPEKEY = 'IMAGETYP'
QUALKEY = 'EXPFLAG'
"""

WFPC2_KEYS = {'OBSKEY': 'IMAGETYP', 'MTKEY': 'MTFLAG', 'SCNKEY': '',
              'FILKEY1': 'FILTNAM1', 'FILKEY2': 'FILTNAM2', 'FILKEY': 'FILTNAM1',
              'APKEY': '', 'TARKEY': 'TARGNAME', 'EXPKEY': 'EXPTIME',
              'FGSKEY': 'FGSLOCK', 'CHINKEY': '', 'DRIZKEY': 'DRIZCORR',
              'TYPEKEY': 'IMAGETYP', 'QUALKEY': 'EXPFLAG'}

DEFAULT_KEYS = {'OBSKEY': 'OBSTYPE', 'MTKEY':' MTFLAG', 'SCNKEY': 'SCAN_TYP',
                'FILKEY1': 'FILTER1', 'FILKEY2': 'FILTER2', 'FILKEY': 'FILTER',
                'APKEY': 'APERTURE', 'TARKEY': 'TARGNAME', 'EXPKEY': 'EXPTIME',
                'FGSKEY': 'FGSLOCK', 'CHINKEY': 'CHINJECT', 'DRIZKEY': 'DRIZCORR',
                'TYPEKEY': 'IMAGETYP', 'QUALKEY': 'EXPFLAG'}
HEADER_KEYS = {'WFPC2': WFPC2_KEYS, 'DEFAULT':DEFAULT_KEYS}

CAL_TARGETS = {'WFPC2': ['INTFLAT', 'UVFLAT', 'VISFLAT', 'KSPOTS',
                         'DARK', 'BIAS', 'EARTH-CALIB'],
               'DEFAULT': ['DARK', 'TUNG', 'BIAS', 'FLAT', 'DEUT', 'EARTH-CAL']
               }

# These definitions are for ACS and WFC3
BAD_DQ_FLAGS = [256,  # full-well saturated pixel
                512,  # bad pixel from reference file
                1024,  # weak charge trap
                2048,  # A-to-D saturated pixel
                4096  # cosmic-ray
]

MIN_LINES = 4  # Minimum number of detected lines for consideration of bad guiding


# Return codes
class Ret_code(Enum):
    """
    Define return status codes for Operations 
    """
    OK = 0
    KEYWORD_UPDATE_PROBLEM = 15
    SBCHRC_DATA = 55 
    NO_VIABLE_DATA = 65

# Annotates level to which image can be aligned according observational parameters
# as described through FITS keywords
class Messages(Enum):
    """
    Define a local classification for OK, Warning, and NoProcess messages
    """

    WARN, NOPROC = -1, -2


def mvm_analyze_wrapper(input_filename, log_level=logutil.logging.DEBUG):
    """
    Thin wrapper for the analyze_data function to return a viability indicator regarding a image for MVM processing.

    Parameters
    ==========
    input_filename : string
        Full fileName of data to be analyzed for viability to be processed as an MVM constituent.

    Returns
    =======
    use_for_mvm : boolean
        Boolean which indicates whether the input image should be used for MVM processing -
        True: use for MVM, False: do NOT use for MVM

    Note: This routine is invoked externally by software used for operations.

    """
    # Set logging level to user-specified level
    log.setLevel(log_level)

    # Invoke the low-level analyze_data routine with type = "MVM"
    filtered_table, _ = analyze_data([input_filename], type = "MVM")

    # There is only one row in this output table
    use_for_mvm = False
    if filtered_table['doProcess'] == 0:
        use_for_mvm = False
        log.warning("Image, {}, cannot be used for MVM processing.  Issue: {}.\n". \
                    format(input_filename, filtered_table['processMsg'][0]))
    else:
        use_for_mvm = True
        log.info("Image, {}, will be used for MVM processing.".format(input_filename))

    return use_for_mvm


def analyze_wrapper(input_file_list, log_level=logutil.logging.DEBUG, use_sbchrc=True, type=""):
    """
    Thin wrapper for the analyze_data function to return a list of viable images.

    Parameters
    ==========
    input_file_list : list
        List containing FLT and/or FLC filenames for all "ipppssoot" input images which comprise
        an associated dataset where 'associated dataset' may be a single image, multiple images,
        an HST association, or a number of HST associations when SVM processing. SVM FLT and/or FLC
        images may also be input for MVM processing.

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the .log file.
        Default value is 20, or 'info'.

    use_sbchrc : bool, optional
        Boolean based upon the environment variable setting of MVM_INCLUDE_SMALL which indicates whether
        or not to generate MVM SkyCell layers from ACS/HRC and ACS/SBC data.  These exposures typically
        only cover a miniscule fraction of a SkyCell and the plate scale of the SkyCell would result in
        a degraded representation of the original ACS/HRC and ACS/SBC data.  A setting of 'False' turns off
        the generation of layers from ACS/HRC and ACS/SBC data, and a setting of 'True' turns the generation on.
        Default = True

    type : string, optional
        String indicating whether this file is for MVM or some other processing.
        If type == "MVM", then Grism/Prism data is ignored.  If type == "" (default) or any
        other string, the Grism/Prism data is considered available for processing unless there is
        some other issue (i.e., exposure time of zero).
        Default = ""

    Returns
    =======
    viable_images_list : list
       List of images which can be used in the drizzle process.
       
    good_index : list
        indices of the viable images in the input_file_list

    return_code : int
        Numeric code indicative of the status of the analysis of the input data.
        These return codes are defined in this module by class Ret_code.

    Note: This routine returns a *list*, as well as a return code, containing only viable images
    instead of a table which provides information, as well as a doProcess bool, regarding each image.
    """
    # Set logging level to user-specified level
    log.setLevel(log_level)
 
    # Analyze the input file list and get the full table assessment
    filtered_table, analyze_data_good_index = analyze_data(input_file_list, type=type)

    # Reduce table to only the data which should be processed (doProcess == 1)
    mask = filtered_table["doProcess"] > 0
    filtered_table = filtered_table[mask]

    # Further reduce table to only the data which is NOT affected by bad guiding
    # This check is only to be done for MVM processing
    if type.upper() == "MVM":
        guide_mask = [not verify_guiding(f) for f in filtered_table["imageName"]]
        filtered_table = filtered_table[guide_mask]

    good_table = None
    good_rows = []
    good_index = []
    process_list = []
    return_value = Ret_code.OK.value

    # MVM processing, but excluding SBC/HRC data
    if use_sbchrc == False and type.upper() == "MVM": 
        # Check the table to determine presence of SBC/HRC data
        if filtered_table:
            for i, old_row in enumerate(filtered_table):
                if old_row["detector"].upper() != "SBC" and old_row["detector"].upper() != "HRC":
                    good_index.append(i)
                    good_rows.append(old_row)

            # The entire filtered_table contains only SBC or HRC data
            if not good_rows:
                log.warning("Only non-viable or SBC/HRC images in the multi-visit table - no processing done.\n")
                return_value = Ret_code.SBCHRC_DATA.value
            # Table contains some non-SBC/non-HRC data for processing
            else:
                good_table = Table(rows=good_rows, names=filtered_table.colnames)

                # Get the list of all "good" files to use for the alignment
                process_list = good_table['imageName'].tolist()
                return_value = Ret_code.OK.value

        # There is already nothing to process based upon the analysis criteria
        else:
            log.warning("No viable images in multi-visit table - no processing done.\n")
            return_value = Ret_code.NO_VIABLE_DATA.value
 
    # SVM processing or MVM processing with SBC/HRC data included in the MVM processing
    else:
        if filtered_table:
            # Get the list of all "good" files to use for the alignment
            process_list = filtered_table['imageName'].tolist()
            good_index = analyze_data_good_index
        else:
            log.warning("No viable images in single/multi-visit table - no processing done.\n")
            return_value = Ret_code.NO_VIABLE_DATA.value
    return process_list, good_index, return_value


def analyze_data(input_file_list, log_level=logutil.logging.DEBUG, type=""):
    """
    Determine if images within the dataset can be aligned

    Parameters
    ==========
    input_file_list : list
        List containing FLT and/or FLC filenames for all "ipppssoot" input images which comprise
        an associated dataset where 'associated dataset' may be a single image, multiple images,
        an HST association, or a number of HST associations when SVM processing. SVM FLT and/or FLC
        images may also be input for MVM processing.

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the .log file.
        Default value is 20, or 'info'.

    type : string
        String indicating whether this file is for MVM or some other processing.
        If type == "MVM", then Grism/Prism data is ignored.  If type == "" (default) or any
        other string, the Grism/Prism data is considered available for processing unless there is
        some other issue (i.e., exposure time of zero).

    Returns
    =======
    output_table : object
        Astropy Table object containing data pertaining to the associated dataset, including
        the do_process bool.  It is intended this table is updated by subsequent functions for
        bookkeeping purposes.
    
    analyze_data_good_index : list
        indices of the viable images in the input_file_list

    Notes
    =====
    The keyword/value pairs below define the "cannot process categories".
    OBSTYPE : is not IMAGING
    MTFLAG : T
    SCAN-TYP : C or D (or !N) (WFC3 only)
    FILTER : G*, PR*,  where G=Grism and PR=Prism
    FILTER1 : G*, PR*, where G=Grism and PR=Prism
    FILTER2 : G*, PR*, where G=Grism and PR=Prism
    TARGNAME : DARK, TUNGSTEN, BIAS, FLAT, EARTH-CALIB, DEUTERIUM (ACS and WFC3)
    TARGNAME : INTFLAT, UVFLAT, VISFLAT, KSPOTS, DARK, BIAS, EARTH-CALIB (WFPC2)
    EXPTIME : 0
    CHINJECT : is not NONE
    DRIZCORR : OMIT, SKIPPED
    EXPFLAG : is not NORMAL


    The keyword/value pairs below define the category which the data can be processed, but
    the results may be compromised
    FGSLOCK : FINE/GYRO, FINE/GY, COARSE, GYROS

    The keyword/value pairs below define the category which the data will be partially processed
    so the Grism/Prism images acquire the same WCS as their corresponding direct images.
    SVM processing only - FILTER or FILTER1 or FILTER2: G*, PR*

    FITS Keywords only for WFC3 data: SCAN_TYP, FILTER, and CHINJECT (UVIS)
    FITS Keywords only for ACS data: FILTER1 and FILTER2
    FITS Keywords only for WFPC2 data: FILTNAM1 and FILTNAM2

    Please be aware of the FITS keyword value NONE vs the Python None.
    """
    # Set logging level to user-specified level
    log.setLevel(log_level)

    analyze_data_good_index = []

    acs_filt_name_list = [DEFAULT_KEYS['FILKEY1'], DEFAULT_KEYS['FILKEY2']]

    # Interpret input filenames and adjust size of column accordingly
    max_name_length = max([len(f) for f in input_file_list])
    name_data_type = 'S{}'.format(max_name_length + 2)  # add a couple of spaces to insure separation of cols

    # Initialize the column entries which will be populated in successive
    # processing steps
    fit_method = ""  # Fit algorithm used for alignment
    catalog = ""     # Astrometric catalog used for alignment
    catalog_sources = 0  # No. of astrometric catalog sources found based on coordinate overlap with image
    found_sources = 0   # No. of sources detected in images
    match_sources = 0   # No. of sources cross matched between astrometric catalog and detected in image
    offset_x = None
    offset_y = None
    rot = None
    scale = None
    rms_x = -1.0
    rms_y = -1.0
    rms_ra = -1.0
    rms_dec = -1.0
    completed = False  # If true, there was no exception and the processing completed all logic
    date_obs = None     # Human readable date
    mjdutc = -1.0      # MJD UTC start of exposure
    fgslock = None
    process_msg = ""
    status = 9999
    compromised = 0
    headerlet_file = ""
    fit_qual = -1

    fit_rms = -1.0
    total_rms = -1.0
    dataset_key = -1.0

    names_array = ('imageName', 'instrument', 'detector', 'filter', 'aperture',
                   'obstype', 'subarray', 'dateObs', 'mjdutc', 'doProcess',
                   'processMsg', 'fit_method', 'catalog', 'foundSources',
                   'catalogSources', 'matchSources', 'offset_x', 'offset_y',
                   'rotation', 'scale', 'rms_x', 'rms_y', 'rms_ra', 'rms_dec',
                   'completed', 'fit_rms', 'total_rms', 'datasetKey', 'status',
                   'fit_qual', 'headerletFile', 'compromised')
    data_type = (name_data_type, 'S20', 'S20', 'S20', 'S20', 'S20', 'b', 'S20', 'f8', 'b',
                 'S30', 'S20', 'S20', 'i4', 'i4', 'i4', 'f8', 'f8', 'f8', 'f8',
                 'f8', 'f8', 'f8', 'f8', 'b', 'f8', 'f8', 'i8', 'i4', 'i4',
                 'S60', 'i4')

    # Create an astropy table
    output_table = Table(names=names_array, dtype=data_type)

    # Loop over the list of images to determine viability for alignment processing
    #
    # Capture the data characteristics before any evaluation so the information is
    # available for the output table regardless of which keyword is used to
    # to determine the data is not viable for alignment.

    for i, input_file in enumerate(input_file_list):
        header_hdu = 0
        header_data = getheader(input_file, header_hdu)

        # Keywords to use potentially for downstream analysis
        instrume = (header_data['INSTRUME']).upper()
        if instrume == 'WFPC2':
            detector = 'PC'
            subarray = False
            hdr_keys = HEADER_KEYS[instrume]
            aperture = 'PC'
            scan_typ = 'N'  # default value of N/A
            mtflag = 'T' if header_data[hdr_keys['MTKEY']] else 'F'
            obstype = 'IMAGING' if (header_data[hdr_keys['OBSKEY']]).upper() == 'EXT' else 'CAL'

        else:
            hdr_keys = HEADER_KEYS['DEFAULT']
            detector = (header_data['DETECTOR']).upper()
            subarray = header_data['SUBARRAY']
            aperture = (header_data[hdr_keys['APKEY']]).upper()
            mtflag = (header_data[hdr_keys['MTKEY']]).upper()
            obstype = (header_data[hdr_keys['OBSKEY']]).upper()

        date_obs = header_data['DATE-OBS']
        mjdutc = header_data['EXPSTART']

        # Obtain keyword values for analysis of viability
        drizcorr = (header_data[hdr_keys['DRIZKEY']]).upper()
        scan_typ = ''
        if instrume == 'WFC3':
            scan_typ = (header_data[hdr_keys['SCNKEY']]).upper()

        sfilter = ''
        if instrume == 'WFC3':
            sfilter = (header_data[hdr_keys['FILKEY']]).upper()
        # Concatenate the two ACS filter names together with an underscore
        # If the filter name is blank, skip it
        if instrume == 'ACS':
            for filtname in acs_filt_name_list:

                # The filter keyword value could be zero or more blank spaces
                # Strip off any leading or trailing blanks
                if header_data[filtname].upper().strip():

                    # If the current filter variable already has some content,
                    # need to append an underscore before adding more text
                    if sfilter:
                        sfilter += '_'
                    sfilter += header_data[filtname].upper().strip()

        # The aperture is only read for informational purposes as it is no
        # longer used for filtering input data.
        targname = (header_data[hdr_keys['TARKEY']]).upper()
        exptime = header_data[hdr_keys['EXPKEY']]
        fgslock = (header_data[hdr_keys['FGSKEY']]).upper()
        imagetype = (header_data[hdr_keys['TYPEKEY']]).upper()
        expflag = (header_data[hdr_keys['QUALKEY']]).upper()

        chinject = 'NONE'
        if instrume == 'WFC3' and detector == 'UVIS':
            chinject = (header_data[hdr_keys['CHINKEY']]).upper()

        # Determine if the image has one of these conditions.  The routine
        # will exit processing upon the first satisfied condition.

        # Check if all science image arrays in the RAW file are filled with zero values
        non_zero_data_in_array = False # start assuming data is zeros
        science_ext_ind_array = count_sci_extensions(input_file, return_ind=True)
        # make sure science extension exists
        if len(science_ext_ind_array)>0:
            for sci_ext_ind in science_ext_ind_array:
                science_data = getdata(input_file, sci_ext_ind)
                # change flag if good data in any science extension array
                if not np.all(science_data==0):
                    non_zero_data_in_array = True
                else:
                    log.warning(
                        f"{input_file} (SCI, {sci_ext_ind}) is all zeros, but processing will continue with the other science extensions."
                    )
    
        else:
            log.warning(f'No science extension in file: {input_file}')

        # Compute if the exposure time is very close to zero as it will be
        # needed when deciding whether or not to use the particular Grism/Prism data
        exposure_time_near_zero = True if math.isclose(exptime, 0.0, abs_tol=1e-5) else False

        no_proc_key = None
        no_proc_value = None
        do_process = True

        # Imaging vs spectroscopic or coronagraphic
        if obstype != 'IMAGING':
            no_proc_key = hdr_keys['OBSKEY']
            no_proc_value = obstype

        # drizzling has been turned off
        elif drizcorr in ['OMIT', 'SKIPPED']:
            no_proc_key = hdr_keys['DRIZKEY']
            no_proc_value = drizcorr

        # Moving target
        elif mtflag == 'T':
            no_proc_key = hdr_keys['MTKEY']
            no_proc_value = mtflag

        # Bostrophidon without or with dwell (WFC3 only)
        elif any([scan_typ == 'C', scan_typ == 'D']):
            no_proc_key = hdr_keys['SCNKEY']
            no_proc_value = scan_typ

        # Calibration target
        elif any(x in targname for x in CAL_TARGETS['DEFAULT'])\
                and instrume != 'WFPC2':
            no_proc_key = hdr_keys['TARKEY']
            no_proc_value = targname

        # WFPC2 sets the imagetyp keyword correctly(?) as something other than EXT for cal observations
        elif any(x in targname for x in CAL_TARGETS['WFPC2'])\
                and instrume == 'WFPC2' and imagetype != 'EXT':
            no_proc_key = hdr_keys['TARKEY']
            no_proc_value = targname

        # Exposure time of effectively zero
        elif math.isclose(exptime, 0.0, abs_tol=1e-5):
            no_proc_key = hdr_keys['EXPKEY']
            no_proc_value = exptime

        # Commanded FGS lock
        elif any(x in fgslock for x in ['GY', 'COARSE']):
            no_proc_key = hdr_keys['FGSKEY']
            no_proc_value = fgslock

        # Charge injection mode
        elif chinject != 'NONE':
            no_proc_key = hdr_keys['CHINKEY']
            no_proc_value = chinject

        # Exposure flag indicates a "issue" happened during an exposure
        # elif (expflag != 'NORMAL'):
        #     no_proc_key = hdr_keys['QUALKEY']
        #     no_proc_value = expflag

        # Ramp filter images should not be processed for MVM products.
        #
        # Filter name which starts with "BLOCK" for internal calibration of SBC
        # The sfilter variable may be the concatenation of two filters (F160_CLEAR)
        #
        # Grism/Prism images are also IMAGING=SPECTROSCOPIC, suppress the "no processing"
        # indicators to allow the Grism/Prism images to be minimally processed for
        # keyword updates.  This was done as a retrofit to allow Grism/Prism images
        # to have the same WCS solution as the direct images in the visit (same detector).
        # The exception to this will be if the Grism/Prism data has a zero exposure time as
        # the WCS will be only "OPUS" under this condition, and this will disrupt processing
        # for all the good data.
        split_sfilter = sfilter.upper().split('_')
        for item in split_sfilter:
            # This is the only circumstance when Grism/Prism data WILL be processed.
            if item.startswith(('G', 'PR')) and not exposure_time_near_zero and type.upper() == "SVM":
                no_proc_key = None
                no_proc_value = None
                log.info("The Grism/Prism data, {}, will be processed.".format(input_file))
            # Grism/Prism WILL NOT be processed primarily if MVM processing or with an exposure time of zero.
            elif item.startswith(('G', 'PR')): 
                if type.upper() == "MVM":
                    no_proc_value += ", Grism/Prism data and MVM processing"
                    log.warning("The Grism/Prism data {} with MVM processing will be ignored.".format(input_file))
                elif exposure_time_near_zero:
                    no_proc_value += ", Grism/Prism data and EXPTIME = 0.0"
                    log.warning("The Grism/Prism data {} with zero exposure time will be ignored.".format(input_file))

            if item.startswith(('BLOCK')):
                no_proc_key = hdr_keys['FILKEY']
                no_proc_value = sfilter

            if item.startswith(('FR')) and type.upper() == "MVM":
                no_proc_key = hdr_keys['FILKEY']
                no_proc_value = "Ramp data and MVM processing"
                log.warning("The Ramp data {} with MVM processing will be ignored.".format(input_file))

        # If no_proc_key is set to a keyword, then this image has been found to not be viable for
        # alignment purposes.
        if no_proc_key is not None:
            if no_proc_key != hdr_keys['FGSKEY']:
                do_process = False
                msg_type = Messages.NOPROC.value
            else:
                msg_type = Messages.WARN.value
                analyze_data_good_index.append(i)

            process_msg = no_proc_key + '=' + str(no_proc_value)

            # Issue message to log file for this data indicating no processing to be done or
            # processing should be allowed, but there may be some issue with the result (e.g.,
            # GYROS mode so some drift)
            generate_msg(input_file, msg_type, no_proc_key, no_proc_value)
        elif non_zero_data_in_array==False:
            do_process=False
            process_msg="SCI data all zeros"
            log.warning(f'Science data for {input_file} filled with zeros. Dataset cannot be aligned.')

        else:
            analyze_data_good_index.append(i)

        # Populate a row of the table
        output_table.add_row([input_file, instrume, detector, sfilter, aperture, obstype, subarray,
                              date_obs, mjdutc, do_process, process_msg, fit_method, catalog,
                              found_sources, catalog_sources, match_sources, offset_x, offset_y,
                              rot, scale, rms_x, rms_y, rms_ra, rms_dec, completed, fit_rms,
                              total_rms, dataset_key, status, fit_qual, headerlet_file,
                              compromised])
        process_msg = ""
    return output_table, analyze_data_good_index


def generate_msg(filename, msg, key, value):
    """ Generate a message for the output log indicating the file/association will not
        be processed as the characteristics of the data are known to be inconsistent
        with alignment.
    """

    log.warning('Dataset ' + filename + ' has (keyword = value) of (' + key + ' = ' + str(value) + ').')
    if msg == Messages.NOPROC.value:
        log.warning('Dataset cannot be aligned.')
    else:
        log.warning('Dataset can be aligned, but the result may be compromised.')


# -----------------------------------------------------------------------------
#  Line detection functions
# -----------------------------------------------------------------------------

def verify_guiding(filename, min_length=33):
    """ Verify whether or not the input image was affected by guiding problems.

    This algorithm evaluates the data from (SCI,1) to try to determine whether
    the image was affected by guiding problems which mimic SCAN mode or GRISM data
    with the trails in an arbitrary angle across the image.

    Parameters
    ==========
    filename : str
        Name of image to be evaluated

    min_length : int, optional
        Minimum length of trails (in pixels) to be detected in image.

    Returns
    ========
    bad_guiding : bool
        Boolean specifying whether or not the image was detected as
        being affected by guiding problems.  Value is True if image
        was affected.

    Note: This function will only be called from analyze_wrapper if the processing
          type is "MVM".  It is deliberately False for other processing (e.g., SVM and
          pipeline).  However, this routine can be invoked directly for pipeline
          processing from runastrodriz.py with the appropriate parameter setting.
    """
    log.info(f"Verifying that {filename} was not affected by guiding problems.")

    hdu = fits.open(filename)
    # Let's start by checking whether the header indicates any problems with
    # the guide stars or tracking.
    # Using .get() insures that this check gets done even if keyword is missing.
    gs_quality = hdu[0].header.get('quality', default="").lower()
    if 'gsfail' in gs_quality or 'tdf-down' in gs_quality:
        hdu.close()
        log.warning(f"Image {filename}'s QUALITY keywords report GUIDING: BAD.")
        return True  # Yes, there was bad guiding...

    # No guide star problems indicated in header, so let's check the
    # data.  There are many instances where the data is compromised
    # despite what values are found in the header.
    data = hdu[("SCI", 1)].data.copy()
    scale_data = hdu[("SCI",1)].header["bunit"].endswith('/S')
    data = np.nan_to_num(data, nan=0.0)  # just to be careful
    if scale_data:
        # Photutils works best in units of DN
        scale_hdr = hdu[0].header if 'exptime' in hdu[0].header else hdu[1].header
        scale_val = scale_hdr['exptime']
        data *= scale_val
    bkg_stats = sigma_clipped_stats(data, maxiters=2)
    bkg_limit = bkg_stats[1] + bkg_stats[2]  # only need a 1-sigma detection limit here...
    log.debug(f"bkg_limit found to be: {bkg_limit:.2f}")

    data -= bkg_limit
    imgarr = np.clip(data, 0, data.max())

    # Build up a mask of all bad pixels and ignore them when looking for
    # sources and linear features
    dqarr = None
    for extn in hdu:
        if 'extname' in extn.header and extn.header['extname'] == 'DQ':
            dqarr = hdu[extn].data.copy()
            break
    if dqarr is not None:
        dqmask = bitfield_to_boolean_mask(dqarr, ignore_flags=BAD_DQ_FLAGS)
    else:
        dqmask = np.ones_like(data)
    # close FITS object (just to be nice to the OS...)
    hdu.close()
    del hdu

    # apply mask now...
    imgarr *= dqmask
    del dqmask  # just to clean up a little

    # Determine rough number of probable sources
    # Trying to ignore small sources (<= 4x4 pixels in size, or npixels < 17)
    #   which are either noise peaks or head-on CRs.
    segm = detect_sources(imgarr, 0, npixels=17)
    if segm is None or segm.nlabels < 2:
        log.debug(f'Did NOT detect enough raw sources in {filename} for guiding verification.')
        return False
    log.debug(f'Detected {segm.nlabels} raw sources in {filename} for guiding verification.')

    src_cat = SourceCatalog(imgarr, segm)
    # Remove likely cosmic-rays based on central_moments classification
    bad_srcs = classify_sources(src_cat, 1.5)
    # Get the label IDs for sources flagged as CRs, IDs are 1-based not 0-based
    pt_srcs = np.where(bad_srcs == 0)[0] + 1
    segm.remove_labels(pt_srcs)
    # If there are no sources left, this is akin to the check above where number
    # of sources < 2, so just return False
    if segm.nlabels == 0:
        log.warning("After removal of suspected cosmic rays, there are no sources remaining in the image.")
        return False
    src_cat = SourceCatalog(imgarr, segm)  # clean up source catalog now...
    num_sources = len(src_cat)

    # trim edges from image to avoid spurious sources
    trim_slice=(slice(2, -2), slice(2, -2))

    # Now determine whether this image was affected by guiding problems
    bad_guiding = lines_in_image(imgarr[trim_slice], num_sources,
                                 min_length=min_length, min_lines=MIN_LINES)
    if bad_guiding:
        log.warning(f"Image {filename}'s GUIDING detected as: BAD.")
    else:
        log.info(f"Image {filename}'s GUIDING detected as: GOOD.")

    return bad_guiding


def detect_lines(image, mask=None, min_length=17):
    """Detect lines in the input image and return list of line parameters """
    lines = {'num': None,
             'startarr': None,
             'endarr': None,
             'angles': None,
             'lengths': None,
             'slopes': None}

    # extract edges from image for faster line detection
    edges = canny(image, sigma=2.5, low_threshold=0, high_threshold=25, mask=mask)

    # Classic straight-line Hough transform
    plines = probabilistic_hough_line(edges, threshold=0,
                                     line_gap=0,
                                     line_length=min_length)
    if len(plines) > 0:
        plines = np.array(plines)
        startarr = plines[:, 0]
        endarr = plines[:, 1]
        rise = endarr[:, 1] - startarr[:, 1]
        run = endarr[:, 0]-startarr[:, 0]
        angles = np.arctan2(rise, run)
        lines['startarr'] = startarr
        lines['endarr'] = endarr
        lines['angles'] = np.rad2deg(angles)
        lines['lengths'] = np.sqrt(rise**2 + run**2)
        lines['slopes'] = np.tan(angles + np.pi/2)
        lines['num'] = len(plines)

    return lines


def lines_in_image(image, num_sources, mask=None, min_length=17, min_lines=4):
    """Determine if image is dominated by linear features

    Parameters
    ----------
    image : ndarray
        Background-subtracted image to detect linear features in

    sensitivity : float, optional
        Increments in degrees for detecting lines in image.

    Returns
    -------
    lines_detected : bool
        Specifies whether or not image is dominated by linear features
    """
    # detect any lines in image
    lines = detect_lines(image, mask=mask, min_length=min_length)
    if lines['num'] is None:
        log.debug(f"No linear features found.")
        return False
    else:
        # If we have a lot of sources, we should have a lot of lines if bad
        # Otherwise, however, min_lines is used to guard against faint fields
        if lines['num'] < max(min_lines, num_sources/10):
            log.debug(f"Only {lines['num']} linear features found.")
            return False

    # perform statistical analysis on detected linear features
    # start by generating a histogram of the angles of all the lines
    # We will ignore lines that are exactly in line with a column (+/- 90)
    # as they are nearly always caused by CTE or saturation bleeding along the columns.
    diff_lines = np.isclose(np.abs(lines['angles']), 90, atol=2.0)
    angles = lines['angles'][~diff_lines]

    angle_bins = np.linspace(-180., 180., 91)
    ahist = np.histogram(angles, bins=angle_bins)

    # if one peak has a more than 10% of all linear features detected,
    # and there are more linear features than lines from saturated columns
    # it is considered as having bad guiding.
    max_angles = ahist[0].max()
    alimit = max(len(angles) / 10.0, diff_lines.sum())

    log.debug(f"Peak number of similar lines: {max_angles} based on a threshold of {alimit}")
    log.debug(f"number of probable sources: {num_sources}")

    # Now check to see if enough detected lines have the same (non-90 deg) orientation
    lines_detected = (max_angles > alimit)
    log.info(f"{max_angles} lines were similar, so linear features were detected.")
    return lines_detected
