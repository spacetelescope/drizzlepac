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
from astropy.io.fits import getheader
from astropy.table import Table
import numpy as np

from stsci.tools import logutil

__taskname__ = 'analyze'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.DEBUG, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)

__all__ = ['analyze_data', 'analyze_wrapper', 'mvm_analyze_wrapper']

# Define global default keyword names for these fields
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

# Annotates level to which image can be aligned according observational parameters
# as described through FITS keywords


class Messages(Enum):
    """
    Define a local classification for OK, Warning, and NoProcess messages
    """

    OK, WARN, NOPROC = 1, -1, -2


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
    """
    # Set logging level to user-specified level
    log.setLevel(log_level)

    # Invoke the low-level analyze_data routine with type = "MVM"
    filtered_table = analyze_data([input_filename], type = "MVM")

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


def analyze_wrapper(input_file_list, log_level=logutil.logging.DEBUG):
    """
    Thin wrapper for the analyze_data function to return a list of viable images.

    Parameters
    ==========
    input_file_list : list
        List containing FLT and/or FLC filenames for all input images which comprise an associated
        dataset where 'associated dataset' may be a single image, multiple images, an HST
        association, or a number of HST associations

    Returns
    =======
    viable_images_list : list
       List of images which can be used in the drizzle process.

    This routine returns a list containing only viable images instead of a table which
    provides information, as well as a doProcess bool, regarding each image.
    """
    # Set logging level to user-specified level
    log.setLevel(log_level)

    process_list = []

    # Analyze the input file list and get the full table assessment
    filtered_table = analyze_data(input_file_list)

    # Extract only the filenames of viable images for processing (i.e., doProcess == 1)
    if filtered_table['doProcess'].sum() == 0:
        log.warning("No viable images in single/multiple visit table - no processing done.\n")
    else:
        # Get the list of all "good" files to use for the alignment
        process_list = filtered_table['imageName'][np.where(filtered_table['doProcess'])]
        process_list = list(process_list)  # Convert process_list from numpy list to regular python list

    return process_list


def analyze_data(input_file_list, log_level=logutil.logging.DEBUG, type=""):
    """
    Determine if images within the dataset can be aligned

    Parameters
    ==========
    input_file_list : list
        List containing FLT and/or FLC filenames for all input images which comprise an associated
        dataset where 'associated dataset' may be a single image, multiple images, an HST
        association, or a number of HST associations

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

    Notes
    =====
    The keyword/value pairs below define the "cannot process categories".
    OBSTYPE : is not IMAGING
    MTFLAG : T
    SCAN-TYP : C or D (or !N)
    FILTER : G*, PR*,  where G=Grism and PR=Prism
    FILTER1 : G*, PR*, where G=Grism and PR=Prism
    FILTER2 : G*, PR*, where G=Grism and PR=Prism
    TARGNAME : DARK, TUNGSTEN, BIAS, FLAT, EARTH-CALIB, DEUTERIUM
    EXPTIME : 0
    CHINJECT : is not NONE
    DRIZCORR : OMIT, SKIPPED

    The keyword/value pairs below define the category which the data can be processed, but
    the results may be compromised
    FGSLOCK : FINE/GYRO, FINE/GY, COARSE, GYROS

    FITS Keywords only for WFC3 data: SCAN_TYP, FILTER, and CHINJECT (UVIS)
    FITS Keywords only for ACS data: FILTER1 and FILTER2

    Please be aware of the FITS keyword value NONE vs the Python None.
    """
    # Set logging level to user-specified level
    log.setLevel(log_level)

    acs_filt_name_list = [FILKEY1, FILKEY2]

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

    for input_file in input_file_list:
        header_hdu = 0
        header_data = getheader(input_file, header_hdu)

        # Keywords to use potentially for downstream analysis
        instrume = (header_data['INSTRUME']).upper()
        detector = (header_data['DETECTOR']).upper()
        subarray = header_data['SUBARRAY']
        date_obs = header_data['DATE-OBS']
        mjdutc = header_data['EXPSTART']

        # Obtain keyword values for analysis of viability
        drizcorr = (header_data[DRIZKEY]).upper()
        obstype = (header_data[OBSKEY]).upper()
        mtflag = (header_data[MTKEY]).upper()
        scan_typ = ''
        if instrume == 'WFC3':
            scan_typ = (header_data[SCNKEY]).upper()

        sfilter = ''
        if instrume == 'WFC3':
            sfilter = (header_data[FILKEY]).upper()
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
        aperture = (header_data[APKEY]).upper()
        targname = (header_data[TARKEY]).upper()
        exptime = header_data[EXPKEY]
        fgslock = (header_data[FGSKEY]).upper()

        chinject = 'NONE'
        if instrume == 'WFC3' and detector == 'UVIS':
            chinject = (header_data[CHINKEY]).upper()

        # Determine if the image has one of these conditions.  The routine
        # will exit processing upon the first satisfied condition.

        # Compute if the exposure time is very close to zero as it will be
        # needed when deciding whether or not to use the particular Grism/Prism data
        is_zero = True if math.isclose(exptime, 0.0, abs_tol=1e-5) else False

        no_proc_key = None
        no_proc_value = None
        do_process = True
        # Imaging vs spectroscopic or coronagraphic
        if obstype != 'IMAGING':
            no_proc_key = OBSKEY
            no_proc_value = obstype

        # drizzling has been turned off
        elif drizcorr in ['OMIT', 'SKIPPED']:
            no_proc_key = DRIZKEY
            no_proc_value = drizcorr

        # Moving target
        elif mtflag == 'T':
            no_proc_key = MTKEY
            no_proc_value = mtflag

        # Bostrophidon without or with dwell (WFC3 only)
        elif any([scan_typ == 'C', scan_typ == 'D']):
            no_proc_key = SCNKEY
            no_proc_value = scan_typ

        # Calibration target
        elif any(x in targname for x in ['DARK', 'TUNG', 'BIAS', 'FLAT', 'DEUT', 'EARTH-CAL']):
            no_proc_key = TARKEY
            no_proc_value = targname

        # Exposure time of effectively zero
        elif math.isclose(exptime, 0.0, abs_tol=1e-5):
            no_proc_key = EXPKEY
            no_proc_value = exptime

        # Commanded FGS lock
        elif any(x in fgslock for x in ['GY', 'COARSE']):
            no_proc_key = FGSKEY
            no_proc_value = fgslock

        # Charge injection mode
        elif chinject != 'NONE':
            no_proc_key = CHINKEY
            no_proc_value = chinject

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
            if item.startswith(('G', 'PR')) and not is_zero and type.upper() != "MVM":
                no_proc_key = None
                no_proc_value = None
                log.info("The Grism/Prism data, {}, will be processed.".format(input_file))
            # Grism/Prism WILL NOT be processed primarily if MVM processing or with an exposure time of zero.
            elif item.startswith(('G', 'PR')): 
                if type.upper() == "MVM":
                    no_proc_value += ", Grism/Prism data and MVM processing"
                    log.warning("The Grism/Prism data {} with MVM processing will be ignored.".format(input_file))
                elif is_zero:
                    no_proc_value += ", Grism/Prism data and EXPTIME = 0.0"
                    log.warning("The Grism/Prism data {} with zero exposure time will be ignored.".format(input_file))
    
            if item.startswith(('BLOCK')):
                no_proc_key = FILKEY
                no_proc_value = sfilter

        # If no_proc_key is set to a keyword, then this image has been found to not be viable for
        # alignment purposes.
        if no_proc_key is not None:
            if no_proc_key != FGSKEY:
                do_process = False
                msg_type = Messages.NOPROC.value
            else:
                msg_type = Messages.WARN.value

            process_msg = no_proc_key + '=' + str(no_proc_value)

            # Issue message to log file for this data indicating no processing to be done or
            # processing should be allowed, but there may be some issue with the result (e.g.,
            # GYROS mode so some drift)
            generate_msg(input_file, msg_type, no_proc_key, no_proc_value)

        # Populate a row of the table
        output_table.add_row([input_file, instrume, detector, sfilter, aperture, obstype, subarray,
                              date_obs, mjdutc, do_process, process_msg, fit_method, catalog,
                              found_sources, catalog_sources, match_sources, offset_x, offset_y,
                              rot, scale, rms_x, rms_y, rms_ra, rms_dec, completed, fit_rms,
                              total_rms, dataset_key, status, fit_qual, headerlet_file,
                              compromised])
        process_msg = ""

    return output_table


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
