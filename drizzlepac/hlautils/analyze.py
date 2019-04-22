""" Utility to analyze an input dataset and determine whether the dataset can be aligned

The function analyze_data opens an input list containing FLT and/or FLC FITS filenames
in order to access the primary header data.  Based upon the values of specific
FITS keywords, the function determines whether or not each file within this dataset
can or should be reconciled against an astrometric catalog and, for multiple images, used
to create a mosaic.
"""
from astropy.io.fits import getheader
from astropy.table import Table
import math
import sys
from enum import Enum

from stsci.tools import logutil

__taskname__ = 'analyze'

log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout)

__all__ = ['analyze_data']

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

# Annotates level to which image can be aligned according observational parameters
# as described through FITS keywords
class Messages(Enum):
    OK, WARN, NOPROC = 1, -1, -2

def analyze_data(input_file_list, **kwargs):
    """
    Determine if images within the dataset can be aligned

    Parameters
    ==========
    input_file_list: list
        List containing FLT and/or FLC filenames for all input images which comprise an associated
        dataset where 'associated dataset' may be a single image, multiple images, an HST association,
        or a number of HST associations

    Returns
    =======
    output_table: object
        Astropy Table object containing data pertaining to the associated dataset, including
        the do_process bool.  It is intended this table is updated by subsequent functions for
        bookkeeping purposes.

    Notes
    =====
    The keyword/value pairs below define the "cannot process categories".
    OBSTYPE : is not IMAGING
    MTFLAG : T
    SCAN-TYP : C or D (or !N)
    FILTER : G*, *POL*, *PRISM*
    FILTER1 : G*, *POL*, *PRISM*
    FILTER2 : G*, *POL*, *PRISM*
    APERTURE : *GRISM*, G*-REF, RAMP, *POL*, *PRISM*
    TARGNAME : DARK, TUNGSTEN, BIAS, FLAT, EARTH-CALIB, DEUTERIUM
    EXPTIME : 0
    CHINJECT : is not NONE

    The keyword/value pairs below define the category which the data can be processed, but
    the results may be compromised
    FGSLOCK : FINE/GYRO, FINE/GY, COARSE, GYROS

    FITS Keywords only for WFC3 data: SCAN_TYP, FILTER, and CHINJECT (UVIS)
    FITS Keywords only for ACS data: FILTER1 and FILTER2

    Please be aware of the FITS keyword value NONE vs the Python None.
    FIX: improve robustness when analyzing filter and aperture names, possibly use PHOTMODE instead
    """

    acs_filt_name_list = [FILKEY1, FILKEY2]

    fit_method = None  # Fit algorithm used for alignment
    catalog = None     # Astrometric catalog used for alignment
    catalog_sources = 0  # Number of astrometric catalog sources determined based upon coordinate overlap with image WCS
    found_sources = 0   # Number of sources detected in images
    match_sources = 0   # Number of sources cross matched between astrometric catalog and detected in image
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
    process_msg = None
    status = 9999
    compromised = 0
    headerlet_file = None
    fit_qual = -1

    fit_rms = -1.0
    total_rms = -1.0
    dataset_key = -1.0

    names_array = ('imageName', 'instrument', 'detector', 'filter', 'aperture', 'obstype', 'subarray',
                  'dateObs', 'mjdutc', 'doProcess', 'processMsg', 'fit_method', 'catalog', 'foundSources',
                  'catalogSources', 'matchSources', 'offset_x', 'offset_y', 'rotation', 'scale', 'rms_x',
                  'rms_y', 'rms_ra', 'rms_dec', 'completed', 'fit_rms', 'total_rms', 'datasetKey', 'status',
                  'fit_qual', 'headerletFile', 'compromised')
    data_type = ('S20', 'S20', 'S20', 'S20', 'S20', 'S20', 'b', 'S20', 'f8', 'b', 'S30', 'S20', 'S20', 'i4',
                'i4', 'i4', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'b', 'f8', 'f8', 'i8', 'i4',
                'i4', 'S30', 'i4')

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
                if len(header_data[filtname].upper().strip()) > 0:

                    # If the current filter variable already has some content,
                    # need to append an underscore before adding more text
                    if len(sfilter) > 0:
                        sfilter += '_'
                    sfilter += header_data[filtname].upper().strip()

        aperture = (header_data[APKEY]).upper()
        targname = (header_data[TARKEY]).upper()
        exptime = header_data[EXPKEY]
        fgslock = (header_data[FGSKEY]).upper()

        chinject = 'NONE'
        if instrume == 'WFC3' and detector == 'UVIS':
            chinject = (header_data[CHINKEY]).upper()

        # Determine if the image has one of these conditions.  The routine
        # will exit processing upon the first satisfied condition.

        no_proc_key = None
        no_proc_value = None
        do_process = True
        # Imaging vs spectroscopic or coronagraphic
        if obstype != 'IMAGING':
            no_proc_key = OBSKEY
            no_proc_value = obstype

        # Moving target
        elif mtflag == 'T':
            no_proc_key = MTKEY
            no_proc_value = mtflag

        # Bostrophidon without or with dwell (WFC3 only)
        elif any([scan_typ == 'C', scan_typ == 'D']):
            no_proc_key = SCNKEY
            no_proc_value = scan_typ

        # Ramp, polarizer, grism, or prism
        elif any(x in aperture for x in ['RAMP', 'POL', 'GRISM', '-REF', 'PRISM']):
            no_proc_key = APKEY
            no_proc_value = aperture

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

        # Filter which does not begin with: 'F'(F###), 'C'(CLEAR), 'N'(N/A), and is not blank
        # The sfilter variable may be the concatenation of two filters (F160_CLEAR)
        elif sfilter and not sfilter.startswith(('F', 'C', 'N')):
            no_proc_key = FILKEY
            no_proc_value = sfilter

        elif '_' in sfilter:
            pos = sfilter.index('_')
            pos += 1

            if sfilter[pos] != 'F' and sfilter[pos] != '' and sfilter[pos] != 'C' and sfilter[pos] != 'N':
                no_proc_key = FILKEY
                no_proc_value = sfilter

        # If no_proc_key is set to a keyword, then this image has been found to not be viable for
        # alignment purposes.
        if (no_proc_key is not None):
            if (no_proc_key != FGSKEY):
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
        output_table.add_row([input_file, instrume, detector, sfilter, aperture, obstype, subarray, date_obs,
                             mjdutc, do_process, process_msg, fit_method, catalog, found_sources,
                             catalog_sources, match_sources, offset_x, offset_y, rot, scale, rms_x, rms_y,
                             rms_ra, rms_dec, completed, fit_rms, total_rms, dataset_key, status, fit_qual,
                             headerlet_file, compromised])
        process_msg = None
    # output_table.pprint(max_width=-1)

    return(output_table)


def generate_msg(filename, msg, key, value):
    """ Generate a message for the output log indicating the file/association will not
        be processed as the characteristics of the data are known to be inconsistent
        with alignment.
    """

    log.info('Dataset ' + filename + ' has (keyword = value) of (' + key + ' = ' + str(value) + ').')
    if msg == Messages.NOPROC.value:
        log.info('Dataset cannot be aligned.')
    else:
        log.info('Dataset can be aligned, but the result may be compromised.')
