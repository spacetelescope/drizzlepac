""" Utility to analyze an input dataset and determine whether the dataset
can be aligned.

The function analyze_data opens an input list containing FLT and/or FLC FITS
filenames in order to access the primary header data.  Based upon the values
of specific FITS keywords, the function determines whether or not each file
within this dataset can or should be reconciled against an astrometric catalog
and, for multiple images, used to create a mosaic.

"""
from astropy.io.fits import getheader
from astropy.table import Table
import math
import sys
from enum import Enum

from stsci.tools import logutil

__taskname__ = 'analyze'

log = logutil.create_logger(
    __name__, level=logutil.logging.INFO, stream=sys.stdout
)

__all__ = ['analyze_data']

# Define global default keyword names for these fields
_OBSKEY = 'OBSTYPE'
_MTKEY = 'MTFLAG'
_SCNKEY = 'SCAN_TYP'
_FILKEY = 'FILTER'
_FILKEY1 = 'FILTER1'
_FILKEY2 = 'FILTER2'
_APKEY = 'APERTURE'
_TARKEY = 'TARGNAME'
_EXPKEY = 'EXPTIME'
_FGSKEY = 'FGSLOCK'
_CHINKEY = 'CHINJECT'


# Annotates level to which image can be aligned according observational
# parameters as described through FITS keywords
class Messages(Enum):
    OK, WARN, NOPROC = 1, -1, -2


def analyze_data(images, **kwargs):
    """
    Determine if images within the CHINKEY dataset can be aligned.

    Parameters
    ----------
    images: list of str
        List containing FLT and/or FLC filenames for all input images which
        comprise an associated dataset where 'associated dataset' may be a
        single image, multiple images, an HST association, or a number of HST
        associations.

    Returns
    -------
    output_table: astropy.table.Table
        `~astropy.table.Table` containing data pertaining to the associated
        dataset, including the ``do_process`` bool. It is intended this table
        is updated by subsequent functions for bookkeeping purposes.

    Notes
    -----
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

    The keyword/value pairs below define the category which the data can be
    processed, but the results may be compromised
    ``FGSLOCK`` : ``FINE/GYRO``, ``FINE/GY``, ``COARSE``, ``GYROS``

    FITS Keywords only for WFC3 data: ``SCAN_TYP``, ``FILTER``, and
    ``CHINJECT`` (UVIS)

    FITS Keywords only for ACS data: ``FILTER1`` and ``FILTER2``

    Please be aware of the FITS keyword value ``'NONE'`` vs the Python `None`.

    **TODO:** improve robustness when analyzing filter and aperture names,
    possibly use ``PHOTMODE`` instead.

    """
    acs_filters = [_FILKEY1, _FILKEY2]

    fit_method = None  # Fit algorithm used for alignment
    catalog = None  # Astrometric catalog used for alignment
    # Number of astrometric catalog sources determined based upon coordinate
    # overlap with image WCS
    catalog_sources = 0
    found_sources = 0  # Number of sources detected in images
    # Number of sources cross matched between astrometric catalog and
    # detected in image:
    match_sources = 0
    offset_x = None
    offset_y = None
    rot = None
    scale = None
    rms_x = -1.0
    rms_y = -1.0
    rms_ra = -1.0
    rms_dec = -1.0
    completed = False
    date_obs = None  # Human readable date
    mjdutc = -1.0  # MJD UTC start of exposure
    fgslock = None
    process_msg = None
    status = 9999
    compromised = 0
    headerlet_file = None
    fit_qual = -1
    fit_rms = -1.0
    total_rms = -1.0
    dataset_key = -1.0

    colnames = (
        'imageName', 'instrument', 'detector', 'filter', 'aperture', 'obstype',
        'subarray', 'dateObs', 'mjdutc', 'doProcess', 'processMsg',
        'fit_method', 'catalog', 'foundSources', 'catalogSources',
        'matchSources', 'offset_x', 'offset_y', 'rotation', 'scale', 'rms_x',
        'rms_y', 'rms_ra', 'rms_dec', 'completed', 'fit_rms', 'total_rms',
        'datasetKey', 'status', 'fit_qual', 'headerletFile', 'compromised'
    )
    data_type = (
        'S20', 'S20', 'S20', 'S20', 'S20', 'S20', 'b', 'S20', 'f8', 'b', 'S30',
        'S20', 'S20', 'i4', 'i4', 'i4', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8',
        'f8', 'f8', 'b', 'f8', 'f8', 'i8', 'i4', 'i4', 'S30', 'i4'
    )

    # Create an astropy table
    output_table = Table(names=colnames, dtype=data_type)

    # Loop over the list of images to determine viability for alignment
    # processing.
    #
    # Capture the data characteristics before any evaluation so the information
    # is available for the output table regardless of which keyword is used to
    # to determine the data is not viable for alignment.
    for input_file in images:
        hdr = getheader(input_file, ext=0)

        # Keywords to use potentially for downstream analysis
        instrume = hdr['INSTRUME'].upper()
        detector = hdr['DETECTOR'].upper()
        subarray = hdr['SUBARRAY']
        date_obs = hdr['DATE-OBS']
        mjdutc = hdr['EXPSTART']

        aperture = hdr[_APKEY].upper()
        targname = hdr[_TARKEY].upper()
        exptime = hdr[_EXPKEY]
        fgslock = hdr[_FGSKEY].upper()

        # Obtain keyword values for analysis of viability
        obstype = hdr[_OBSKEY].upper()
        mtflag = hdr[_MTKEY].upper()

        scan_typ = sfilter = ''
        chinject = 'NONE'
        if instrume == 'WFC3':
            scan_typ = hdr[_SCNKEY].upper()
            sfilter = hdr[_FILKEY].upper()
            if detector == 'UVIS':
                chinject = hdr[_CHINKEY].upper()

        elif instrume == 'ACS':
            # Concatenate the two ACS filter names together with an underscore
            # If the filter name is blank, skip it
            sfilter = '_'.join(
                [hdr[f].strip().upper() for f in acs_filters if hdr[f].strip()]
            )

        # Determine if the image has one of these conditions.  The routine
        # will exit processing upon the first satisfied condition.
        no_proc_key = None
        no_proc_value = None
        do_process = True

        if obstype != 'IMAGING':
            # Imaging vs spectroscopic or coronagraphic
            no_proc_key = _OBSKEY
            no_proc_value = obstype

        elif mtflag == 'T':
            # Moving target
            no_proc_key = _MTKEY
            no_proc_value = mtflag

        elif any([scan_typ == 'C', scan_typ == 'D']):
            # Bostrophidon without or with dwell (WFC3 only)
            no_proc_key = _SCNKEY
            no_proc_value = scan_typ

        elif any(x in aperture for x in ['RAMP', 'POL', 'GRISM', '-REF',
                                         'PRISM']):
            # Ramp, polarizer, grism, or prism
            no_proc_key = _APKEY
            no_proc_value = aperture

        elif any(x in targname for x in ['DARK', 'TUNG', 'BIAS', 'FLAT',
                                         'DEUT', 'EARTH-CAL']):
            # Calibration target
            no_proc_key = _TARKEY
            no_proc_value = targname

        elif math.isclose(exptime, 0.0, abs_tol=1e-5):
            # Exposure time of effectively zero
            no_proc_key = _EXPKEY
            no_proc_value = exptime

        elif any(x in fgslock for x in ['GY', 'COARSE']):
            # Commanded FGS lock
            no_proc_key = _FGSKEY
            no_proc_value = fgslock

        elif chinject != 'NONE':
            # Charge injection mode
            no_proc_key = _CHINKEY
            no_proc_value = chinject

        elif sfilter and not sfilter.startswith(('F', 'C', 'N')):
            # Filter which does not begin with: 'F'(F###), 'C'(CLEAR),
            # 'N'(N/A), and is not blank. The sfilter variable may be the
            # concatenation of two filters (F160_CLEAR):
            no_proc_key = _FILKEY
            no_proc_value = sfilter

        elif '_' in sfilter:
            k = sfilter[sfilter.index('_') + 1]
            if k not in ['F', '', 'C', 'N']:
                no_proc_key = _FILKEY
                no_proc_value = sfilter

        # If no_proc_key is set to a keyword, then this image has been found
        # to not be viable for alignment purposes.
        if no_proc_key is not None:
            if no_proc_key != _FGSKEY:
                do_process = False
                msg_type = Messages.NOPROC.value
            else:
                msg_type = Messages.WARN.value

            process_msg = no_proc_key + '=' + str(no_proc_value)

            # Issue message to log file for this data indicating no processing
            # to be done or processing should be allowed, but there may be
            # some issue with the result (e.g., GYROS mode so some drift):
            generate_msg(input_file, msg_type, no_proc_key, no_proc_value)

        # Populate a row of the table
        output_table.add_row([
            input_file, instrume, detector, sfilter, aperture, obstype,
            subarray, date_obs, mjdutc, do_process, process_msg, fit_method,
            catalog, found_sources, catalog_sources, match_sources,
            offset_x, offset_y, rot, scale, rms_x, rms_y, rms_ra, rms_dec,
            completed, fit_rms, total_rms, dataset_key, status, fit_qual,
            headerlet_file, compromised
        ])
        process_msg = None

    return output_table


def generate_msg(filename, msg, key, value):
    """ Generate a message for the output log indicating the file/association
    will not be processed as the characteristics of the data are known to be
    inconsistent with alignment.

    """
    log.info('Dataset {:s} has (keyword = value) of ({} = {}).'
             .format(filename, key, value))

    if msg == Messages.NOPROC.value:
        log.info('Dataset cannot be aligned.')
    else:
        log.info('Dataset can be aligned, but the result may be compromised.')
