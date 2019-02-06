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
import pdb
__all__ = ['analyze_data']

from enum import Enum
 
# Annotates level to which image can be aligned according observational parameters
# as described through FITS keywords
class Messages(Enum):
    OK, WARN, NOPROC = 1, -1, -2
 
def analyze_data(inputFileList, **kwargs):
    """
    Determine if images within the dataset can be aligned
   
    Parameters
    ==========
    inputFileList: list
        List containing FLT and/or FLC filenames for all input images which comprise an associated 
        dataset where 'associated dataset' may be a single image, multiple images, an HST association, 
        or a number of HST associations

    Returns
    =======
    outputTable: object
        Astropy Table object containing data pertaining to the associated dataset, including 
        the doProcess bool.  It is intended this table is updated by subsequent functions for
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
    OBSKEY = 'OBSTYPE'
    MTKEY  = 'MTFLAG'
    SCNKEY = 'SCAN_TYP'
    FILKEY = 'FILTER'
    FILKEY1 = 'FILTER1'
    FILKEY2 = 'FILTER2'
    APKEY  = 'APERTURE'
    TARKEY = 'TARGNAME'
    EXPKEY = 'EXPTIME'
    FGSKEY = 'FGSLOCK'
    CHINKEY = 'CHINJECT'

    acsFiltNameList = [FILKEY1, FILKEY2]

    catalog = None     # Astrometric catalog used for alignment
    catalogSources = 0 # Number of astrometric catalog sources determined based upon coordinate overlap with image WCS
    foundSources = 0   # Number of sources detected in images
    matchSources = 0   # Number of sources cross matched between astrometric catalog and detected in image
    rms_x = -1.0
    rms_y = -1.0
    rms_ra = -1.0
    rms_dec = -1.0
    chisq_x = -1.0
    chisq_y = -1.0
    completed = False # If true, there was no exception and the processing completed all logic
    dateObs = None     # Human readable date
    mjdutc = -1.0      # MJD UTC start of exposure
    fgslock = None
    processMsg = None
    status = 9999
    compromised = 0
    headerletFile = None

    fit_rms = -1.0
    total_rms = -1.0
    datasetKey = -1.0

    namesArray = ('imageName', 'instrument', 'detector', 'filter', 'aperture', 'obstype',
            'subarray', 'dateObs', 'mjdutc', 'doProcess', 'processMsg', 'catalog', 'foundSources',
            'catalogSources','matchSources', 'rms_x', 'rms_y', 'rms_ra', 'rms_dec', 'completed',
            'fit_rms', 'total_rms', 'datasetKey', 'status', 'headerletFile')
    dataType = ('S20', 'S20', 'S20', 'S20', 'S20', 'S20', 'b', 'S20', 'f8', 'b', 'S30',
            'S20', 'i4', 'i4', 'i4', 'f8', 'f8', 'f8', 'f8', 'b', 'f8', 'f8', 'i8', 'i4', 'S30')

    # Create an astropy table
    outputTable = Table(names=namesArray,dtype=dataType)

    # Loop over the list of images to determine viability for alignment processing
    #
    # Capture the data characteristics before any evaluation so the information is
    # available for the output table regardless of which keyword is used to 
    # to determine the data is not viable for alignment.

    for inputFile in inputFileList:

        header_hdu  = 0
        header_data = getheader(inputFile, header_hdu)

        # Keywords to use potentially for downstream analysis
        instrume = (header_data['INSTRUME']).upper()
        detector = (header_data['DETECTOR']).upper()
        subarray = header_data['SUBARRAY']
        dateObs  = header_data['DATE-OBS']
        mjdutc  = header_data['EXPSTART']

        # Obtain keyword values for analysis of viability
        obstype  = (header_data[OBSKEY]).upper()
        mtflag   = (header_data[MTKEY]).upper()
        
        scan_typ = ''
        if instrume == 'WFC3':
            scan_typ = (header_data[SCNKEY]).upper()

        sfilter = ''
        if instrume == 'WFC3':
            sfilter  = (header_data[FILKEY]).upper()
        # Concatenate the two ACS filter names together with an underscore
        # If the filter name is blank, skip it
        if instrume == 'ACS':
            for filtname in acsFiltNameList:

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

        noProcKey   = None
        noProcValue = None
        doProcess = True
        # Imaging vs spectroscopic or coronagraphic
        if obstype != 'IMAGING':
            noProcKey   = OBSKEY
            noProcValue = obstype 

        # Moving target
        elif mtflag == 'T':
            noProcKey   = MTKEY
            noProcValue = mtflag 

        # Bostrophidon without or with dwell (WFC3 only)
        elif any ([scan_typ == 'C', scan_typ == 'D']):
            noProcKey   = SCNKEY
            noProcValue = scan_typ

        # Filter which does not begin with: 'F'(F###), 'C'(CLEAR), 'N'(N/A), and is not blank
        # The sfilter variable may be the concatenation of two filters (F160_CLEAR)
        elif sfilter[0] != 'F' and sfilter[0] != '' and sfilter[0] != 'C' and sfilter[0] != 'N': 
            noProcKey   = FILKEY
            noProcValue = sfilter

        elif '_' in sfilter:
            pos = sfilter.index('_')
            pos += 1

            if sfilter[pos] != 'F' and sfilter[pos] != '' and sfilter[pos] != 'C' and sfilter[pos] != 'N': 
                noProcKey   = FILKEY
                noProcValue = sfilter

        # Ramp, polarizer, grism, or prism 
        elif any (x in aperture for x in ['RAMP', 'POL', 'GRISM', '-REF', 'PRISM']):
            noProcKey   = APKEY
            noProcValue = aperture 

        # Calibration target
        elif any (x in targname for x in ['DARK', 'TUNG', 'BIAS', 'FLAT', 'DEUT', 'EARTH-CAL']):
            noProcKey   = TARKEY
            noProcValue = targname

        # Exposure time of effectively zero
        elif math.isclose(exptime, 0.0, abs_tol=1e-5):
            noProcKey   = EXPKEY
            noProcValue = exptime 

        # Commanded FGS lock
        elif any (x in fgslock for x in ['GY', 'COARSE']):
            noProcKey   = FGSKEY
            noProcValue = fgslock

        # Charge injection mode
        elif chinject != 'NONE':
            noProcKey   = CHINKEY
            noProcValue = chinject

        # If noProcKey is set to a keyword, then this image has been found to not be viable for
        # alignment purposes.
        if (noProcKey is not None):
            if (noProcKey != FGSKEY):
                doProcess = False
                msgType = Messages.NOPROC.value
            else:
                msgType = Messages.WARN.value

            processMsg = noProcKey + '=' + str(noProcValue)

            # Issue message to log file for this data indicating no processing to be done or 
            # processing should be allowed, but there may be some issue with the result (e.g., 
            # GYROS mode so some drift)
            generate_msg(inputFile, msgType, noProcKey, noProcValue)

        # Populate a row of the table
        outputTable.add_row([inputFile, instrume, detector, sfilter, aperture, obstype,
                             subarray, dateObs, mjdutc, doProcess, processMsg, catalog, 
                             foundSources, catalogSources, matchSources, rms_x, rms_y, 
                             rms_ra, rms_dec, completed, fit_rms, total_rms, datasetKey,
                             status, headerletFile])
    #outputTable.pprint(max_width=-1)

    return(outputTable)


def generate_msg(filename, msg, key, value):
    """ Generate a message for the output log indicating the file/association will not
        be processed as the characteristics of the data are known to be inconsistent
        with alignment.
    """

    print('\nDataset ' + filename + ' has (keyword = value) of (' + key + ' = ' + str(value) + ').')
    print(msg)
    if msg == Messages.NOPROC.value:
        print('Dataset cannot be aligned.\n')
    else:
        print('Dataset can be aligned, but the result may be compromised.')

