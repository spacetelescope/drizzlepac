"""Utilities to support creation of Hubble Advanced Pipeline(HAP) products.

"""
import sys
import os
import shutil

import numpy as np

from astropy.io import fits as fits
from astropy.io.fits import Column
from astropy.time import Time
from stsci.tools import logutil
from stsci.tools.fileutil import countExtn
from stwcs import wcsutil

from .cell_utils import SkyFootprint

LEVEL_DEFS = {1: 'single exposure product', 2: 'filter product', 3: 'total detection product'}
HAPCOLNAME = 'HAPEXPNAME'
PHOT_KEYWORDS = ['PHOTMODE', 'PHOTFLAM', 'PHOTFNU', 'PHOTZPT', 'PHOTPLAM', 'PHOTBW']

__taskname__ = 'processing_utils'

log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout)

def get_rules_file(product, rules_type="", rules_root=None):
    """Copies default HAP rules file to local directory.

    This function enforces the naming convention for rules files
    provided as part of this package; namely,

        <instrument name>_{<rules type>_}_header_hap.rules

    where <instrument name> is the lower-case value of the INSTRUME keyword
          for the exposure, and, optionally, <rules type> is the lower-case
          value of the `rules_type` parameter.

    The default rules file will be found in the package's installation directory,
    renamed for the input exposure and copied into the local working directory.
    This will allow AstroDrizzle to find the rules file when drizzling the exposure.

    Parameters
    ----------
    product : str
        Filename of the input exposure that the rules file applies to

    rules_type : str, optional
        Specifies the type of processing being performed on the input
        exposure.  Valid values: blank/empty string (''), 'SVM' or 'svm'
        for SVM processing (default) and 'MVM' or 'mvm' for MVM processing.

    rules_root : str, optional
        Name of output product to use as the rootname for the output rules file.
        Specifying this filename indicates that the default rules for combining
        multiple inputs should be used.
        If None, output rules file will be derived from the product filename and
        indicates that the rules file should be specific to a single input.

    Returns
    -------
    new_rules_name : str
        filename of the new rules file copied into the current working directory
    """
    # Interpret rules_type to guard against a variety of possible inputs
    # and to insure all values are converted to lower-case for use in
    # filenames
    rules_type = "" if rules_type == None else rules_type.strip(' ').lower()
    hdu, closefits = _process_input(product)
    phdu = hdu[0].header
    instrument = phdu['instrume']

    # Append '_single' to rules_type if single image product
    if rules_root is None:
        if rules_type == "":
            rules_type = 'single'
        else:
            rules_type = '_'.join([rules_type, 'single'])
        rootname = '_'.join(product.split("_")[:-1])
    else:
        rootname = '_'.join(rules_root.split("_")[:-1])


    # Create rules name prefix here
    # The trailing rstrip guards against a blank rules_type, which
    # is the default for SVM processing.
    rules_prefix = '_'.join([instrument.lower(), rules_type]).rstrip('_')
    code_dir = os.path.abspath(__file__)
    base_dir = os.path.dirname(os.path.dirname(code_dir))
    def_rules_name = "{}_header_hap.rules".format(rules_prefix)
    new_rules_name = "{}_header_hap.rules".format(rootname)
    rules_filename = os.path.join(base_dir, 'pars', def_rules_name)
    new_rules_filename = os.path.join(os.getcwd(), new_rules_name)
    log.debug(f'Copying \n\t{rules_filename} \nto \n\t{new_rules_filename}')
    if new_rules_name not in os.listdir('.'):
        shutil.copy(rules_filename, new_rules_filename)

    return new_rules_name


def build_logname(input_filename, process_type='svm'):
    """Build the log filename based on the input filename"""

    suffix = '.{}'.format(input_filename.split('.')[-1])
    if isinstance(input_filename, str):  # input file is a poller file -- easy case
        logname = input_filename.replace(suffix, '.log')
    else:
        logname = '{}_process.log'.format(process_type)

    return logname


def refine_product_headers(product, total_obj_list):
    """Refines output product headers to include values not available to AstroDrizzle.

    A few header keywords need to have values computed to reflect the type of product being
    generated, which in some cases can only be done using information passed in from the calling
    routine.  This function insures that all header keywords have been populated with values
    appropriate for the type of product being processed.

    Parameters
    -----------
    product : str or object
        Filename or HDUList object for product to be updated

    total_obj_list: list
        List of TotalProduct objects which are composed of Filter and Exposure
        Product objects

    """
    hdu, closefits = _process_input(product)
    phdu = hdu[0].header
    # Insure rootname and filename keywords matches actual filename
    phdu['rootname'] = '_'.join(product.split('_')[:-1])
    phdu['filename'] = product

    # Determine level of the product
    # Get the level directly from the product instance
    level = None
    for tdp in total_obj_list:
        member = tdp.find_member(product)
        if member is None:
            continue
        level = member.haplevel
        if level:
            break
    if level is None:
        level = 1

    # Update PINAME keyword
    phdu['piname'] = phdu['pr_inv_l']

    # Start by updating the S_REGION keyword.
    compute_sregion(hdu)

    # Compute numexp as number of exposures NOT chips
    input_exposures = list(set([kw[1].split('[')[0] for kw in phdu['d*data'].items()]))
    if level == 1:
        ipppssoots = [fits.getval(fname, 'rootname') for fname in input_exposures]
        phdu['ipppssoo'] = ';'.join(ipppssoots)
    phdu['numexp'] = len(input_exposures)

    # Convert dates to ISO format
    phdu['date-beg'] = (Time(phdu['expstart'], format='mjd').iso, "Starting Date and Time")
    phdu['date-end'] = (Time(phdu['expend'], format='mjd').iso, "Ending Date and Time")

    phdu['equinox'] = hdu[('sci', 1)].header['equinox'] if 'equinox' in hdu[('sci', 1)].header else 2000.0

    # Re-format ACS filter specification
    if phdu['instrume'] == 'ACS':
        phdu['filter'] = get_acs_filters(hdu, delimiter=';')

    # Insure PHOT* keywords are always in SCI extension
    for pkw in PHOT_KEYWORDS:
        if pkw in phdu:
            hdu[('sci', 1)].header[pkw] = (phdu[pkw], phdu.cards[pkw].comment)
            del phdu[pkw]

    # Apply any additional inputs to drizzle product header
    if level:
        hdu[0].header['haplevel'] = (level, "Classification level of this product")

        # Reset filter specification for total detection images which combine filters
        if 'total' in phdu['rootname']:
            phdu['filter'] = 'detection'

        # Build HAP table
        # if 'total' in product: level = 3
        update_hdrtab(hdu, level, total_obj_list, input_exposures)

    # close file if opened by this function
    if closefits:
        hdu.close()

def get_acs_filters(image, delimiter=';', all=False):
    hdu, closefits = _process_input(image)
    filters = [kw[1] for kw in hdu[0].header['filter?'].items()]
    acs_filters = []
    for f in filters:
        if ('clear' not in f.lower() and not all) or all:
            acs_filters.append(f)

    if not acs_filters:
        acs_filters = ['clear']
    acs_filters = delimiter.join(acs_filters)

    return acs_filters


def update_hdrtab(image, level, total_obj_list, input_exposures):
    """Build HAP entry table extension for product"""
    # Convert input_exposure filenames into HAP product filenames
    name_col = []
    orig_tab = image['hdrtab'].data
    # get the name of the product so it can be selected from
    # the total_obj_list for updating
    update_filename = image[0].header['filename']
    for tot_obj in total_obj_list:
        # Get the HAPProduct object for the input image to be updated
        # The '.find_member()' method looks for exposure, filter and
        # total level product.
        img_obj = tot_obj.find_member(update_filename)
        if img_obj is None:
            # Didn't find the input image in this total_obj instance,
            # try another...
            continue
        # if tot_obj.drizzle_filename != update_filename:
        #     continue
        # Only for the total_obj_list entry that matches the input image
        # should we build the list of new rootnames
        for row in orig_tab:
            rootname = str(row['rootname'])

            # The rootname is ipppssoot, but the expname is only contains ipppssoo,
            # so remove the last character for the comparisons
            rootname = rootname[0:-1]

            for expname in input_exposures:
                if rootname in expname:
                    # Convert input exposure names into HAP names
                    for exposure in tot_obj.edp_list:
                        if rootname in exposure.full_filename:
                            name_col.append(exposure.product_basename)
                            break

    hdrtab_cols = orig_tab.columns
    if name_col:
        # define new column with HAP expname
        max_len = min(max([len(name) for name in name_col]), 51)
        hapcol = Column(array=np.array(name_col, dtype=np.str), name=HAPCOLNAME, format='{}A'.format(max_len + 4))
        newcol = fits.ColDefs([hapcol])
        hdrtab_cols += newcol

    # define new extension
    haphdu = fits.BinTableHDU.from_columns(hdrtab_cols)
    haphdu.header['extname'] = 'HDRTAB'
    haphdu.header['extver'] = 1
    # remove old extension
    del image['hdrtab']
    # replace with new extension
    image.append(haphdu)


def compute_sregion(image, extname='SCI'):
    """Compute the S_REGION keyword for a given WCS.

    Parameters
    -----------
    image : Astropy io.fits  HDUList object
        Image to update with the S_REGION keyword in each of the SCI extensions.

    extname : str, optional
        EXTNAME value for extension containing the WCS(s) to be updated
    """
    # This function could, conceivably, be called directly...
    hdu, closefits = _process_input(image)

    # Find all extensions to be updated
    numext = countExtn(hdu, extname=extname)

    for extnum in range(1, numext + 1):
        sciext = (extname, extnum)
        if 'd001data' not in hdu[0].header:
            sregion_str = 'POLYGON ICRS '
            # Working with FLT/FLC file, so simply use
            #  the array corners directly
            extwcs = wcsutil.HSTWCS(hdu, ext=sciext)
            footprint = extwcs.calc_footprint(center=True)
            for corner in footprint:
                sregion_str += '{} {} '.format(corner[0], corner[1])
        else:
            if hdu[(extname, extnum)].data.min() == 0 and hdu[(extname, extnum)].data.max() == 0:
                continue
            # Working with a drizzled image, so we need to
            # get all the corners from each of the input files
            footprint = find_footprint(hdu, extname=extname, extnum=extnum)
            sregion_str = ''
            for region in footprint.corners:
                # S_REGION string should contain a separate POLYGON
                # for each region or chip in the SCI array
                sregion_str += 'POLYGON ICRS '
                for corner in region:
                    sregion_str += '{} {} '.format(corner[0], corner[1])

        hdu[sciext].header['s_region'] = sregion_str

    # close file if opened by this functions
    if closefits:
        hdu.close()


def find_footprint(hdu, extname='SCI', extnum=1):
    """Extract the footprints from each input file

    Determine the composite of all the input chip's corners
    as the footprint for the entire set of overlapping images
    that went into creating this drizzled image.

    Parameters
    ===========
    hdu : str or `fits.HDUList`
        Filename or HDUList for a drizzled image

    extname : str, optional
        Name of science array extension (extname value)

    extnum : int, optional
        EXTVER value for science array extension

    Returns
    ========
    footprint : ndarray
        Array of RA/Dec for the 4 corners that comprise the
        footprint on the sky for this mosaic.  Values were
        determined from northern-most corner counter-clockwise
        to the rest.

    """
    # extract WCS from this product
    meta_wcs = wcsutil.HSTWCS(hdu, ext=(extname, extnum))
    # create SkyFootprint object for all input_files to determine footprint
    footprint = SkyFootprint(meta_wcs=meta_wcs)
    # create mask of all input chips as they overlap on the product WCS
    footprint.extract_mask(hdu.filename())
    # Now, find the corners from this mask
    footprint.find_corners()

    return footprint


def interpret_sregion(image, extname='SCI'):
    """Interpret the S_REGION keyword as a list of RA/Dec points"""
    # This function could, conceivably, be called directly...
    hdu, closefits = _process_input(image)

    # Find all extensions to be updated
    numext = countExtn(hdu, extname=extname)
    sregions = []
    for extnum in range(1, numext + 1):
        sregions.append(fits.getval(image, 's_region', ext=(extname, extnum)))

    coords = []
    for region_str in sregions:
        regions = region_str.split('POLYGON ICRS ')
        regions = [val for val in regions if val.strip(' ') != '']
        for region in regions:
            radec_str = np.array(region.strip(' ').split(' '), dtype=np.float64)
            coords.append(radec_str.reshape((radec_str.shape[0] // 2, 2)))

    return coords


def _process_input(input):
    """Verify that input is an Astropy HDUList object opened in 'update' mode

    Parameters
    ----------
    input : str or object
        Filename of input image or HDUList object for image to be processed

    Returns
    --------
    hdu : object
        Astropy HDUList object of input opened in

    closefits : bool
        Boolean which indicates whether input should be closed when processing has been completed.

    """
    # Process input product
    if isinstance(input, str):
        hdu = fits.open(input, mode='update')
        closefits = True
    else:
        hdu = input
        closefits = False

    # Insure that input has been opened in update mode to allow for new values to be saved to headers
    if hdu._file.mode != 'update':
        filename = hdu.filename()
        if filename:
            hdu.close()
            hdu = fits.open(filename, mode='update')
            closefits = True
        else:
            log.error("Input object could not be opened in 'update' mode.")
            raise ValueError

    return hdu, closefits


def append_trl_file(trlfile, drizfile, clean=True):
    """ Append log file to already existing log or trailer file.

    Parameters
    -----------
    clean : bool
        Remove the `drizfile` or not when finished appending it to `trlfile`
    """
    if not os.path.exists(drizfile):
        return
    # Open already existing trailer file for appending
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

    if clean:
        # Now, clean up astrodrizzle trailer file
        os.remove(drizfile)


def make_section_str(str="", width=60, symbol='-'):
    """Generate string for starting/ending sections of log messages"""
    strlen = len(str)
    dash_len = width - strlen if strlen < width else 0
    side_dash = symbol * (dash_len // 2)
    section_str = "{}{}{}".format(side_dash, str, side_dash)
    return section_str


def find_flt_keyword (hdu, keyword, extname="SCI", extver=1):
    """Look for the FITS keyword in the Primary and specified extension header

    This routine deliberately does not use astropy.io.fits getval() as it
    should not be used in application code due to its self-documented inefficiency.
    This routine is only a STOPGAP to be used as a backstop to handle the
    issue when refine_product_headers() does not successfully update the drizzle
    product headers.

    Parameters
    ----------
    hdu : HDUList
    List of Header Data Unit objects

    keyword : str
    Keyword of interest

    extname : str
    Named extension to check in case the keyword is not found in the Primary
    e.g., (extname, extver) ("SCI", 2)

    extver : int
    Specified extver of an extname to check in case the keyword is not found in the Primary
    e.g., (extname, extver) ("SCI", 2)

    Returns
    -------
    value : float
    """

    # Look for the keyword in the Primary header
    try:
        value = hdu[0].header[keyword]
        log.info("Read keyword {} from the Primary header.".format(keyword))
    # Oops.  Look for the keyword in the extension
    except KeyError:
        try:
            sciext = (extname, extver)
            value = hdu[sciext].header[keyword]
            log.info("Read keyword from the {} header.".format(sciext))
        except KeyError as err:
            log.error("Keyword {} does not exist in the Primary or {} extension.".format(keyword, sciext))
            raise

    # Now ensure the returned variable is of the proper type
    try:
        value = float(value)
        log.info("Keyword {}: {}.".format(keyword, value))
    except (ValueError, TypeError) as err:
            log.error("The value of keyword, {}, cannot be cast as a floating point value.".format(keyword))
            raise

    return value


def _find_open_files():
    """Useful utility function to identify any open file handles during processing."""
    import psutil

    cwd = os.getcwd()
    log.info("Checking for open files in {}".format(cwd))

    open_files = []
    for proc in psutil.process_iter():
        try:
            for pfile, pid in proc.open_files():
                if os.path.dirname(pfile) == cwd:
                    open_files.append(pfile)
        except psutil.AccessDenied:
            continue
        except FileNotFoundError:
            continue

    return open_files
