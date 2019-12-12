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


LEVEL_DEFS = {1: 'single exposure product', 2: 'filter product', 3: 'total detection product'}
HAPCOLNAME = 'HAPEXPNAME'
PHOT_KEYWORDS = ['PHOTMODE', 'PHOTFLAM', 'PHOTFNU', 'PHOTZPT', 'PHOTPLAM', 'PHOTBW']

__taskname__ = 'processing_utils'

log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout)

def get_rules_file(product):
    """Copies default HLA rules file to local directory."""
    hdu, closefits = _process_input(product)
    phdu = hdu[0].header
    instrument = phdu['instrume']
    code_dir = os.path.abspath(__file__)
    base_dir = os.path.dirname(os.path.dirname(code_dir))
    rules_name = "{}_header_hla.rules".format(instrument.lower())
    rules_filename = os.path.join(base_dir, 'pars', rules_name)
    if rules_name not in os.listdir('.'):
        shutil.copy(rules_filename, os.getcwd())

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
    level = 1 if len(phdu['rootname'].split('_')[-1]) > 6 else 2

    # Update PINAME keyword
    phdu['piname'] = phdu['pr_inv_l']

    # Start by updating the S_REGION keyword.
    compute_sregion(hdu)

    # Compute numexp as number of exposures NOT chips
    input_exposures = list(set([kw[1].split('[')[0] for kw in phdu['d*data'].items()]))
    if level == 1:
        ipppssoots = [fname.split('_')[0] for fname in input_exposures]
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

    for row in orig_tab:
        rootname = str(row['rootname'])
        for expname in input_exposures:
            if rootname in expname:
                if level == 1:
                    # Intrepret inputs as exposures (FLT/FLC) filename not HAP names
                    name_col.append(expname)
                else:
                    # Convert input exposure names into HAP names
                    foundit = False
                    for tot_obj in total_obj_list:
                        for exposure in tot_obj.edp_list:
                            if rootname in exposure.full_filename:
                                name_col.append(exposure.drizzle_filename)
                                foundit = True
                                break

    # define new column with HAP expname
    max_len = min(max([len(name) for name in name_col]), 51)
    hapcol = Column(array=np.array(name_col, dtype=np.str), name=HAPCOLNAME, format='{}A'.format(max_len + 4))
    newcol = fits.ColDefs([hapcol])

    # define new extension
    haphdu = fits.BinTableHDU.from_columns(orig_tab.columns + newcol)
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
        sregion_str = 'POLYGON ICRS '
        sciext = (extname, extnum)
        extwcs = wcsutil.HSTWCS(hdu, ext=sciext)
        footprint = extwcs.calc_footprint(center=True)
        for corner in footprint:
            sregion_str += '{} {} '.format(corner[0], corner[1])
        hdu[sciext].header['s_region'] = sregion_str

    # close file if opened by this functions
    if closefits:
        hdu.close()

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
