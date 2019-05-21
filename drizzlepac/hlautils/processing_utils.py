"""Utilities to support creation of Hubble Advanced Pipeline(HAP) products.

"""
import sys

from astropy.io import fits as fits
from astropy.time import Time
from stsci.tools import logutil
from stsci.tools.fileutil import countExtn
from stwcs import wcsutil


__taskname__ = 'processing_utils'

log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout)

def refine_product_headers(product, level=None):
    """Refines output product headers to include values not available to AstroDrizzle.

    A few header keywords need to have values computed to reflect the type of product being
    generated, which in some cases can only be done using information passed in from the calling
    routine.  This function insures that all header keywords have been populated with values
    appropriate for the type of product being processed.

    Parameters
    -----------
    product : str or object
        Filename or HDUList object for product to be updated

    level : int, optional
        If defined, will add the 'LEVEL' keyword to the header of the product as defined by the calling
        routine. For example, 'level=2' to add/update the 'LEVEL' keyword for a level 2 (filter image)
        product.
    """

    hdu, closefits = _process_input(product)
    phdu = hdu[0].header

    # Start by updating the S_REGION keyword.
    compute_sregion(hdu)

    # Compute numexp as number of exposures NOT chips
    input_exposures = set([kw[1].split('[')[0] for kw in phdu['d*data'].items()])
    phdu['numexp'] = len(input_exposures)

    # Convert dates to ISO format
    phdu['date-beg'] = (Time(phdu['expstart'], format='mjd').iso, "Starting Date and Time")
    phdu['date-end'] = (Time(phdu['expend'], format='mjd').iso, "Ending Date and Time")

    # Re-format ACS filter specification
    if phdu['instrume'] == 'ACS':
        acs_filters = [kw[1] for kw in phdu['filter?'].items()]
        acs_filters = ';'.join(acs_filters)
        phdu['filter'] = acs_filters

    # Apply any additional inputs to drizzle product header
    if level:
        hdu[0].header['level'] = (level, "Classification level of this product")

        # Reset filter specification for total detection images which combine filters
        if level > 2:
            phdu['filter'] = 'detection'

    # close file if opened by this function
    if closefits:
        hdu.close()


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
