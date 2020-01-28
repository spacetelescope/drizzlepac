"""read_hla_image.py: read an HLA image (from a file or from HLA webservice) and correct WCS

R. White, 2019 June 20
"""
import os
from astropy.io import fits
from drizzlepac.devutils.comparison_tools.read_hla import omegaxyz


def read_hla_image(dataset, applyomega=True, verbose=False, url="https://hla.stsci.edu/cgi-bin/getdata.cgi"):

    """Return a astropy.io.fits HDUlist for the given HLA dataset.
    

    If applyomega is true, the HSC astrometry correction is applied to the WCS in the header.

    Parameters
    ----------
    dataset : str
        image to be read in. If dataset is a filename, the file is opened; if dataset is a
        HLA datasetname, the file is retrieved from the HLA web server.

    applyomega : bool, optional
        If true, applies the HSCv3 astrometric correction to the RA and Dec columns Default value = True

    verbose : bool, optional
        If true, display detailed information about the results. Default value = False

    url : str, optional
        web service URL. Default value = 'https://hla.stsci.edu/cgi-bin/getdata.cgi'

    Returns
    -------
    fh : astropy.io.fits.hdu.hdulist.HDUList object
        fits file contents
    """

    if not os.path.exists(dataset):
        dataloc = "{}?dataset={}".format(url, dataset)
    else:
        dataloc = dataset
    if applyomega:
        fh = omegaxyz.updatefits(dataloc, verbose=verbose)
    else:
        fh = fits.open(dataloc)
    return fh


def write_hla_image(dataset, outfile, applyomega=True, overwrite=False, verbose=False,
                    url="https://hla.stsci.edu/cgi-bin/getdata.cgi"):

    """Write FITS file for the HLA dataset to outfile.
    
    If dataset is a filename, the file is opened; if dataset is a
    HLA datasetname, the file is retrieved from the HLA web server.
    If applyomega is true, the HSC astrometry correction is applied to the WCS in the header.

    Parameters
    ----------
    dataset : str
        dataset name. If dataset is a filename, the file is opened; if dataset is a
        HLA datasetname, the file is retrieved from the HLA web server.

    outfile : str
        output file name.

    applyomega : bool, optional
        If true, applies the HSCv3 astrometric correction to the RA and Dec columns Default value = True

    overwrite : bool, optional
        If a file with the same name as input argument outfile, overwrite? Default value = False

    verbose : bool, optional
        If true, display detailed information about the results. Default value = False

    url : str, optional
        web service URL. Default value = 'https://hla.stsci.edu/cgi-bin/getdata.cgi'

    Returns
    -------
    """

    if not os.path.exists(dataset):
        dataloc = "{}?dataset={}".format(url, dataset)
    else:
        dataloc = dataset
    if applyomega:
        omegaxyz.updatefits(dataloc, outfile, verbose=verbose, overwrite=overwrite)
    else:
        fh = fits.open(dataloc)
        fh.writeto(outfile, overwrite=overwrite)


if __name__ == "__main__":
    dataset = 'hst_10188_10_acs_wfc_f814w'
    fh = read_hla_image(dataset, verbose=True)
    print("Retrieved FITS HDUlist for {} with {} extensions".format(dataset, len(fh)))
    outfile = dataset+'_corr.fits'
    write_hla_image(dataset, outfile, verbose=True, overwrite=True)
    print("Wrote {} to {}".format(dataset, outfile))

