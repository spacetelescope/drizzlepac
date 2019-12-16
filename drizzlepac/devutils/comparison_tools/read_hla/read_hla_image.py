"""read_hla_image.py: read an HLA image (from a file or from HLA webservice) and correct WCS

R. White, 2019 June 20
"""

from __future__ import print_function
import omegaxyz, os
from astropy.io import fits

def read_hla_image(dataset, applyomega=True, verbose=False,
        url="https://hla.stsci.edu/cgi-bin/getdata.cgi"):

    """Return a astropy.io.fits HDUlist for the given HLA dataset.
    
    If dataset is a filename, the file is opened; if dataset is a
    HLA datasetname, the file is retrieved from the HLA web server.
    If applyomega is true, the HSC astrometry correction is applied to the WCS in the header.
    """

    if not os.path.exists(dataset):
        dataloc = "{}?dataset={}".format(url,dataset)
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
    """

    if not os.path.exists(dataset):
        dataloc = "{}?dataset={}".format(url,dataset)
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
    print("Retrieved FITS HDUlist for {} with {} extensions".format(dataset,len(fh)))
    outfile = dataset+'_corr.fits'
    write_hla_image(dataset, outfile, verbose=True, overwrite=True)
    print("Wrote {} to {}".format(dataset,outfile))

