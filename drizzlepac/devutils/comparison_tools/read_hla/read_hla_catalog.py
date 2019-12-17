"""read_hla_catalog.py: read an HLA catalog (from a file or from HLA webservice) and correct astrometry

R. White, 2019 June 12
"""

from __future__ import print_function
from . import omegaxyz
import os
from astropy.table import Table

formatlist = [ "text", "votable", "tsv", "csv", "html", "kml", "json" ]
catlist = [ "sex", "dao" ]

def read_hla_catalog(dataset, cattype="sex", catformat="csv", trim=True, multiwave=False,
        applyomega=True, verbose=False, url="https://hla.stsci.edu/HLA/Catalogs/HLAcat.aspx"):

    """Return a astropy table with the given HLA catalog.
    
    Parameters:
    dataset     HLA dataset name (e.g. "hst_10188_10_acs_wfc_f814w")
    cattype     "sex" or "dao"
    catformat   One of the formats in formatlist
    trim        If false, includes sources flagged as bad
    multiwave   Return multiwavelength version of catalog (fewer columns with multiple filter info)
    verbose     If true, print information about the results
    applyomega  If true, applies the HSCv3 astrometric correction to the RA and Dec columns
    url         URL for the HLA catalog access (default should be fine)

    Returns an astropy.table with the catalog.  
    """

    cattype = cattype.lower()
    if cattype not in catlist:
        raise ValueError("cattype '{}' must be one of {}".format(cattype, " ".join(catlist)))
    if catformat.lower() not in formatlist:
        raise ValueError("catformat '{}' must be one of {}".format(catformat, " ".join(formatlist)))
    # get URL for this catalog
    params = ['catalog='+cattype, 'format='+catformat]
    if trim:
        suffix = "{}phot_trm.cat".format(cattype)
    else:
        suffix = "{}phot.cat".format(cattype)
        params.append('ignoreflag=T')
    if multiwave:
        params.append('multiwave=T')
        f = '_'.join(dataset.split('_')[:5])
        dataset = '{}_total'.format(f)
        catname = '{}_multiwave'.format(f)
    else:
        catname = dataset
    params.append('image={}'.format(dataset))
    params.append('filename={}_{}'.format(catname,suffix))
    caturl = "{}?{}".format(url, "&".join(params))
    if catformat == "votable":
        cat = Table.read(caturl,format=catformat)
    else:
        # csv (and other?) formats start with a comment
        cat = Table.read(caturl,format=catformat,comment='#')
    if verbose:
        if multiwave:
            print("Retrieved multiwave {} catalog for {} with {} rows".format(cattype,dataset,len(cat)))
        else:
            print("Retrieved {} catalog for {} with {} rows".format(cattype,dataset,len(cat)))
        ncols = 10
        for i in range(0,len(cat.colnames),ncols):
            cc = cat.colnames[i:i+ncols]
            print(("{:15} "*len(cc)).format(*cc))
    if applyomega:
        omega = omegaxyz.getomegaxyz(dataset)
        if omega == (0.0,0.0,0.0):
            if verbose:
                print("Omega correction is zero for this visit, no correction applied")
        else:
            if verbose:
                print("Applying omega correction to astrometry")
            # names of RA and Dec columns varies for different catalogs
            racol, deccol = get_radec_cols(cat)
            newra, newdec = omegaxyz.applyomegacat(cat[racol], cat[deccol], omega)
            cat[racol] = newra
            cat[deccol] = newdec
    return cat


def get_radec_cols(cat):

    """Return names of the RA and Dec columns for this catalog."""

    # create a column name dictionary
    d = dict([(x,1) for x in cat.colnames])
    for v in ["ra", "RA", "ALPHA_J2000", "alpha_j2000"]:
        if v in d:
            racol = v
            break
    else:
        raise ValueError("RA column not found in", " ".join(cat.colnames))
    for v in ["dec", "DEC", "Dec", "DELTA_J2000", "delta_j2000"]:
        if v in d:
            deccol = v
            break
    else:
        raise ValueError("RA column not found in", " ".join(cat.colnames))
    return racol, deccol


if __name__ == "__main__":
    # run through all combinations of catalog types and corrections
    dataset = 'hst_10188_10_acs_wfc_f814w'
    for cattype in ("sex", "dao"):
        for multiwave in (False, True):
            for applyomega in (False, True):
                cat = read_hla_catalog(dataset, cattype=cattype, applyomega=applyomega, multiwave=multiwave, verbose=True)
                cat[:2].pprint()
