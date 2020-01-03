"""Functions to get and apply HSC infinitesimal rotations to correct small shifts and rotations in spherical coordinates

Uses formulae from `Budavari & Lubow (2012, ApJ, 761, 188) <http://adsabs.harvard.edu/abs/2012ApJ...761..188B>`_

vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4 ai :

2019 June 12, Rick White
"""
import math
import os
import pdb
import sys
import numpy as np
from astropy.io import fits
import requests

# cache of omega vectors that have already been retrieved
getomegaxyz_cache = {}


def getomegaxyz(dataset, service="https://hla.stsci.edu/cgi-bin/getomega.cgi"):
    """Return tuple (omegax, omegay, omegaz) with rotation vector for an image. (0, 0, 0) will be returned if the image
    does not have a rotation correction.

    Parameters
    ----------
    dataset : str
        dataset name

    service : str, optional
        URL of web service used to get omega values. Default value = https://hla.stsci.edu/cgi-bin/getomega.cgi

    Returns
    -------
    omega : tuple
        3-element tuple containing omega X, Y, and Z components. (0, 0, 0) if the image does not have a rotation
        correction
    """

    if not dataset:
        return (0.0, 0.0, 0.0)
    dataset = getdataset(dataset.lower())
    if dataset in getomegaxyz_cache:
        return getomegaxyz_cache[dataset]
    r = requests.get(service, params=dict(image=dataset))
    result = r.json()
    omega = tuple(result.get("omega", (0.0, 0.0, 0.0)))
    getomegaxyz_cache[dataset] = omega
    return omega


getwcs_cache = {}
def getwcs(dataset, applyomega=True, service="https://hla.stsci.edu/cgi-bin/fitscut.cgi"):
    """Return dictionary with WCS information for an image.

    Parameters
    ----------
    dataset : str
        dataset name

    applyomega : bool, optional
        Apply HSC correction? Default value = True

    service : str, optional
        web service URL. Default value = https://hla.stsci.edu/cgi-bin/fitscut.cgi

    Returns
    -------
    updated WCS values
    """

    if not dataset:
        raise ValueError("Undefined dataset '{}'".format(dataset))
    key = (dataset, bool(applyomega))
    if key in getwcs_cache:
        return getwcs_cache[key]
    r = requests.get(service, params=dict(red=dataset, getwcs=1, applyomega=applyomega))
    result = r.json()
    getwcs_cache[key] = result
    return result


def crossproduct(a, b):
    """Return cross product (a X b)

    Parameters
    ----------
    a : 3-element list of floats
        first set of values in cross-product calculation

    b : 3-element list of floats
        second set of values in cross-product calculation

    Returns
    -------
    c :  3-element list of floats
        result of cross-product calculation
    """
    c = [0]*3
    c[0] = a[1]*b[2] - a[2]*b[1]
    c[1] = a[2]*b[0] - a[0]*b[2]
    c[2] = a[0]*b[1] - a[1]*b[0]
    return c


def dotproduct(a, b):
    """Return a . b

    Parameters
    ----------
    a : 3-element list of floats
        first set of values in dot-product calculation

    b : 3-element list of floats
        second set of values in dot-product calculation

    Returns
    -------
    c :  3-element list of floats
        result of dot-product calculation
    """
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]


def getdataset(fitsfile):
    """Extract visit identifier from FITS filename
    
    Returns None if this does not look like an HLA dataset name
    This works on either a plain dataset name (hst_10188_10_acs_wfc_f814) or a
    getdata.cgi URL (http://hla.stsci.edu/getdata.cgi?dataset=hst_10888_10_acs_wfc_f814w)

    Parameters
    ----------
    fitsfile : str
        fits file name

    Returns
    -------
    vis_id : str
        visit identifier value. Returns logical 'None' if this does not look like an HLA dataset name
    """

    f = os.path.split(fitsfile)[-1].split('_')
    if len(f) < 5 or f[0].lower() != 'hst':
        # perhaps this is a getdata URL?
        ff = fitsfile.split('dataset=')
        f = os.path.split(ff[-1])[-1].split('_')
        if len(f) < 5 or f[0].lower() != 'hst':
            return None
    if f[3].lower() == 'wfpc2':
        # look for shift for either WFPC2 or PC using WFPC2 catalogs
        f[4] = 'wfpc2'
    vis_id = '_'.join(f[:5])
    return vis_id


def applyomegawcs(filename, crval, cdmatrix, omega=None):

    """Get the omega values for this dataset and return the updated (crval, cdmatrix) values

    Parameters
    ----------
    filename : str
        dataset filename

    crval : list of floats
        Right Ascension and Declination position at the reference pixel

    cdmatrix : list of floats
        Description of linear distortions: plate scale, rotation, and skew of the image

    omega : tuple, optional. Default value = None
        3-element tuple containing omega X, Y, and Z components.

    Returns
    -------
    Updated versions of input variables crval and omega
    """
    if (not crval) or (not cdmatrix):
        return (crval, cdmatrix)

    if not omega:
        # parse the filename and extract the visit identification
        dataset = getdataset(filename)
        if not dataset:
            return (crval, cdmatrix)
        omega = getomegaxyz(dataset)

    if omega[0] != 0 or omega[1] != 0 or omega[2] != 0:
        d2r = math.pi/180
        # apply rotation to the tangent point
        ra0 = crval[0]*d2r
        dec0 = crval[1]*d2r
        cdec0 = math.cos(dec0)
        p0 = [cdec0*math.cos(ra0), cdec0*math.sin(ra0), math.sin(dec0)]
        dp0 = crossproduct(omega, p0)
        decnew = math.asin(p0[2]+dp0[2])
        ranew = math.atan2(p0[1]+dp0[1], p0[0]+dp0[0])
        crval[0] = ranew/d2r
        crval[1] = decnew/d2r
        # compute angle of rotation
        # the 2 terms are the rotation from omega and from the shift of the north
        # vector at the new reference position
        rot = math.atan(dotproduct(p0, omega)) + math.asin(math.sin(decnew)*math.sin(ra0-ranew))
        cth = math.cos(rot)
        sth = math.sin(rot)
        cd = cdmatrix[:]
        cdmatrix[0] = cth*cd[0] - sth*cd[2]
        cdmatrix[1] = cth*cd[1] - sth*cd[3]
        cdmatrix[2] = sth*cd[0] + cth*cd[2]
        cdmatrix[3] = sth*cd[1] + cth*cd[3]

    return (crval, cdmatrix)


def getdeltas(filename, crval, cdmatrix, omega=None):
    """Get the omega values for this dataset and return the shifts and rotation.

    Similar to applyomegaxyz but returns deltas instead of offsets

    Parameters
    ----------
    filename : str
        dataset filename

    crval : list of floats
        Right Ascension and Declination position at the reference pixel

    cdmatrix : list of floats
        Description of linear distortions: plate scale, rotation, and skew of the image

    omega : tuple, optional. Default value = None
        3-element tuple containing omega X, Y, and Z components.

    Returns
    -------
    dra : float
        delta ra in arcseconds

    ddec : float
        delta dec in arcseconds

    rot : float
        delta rotation in degrees
    """

    dx = dy = rot = 0.0
    if (not crval) or (not cdmatrix):
        return (dx, dy, rot)

    if not omega:
        # parse the filename and extract the visit identification
        dataset = getdataset(filename)
        if not dataset:
            return (dx, dy, rot)
        omega = getomegaxyz(dataset)

    if omega[0] != 0 or omega[1] != 0 or omega[2] != 0:
        # apply rotation to the tangent point
        ra0 = math.radians(crval[0])
        dec0 = math.radians(crval[1])
        cdec0 = math.cos(dec0)
        p0 = [cdec0*math.cos(ra0), cdec0*math.sin(ra0), math.sin(dec0)]
        dp0 = crossproduct(omega, p0)
        decnew = math.asin(p0[2]+dp0[2])
        ranew = math.atan2(p0[1]+dp0[1], p0[0]+dp0[0])
        dx = math.degrees((ranew - ra0))
        if dx > 180:
            dx = dx-360
        elif dx < -180:
            dx = dx+360
        dx = dx*math.cos(dec0) * 3600
        dy = math.degrees(((decnew-dec0)*3600))
        # compute angle of rotation
        # the 2 terms are the rotation from omega and from the shift of the north
        # vector at the new reference position
        rot = math.degrees((math.atan(dotproduct(p0, omega)) + math.asin(math.sin(decnew)*math.sin(ra0-ranew))))

    return (dx, dy, rot)


def updatefits(infile, outfile=None, dataset=None, omega=None, wcsname='HLA_HSC', verbose=False, overwrite=False):
    """Read input FITS file infile, update astrometry keyword, write outfile

    Default is to assume infile gives the dataset name.

    Parameters
    ----------
    infile : str
        fits image filename

    outfile : str, optional
        optional output filename. If outfile is omitted (or None), the data is not written but the function returns the
        updated astropy.io.fits object. Default value = None

    dataset : str, optional
        dataset filename. If it is specified then it is used as the dataset name for the wcs query. Default value = None

    omega : tuple, optional. Default value = None
        3-element tuple containing omega X, Y, and Z components.

    wcsname : str, optional
        WCS name. Default value = 'HLA_HSC'

    verbose : bool, optional
        display extra information? Default value = False

    overwrite : bool, optional
        if outfile is specified and a file with the same name already exists, overwrite? Default value = False

    Returns
    -------
    pout : astropy.io.fits object
        Updated fits data object
    """

    if not omega:
        if not dataset:
            # parse the filename and extract the visit identification
            dataset = getdataset(infile)
            if not dataset:
                raise ValueError("Unable to determine dataset for file")
        omega = getomegaxyz(dataset)
    if verbose:
        print("omega=", omega)
    nonzero = omega[0] != 0 or omega[1] != 0 or omega[2] != 0

    if verbose:
        print("reading", infile)
    pin = fits.open(infile)
    if outfile:
        if os.path.exists(outfile):
            if overwrite:
                if verbose:
                    print("Replacing existing file", outfile)
                os.remove(outfile)
            else:
                raise ValueError("Output file {} exists; specify overwrite=True to replace it".format(outfile))
        if verbose:
            print("creating", outfile)
        pout = fits.open(outfile, mode='append')
    else:
        if verbose:
            print("creating HDUlist for output")
        pout = fits.HDUList()
    for i, hdu in enumerate(pin):
        if nonzero:
            # update WCS keywords if present
            try:
                # skip update if the wcsname indicates correction has already been made
                oldwcsname = hdu.header.get("wcsname", "")
                if not oldwcsname.upper().find(wcsname.upper()) >= 0:
                    crval = [hdu.header['crval1'], hdu.header['crval2']]
                    try:
                        cdmatrix = [hdu.header['cd1_1'], hdu.header['cd1_2'], hdu.header['cd2_1'], hdu.header['cd2_2']]
                    except KeyError:
                        # try computing CD matrix from CROTA2 and CDELT
                        # delete those keywords if present and switch to CD matrix
                        cdelt1 = hdu.header['cdelt1']
                        cdelt2 = hdu.header['cdelt2']
                        try:
                            crota2 = hdu.header['crota2'] * math.pi / 180
                            del hdu.header['crota2']
                        except KeyError:
                            crota2 = 0.0
                        sinrot = math.sin(crota2)
                        cosrot = math.cos(crota2)
                        cdmatrix = [cdelt1*cosrot, -cdelt2*sinrot, cdelt1*sinrot, cdelt2*cosrot]
                        del hdu.header['cdelt1']
                        del hdu.header['cdelt2']
                    crval, cdmatrix = applyomegawcs(None, crval, cdmatrix, omega=omega)
                    if crval[0] < 0:
                        crval[0] = crval[0] + 360
                    elif crval[0] > 360:
                        crval[0] = crval[0] - 360
                    hdu.header['crval1'] = crval[0]
                    hdu.header['crval2'] = crval[1]
                    hdu.header['cd1_1'] = cdmatrix[0]
                    hdu.header['cd1_2'] = cdmatrix[1]
                    hdu.header['cd2_1'] = cdmatrix[2]
                    hdu.header['cd2_2'] = cdmatrix[3]
                    hdu.header['wcsname'] = wcsname
                    if verbose:
                        print("updated WCS in extension {} ({})".format(i, hdu.header.get('extname', 'primary')))
                else:
                    if verbose:
                        print("wcsname in extension {} = {}, no update".format(i, hdu.header['wcsname']))
            except KeyError:
                # OK, no WCS
                if verbose:
                    print("no WCS found in extension {} ({})".format(i, hdu.header.get('extname', 'primary')))
        pout.append(hdu)
        if outfile:
            pout.flush()
        # control memory usage
        hdu.data = None
        pout[i].data = None
    pin.close()
    if outfile:
        pout.close()
        if verbose:
            print("wrote", outfile)
    else:
        return pout


def applyomegacat(rain, decin, omega, radians=False):
    """Apply the 3-element infinitesimal rotation vector to a set of RA and Dec positions.
    
    Usage: raout, decout = applyomegacat(rain, decin, omega)

    Parameters
    ----------
    rain : numpy.ndarray
        Input right ascension positions in degrees (or radians)

    decin : numpy.ndarray
        Input declination positions in degrees (or radians)

    omega : numpy.ndarray
        3-element infinitesimal rotation vector

    radians : bool, optional
        If True, ra/dec values are in radians instead of degrees. Default value = False

    Returns
    -------
    raout : float
        RA output position in degrees (or radians if input argument radians = True)

    decout : float
        Dec output position in degrees (or radians if input argument radians = True)
    """
    xyz = radec2xyz(rain, decin, radians=radians)
    xyz += np.cross(omega, xyz)
    raout, decout = xyz2radec(xyz, radians=radians)
    return raout, decout


def radec2xyz(ra, dec, radians=False):

    """Convert RA, Dec to Cartesian (x, y, z) coordinates

    Usage: xyz = radec2xyz(ra, dec)
    
    Important Notes:
    
    - inputs *ra* and *dec* must match in shape

    Parameters
    ----------
    ra : numpy.ndarray
        Input right ascension positions in degrees

    dec : numpy.ndarray
        Input declination positions in degrees

    radians : bool, optional
        If True, ra/dec values are in radians instead of degrees. Default value = False

    Returns
    -------
    c : numpy.ndarray
        [\*,3] array with normalized cartesian coordinates
    """

    ra = np.asarray(ra)
    dec = np.asarray(dec)
    s = ra.shape
    if s != dec.shape:
        raise ValueError("ra, dec must be same-shape arrays")
    if not radians:
        dtor = np.pi/180
        ra = ra * dtor
        dec = dec * dtor
    c = np.empty(s + (3, ), dtype=float)
    cdec = np.cos(dec)
    c[:, 0] = np.cos(ra)*cdec
    c[:, 1] = np.sin(ra)*cdec
    c[:, 2] = np.sin(dec)
    return c


def xyz2radec(xyz, radians=False):

    """
    Convert Cartesian (x, y, z) coordinates to RA, Dec

    Usage: ra, dec = xyz2radec(xyz)

    Parameters
    ----------
    xyz : numpy.ndarray
        [\*,3] array with normalized cartesian coordinates. May be multi-dimensional but last dimension must be 3,
        e.g., shape = (10,10,3).

    radians : bool, optional
        If True, ra/dec values are in radians instead of degrees. Default value = False

    Returns
    -------
    ra : numpy.ndarray
        right ascension values with normalized cartesian coordinates in degrees

    dec : numpy.ndarray
        declination values with normalized cartesian coordinates in degrees
    """

    xyz = np.asarray(xyz)
    s = xyz.shape
    if s[-1] != 3:
        raise ValueError('xyz last dimension must be 3')
    # reshape to be a 2-D array [n,3]
    n = xyz.size//3
    c = np.reshape(xyz, (n, 3))
    # normalize to unity (for safety)
    norm = np.sqrt((c**2).sum(axis=-1))
    c = c/norm[:, None]

    dec = np.arcsin(c[:, 2])
    ra = np.arctan2(c[:, 1], c[:, 0])
    # force into range 0 to 2*pi
    w = np.where(ra < 0)
    ra[w] = ra[w] + 2*np.pi
    if not radians:
        # convert to degrees
        radeg = 180/np.pi
        ra = ra * radeg
        dec = dec * radeg
    # reshape using original dimensions of xyz
    ra = ra.reshape(s[:-1])
    dec = dec.reshape(s[:-1])
    return (ra, dec)


if __name__ == "__main__":
    dataset = 'HST_05397_2V_WFPC2_WFPC2'
    v = getomegaxyz(dataset)
    print(dataset, v)
    dataset = 'hst_12109_17_wfc3_ir'
    v = getomegaxyz(dataset)
    print(dataset, v)
    dataset = 'hst_10188_10_acs_wfc_f814w'
    v = getomegaxyz(dataset)
    print(dataset, v)
    v = getwcs(dataset)
    print(dataset, v)
    v = getwcs(dataset, applyomega=False)
    print(dataset, v)
    # filename = '/ifs/public/hst/hla/acs/V10.0/10188/10188_10/hst_10188_10_acs_wfc_f814w_drz.fits'
    # hdu = fits.open(filename)[1]
    # crval = [hdu.header['crval1'], hdu.header['crval2']]
    # cdmatrix = [hdu.header['cd1_1'], hdu.header['cd1_2'], hdu.header['cd2_1'], hdu.header['cd2_2']]
    # v = getdeltas(dataset, crval, cdmatrix)
    # print(dataset, v)
