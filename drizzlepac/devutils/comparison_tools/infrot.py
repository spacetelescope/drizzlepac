# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4 ai :
"""
Functions to compute and apply infinitesimal rotations to correct small shifts
and rotations in spherical coordinates

Uses formulae from `Budavari & Lubow (2012, ApJ, 761, 188) <http://adsabs.harvard.edu/abs/2012ApJ...761..188B>`_

Created by R. White from IDL versions, 2017 January 12

Added getxyshiftrot function for X,Y pixel coords, 2017 June 15

Path
----
HLApipeline/regression_testing/infrot.py

Dependencies
------------
None.

Inputs
------
None.

Classes and Functions
---------------------
"""

import numpy as np

def getxyshiftrot(x1,y1,x2,y2,xref=None,yref=None):

    """
    Compute xshift, yshift, rotation to shift/rotate positions in
    first set of arrays to match those in second set.  Positions are in pixels for
    an assumed small field.

    Calls getinfrot to do the fit.
    
    Usage: xshift, yshift, rotation = getxyshiftrot(x1,y1,x2,y2)

    Important notes:
    
    - Input arrays x1, y1, x2, and y2 positions are in pixels.
    - Position lists must be matched with the same number of entries.
    
    :param x1: Input array #1, X components.
    :param y1: Input array #1, Y components.
    :param x2: Input array #2, X components.
    :param y2: Input array #2, Y components.
    :param xref: Position of X reference pixel. If not specified, default value = mean(x1)
    :param yref: Position of y reference pixel. If not specified, default value = mean(y1)
    :type x1: numpy.ndarray
    :type y1: numpy.ndarray
    :type x2: numpy.ndarray
    :type y2: numpy.ndarray
    :type xref: float
    :type yref: float
    :return: X and Y shifts in pixels and rotation angle in degrees
    """

    x1 = np.asarray(x1)
    y1 = np.asarray(y1)
    x2 = np.asarray(x2)
    y2 = np.asarray(y2)

    if x1.ndim != 1 or y1.ndim != 1 or x2.ndim != 1 or y2.ndim != 1:
        raise ValueError("x1,y1,x2,y2 must be 1-D arrays")
    n = x1.size
    if y1.size != n or x2.size != n or y2.size != n:
        raise ValueError("x1,y1,x2,y2 must be matching-length 1-D arrays")

    # scale values to fake RA & Dec with range 0.05 degrees (approximate HST camera field)
    if xref is None: xref = x1.mean()
    if yref is None: yref = y1.mean()
    pxsize = x1.max() - x1.min()
    pysize = y1.max() - y1.min()
    pscale = 0.05 / max(pxsize,pysize)
    ra1 = (xref - x1) * pscale
    ra2 = (xref - x2) * pscale
    dec1 = (y1 - yref) * pscale
    dec2 = (y2 - yref) * pscale
    omega = getinfrot(ra1,dec1,ra2,dec2)
    # convert to degrees
    omega = omega * (180/np.pi)
    xshift = -omega[2]/pscale
    yshift = -omega[1]/pscale
    rotation = omega[0]
    return (xshift, yshift, rotation)


def getinfrot(ra1, dec1, ra2, dec2, radians=False):

    """
    Compute 3-element infinitesimal rotation vector to shift/rotate positions in
    first set of arrays to match those in second set.
    
    Uses formulae from `Budavari & Lubow (2012, ApJ, 761, 188) <http://adsabs.harvard.edu/abs/2012ApJ...761..188B>`_
    
    Usage: omega = getinfrot(ra1,dec1,ra2,dec2)
    
    Important Notes:
    
    - Input arrays ra1, dec1, ra2, and dec2 positions are in degrees
    - Position lists must be matched with the same number of entries.
    
    :param ra1: Input array #1, right ascension components. 
    :param dec1: Input array #1, declination components. 
    :param ra2: Input array #2, right ascension components. 
    :param dec2: Input array #2, declination components. 
    :param radians: If True, ra/dec values are in radians instead of degrees. Default value = False
    :type ra1: numpy.ndarray
    :type dec1: numpy.ndarray 
    :type ra2: numpy.ndarray
    :type dec2: numpy.ndarray
    :type radians: Boolean
    :return: 3-element infinitesimal rotation vector omega such that r1(shifted) = r1 + cross(omega,r1)
    """

    ra1 = np.asarray(ra1)
    dec1 = np.asarray(dec1)
    ra2 = np.asarray(ra2)
    dec2 = np.asarray(dec2)
    if ra1.ndim != 1 or dec1.ndim != 1 or ra2.ndim != 1 or dec2.ndim != 1:
        raise ValueError("ra1,dec1,ra2,dec2 must be 1-D arrays")
    n = ra1.size
    if dec1.size != n or ra2.size != n or dec2.size != n:
        raise ValueError("ra1,dec1,ra2,dec2 must be matching-length 1-D arrays")

    c = radec2xyz(ra2,dec2,radians=radians)
    r = radec2xyz(ra1,dec1,radians=radians)

    b = np.cross(r,c)
    a = -vdyadicp(r,r)
    if b.ndim == 2:
        b = b.sum(axis=0)
        a = a.sum(axis=0)
    np.fill_diagonal(a, a.diagonal() + n)
    return np.dot(np.linalg.inv(a), b)


def applyomega(rain, decin, omega, radians=False):

    """
    Apply the 3-element infinitesimal rotation vector computed by getinfrot()
    to a set of RA and Dec positions.
    
    Usage: raout, decout = applyomega(rain, decin, omega)
   
    :param rain: Input right ascension positions in degrees (or radians)
    :param decin: Input declination positions in degrees (or radians)
    :param omega: 3-element infinitesimal rotation vector
    :param radians: If True, ra/dec values are in radians instead of degrees. Default value = False
    :type rain: numpy.ndarray
    :type decin: numpy.ndarray
    :type omega: numpy.ndarray
    :type radians: Boolean
    :return: the RA and Dec output positions in degrees (or radians)
    """

    xyz = radec2xyz(rain,decin,radians=radians)
    xyz += np.cross(omega,xyz)
    raout, decout = xyz2radec(xyz, radians=radians)
    return raout, decout

def radec2xyz(ra, dec, radians=False):

    """
    Convert RA, Dec to Cartesian (x, y, z) coordinates

    Usage: xyz = radec2xyz(ra, dec)
    
    Important Notes:
    
    - inputs *ra* and *dec* must match in shape
    
    :param ra: Input right ascension positions in degrees
    :param dec: Input declination positions in degrees 
    :param radians: If True, ra/dec values are in radians instead of degrees. Default value = False
    :type ra: numpy.ndarray 
    :type dec: numpy.ndarray 
    :type radians: Boolean
    :return: [\*,3] array with normalized cartesian coordinates
    """

    ra = np.asarray(ra)
    dec = np.asarray(dec)
    s = ra.shape
    if s != dec.shape:
        raise ValueError("ra,dec must be same-shape arrays")
    if not radians:
        dtor =  np.pi/180
        ra = ra * dtor
        dec = dec * dtor
    c = np.empty(s + (3,), dtype=float)
    cdec = np.cos(dec)
    c[:,0] = np.cos(ra)*cdec
    c[:,1] = np.sin(ra)*cdec
    c[:,2] = np.sin(dec)
    return c


def xyz2radec(xyz, radians=False):

    """
    Convert Cartesian (x, y, z) coordinates to RA, Dec

    Usage: ra, dec = xyz2radec(xyz)
   
    :param xyz: [\*,3] array with normalized cartesian coordinates.  May be multi-dimensional but last dimension must be 3, e.g., shape = (10,10,3).
    :param radians: If True, ra/dec values are returned in radians instead of degrees. Default value = False
    :type xyz: numpy.ndarray
    :type radians: Boolean
    :return: ra and dec arrays with normalized cartesian coordinates in degrees
    """

    xyz = np.asarray(xyz)
    s = xyz.shape
    if s[-1] != 3:
        raise ValueError('xyz last dimension must be 3')
    # reshape to be a 2-D array [n,3]
    n = xyz.size/3
    c = np.reshape(xyz,(n,3))
    # normalize to unity (for safety)
    norm = np.sqrt((c**2).sum(axis=-1))
    c = c/norm[:,None]

    dec = np.arcsin(c[:,2])
    ra = np.arctan2(c[:,1],c[:,0])
    # force into range 0 to 2*pi
    w = np.where(ra < 0)
    ra[w] = ra[w] + 2*np.pi
    if not radians:
        # convert to degrees
        radeg =  180/np.pi
        ra = ra * radeg
        dec = dec * radeg
    # reshape using original dimensions of xyz
    ra = ra.reshape(s[:-1])
    dec = dec.reshape(s[:-1])
    return (ra, dec)

def vdyadicp(a, b):

    """
    Return dyadic product of arrays a, b
 
    :param a: [\*,3] input array.  May be multidimensional (e.g. 5,4,3) as long as last dimension is 3.
    :param b: [\*,3] input array.  May be multidimensional (e.g. 5,4,3) as long as last dimension is 3.
    :type a: numpy.ndarray
    :type b: numpy.ndarray
    :return: [\*, 3, 3] array of outer products
    """

    a = np.asarray(a)
    b = np.asarray(b)
    s = a.shape
    if s[-1] != 3:
        raise ValueError('a,b last dimension must be 3')
    if s != b.shape:
        raise ValueError('a,b must have same shape [...,3]')
    return np.einsum('...i,...j->...ij', a, b)


if __name__ == "__main__":
    print('='*40)
    ra1 = np.array([1.0, 1.1, 1.2, 1.3])
    dec1 = 45.0 + np.zeros(ra1.size,dtype=float)
    ra2 = ra1
    dec2 = dec1 + (ra1-ra1[0])*1.e-3
    omega = getinfrot(ra1, dec1, ra2, dec2)
    print('ra1  ',ra1)
    print('dec1 ',dec1)
    print('omega',omega)
    print('ra2  ',ra2)
    print('dec2 ',dec2)
    print('rotated values')
    ra1r, dec1r = applyomega(ra1,dec1,omega)
    print('ra1r ',ra1r)
    print('dec1r',dec1r)
    print('check')
    print('ra diff ', ra1r-ra2)
    print('dec diff', dec1r-dec2)

    print('='*40)
    ra1 = np.array([1.0, 1.1, 1.2, 1.3])
    dec1 = 45.0 + np.zeros(ra1.size,dtype=float)
    ra2 = ra1
    dec2 = dec1 + 1.e-3
    omega = getinfrot(ra1, dec1, ra2, dec2)
    print('ra1  ',ra1)
    print('dec1 ',dec1)
    print('omega',omega)
    print('ra2  ',ra2)
    print('dec2 ',dec2)
    print('rotated values')
    ra1r, dec1r = applyomega(ra1,dec1,omega)
    print('ra1r ',ra1r)
    print('dec1r',dec1r)
    print('ra diff ', ra1r-ra2)
    print('dec diff', dec1r-dec2)
