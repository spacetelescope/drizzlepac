"""
    ##################
    # DEVELOPMENT NOTE:
    #
    # This code needs to be refactored into a class for computing
    #   and applying the fit.
    #
    ##################

:Authors: Warren Hack, Mihai Cara

:License: :doc:`LICENSE`

"""
import numpy as np
from numpy import linalg as npla
from stsci.tools import logutil


# This is specifically NOT intended to match the package-wide version information.
__version__ = '0.4.0'
__version_date__ = '10-Oct-2014'


log = logutil.create_logger(__name__, level=logutil.logging.NOTSET)

if hasattr(np, 'float128'):
    ndfloat128 = np.float128
elif hasattr(np, 'float96'):
    ndfloat128 = np.float96
else:
    ndfloat128 = np.float64


def iter_fit_shifts(xy,uv,nclip=3,sigma=3.0):
    """ Perform an iterative-fit with 'nclip' iterations
    """
    fit = fit_shifts(xy,uv)
    if nclip is None: nclip = 0
    # define index to initially include all points
    for n in range(nclip):
        resids = compute_resids(xy,uv,fit)
        resids1d = np.sqrt(np.power(resids[:,0],2)+np.power(resids[:,1],2))
        sig = resids1d.std()
        # redefine what pixels will be included in next iteration
        goodpix = resids1d < sigma*sig
        xy = xy[goodpix]
        uv = uv[goodpix]
        fit = fit_shifts(xy,uv)

    fit['img_coords'] = xy
    fit['ref_coords'] = uv

    return fit


def iter_fit_arrays(xy,uv,nclip=3,sigma=3.0):
    """ Perform an iterative-fit with 'nclip' iterations
    """
    fit = fit_arrays(xy,uv)

    if nclip is None: nclip = 0
    # define index to initially include all points
    for n in range(nclip):
        resids = compute_resids(xy,uv,fit)
        resids1d = np.sqrt(np.power(resids[:,0],2)+np.power(resids[:,1],2))
        sig = resids1d.std()
        # redefine what pixels will be included in next iteration
        goodpix = resids1d < sigma*sig
        xy = xy[goodpix]
        uv = uv[goodpix]
        fit = fit_arrays(xy,uv)

    fit['img_coords'] = xy
    fit['ref_coords'] = uv
    return fit


def iter_fit_all(xy,uv,xyindx,uvindx,
                    xyorig=None,uvorig=None,
                    mode='rscale',nclip=3,sigma=3.0,minobj=3,
                    center=None,verbose=False):

    if not isinstance(xy,np.ndarray):
        # cast input list as numpy ndarray for fitting
        xy = np.array(xy)
    if not isinstance(uv,np.ndarray):
        # cast input list as numpy ndarray for fitting
        uv = np.array(uv)

    if xy.shape[0] < nclip:
        log.warning('The number of sources for the fit < number of clipping iterations.',
            '    Resetting number of clipping iterations to 0.')
        nclip=0

    if center is None:
        xcen = uv[:,0].mean()
        ycen = uv[:,1].mean()
        center = [xcen,ycen]
    xy -= center
    uv -= center

    fit = fit_all(xy, uv, mode=mode, center=center, verbose=verbose)
    npts = xy.shape[0]
    npts0 = 0
    if nclip is None: nclip = 0
    # define index to initially include all points
    for n in range(nclip):
        if 'resids' in fit:
            resids = fit['resids']
        else:
            resids = compute_resids(xy, uv, fit)

        # redefine what pixels will be included in next iteration
        whtfrac = npts/(npts-npts0-1.0)
        cutx = sigma*(fit['rms'][0]*whtfrac)
        cuty = sigma*(fit['rms'][1]*whtfrac)

        goodx = (np.abs(resids[:,0]) < cutx)
        goody = (np.abs(resids[:,1]) < cuty)
        goodpix = np.bitwise_and(goodx,goody)

        if np.where(goodpix == True)[0].shape[0] > 2:
            npts0 = npts - goodpix.shape[0]
            xy = xy[goodpix]
            uv = uv[goodpix]
            xyindx = xyindx[goodpix]
            uvindx = uvindx[goodpix]
            if xyorig is not None:
                xyorig = xyorig[goodpix]
            if uvorig is not None:
                uvorig = uvorig[goodpix]
            fit = fit_all(xy, uv, mode=mode, center=center, verbose=False)
            del goodpix,goodx,goody
        else:
            break

    fit['img_coords'] = xy
    fit['ref_coords'] = uv
    fit['img_indx'] = xyindx
    fit['ref_indx'] = uvindx
    fit['img_orig_xy'] = xyorig
    fit['ref_orig_xy'] = uvorig
    fit['fit_xy'] = np.dot(xy - fit['offset'],
                           np.linalg.inv(fit['fit_matrix'])) + center

    return fit


def fit_all(xy,uv,mode='rscale',center=None,verbose=True):
    """ Performs an 'rscale' fit between matched lists of pixel positions xy and uv"""
    if mode not in ['general', 'shift', 'rscale']:
        mode = 'rscale'
    if not isinstance(xy,np.ndarray):
        # cast input list as numpy ndarray for fitting
        xy = np.array(xy)
    if not isinstance(uv,np.ndarray):
        # cast input list as numpy ndarray for fitting
        uv = np.array(uv)

    if mode == 'shift':
        logstr = 'Performing "shift" fit'
        if verbose:
            print(logstr)
        else:
            log.info(logstr)
        result = fit_shifts(xy, uv)

    elif mode == 'general':
        logstr = 'Performing "general" fit'
        if verbose:
            print(logstr)
        else:
            log.info(logstr)
        result = fit_general(xy, uv)

    else:
        logstr = 'Performing "rscale" fit'
        if verbose:
            print(logstr)
        else:
            log.info(logstr)
        result = geomap_rscale(xy, uv, center=center)

    return result


def fit_shifts(xy, uv):
    """ Performs a simple fit for the shift only between
        matched lists of positions 'xy' and 'uv'.

        Output: (same as for fit_arrays)
        =================================
        DEVELOPMENT NOTE:
            Checks need to be put in place to verify that
            enough objects are available for a fit.
        =================================
    """
    diff_pts = xy - uv
    Pcoeffs = np.array([1.0,0.0,diff_pts[:,0].mean(dtype=np.float64)])
    Qcoeffs = np.array([0.0,1.0,diff_pts[:,1].mean(dtype=np.float64)])

    fit = build_fit(Pcoeffs, Qcoeffs, 'shift')
    resids = diff_pts - fit['offset']
    fit['resids'] = resids
    fit['rms'] = resids.std(axis=0)
    fit['rmse'] = float(np.sqrt(np.mean(2 * resids**2)))
    fit['mae'] = float(np.mean(np.linalg.norm(resids, axis=1)))

    return fit


def fit_general(xy, uv):
    """ Performs a simple fit for the shift only between
        matched lists of positions 'xy' and 'uv'.

        Output: (same as for fit_arrays)
        =================================
        DEVELOPMENT NOTE:
            Checks need to be put in place to verify that
            enough objects are available for a fit.
        =================================
    """
    # Set up products used for computing the fit
    gxy = uv.astype(ndfloat128)
    guv = xy.astype(ndfloat128)
    Sx = gxy[:,0].sum()
    Sy = gxy[:,1].sum()
    Su = guv[:,0].sum()
    Sv = guv[:,1].sum()

    Sux = np.dot(guv[:,0], gxy[:,0])
    Svx = np.dot(guv[:,1], gxy[:,0])
    Suy = np.dot(guv[:,0], gxy[:,1])
    Svy = np.dot(guv[:,1], gxy[:,1])
    Sxx = np.dot(gxy[:,0], gxy[:,0])
    Syy = np.dot(gxy[:,1], gxy[:,1])
    Sxy = np.dot(gxy[:,0], gxy[:,1])

    n = len(xy[:,0])
    M = np.array([[Sx, Sy, n], [Sxx, Sxy, Sx], [Sxy, Syy, Sy]])
    U = np.array([Su, Sux, Suy])
    V = np.array([Sv, Svx, Svy])

    # The fit solutioN...
    # where
    #   u = P0 + P1*x + P2*y
    #   v = Q0 + Q1*x + Q2*y
    #
    try:
        invM = np.linalg.inv(M.astype(np.float64))
    except np.linalg.LinAlgError:
        raise SingularMatrixError(
            "Singular matrix: suspected colinear points."
        )
    P = np.dot(invM, U).astype(np.float64)
    Q = np.dot(invM, V).astype(np.float64)
    if not (np.all(np.isfinite(P)) and np.all(np.isfinite(Q))):
        raise ArithmeticError('Singular matrix.')

    # Return the shift, rotation, and scale changes
    result = build_fit(P, Q, 'general')
    resids = xy - np.dot(uv, result['fit_matrix']) - result['offset']
    result['rms'] = resids.std(axis=0)
    result['resids'] = resids
    result['rmse'] = float(np.sqrt(np.mean(2 * resids**2)))
    result['mae'] = float(np.mean(np.linalg.norm(resids, axis=1)))

    return result


def fit_arrays(uv, xy):
    """ Performs a generalized fit between matched lists of positions
        given by the 2 column arrays xy and uv.

        This function fits for translation, rotation, and scale changes
        between 'xy' and 'uv', allowing for different scales and
        orientations for X and Y axes.

        =================================
        DEVELOPMENT NOTE:
            Checks need to be put in place to verify that
            enough objects are available for a fit.
        =================================

        Output:
           (Xo,Yo),Rot,(Scale,Sx,Sy)
           where
                Xo,Yo:  offset,
                Rot:    rotation,
                Scale:  average scale change, and
                Sx,Sy:  scale changes in X and Y separately.

        Algorithm and nomenclature provided by: Colin Cox (11 Nov 2004)
    """

    if not isinstance(xy,np.ndarray):
        # cast input list as numpy ndarray for fitting
        xy = np.array(xy)
    if not isinstance(uv,np.ndarray):
        # cast input list as numpy ndarray for fitting
        uv = np.array(uv)

    # Set up products used for computing the fit
    Sx = xy[:,0].sum()
    Sy = xy[:,1].sum()
    Su = uv[:,0].sum()
    Sv = uv[:,1].sum()

    Sux = np.dot(uv[:,0], xy[:,0])
    Svx = np.dot(uv[:,1], xy[:,0])
    Suy = np.dot(uv[:,0], xy[:,1])
    Svy = np.dot(uv[:,1], xy[:,1])
    Sxx = np.dot(xy[:,0], xy[:,0])
    Syy = np.dot(xy[:,1], xy[:,1])
    Sxy = np.dot(xy[:,0], xy[:,1])

    n = len(xy[:,0])
    M = np.array([[Sx, Sy, n], [Sxx, Sxy, Sx], [Sxy, Syy, Sy]])
    U = np.array([Su, Sux, Suy])
    V = np.array([Sv, Svx, Svy])

    # The fit solutioN...
    # where
    #   u = P0 + P1*x + P2*y
    #   v = Q0 + Q1*x + Q2*y
    #
    try:
        invM = np.linalg.inv(M.astype(np.float64))
    except np.linalg.LinAlgError:
        raise SingularMatrixError(
            "Singular matrix: suspected colinear points."
        )
    P = np.dot(invM, U).astype(np.float64)
    Q = np.dot(invM, V).astype(np.float64)
    if not (np.all(np.isfinite(P)) and np.all(np.isfinite(Q))):
        raise ArithmeticError('Singular matrix.')

    # Return the shift, rotation, and scale changes
    return build_fit(P, Q, 'general')


def build_fit(P, Q, fitgeom):
    # Build fit matrix:
    fit = np.dstack((P[:2],Q[:2]))[0]

    # determinant of the transformation
    det = P[0]*Q[1] - P[1]*Q[0]
    sdet = np.sign(det)
    proper = sdet >= 0

    # Create a working copy (no reflections) for computing transformation
    # parameters (scale, rotation angle, skew):
    wfit = fit.copy()

    # Default skew:
    skew = 0.0

    if fitgeom == 'shift':
        return { 'offset' : (P[2], Q[2]),
                 'fit_matrix' : fit,
                 'rot' : 0.0,
                 'rotxy': (0.0, 0.0, 0.0, skew),
                 'scale' : (1.0, 1.0, 1.0),
                 'coeffs' : (P, Q),
                 'skew' : skew,
                 'proper' : proper,
                 'fitgeom' : fitgeom }

    # Compute average scale:
    s = np.sqrt(np.abs(det))
    # Compute scales for each axis:
    if fitgeom == 'general':
        sx = np.sqrt(P[0]**2 + Q[0]**2)
        sy = np.sqrt(P[1]**2 + Q[1]**2)
    else:
        sx = s
        sy = s

    # Remove scale from the transformation matrix:
    wfit[0,:] /= sx
    wfit[1,:] /= sy

    # Compute rotation angle as if we have a proper rotation.
    # This will also act as *some sort* of "average rotation" even for
    # transformations with different rot_x and rot_y:
    prop_rot = np.rad2deg(np.arctan2(wfit[1,0] - sdet*wfit[0,1],
                                     wfit[0,0] + sdet*wfit[1,1])) % 360.0

    if proper and fitgeom == 'rscale':
        rotx = prop_rot
        roty = prop_rot
        rot = prop_rot
        skew = 0.0
    else:
        rotx = np.rad2deg(np.arctan2(-wfit[0,1], wfit[0,0])) % 360.0
        roty = np.rad2deg(np.arctan2(wfit[1,0], wfit[1,1])) % 360.0
        rot = 0.5 * (rotx + roty)
        skew = roty - rotx

    return { 'offset' : (P[2], Q[2]),
             'fit_matrix' : fit,
             'rot' : prop_rot,
             'rotxy': (rotx, roty, rot, skew),
             'scale' : (s, sx, sy),
             'coeffs' : (P, Q),
             'skew' : skew,
             'proper' : proper,
             'fitgeom' : fitgeom }


def buildFitMatrix(rot, scale=1):
    if hasattr(rot, '__iter__'):
        rx = rot[0]
        ry = rot[1]
    else:
        rx = float(rot)
        ry = rx
    if hasattr(scale, '__iter__'):
        sx = scale[0]
        sy = scale[1]
    else:
        sx = float(scale)
        sy = sx
    m = np.array(
        [
            [ sx*np.cos(np.deg2rad(rx)), -sx*np.sin(np.deg2rad(rx)) ],
            [ sy*np.sin(np.deg2rad(ry)),  sy*np.cos(np.deg2rad(ry)) ]
        ]
    )
    return m


def apply_old_coeffs(xy,coeffs):
    """ Apply the offset/shift/rot values from a linear fit
        to an array of x,y positions.
    """
    _theta = np.deg2rad(coeffs[1])
    _mrot = np.zeros(shape=(2,2),dtype=np.float64)
    _mrot[0] = (np.cos(_theta),np.sin(_theta))
    _mrot[1] = (-np.sin(_theta),np.cos(_theta))

    new_pos = (np.dot(xy,_mrot)/coeffs[2][0]) + coeffs[0]

    return new_pos


def apply_fit(xy,coeffs):
    """ Apply the coefficients from a linear fit to
        an array of x,y positions.

        The coeffs come from the 'coeffs' member of the
        'fit_arrays()' output.
    """
    x_new = coeffs[0][2] + coeffs[0][0]*xy[:,0] + coeffs[0][1]*xy[:,1]
    y_new = coeffs[1][2] + coeffs[1][0]*xy[:,0] + coeffs[1][1]*xy[:,1]

    return x_new,y_new


def compute_resids(xy,uv,fit):
    """ Compute the residuals based on fit and input arrays to the fit
    """
    print('FIT coeffs: ',fit['coeffs'])
    xn,yn = apply_fit(uv,fit['coeffs'])
    resids = xy - np.transpose([xn,yn])
    return resids


##### My interpretation of geomap 'rscale' fitting based on 'lib/geofit.x'
def geomap_rscale(xyin,xyref,center=None):
    """
    Set up the products used for computing the fit derived using the code from
    lib/geofit.x for the function 'geo_fmagnify()'. Comparisons with results from
    geomap (no additional clipping) were made and produced the same results
    out to 5 decimal places.

    Output
    ------
    fit: dict
        Dictionary containing full solution for fit.
    """
    if center is not None:
        xcen = center[0]
        ycen = center[1]
    else:
        xcen = xyref[:,0].mean()
        ycen = xyref[:,1].mean()

    dx = xyref[:,0].astype(ndfloat128)
    dy = xyref[:,1].astype(ndfloat128)
    du = xyin[:,0].astype(ndfloat128)
    dv = xyin[:,1].astype(ndfloat128)

    n = xyref.shape[0]
    Sx = dx.sum()
    Sy = dy.sum()
    Su = du.sum()
    Sv = dv.sum()
    xr0 = Sx/n
    yr0 = Sy/n
    xi0 = Su/n
    yi0 = Sv/n
    Sxrxr = np.power((dx-xr0),2).sum()
    Syryr = np.power((dy-yr0),2).sum()
    Syrxi = ((dy-yr0)*(du-xi0)).sum()
    Sxryi = ((dx-xr0)*(dv-yi0)).sum()
    Sxrxi = ((dx-xr0)*(du-xi0)).sum()
    Syryi = ((dy-yr0)*(dv-yi0)).sum()

    rot_num = Sxrxi * Syryi
    rot_denom = Syrxi * Sxryi
    if rot_num == rot_denom: det = 0.0
    else: det = rot_num - rot_denom
    if (det < 0):
        rot_num = Syrxi + Sxryi
        rot_denom = Sxrxi - Syryi
    else:
        rot_num = Syrxi - Sxryi
        rot_denom = Sxrxi + Syryi
    if rot_num == rot_denom: theta = 0.0
    else:
        theta = np.rad2deg(np.arctan2(rot_num,rot_denom))
        if theta < 0:
            theta += 360.0

    ctheta = np.cos(np.deg2rad(theta))
    stheta = np.sin(np.deg2rad(theta))
    s_num = rot_denom*ctheta + rot_num*stheta
    s_denom = Sxrxr + Syryr
    if s_denom < 0:
        mag = 1.0
    else:
        mag = s_num/s_denom

    if det < 0:
        # "flip" y-axis (reflection about x-axis *after* rotation)
        # NOTE: keep in mind that 'fit_matrix'
        #       is the transposed rotation matrix.
        sthetax = -mag*stheta
        cthetay = -mag*ctheta
    else:
        sthetax = mag*stheta
        cthetay = mag*ctheta
    cthetax = mag*ctheta
    sthetay = mag*stheta

    sdet = np.sign(det)
    xshift = (xi0 - (xr0*cthetax + sdet*yr0*sthetax)).astype(np.float64)
    yshift = (yi0 - (-sdet*xr0*sthetay + yr0*cthetay)).astype(np.float64)

    P = np.array([ cthetax, sthetay, xshift],dtype=np.float64)
    Q = np.array([ -sthetax, cthetay, yshift],dtype=np.float64)

    # Return the shift, rotation, and scale changes
    result = build_fit(P, Q, fitgeom='rscale')
    resids = xyin - np.dot((xyref), result['fit_matrix']) - result['offset']
    result['rms'] = resids.std(axis=0)
    result['resids'] = resids
    result['rmse'] = float(np.sqrt(np.mean(2 * resids**2)))
    result['mae'] = float(np.mean(np.linalg.norm(resids, axis=1)))

    return result
