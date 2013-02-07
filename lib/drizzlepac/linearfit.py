import numpy as np
from numpy import linalg as npla
from stsci.tools import logutil

"""
    ##################
    # DEVELOPMENT NOTE:
    #
    # This code needs to be refactored into a class for computing 
    #   and applying the fit. 
    #
    ##################
    
"""

# This is specifically NOT intended to match the package-wide version information.
__version__ = '0.3.1'
__vdate__ = '22-Dec-2005'


log = logutil.create_logger(__name__)

def RADTODEG(rad):
    return (rad * 180. / np.pi)

def DEGTORAD(deg):
    return (deg * np.pi / 180.)

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

    fit = fit_all(xy,uv,mode=mode,center=center,verbose=verbose)
    npts = xy.shape[0]
    npts0 = 0
    if nclip is None: nclip = 0
    # define index to initially include all points
    for n in range(nclip):
        if 'resids' in fit: 
            resids = fit['resids']
        else:
            resids = compute_resids(xy,uv,fit)

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
            fit = fit_all(xy,uv,mode=mode,center=center,verbose=False)
            del goodpix,goodx,goody
        else:
            break
    
    fit['img_coords'] = xy
    fit['ref_coords'] = uv
    fit['img_indx'] = xyindx
    fit['ref_indx'] = uvindx
    fit['img_orig_xy'] = xyorig
    fit['ref_orig_xy'] = uvorig
    
    return fit

def fit_all(xy,uv,mode='rscale',center=None,verbose=True):
    """ Performs an 'rscale' fit between matched lists of pixel positions xy and uv"""    
    if mode not in ['general','shift','rscale']: 
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
            print logstr
        else:
            log.info(logstr)
        result = fit_shifts(xy,uv)
    elif mode == 'general':
        logstr = 'Performing "general" fit'
        if verbose:
            print logstr
        else:
            log.info(logstr)
        # Set up products used for computing the fit
        gxy = uv    
        guv = xy
        Sx = gxy[:,0].sum()
        Sy = gxy[:,1].sum()
        Su = guv[:,0].sum()
        Sv = guv[:,1].sum()
        
        Sux = np.dot(guv[:,0],gxy[:,0])
        Svx = np.dot(guv[:,1],gxy[:,0])
        Suy = np.dot(guv[:,0],gxy[:,1])
        Svy = np.dot(guv[:,1],gxy[:,1])
        Sxx = np.dot(gxy[:,0],gxy[:,0])
        Syy = np.dot(gxy[:,1],gxy[:,1])
        Sxy = np.dot(gxy[:,0],gxy[:,1])
        Syx = np.dot(gxy[:,1],gxy[:,0])
        
        n = len(xy[:,0])
        M = np.array([[Sx, Sy, n], [Sxx, Sxy, Sx], [Syx, Syy, Sy]])
        U = np.array([Su,Sux,Suy])
        V = np.array([Sv,Svx,Svy])
        
        # The fit solutioN...
        # where 
        #   u = P0 + P1*x + P2*y
        #   v = Q0 + Q1*x + Q2*y
        #
        P = np.dot(npla.inv(M),U)
        Q = np.dot(npla.inv(M),V)
        
        # Return the shift, rotation, and scale changes
        result = build_fit(P,Q) 
        result['fit_matrix'] = np.array([[P[0],Q[0]],[P[1],Q[1]]])
        resids = guv - np.dot((gxy),result['fit_matrix']) - result['offset']
        rms = [resids[:,0].std(),resids[:,1].std()]
        result['rms'] = rms
        result['resids'] = resids
    else:
        logstr = 'Performing "rscale" fit'
        if verbose:
            print logstr
        else:
            log.info(logstr)
        result = geomap_rscale(xy,uv,center=center)

    return result
def fit_shifts(xy,uv):
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
    Pcoeffs = np.array([1.0,0.0,diff_pts[:,0].mean()])
    Qcoeffs = np.array([0.0,1.0,diff_pts[:,1].mean()])
    
    fit = build_fit(Pcoeffs,Qcoeffs)
    resids = diff_pts - fit['offset']
    rms = [resids[:,0].std(),resids[:,1].std()]
    fit['resids'] = resids
    fit['rms'] = rms
    fit['fit_matrix'] = np.array([[1.0,0.0],[0.0,1.0]])

    return fit

def fit_arrays(uv,xy):
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
    
    Sux = np.dot(uv[:,0],xy[:,0])
    Svx = np.dot(uv[:,1],xy[:,0])
    Suy = np.dot(uv[:,0],xy[:,1])
    Svy = np.dot(uv[:,1],xy[:,1])
    Sxx = np.dot(xy[:,0],xy[:,0])
    Syy = np.dot(xy[:,1],xy[:,1])
    Sxy = np.dot(xy[:,0],xy[:,1])
    
    n = len(xy[:,0])
    M = np.array([[Sx, Sy, n], [Sxx, Sxy, Sx], [Sxy, Syy, Sy]])
    U = np.array([Su,Sux,Suy])
    V = np.array([Sv,Svx,Svy])
    
    # The fit solutioN...
    # where 
    #   u = P0 + P1*x + P2*y
    #   v = Q0 + Q1*x + Q2*y
    #
    P = np.dot(npla.inv(M),U)
    Q = np.dot(npla.inv(M),V)
    #P = N.array([-0.434589, -0.893084, 285.420816])
    #Q = N.array([0.907435, -0.433864, 45.553862])
    
    # Return the shift, rotation, and scale changes
    return build_fit(P,Q)
    
def build_fit(P,Q):

    # Extract the results from P and Q
    det = P[0]*Q[1] - P[1]*Q[0]
    if det > 0:
        p = 1
    else:
        p = -1
    
    theta = np.arctan2(P[1] - p*Q[0], p*P[0] + Q[1]) 
    theta_deg = RADTODEG(theta) % 360.0
    
    avg_scale = (((p*P[0]+Q[1])*np.cos(theta)) + ((P[1] - p*Q[0])*np.sin(theta)) )/2
    alpha = np.arcsin( (-p*P[0]*np.sin(theta)) - (p*Q[0]*np.cos(theta)))/(2*avg_scale)
    d = ( ((p*P[0] - Q[1])*np.cos(theta)) - ((P[1]+p*Q[0])*np.sin(theta)))/(2*np.cos(alpha))
    
    scale_x = avg_scale + d
    scale_y = avg_scale - d

    return {'offset':(P[2],Q[2]),'fit_matrix':np.array([[P[0],Q[0]],[P[1],Q[1]]]),'rot':theta_deg,'scale':(avg_scale,scale_x,scale_y),'coeffs':(P,Q)}

def apply_old_coeffs(xy,coeffs):
    """ Apply the offset/shift/rot values from a linear fit 
        to an array of x,y positions.
    """
    _theta = DEGTORAD(coeffs[1])
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
    print 'FIT coeffs: ',fit['coeffs']
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

    dx = xyref[:,0] 
    dy = xyref[:,1] 
    du = xyin[:,0] 
    dv = xyin[:,1]
    
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
        rot_denom = -Sxrxi + Syryi
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
        cthetax = -mag*ctheta
        sthetay = -mag*stheta
    else:
        cthetax = mag*ctheta
        sthetay = mag*stheta
    sthetax = mag*stheta
    cthetay = mag*ctheta
    
    xshift = xi0 - (xr0*cthetax + yr0*sthetax)
    yshift = yi0 - (-xr0*sthetay + yr0*cthetay)
    rotmat = np.array([[cthetax,-sthetax],[sthetay,cthetay]])
    
    P = np.array([cthetax,-sthetax,xshift])
    Q = np.array([sthetay,cthetay,yshift])
    resids = xyin - np.dot((xyref),rotmat) - [xshift,yshift]
    rms = [resids[:,0].std(),resids[:,1].std()]
    
    rscale_fit = {'offset':(xshift,yshift),'rot':theta,'scale':(mag,mag,mag),'fit_matrix':rotmat,'coeffs':(P,Q),'resids':resids,'rms':rms}
    return rscale_fit
    
