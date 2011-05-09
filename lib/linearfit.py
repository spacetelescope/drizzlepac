import numpy as np
from numpy import linalg as npla

"""
    ##################
    # DEVELOPMENT NOTE:
    #
    # This code needs to be refactored into a class for computing 
    #   and applying the fit. 
    #
    ##################
    
"""

__version__ = '0.3.1 (22-Dec-2005)'

def RADTODEG(rad):
    return (rad * 180. / np.pi)

def DEGTORAD(deg):
    return (deg * np.pi / 180.)

def iter_fits_shifts(xy,uv,nclip=1,sigma=3.0):
    """ Perform an iterative-fit with 'nclip' iterations 
    """
    pass
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

    return build_fit(Pcoeffs,Qcoeffs)

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

    return {'offset':(P[2],Q[2]),'rot':theta_deg,'scale':(avg_scale,scale_x,scale_y),'coeffs':(P,Q)}

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
    xn,yn = apply_fit(uv,fit['coeffs'])
    resids = xy - np.transpose([xn,yn])
    return resids