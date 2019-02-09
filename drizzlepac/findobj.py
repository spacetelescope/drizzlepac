"""
A suite of functions for finding sources in images.

:Authors: Warren Hack, Mihai Cara

:License: :doc:`LICENSE`

"""
import sys
import math

import numpy as np
from scipy import signal, ndimage

import stsci.imagestats as imagestats

from . import cdriz

__all__ = ['gaussian1', 'gausspars', 'gaussian', 'moments', 'errfunc',
           'findstars', 'apply_nsigma_separation', 'xy_round',
           'precompute_sharp_round', 'sharp_round', 'roundness', 'immoments',
           'nmoment', 'centroid', 'cmoment', 'central_moments', 'covmat',
           'help', 'getHelpAsString']


#def gaussian(amplitude, xcen, ycen, xsigma, ysigma):
#from numpy import *

FWHM2SIG = 2*np.sqrt(2*np.log(2))

#def gaussian1(height, x0, y0, fwhm, nsigma=1.5, ratio=1., theta=0.0):
def gaussian1(height, x0, y0, a, b, c):
    """
    height - the amplitude of the gaussian
    x0, y0, - center of the gaussian
    a, b, c - ellipse parameters (coefficients in the quadratic form)
    """

    return lambda x, y: height * np.exp(-0.5* (a*(x-x0)**2 + b*(x-x0)*(y-y0) + c*(y-y0)**2))


def gausspars(fwhm, nsigma=1.5, ratio=1, theta=0.):
    """
        height - the amplitude of the gaussian
        x0, y0, - center of the gaussian
        fwhm - full width at half maximum of the observation
        nsigma - cut the gaussian at nsigma
        ratio = ratio of xsigma/ysigma
        theta - angle of position angle of the major axis measured
        counter-clockwise from the x axis

        Returns dimensions nx and ny of the elliptical kernel as well as the
        ellipse parameters a, b, c, and f when defining an ellipse through the
        quadratic form: a*(x-x0)^2+b(x-x0)*(y-y0)+c*(y-y0)^2 <= 2*f
    """

    xsigma = fwhm / FWHM2SIG
    ysigma = ratio * xsigma

    f = nsigma**2/2.

    theta = np.deg2rad(theta)
    cost  = np.cos(theta)
    sint  = np.sin(theta)

    if ratio == 0: # 1D Gaussian

        if theta == 0 or theta == 180:
            a = 1/xsigma**2
            b = 0.0
            c = 0.0
        elif theta  == 90:
            a = 0.0
            b = 0.0
            c = 1/xsigma**2
        else:
            print('Unable to construct 1D Gaussian with these parameters\n')
            raise ValueError

        nx = 2 * int(max(2, (xsigma*nsigma*np.abs(cost))))+1
        ny = 2 * int(max(2, (xsigma*nsigma*np.abs(sint))))+1

    else: #2D gaussian

        xsigma2 = xsigma * xsigma
        ysigma2 = ysigma * ysigma

        a = cost**2/xsigma2 + sint**2/ysigma2
        b = 2 * cost * sint *(1.0/xsigma2-1.0/ysigma2)
        c = sint**2/xsigma2 + cost**2/ysigma2

        d = b**2 - 4*a*c # discriminant

        #        nx = int(2*max(2, math.sqrt(-8*c*f/d)))+1
        #        ny = int(2*max(2, math.sqrt(-8*a*f/d)))+1
        nx = 2 * int(2*max(1, nsigma*math.sqrt(-c/d)))+1
        ny = 2 * int(2*max(1, nsigma*math.sqrt(-a/d)))+1

    return nx, ny, a, b, c, f


def gaussian(height, center_x, center_y, width_x, width_y):
    #Returns a gaussian function with the given parameters
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)


def moments(data,cntr):
    """
    Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments.

    """
    total = data.sum()
    #X, Y = np.indices(data.shape)
    #x = (X*data).sum()/total
    #y = (Y*data).sum()/total
    x,y = cntr
    xi = int(x)
    yi = int(y)
    if xi < 0 or xi >= data.shape[1] or yi < 0 or yi >= data.shape[0]:
        raise ValueError
    col = data[:, xi]
    width_x = np.sqrt(abs(((np.arange(col.size)-y)**2*col).sum()/col.sum()))
    row = data[yi, :]
    width_y = np.sqrt(abs(((np.arange(row.size)-x)**2*row).sum()/row.sum()))
    height = data.max()
    return height, x, y, width_x, width_y

def errfunc(p, *args):

    func = gaussian1(*p)
    ret =np.ravel(func(*args[1:]) - args[0])
    return ret

def findstars(jdata, fwhm, threshold, skymode,
              peakmin=None, peakmax=None, fluxmin=None, fluxmax=None,
              nsigma=1.5, ratio=1.0, theta=0.0,
              use_sharp_round=False,mask=None,
              sharplo=0.2,sharphi=1.0,roundlo=-1.0,roundhi=1.0):

    # store input image size:
    (img_ny, img_nx) = jdata.shape

    # Define convolution inputs
    nx, ny, a, b, c, f = gausspars(fwhm, nsigma=nsigma, ratio= ratio, theta=theta)

    xc = nx//2
    yc = ny//2

    yin, xin = np.mgrid[0:ny, 0:nx]
    kernel = gaussian1(1.0, xc, yc, a, b, c)(xin,yin)

    # define size of extraction box for each source based on kernel size
    grx = xc
    gry = yc

    # DAOFIND STYLE KERNEL "SHAPE"
    rmat    = np.sqrt((xin-xc)**2 + (yin-yc)**2)
    rmatell = a*(xin-xc)**2 + b*(xin-xc)*(yin-yc) + c*(yin-yc)**2
    xyrmask = np.where((rmatell <= 2*f) | (rmat <= 2.001),1,0).astype(np.int16)
    # Previous *style* computation for kernel "shape":
    #xyrmask = np.where(rmat <= max(grx,gry),1,0).astype(np.int16)

    npts = xyrmask.sum()

    rmask = kernel*xyrmask
    denom = (rmask*rmask).sum() - rmask.sum()**2/npts
    nkern = (rmask - (rmask.sum()/npts))/denom # normalize kernel to preserve
                                               # fluxes for thresholds
    nkern *= xyrmask

    # initialize values used for getting source centers
    relerr = 1./((rmask**2).sum() - (rmask.sum()**2/xyrmask.sum()))

    xsigsq = (fwhm / FWHM2SIG)**2
    ysigsq = (ratio**2) * xsigsq

    # convolve image with gaussian kernel
    convdata = signal.convolve2d(jdata, nkern, boundary='symm', mode='same').astype(np.float32)

    # clip image to create regions around each source for segmentation
    if mask is None:
        tdata=np.where(convdata > threshold, convdata, 0)
    else:
        tdata=np.where((convdata > threshold) & mask, convdata, 0)

    # segment image and find sources
    s = ndimage.morphology.generate_binary_structure(2, 2)
    ldata, nobj = ndimage.label(tdata, structure=s)
    fobjects = ndimage.find_objects(ldata)

    fluxes = []
    fitind = []
    if nobj < 2:
        print('No objects found for this image. Please check value of "threshold".')
        return fitind,fluxes

    # determine center of each source, while removing spurious sources or
    # applying limits defined by the user
    ninit  = 0
    ninit2 = 0

    s2m, s4m = precompute_sharp_round(nx, ny, xc, yc)

    satur  = False # Default assumption if use_sharp_round=False
    sharp  = None
    round1 = None
    round2 = None

    for ss,n in zip(fobjects,range(len(fobjects))):
        ssx = ss[1].stop - ss[1].start
        ssy = ss[0].stop - ss[0].start
        if ssx >= tdata.shape[1]-1 or ssy >= tdata.shape[0]-1:
            continue

        yr0 = ss[0].start - gry
        yr1 = ss[0].stop  + gry + 1
        if yr0 <= 0 or yr1 >= img_ny: continue # ignore sources within ny//2 of edge

        xr0 = ss[1].start - grx
        xr1 = ss[1].stop  + grx + 1
        if xr0 <= 0 or xr1 >= img_nx: continue # ignore sources within nx//2 of edge

        ssnew = (slice(yr0,yr1),slice(xr0,xr1))
        region = tdata[ssnew]

        cntr = centroid(region)

        # Define region centered on max value in object (slice)
        # This region will be bounds-checked to insure that it only accesses
        # a valid section of the image (not off the edge)
        maxpos = (int(cntr[1]+0.5)+ssnew[0].start,int(cntr[0]+0.5)+ssnew[1].start)
        yr0 = maxpos[0] - gry
        yr1 = maxpos[0] + gry + 1
        if yr0 < 0 or yr1 > img_ny:
            continue
        xr0 = maxpos[1] - grx
        xr1 = maxpos[1] + grx + 1
        if xr0 < 0 or xr1 > img_nx:
            continue

        # Simple Centroid on the region from the input image
        jregion = jdata[yr0:yr1,xr0:xr1]
        src_flux = jregion.sum()
        src_peak = jregion.max()

        if (peakmax is not None and src_peak >= peakmax):
            continue
        if (peakmin is not None and src_peak <= peakmin):
            continue

        if fluxmin and src_flux <= fluxmin:
            continue
        if fluxmax and src_flux >= fluxmax:
            continue

        datamin = jregion.min()
        datamax = jregion.max()

        if use_sharp_round:
            # Compute sharpness and first estimate of roundness:
            dregion = convdata[yr0:yr1,xr0:xr1]
            satur, round1, sharp = \
                sharp_round(jregion, dregion, xyrmask, xc, yc,
                            s2m, s4m, nx, ny, datamin, datamax)
            # Filter sources:
            if sharp is None or (sharp < sharplo or sharp > sharphi):
                continue
            if round1 is None or (round1 < roundlo or round1 > roundhi):
                continue

        px, py, round2 = xy_round(jregion, grx, gry, skymode,
                                  kernel, xsigsq, ysigsq, datamin, datamax)

        # Filter sources:
        if px is None:
            continue

        if use_sharp_round and not satur and \
           (round2 is None or round2 < roundlo or round2 > roundhi):
            continue

        fitind.append((px + xr0, py + yr0, sharp, round1, round2))
        # compute a source flux value
        fluxes.append(src_flux)

    fitindc, fluxesc = apply_nsigma_separation(fitind, fluxes, fwhm*nsigma / 2)

    return fitindc, fluxesc


def apply_nsigma_separation(fitind,fluxes,separation,niter=10):
    """
    Remove sources which are within nsigma*fwhm/2 pixels of each other, leaving
    only a single valid source in that region.

    This algorithm only works for sources which end up sequentially next to each other
    based on Y position and removes enough duplicates to make the final source list more
    managable.  It sorts the positions by Y value in order to group those at the
    same positions as much as possible.
    """
    for n in range(niter):
        if len(fitind) < 1:
            break
        fitarr = np.array(fitind,np.float32)
        fluxarr = np.array(fluxes,np.float32)
        inpind = np.argsort(fitarr[:,1])
        npind = fitarr[inpind]
        fluxind = fluxarr[inpind]
        fitind = npind.tolist()
        fluxes = fluxind.tolist()
        dx = npind[1:,0] - npind[:-1,0]
        dy = npind[1:,1] - npind[:-1,1]
        dr = np.sqrt(np.power(dx,2)+np.power(dy,2))
        nsame = np.where(dr <= separation)[0]
        if nsame.shape[0] > 0:
            for ind in nsame[-1::-1]:
                #continue # <- turn off filtering by source separation
                del fitind[ind]
                del fluxes[ind]
        else:
            break
    return fitind,fluxes

def xy_round(data,x0,y0,skymode,ker2d,xsigsq,ysigsq,datamin=None,datamax=None):
    """ Compute center of source
    Original code from IRAF.noao.digiphot.daofind.apfind ap_xy_round()
    """
    nyk,nxk = ker2d.shape
    if datamin is None:
        datamin = data.min()
    if datamax is None:
        datamax = data.max()

    # call C function for speed now...
    xy_val = cdriz.arrxyround(data,x0,y0,skymode,ker2d,xsigsq,ysigsq,datamin,datamax)
    if xy_val is None:
        x = None
        y = None
        round = None
    else:
        x = xy_val[0]
        y = xy_val[1]
        round = xy_val[2]

    return x,y,round


def precompute_sharp_round(nxk, nyk, xc, yc):
    """
    Pre-computes mask arrays to be used by the 'sharp_round' function
    for roundness computations based on two- and four-fold symmetries.
    """

    # Create arrays for the two- and four-fold symmetry computations:
    s4m = np.ones((nyk,nxk),dtype=np.int16)
    s4m[yc, xc] = 0

    s2m = np.ones((nyk,nxk),dtype=np.int16)
    s2m[yc, xc] = 0
    s2m[yc:nyk, 0:xc]     = -1;
    s2m[0:yc+1, xc+1:nxk] = -1;

    return s2m, s4m


def sharp_round(data, density, kskip, xc, yc, s2m, s4m, nxk, nyk,
                datamin, datamax):
    """
    sharp_round -- Compute first estimate of the roundness and sharpness of the
    detected objects.

    A Python translation of the AP_SHARP_ROUND IRAF/DAOFIND function.
    """

    # Compute the first estimate of roundness:
    sum2 = np.sum(s2m*density)
    sum4 = np.sum(s4m*abs(density))
    if sum2 == 0.0:
        round = 0.0
    elif sum4 <= 0.0: # eps?
        round = None
    else:
        round = 2.0 * sum2 / sum4

    # Eliminate the sharpness test if the central pixel is bad:
    mid_data_pix = data[yc, xc]
    mid_dens_pix = density[yc, xc]
    if mid_data_pix > datamax:
        return True, round, None
    if mid_data_pix < datamin:
        return False, round, None

    ########################
    # Sharpness statistics:

    satur = np.max(kskip*data) > datamax

    # Exclude pixels (create a mask) outside the [datamin, datamax] range:
    uskip = np.where((data >= datamin) & (data <= datamax), 1, 0)
    # Update the mask with the "skipped" values from the convolution kernel:
    uskip *= kskip
    # Also, exclude central pixel:
    uskip[yc, xc] = 0

    npixels = np.sum(uskip)
    if (npixels < 1 or mid_dens_pix <= 0.0):
        return satur, round, None

    sharp = (mid_data_pix - np.sum(uskip*data)/npixels) / mid_dens_pix
    #sharp = (mid_data_pix - np.mean(uskip*data)) / mid_dens_pix

    return satur, round, sharp


def roundness(im):
    """
    from astropy.io import fits as pyfits
    data=pyfits.getdata('j94f05bgq_flt.fits',ext=1)
    star0=data[403:412,423:432]
    star=data[396:432,3522:3558]
    In [53]: findobj.roundness(star0)
    Out[53]: 0.99401955054989544
    In [54]: findobj.roundness(star)
    Out[54]: 0.83091919980660645

    """
    perimeter = im.shape[0]*2 +im.shape[1]*2 -4
    area = im.size
    return 4*np.pi*area/perimeter**2

def immoments(im, p,q):
    x = list(range(im.shape[1]))
    y = list(range(im.shape[0]))
    #coord=np.array([x.flatten(),y.flatten()]).T
    """
    moment = 0
    momentx = 0
    for i in x.flatten():
        moment+=momentx
        sumx=0
        for j in y.flatten():
            sumx+=i**0*j**0*star0[i,j]
    """
    moment = np.sum([i**p*j**q*im[i,j] for j in x for i in y], dtype=np.float64)
    return moment
#ss=[i**0*j**0*list(star0[i,j].flatten()) for i in list(x.flatten()) for j in list(y.flatten())]
def nmoment(im,p,q):
    m = immoments(im,p,q)
    nmoment = m/np.sum(im, dtype=np.float64)
    return nmoment

def centroid(im):
    """
    Computes the centroid of an image using the image moments:

    centroid = {m10/m00, m01/m00}

    These calls point to Python version of moments function
    m00 = immoments(im,0,0)
    m10 = immoments(im, 1,0)
    m01 = immoments(im,0,1)

    """
    # These calls point to Python version of moments function
    m00 = cdriz.arrmoments(im,0,0)
    m10 = cdriz.arrmoments(im, 1,0)
    m01 = cdriz.arrmoments(im,0,1)

    ycen = m10 / m00
    xcen = m01 / m00
    return xcen, ycen


def cmoment(im,p,q):
    xcen,ycen = centroid(im)
    #x,y=np.meshgrid(range(403,412),range(423,432))
    x = list(range(im.shape[1]))
    y = list(range(im.shape[0]))
    mu = np.sum([(i-xcen)**p * (j-ycen)**q * im[i,j] for i in y for j in x],
                dtype=np.float64)
    return mu

def central_moments(im):
    xcen,ycen = centroid(im)
    mu00 = cmoment(im,p=0,q=0)
    mu01 = 0.
    mu10 = 0.
    mu11 = immoments(im,1,1) - xcen * immoments(im,0,1)
    mu20 = immoments(im,2,0) - xcen * immoments(im,1,0)
    mu02 = immoments(im,0,2) - ycen*immoments(im,0,1)
    mu21 = immoments(im,2,1) - 2*xcen*immoments(im,1,1) - ycen*immoments(im,2,0) + \
         2*xcen**2*immoments(im,0,1)
    mu12 = immoments(im,1,2) - 2*ycen*immoments(im,1,1) - xcen*immoments(im,0,2) + \
         2*ycen**2*immoments(im,1,0)
    mu30 = immoments(im,3,0) - 3*xcen*immoments(im,2,0) + 2*xcen**2*immoments(im,1,0)
    mu03 = immoments(im,0,3) - 3*ycen*immoments(im,0,2) + 2*ycen**2*immoments(im,0,1)
    cmoments = {'mu00': mu00,
                'mu01': mu01,
                'mu10': mu10,
                'mu11': mu11,
                'mu20': mu20,
                'mu02': mu02,
                'mu21': mu21,
                'mu12': mu12,
                'mu30': mu30,
                'mu03': mu03
                }
    return cmoments

def covmat(im):
    cmoments = central_moments(im)
    nmu20 = cmoments['mu20'] / cmoments['mu00']
    nmu02 = cmoments['mu02'] / cmoments['mu00']
    nmu11 = cmoments['mu11'] / cmoments['mu00']
    covmat = np.array([[nmu20, nmu11],[nmu11,nmu02]])
    return covmat
