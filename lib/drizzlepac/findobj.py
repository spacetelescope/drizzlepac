
from __future__ import division
import sys

import math
import numpy as np
from stsci import convolve
from stsci import ndimage as ndim

import stsci.imagestats as imagestats
from . import cdriz

#def gaussian(amplitude, xcen, ycen, xsigma, ysigma):
#from numpy import *

fwhm2sig = 2*np.sqrt(2*np.log(2))

#def gaussian1(height, x0, y0, a, b, c):
def gaussian1(height, x0, y0, fwhm, nsigma=1.5, ratio=1., theta=0.0):
    """
    height - the amplitude of the gaussian
    x0, y0, - center of the gaussian
    fwhm - full width at half maximum of the observation
    nsigma - cut the gaussian at nsigma
    ratio = ratio of xsigma/ysigma
    theta - angle of position angle of the major axis measured
    counter-clockwise from the x axis
    """

    xsigma = (1./fwhm2sig) * fwhm # computes a kernel that matches daofind
    ysigma = ratio * xsigma
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
            print 'Unable to construct 1D Gaussian with these parameters\n'
            raise ValueError
    else: #2D gaussian
        a = (math.cos(theta)**2/xsigma**2 + math.sin(theta)**2/ysigma**2)
        b = 2 * math.cos(theta) * math.sin(theta) *(1/xsigma**2-1./ysigma**2)
        c = (math.sin(theta)**2/xsigma**2 + math.cos(theta)**2/ysigma**2)

    discrim = b**2 - 4*a*c
    f = nsigma**2/2.
    #nx = int(2*max(2, math.sqrt(-8*c*f/discrim)))
    #ny = int(2*max(2, math.sqrt(-8*a*f/discrim)))

    return lambda x, y: height * np.exp(-0.5* (a*(x-x0)**2 + 2*b*(x-x0)*(y-y0) + c*(y-y0)**2))


def gausspars(fwhm, nsigma=1.5, ratio=1, theta=0.):
    xsigma = (1/fwhm2sig) * fwhm
    ysigma = ratio * xsigma
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
            print 'Unable to construct 1D Gaussian with these parameters\n'
            raise ValueError
        theta = np.deg2rad(theta)
        nx = 2 * int(max(2, (xsigma*nsigma*np.abs(np.cos(theta)))))+1
        ny = 2 * int(max(2, (xsigma*nsigma*np.abs(np.sin(theta)))))+1
    else: #2D gaussian
        theta = np.deg2rad(theta)
        a = math.cos(theta)**2/(xsigma**2) + math.sin(theta)**2/(ysigma**2)
        b = 2 * math.cos(theta) * math.sin(theta) *((1.0/xsigma**2)-(1./ysigma**2))
        c = math.sin(theta)**2/xsigma**2 + math.cos(theta)**2/ysigma**2
        discrim = b**2 - 4*a*c
        f = nsigma**2/2.
        nx = int(2*max(2, math.sqrt(-8*c*f/discrim)))+1
        ny = int(2*max(2, math.sqrt(-8*a*f/discrim)))+1

    return nx, ny

def gaussian(height, center_x, center_y, width_x, width_y):
    #Returns a gaussian function with the given parameters
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)


def moments(data,cntr):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
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
                ratio=1, nsigma=1.5, theta=0.):

    # Define convolution inputs
    nx, ny = gausspars(fwhm, nsigma=nsigma, ratio= ratio, theta=theta)
    xin, yin = np.mgrid[0:nx, 0:ny]
    kernel = gaussian1(1.0, nx//2, ny//2, fwhm)(xin,yin)  #+np.random.random(xin.shape)
    # define size of extraction box for each source based on kernel size
    gradius = nx//2

    rmat = np.sqrt((xin-nx//2)**2 + (yin-ny//2)**2)
    xyrmask = np.where(rmat <= gradius,1,0).astype(np.int16)
    npts = xyrmask.sum()

    rmask = kernel*xyrmask
    denom = (rmask*rmask).sum() - rmask.sum()**2/npts
    nkern = (rmask - (rmask.sum()/npts))/denom # normalize kernel to preserve fluxes for thresholds
    nkern *= xyrmask

    # initialize values used for getting source centers
    relerr = 1./((rmask**2).sum() - (rmask.sum()**2/xyrmask.sum()))

    xsigsq = (fwhm/fwhm2sig)**2
    ysigsq = (fwhm/fwhm2sig)**2

    # convolve image with gaussian kernel
    convdata = convolve.convolve2d(jdata, nkern).astype(np.float32)

    # clip image to create regions around each source for segmentation
    #tdata=np.where(convdata > skymode*2.0, convdata, 0)
    tdata=np.where(convdata > threshold, convdata, 0)
    # segment image and find sources
    s = ndim.generate_binary_structure(2,2)
    ldata,nobj=ndim.label(tdata,structure=s)
    fobjects = ndim.find_objects(ldata)
    #print 'Number of potential sources: ',nobj

    fluxes = []
    fitind = []
    if nobj < 2:
        print 'No objects found for this image. Please check value of "threshold".'
        return fitind,fluxes

    # determine center of each source, while removing spurious sources or
    # applying limits defined by the user
    ninit = 0
    ninit2 = 0
    minxy = gradius * 2 + 1

    for ss,n in zip(fobjects,range(len(fobjects))):
        ssx = ss[1].stop - ss[1].start
        ssy = ss[0].stop - ss[0].start 
        if ssx >= tdata.shape[1]-1 or ssy >= tdata.shape[0]-1:
            continue
        yr0 = ss[0].start-gradius
        yr1 = ss[0].stop+gradius+1
        if yr0 <= 0 or yr1 >= jdata.shape[0]: continue # ignore sources within ny//2 of edge
        #if yr0 <= 0: yr0 = 0
        #if yr1 >= jdata.shape[0]: yr1 = jdata.shape[0]

        xr0 = ss[1].start - gradius
        xr1 = ss[1].stop + gradius+1
        if xr0 <= 0 or xr1 >= jdata.shape[1]: continue # ignore sources within nx//2 of edge
        #if xr0 <= 0: xr0 = 0
        #if xr1 >= jdata.shape[1]: xr1 = jdata.shape[1]

        ssnew = (slice(yr0,yr1),slice(xr0,xr1))
        region = tdata[ssnew]

        #if region.shape[0] < minxy or region.shape[1] < minxy:
        #    continue

        cntr = centroid(region)

        # Define region centered on max value in object (slice)
        # This region will be bounds-checked to insure that it only accesses
        # a valid section of the image (not off the edge)
        maxpos = (int(cntr[1]+0.5)+ssnew[0].start,int(cntr[0]+0.5)+ssnew[1].start)
        yr0 = maxpos[0]-gradius
        yr1 = maxpos[0]+gradius+1
        if yr0 < 0 or yr1 > jdata.shape[0]:
            continue
        xr0 = maxpos[1] - gradius
        xr1 = maxpos[1] + gradius+1
        if xr0 < 0 or xr1 > jdata.shape[1]:
            continue

        #ninit += 1
        # Simple Centroid on the region from the convoluted image
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

        #ninit2 += 1
        px,py,pround = xy_round(jregion,gradius,gradius,skymode,
                kernel,xsigsq,ysigsq)

        if px is None:
            continue

        fitind.append((px+xr0,py+yr0))
        # compute a source flux value
        fluxes.append(jregion.sum())

    fitindc,fluxesc = apply_nsigma_separation(fitind,fluxes,fwhm*nsigma/2)

    #print 'ninit: ',ninit,'   ninit2: ',ninit2,' final n: ',len(fitind)
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

def roundness(im):
    """
    import pyfits
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
    x = range(im.shape[1])
    y = range(im.shape[0])
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
    """
    """
    # These calls point to Python version of moments function
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
    x = range(im.shape[1])
    y = range(im.shape[0])
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
