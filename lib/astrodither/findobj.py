
from __future__ import division 
import sys

import math
import numpy as np
import scipy as sp
from stsci import convolve
from stsci import ndimage as ndim

import stsci.imagestats as imagestats

#def gaussian(amplitude, xcen, ycen, xsigma, ysigma):
#from numpy import *
#from scipy import optimize

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
    
    xsigma = 0.42466 * fwhm * np.sqrt(2) # computes a kernel that matches daofind
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
    
    return lambda x, y: height * np.exp(- (a*(x-x0)**2 + 2*b*(x-x0)*(y-y0) + c*(y-y0)**2))
    

def gausspars(fwhm, nsigma=1.5, ratio=1, theta=0.):
    xsigma = 0.42466 * fwhm
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
        a = math.cos(theta)**2/xsigma**2 + math.sin(theta)**2/ysigma**2
        b = 2 * math.cos(theta) * math.sin(theta) *(1/xsigma**2-1./ysigma**2)
        c = math.sin(theta)**2/xsigma**2 + math.cos(theta)**2/ysigma**2
        
    discrim = b**2 - 4*a*c
    f = nsigma**2/2.
    nx = int(2*max(2, math.sqrt(-8*c*f/discrim)))
    ny = int(2*max(2, math.sqrt(-8*a*f/discrim)))

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

def fitgaussian(data,cntr):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit
    """
    from scipy import optimize
    if data.shape[1] < 3 or data.shape[0] < 3:
        raise ValueError

    params = moments(data,cntr)
    #params = gausspars(data, fwhm) 
    #errorfunction = lambda p: np.ravel(gaussian1(*p)(*np.indices(data.shape)) -
    #                             data)
    y,x = np.indices(data.shape)

    p = optimize.leastsq(errfunc, params, args=(data, x,y), maxfev=10)
    return p[0]

def findstars(jdata, fwhm, threshold, skymode, datamax=None, ratio=1, nsigma=1.5, theta=0.):
    # Define convolution inputs
    nx, ny = gausspars(fwhm, nsigma=nsigma, ratio= ratio, theta=theta)
    xin, yin = np.mgrid[0:nx+1, 0:ny+1]
    gsigma = 2*np.sqrt(2*np.log(2))*fwhm
    kernel = gaussian1(threshold, nx/2., ny/2., fwhm)(xin,yin)  #+np.random.random(xin.shape)
    kernel /= kernel.sum() # normalize kernel to preserve fluxes for thresholds

    # initialize values used for getting source centers
    nkern = kernel.copy()/kernel.max() 
    kern_f = (nsigma**2)/2.0
    xmat = xin - (nx/2)
    ymat = yin - (ny/2)
    rmat = np.sqrt(xmat**2+ymat**2)
    xyrmask = np.where(rmat <=2.001,1,0).astype(np.int16)
    rmask = nkern*xyrmask
    relerr = 1./((rmask**2).sum() - (rmask.sum()**2/xyrmask.sum()))

    xsigsq = (fwhm/fwhm2sig)**2
    ysigsq = (fwhm/fwhm2sig)**2

    # convolve image with gaussian kernel
    convdata = convolve.convolve2d(jdata, kernel)

    # clip image to create regions around each source for segmentation
    tdata=np.where(convdata > skymode*2.0, convdata, 0)
    # segment image and find sources
    ldata,nobj=ndim.label(tdata)
    fobjects = ndim.find_objects(ldata)

    # define size of extraction box for each source based on kernel size
    gradius = kernel.shape[0]//2

    fluxes = []
    fitind = []
    if nobj < 2:
        print 'No objects found for this image. Please check value of "threshold".'
        return fitind,lind
    #print 'Initially identified ',nobj,' possible sources for threshold = ',threshold
    # determine center of each source, while removing spurious sources or
    # applying limits defined by the user 
    ninit = 0
    ninit2 = 0
    for ss in fobjects:
        region = tdata[ss]
        # skip sources which have pixel values greater than the user set limit
        if (region.max() < (threshold/relerr)):
            continue
        
        cntr = centroid(region)

        # Define region centered on max value in object (slice)
        # This region will be bounds-checked to insure that it only accesses
        # a valid section of the image (not off the edge)
        maxpos = (int(cntr[1]+0.5)+ss[0].start,int(cntr[0]+0.5)+ss[1].start)
        yr0 = maxpos[0]-gradius
        yr1 = maxpos[0]+gradius+1
        if yr0 < 0 or yr1 > jdata.shape[0]: 
            continue
        xr0 = maxpos[1] - gradius
        xr1 = maxpos[1] + gradius+1
        if xr0 < 0 or xr1 > jdata.shape[1]:
            continue
        
        ninit += 1
        # Simple Centroid on the region from the convoluted image 
        jregion = jdata[yr0:yr1,xr0:xr1]
        region = tdata[yr0:yr1,xr0:xr1]
        
        if (datamax is not None and region.max() >= datamax):
            continue
        ninit2 += 1

        """
        # Refine position using gaussian fitting based on scipy 
        try:
            p = fitgaussian(jregion,cntr)
            if np.isnan(p[1]) or np.isnan(p[2]):
                continue
            fitind.append((p[1]+xr0, p[2]+yr0))
        except ValueError:
            pass
        """
        px,py,pround = xy_round(jregion,gradius,gradius,skymode,
                nkern,xsigsq,ysigsq,datamax=datamax)
        if px is None:
            continue
        fitind.append((px+xr0,py+yr0))
        # compute a source flux value
        fluxes.append(jregion.sum())
    #print 'ninit: ',ninit,'   ninit2: ',ninit2,' final n: ',len(fitind)
    return fitind, fluxes
    
def xy_round(data,x0,y0,skymode,ker2d,xsigsq,ysigsq,datamin=None,datamax=None):
    """ Compute center of source
    Original code from IRAF.noao.digiphot.daofind.apfind ap_xy_round()
    """
    nyk,nxk = ker2d.shape
    if datamin is None:
        datamin = data.min()
    if datamax is None:
        datamax = data.max()

    # These are all 0-based indices
    xhalf = (nxk / 2.0) - 0.5
    yhalf = (nyk / 2.0) - 0.5
    xmiddle = nxk // 2
    ymiddle = nyk // 2

    # Initialize the x fit.
    sumgd = 0.0
    sumgsq = 0.0
    sumg = 0.0
    sumd = 0.0
    sumdx = 0.0
    sdgdx = 0.0
    sdgdxsq = 0.0
    sddgdx = 0.0
    sgdgdx = 0.0
    p = 0.0
    n = 0

    # Compute the sums required for the x fit.
    for k in range(nxk):
        sg = 0.0
        sd = 0.0
        for j in range(nyk): 
            wt = float(ymiddle+1 - abs (j - ymiddle))
            pixval = data[y0-ymiddle+j,x0-xmiddle+k]
            if (pixval < datamin or pixval > datamax):
                sg=None
                break
            
            sd = sd + (pixval - skymode) * wt
            sg = sg + ker2d[j,k] * wt

        if (sg <= 0.0):
            sg = None
            break
        dxk = xmiddle-k
        wt = float(xmiddle+1 - abs (dxk))  
        sumgd += wt * sg * sd
        sumgsq += wt * sg ** 2
        sumg += wt * sg
        sumd += wt * sd
        sumdx += wt * sd * dxk

        p += wt
        n += 1
        dgdx = sg * dxk
        sdgdxsq += wt * dgdx ** 2
        sdgdx += wt * dgdx
        sddgdx += wt * sd * dgdx
        sgdgdx += wt * sg * dgdx

    # Need at least three points to estimate the x height, position
    # and local sky brightness of the star.

    if sg is None or (n <= 2 or p <= 0.0) :
        return None,None,None

    # Solve for the height of the best-fitting gaussian to the
    # xmarginal. Reject the star if the height is non-positive.

    hx = sumgsq - (sumg ** 2) / p
    if (hx <= 0.0):
        return None,None,None

    hx = (sumgd - sumg * sumd / p) / hx
    if (hx <= 0.0):
        return None,None,None

    # Solve for the new x centroid.
    skylvl = (sumd - hx * sumg) / p
    dx = (sgdgdx - (sddgdx - sdgdx * (hx * sumg + skylvl * p))) / \
        (hx * sdgdxsq / xsigsq)
    if (abs (dx) > xhalf):
        if (sumd == 0.0):
            dx = 0.0
        else:
            dx = sumdx / sumd 
        if (abs (dx) > xhalf):
            dx = 0.0
            
    x = int(x0)+1 + dx 

    # Initialize y fit.
    sumgd = 0.0
    sumgsq = 0.0
    sumg = 0.0
    sumd = 0.0
    sumdx = 0.0
    sdgdx = 0.0
    sdgdxsq = 0.0
    sddgdx = 0.0
    sgdgdx = 0.0
    p = 0.0
    n = 0

    for j in range(nyk):
        sg = 0.0
        sd = 0.0
        for k in range(nxk) :
            wt = float(xmiddle+1 - abs (k - xmiddle))
            pixval = data[y0-ymiddle+j,x0-xmiddle+k]
            if (pixval < datamin or pixval > datamax):
                sg = None
                break
            sd += (pixval - skymode) * wt
            sg += ker2d[j,k] * wt
        
        if (sg <= 0.0):
            sg = None
            break
        dyj = ymiddle-j
        wt = float(ymiddle+1 - abs (j - ymiddle))

        sumgd += wt * sg * sd
        sumgsq += wt * sg ** 2
        sumg += wt * sg
        sumd += wt * sd   
        sumdx += wt * sd * dyj
        p = p + wt
        n = n + 1
        dgdx = sg * dyj
        sdgdx += wt * dgdx
        sdgdxsq += wt * dgdx ** 2
        sddgdx += wt * sd * dgdx
        sgdgdx += wt * sg * dgdx

    # Need at least three points to estimate the y height, position
    # and local sky brightness of the star.

    if sg is None or (n <= 2 or p <= 0.0):
        return None,None,None

    # Solve for the height of the best-fitting gaussian to the
    # y marginal. Reject the star if the height is non-positive.

    hy = sumgsq - (sumg ** 2) / p 
    if (hy <= 0.0):
        return None,None,None

    hy = (sumgd - sumg * sumd / p) / (sumgsq - (sumg ** 2) / p)
    if (hy <= 0.0):
        return None,None,None
    
    # Solve for the new x centroid.
    skylvl = (sumd - hy * sumg) / p
    dy = (sgdgdx - (sddgdx - sdgdx * (hy * sumg + skylvl * p))) / \
        (hy * sdgdxsq / ysigsq)
    if (abs (dy) > yhalf):
        if (sumd == 0.0):
            dy = 0.0
        else :
            dy = sumdx / sumd 
            
        if (abs (dy) > yhalf):
            dy = 0.0
    
    y = int(y0)+1 + dy
    round = 2.0 * (hx - hy) / (hx + hy)
    return x,y,round

def roundness(im):
    """
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
    m00 = immoments(im,0,0)
    m10 = immoments(im, 1,0)
    m01 = immoments(im,0,1)
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

"""
def central_moments(im, xcen, ycen, p, q):
    moment = 0.

    moment += 
    for (y = 0, y < ny; y++) {;
        for (x = 0; x < nx; x++){
            moment += (pow((x-xcen),p) * pow((y-ycen),q) * fimage[x + y*nx]);
        }
    }

    return moment;
""" 
#def run(jdata, fwhm=2.5, ratio=1, nsigma=1.5, theta=0., threshold=4., skysigma=
"""
from pylab import *
# Create the gaussian data
Xin, Yin = mgrid[0:201, 0:201]
data = gaussian(3, 100, 100, 20, 40)(Xin, Yin) + np.random.random(Xin.shape)
data1 = gaussian1(5, 100, 200, .2)(Xin, Yin) + np.random.random(Xin.shape)
matshow(data1, cmap=cm.gist_earth_r)

params = fitgaussian(data1)
fit = gaussian(*params)

contour(fit(*indices(data.shape)), cmap=cm.copper)
ax = gca()
(height, x, y, width_x, width_y) = params
"""
#text(0.95, 0.05, 
"""
#x : %.1f
#y : %.1f
#width_x : %.1f
#width_y : %.1f
#%(x, y, width_x, width_y),
#        fontsize=16, horizontalalignment='right',
#        verticalalignment='bottom', transform=ax.transAxes)
"""
'''
histograms
-ACS:
h=im.histogram1d((data.flatten()).astype(np.int), 1900, 500, data.min())
x=np.arange(h.minValue,h.maxValue, 500)
plt.plot(x,np.log10(h.histogram))
threshold=25000
objacs=findobj.findstars(data,2.5,60000)
objacsar=np.array(objacs)
plt.imshow(data, vmin=-1100, vmax=1000, aspect='auto')
plt.scatter(objacsar[:,1],objacsar[:,0], marker='+', color='r')


WFPC2 - chip4 - almost no stars
hu4=im.histogram1d((u4.flatten()).astype(np.int), 512, 100, u4.min())
xu4=np.arange(hu4.minValue,hu4.maxValue, 100)
plt.plot(xu4,np.log10(hu4.histogram))

threshold=450 (500)
'''

'''
testing detection
ACS:
objacs=findobj.findstars(data,2.5,25000)
obj=np.array(objacs)
objacsar = np.roll(obj,1)
w1=wcsutil.HSTWCS('j94f05bgq_flt.fits', ext=1)
w4=wcsutil.HSTWCS('j94f05bgq_flt.fits', ext=4)
outwcs=utils.output_wcs([w1,w4])
sky=w1.all_pix2sky(objacsar+1,1)
tan=outwcs.wcs_sky2pix(sky,1)
datasingle=pyfits.getdata('j94f05bgq_single_sci.fits')
plt.imshow(datasingle, vmin=-11, vmax=150, aspect='auto')
plt.scatter(tan[:,0],tan[:,1], marker='+', color='r')

f=open('tvmarkobj.txt', 'w')
for i in tan:
    f.writelines([str(i[0])," ", str(i[1]), '\n'])
    
f.close()
'''

'''
compare daofind
f=open('j94f05bgq_flt.fits.coo.1')
lines=f.readlines()
f.close()
nlines=[l.strip().split() for l in lines[43:]]
nlinesarr=np.array(nlines)
xdao=nlinesarr[:,0].astype(np.float)
ydao=nlinesarr[:,1].astype(np.float)
skydao=w1.all_pix2sky(xdao,ydao,1)
tandao=outwcs.wcs_sky2pix(np.array(skydao).T,1)
plt.scatter(tan[:,0],tan[:,1], marker='+', color='r')
'''

"""
#imports
import numpy as np
import pyfits
import findobj
from stwcs.distortion import utils
from stwcs import wcsutil
from matplotlib import pyplot as plt

"""
