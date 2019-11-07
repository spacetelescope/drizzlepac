#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4 ai :
"""
A set of routines used in matching sources by offset histograms

Path
----
HLApipeline/regression_testing/starmatch_hist.py

Dependencies
------------
* HLApipeline/regression_testing/infrot.py

Inputs
------
None.

Classes and Functions
---------------------
"""

import os,sys,numpy
import scipy.special, scipy.signal
from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.morphology import generate_binary_structure
import astropy.io.fits as pyfits
from astropy.table import Table

from drizzlepac import util
from drizzlepac.devutils.comparison_tools import infrot
from stsci.tools import logutil

log = logutil.create_logger('starmatch_hist', level=logutil.logging.INFO, stream=sys.stdout)

msgunit = sys.stdout
@util.with_logging
def run(source_list_dict, minimum_match=10, xref=0.0, yref=0.0, postarg=None):
    """
    Match source lists in x,y coordinates allowing for possible shift & rotation

    :param source_list_dict: dictionary indexed by coo filename with number of sources in each file as value
    :param minimum_match: minimum number of matches in peak. Default value = '10'.
    :param xref: X reference pixel in image (default 0, which is a very bad value if image rotates) 
    :param yref: Y reference pixel in image (default 0, which is a very bad value if image rotates) 
    :param postarg: dictionary indexed by coo filename with (dx, dy) in pixels for postarg (default is to read FITS headers)
    :type source_list_dict: dictionary
    :type minimum_match: integer
    :type xref: float
    :type yref: float
    :type postarg: dictionary
    :returns: dictionary of lists of matching sourcelist indicies indexed by coo filename.
    """
    out_dict={}
    ### assumes list with greatest number of sources should be the reference file
    (coo_ref,number)=max(iter(source_list_dict.items()), key=lambda x:x[1])
    if postarg is None:
        postarg = getpostarg(source_list_dict)
    for coo_img in list(source_list_dict.keys()):
        matched_flag=True
        if coo_img == coo_ref:
            matched_list=['#This is reference']
        else:
            matched_list,matching_lines_ref,matching_lines_img=match_cat_histogram(coo_ref,coo_img,maxdiff=50, step=1, verbose=True, extra_verbose=True,
                    minimum_match=minimum_match, xref=xref, yref=yref,postarg_ref=postarg[coo_ref], postarg_img=postarg[coo_img])
            if not matched_list:
                matched_list=['#Not a thing found']
            out_dict[coo_ref] = matching_lines_ref
            out_dict[coo_img] = matching_lines_img
        if len(matched_list) < minimum_match:
            matched_flag = False
            if coo_img == coo_ref:
                log.info("\n %s is reference image, starmatch_hist leaving WCS alone\n" %(coo_img,))
            else:
                log.info("\n %s -> %s starmatch_hist failed! Leaving WCS alone\n" %(coo_img,coo_ref))
    return(out_dict)


def read_cat_file(catfile):
    """
    read in specified catalog file or daophot output file and force result to be a 2-D array if it is not empty
    
    :param catfile: catalog file to read in
    :type catfile: string
    :returns: Catalog file information (reshaped if need be)
    """
    if catfile.endswith("point-cat.ecsv"):
        ecsvData = Table.read(catfile, format='ascii.ecsv') #now will read in and process HAP catalog data
        data = numpy.stack((ecsvData["X-Center"].data, ecsvData["Y-Center"].data), axis=-1)
    elif catfile.endswith("segment-cat.ecsv"):
        ecsvData = Table.read(catfile, format='ascii.ecsv') #now will read in and process HAP catalog data
        data = numpy.stack((ecsvData["X-Centroid"].data, ecsvData["Y-Centroid"].data), axis=-1)

    else:
        try:
            daoData = Table.read(catfile, format='ascii.daophot') #now will read in and process daophot data
            data = numpy.stack((daoData["XCENTER"].data, daoData["YCENTER"].data), axis=-1)
        except:
            try:
                data = numpy.loadtxt(catfile, comments='#', skiprows=0,usecols=(0,1)) #orig. loadtxt call for .coo files
            except:
                data = numpy.loadtxt(catfile, comments='#', usecols=(0, 1), delimiter=',', skiprows=1) #new loadtxt call for daophot.txt, sexphot.txt sourcelists
    # force result to be 2-D if it is not empty
    if data.size == 2:
        data = data.reshape(1, data.size)
    return data


def getpostarg(source_list_dict):
    """
    Read FITS files used for the source list to get POSTARG info

    :param source_list_dict: dictionary
    :returns: dictionary of POSTARG info, keyed by image name
    """
    rv = {}
    for coo_img in source_list_dict:
        fitsfile = coo_img.replace(coo_img.split("_")[-1],"drz.fits")
        if not os.path.exists(fitsfile):
            fitsfile.replace("drz.fits","drc.fits")
        if (os.path.exists(fitsfile) and fitsfile.endswith(".fits")):
            try:
                fh = pyfits.open(fitsfile)
                hdr = fh[0].header
                # assume image has north up (which is true for HLA images)
                pixsize = abs(hdr['cd1_1'])*3600.0 #TODO: figure out how to get "cd1_1" (or eqivlent) value from final drizzle-combined (e.g. hst_10265_01_acs_wfc_f606w_drz.fits or hst_10265_01_acs_wfc_total_drz.fits) images
                postarg1 = hdr['postarg1']/pixsize
                postarg2 = hdr['postarg2']/pixsize
                pa_aper = hdr['pa_aper']*numpy.pi/180 #TODO: figure out how to get "pa_aper" (or eqivlent) value from final drizzle-combined (e.g. hst_10265_01_acs_wfc_f606w_drz.fits or hst_10265_01_acs_wfc_total_drz.fits) images
                fh.close()
                cospa = numpy.cos(pa_aper)
                sinpa = numpy.sin(pa_aper)
                dx = cospa*postarg1-sinpa*postarg2
                dy = sinpa*postarg1+cospa*postarg2
                rv[coo_img] = (dx, dy)
            except:
                log.info("Warning: unable to fetch POSTARG info for image {}. Using POSTARG=0".format(fitsfile))
                rv[coo_img] = (0.0, 0.0)
        else:
            log.info('FITS file %s not found, using POSTARG=0' % fitsfile)
            rv[coo_img] = (0.0,0.0)
    return rv


def match_cat_histogram(catHH, catWW, maxdiff=50, step=1, verbose = False, extra_verbose = True, minimum_match=10, xref=0.0, yref=0.0,
        postarg_ref=(0.0,0.0), postarg_img=(0.0,0.0)):
    """
    Procedure to find a coordinate match between two images
    in the presence of substantial CR contamination
    
    Modified to work on pixel coordinates rather than RA, dec.
    Assumes rectified tangent plane projection.  Note that
    the step and tolerance (maxdiff) need to be in pixels.
    
    Default operation assumes that there is not a large rotation
    (a small rotation is acceptable and is included in the matching)
    and the shift is limited to less than maxdiff pixels in each
    direction.
    
    Procedure includes the following steps:
    
    #. Compute two-dimensional histogram of source-to-source distances
    #. Find best peak in histogram, determine significance
    #. Match sources on basis of peak
    #. Refine relative offset
    #. Return matched column list of coordinates in the sense: ``X_ref Y_ref    X_transform Y_transform`` (ref and transform can be switched, in this order just to be consistent with XYXYMATCH output, hence right input for GEOTRAN) 
    
    Tested.
    
    :param catHH: reference catalog name
    :param catWW: coordinate catalog name
    :param maxdiff: maximum difference value to use when computing histograms in *findOffset()*. Default value = 50.
    :param step: max separation used by *catMatch()* (in pixels) to match sources in the catalogs. Default value = 1.
    :param verbose: Verbose output (True/False)? Default value = False
    :param extra_verbose: Even more verbose output (True/False)? Default value = False
    :param minimum_match: Minimum number of matches. Default value = 10.
    :param xref: x-axis pixel offset to place reference pixel at zero. Default value = 0.0
    :param yref: y-axis pixel offset to place reference pixel at zero. Default value = 0.0
    :param postarg_ref: tuple (dx,dy) with postarg pixel offsets for reference image
    :param postarg_img: tuple (dx,dy) with postarg pixel offsets for coordinate image
    :type catHH: string
    :type catWW: string
    :type maxdiff: integer
    :type step: integer
    :type verbose: Boolean
    :type extra_verbose: Boolean
    :type minimum_match: integer
    :type xref: float
    :type yref: float
    :type postarg_ref: tuple (float, float)
    :type postarg_img: tuple (float, float)
    :returns: matched column list of coordinates, list of catHH lines that match catWW lines, and list of catWW lines that match catHH lines
    """
    if extra_verbose:
        log.info("Matching %s --> %s using the Stefano method."%(catWW.split("/")[-1], catHH.split("/")[-1]))
    # Read the catalogs - Column 0 and 1 are the X and Y respectively
    chh = read_cat_file(catHH)
    if chh.size == 0:
        log.info("Empty reference file {}".format(catHH))
        return []
    x1 = chh[:,0]
    y1 = chh[:,1]
    cww = read_cat_file(catWW)
    if cww.size == 0:
        log.info("Empty coordinate file {}".format(catWW))
        return []
    x2 = cww[:,0]
    y2 = cww[:,1]

    # offset to put reference pixel at zero
    # this is what makes geomap work when there is rotation
    x1 -= xref
    y1 -= yref
    x2 -= xref
    y2 -= yref

    # pdb.set_trace() #- for de-bugging
    rv = findOffsetAndRotation(x1,y1, x2,y2, maxdiff, postarg_ref, postarg_img, verbose = verbose)
    if rv is None:
        if verbose: log.info("No significant peak found, skipping fine match")
        return ["#No significant peak found"],[],[]
    base_offset, rotangle = rv

    if verbose:
        log.info("base_offset, rotangle {} {}".format(base_offset, rotangle))
    # Now match sources within bin; use radius = twice the bin size
    index = catMatch(x1,y1, x2,y2, 2.*step, base_offset, rotangle)
    # index is an array of length = n(cat1) with is -1 if the ith source
    # has no match, and points to the matching source in cat2 if there
    # is a match.  Thus the matched positions are:
    ww = numpy.where(index > -1)[0]
    nsub = numpy.size(ww)
    subx1 = x1[ww]
    suby1 = y1[ww]
    subx2 = x2[index[ww]]
    suby2 = y2[index[ww]]
    if verbose: log.info("{}".format(nsub))

    # iterate to tighten up match if there is a tight cluster
    # this is not useful if there is a rotation, so skip it 
    wsub = numpy.arange(nsub)
    if rotangle == 0:
        for iter in range(0,5):
            dx = subx1-subx2
            dy = suby1-suby2
            meandx = numpy.mean(dx[wsub])
            meandy = numpy.mean(dy[wsub])
            dist = numpy.sqrt((dx-meandx)**2+(dy-meandy)**2)
            rmsdist = numpy.sqrt(numpy.mean(dist[wsub]**2))
            wsub = numpy.where(dist <= 2.*rmsdist)[0]
            if verbose: log.info("{} {} {} {}".format(meandx, meandy, rmsdist, wsub.size))

    nsub = wsub.size
    output_array = numpy.empty((nsub,4), dtype=float)
    output_array[:,0] = subx1[wsub]
    output_array[:,1] = suby1[wsub]
    output_array[:,2] = subx2[wsub]
    output_array[:,3] = suby2[wsub]

    matching_lines_HH = ww[wsub]
    matching_lines_WW = index[ww][wsub]
    # print "Sourcelist Matching Results"
    # print "Reference sourcelist:  {} of {} total sources matched ({} %)".format(len(matching_lines_HH),len(x1),100.0*(float(len(matching_lines_HH))/float(len(x1))))
    # print "Comparison sourcelist: {} of {} total sources matched ({} %)".format(len(matching_lines_WW),len(x2),100.0*(float(len(matching_lines_WW))/float(len(x2))))
    if nsub >= minimum_match:
        # compute the shift and rotation
        xshift, yshift, rotation = infrot.getxyshiftrot(output_array[:, 0], output_array[:, 1], output_array[:, 2], output_array[:, 3], xref=0.0, yref=0.0)
        log.info("infrot xshift yshift rotation {} {} {}".format(xshift, yshift, rotation))

    output_array = output_array.astype(str)
    column_list = [None]*nsub
    for i, v in enumerate(output_array):
        column_list[i] = ' '.join(v)
    return column_list,matching_lines_HH,matching_lines_WW


def findOffsetAndRotation(x1,y1, x2, y2, maxdiff, postarg1, postarg2, verbose = False, oversample=5):
    """
    Determine the offset of the positions in x1,y1 and x2,y2 using the
    histogram method

    This version uses an algorithm based on my xymatch function
    (as in catMatch).  This uses a binning size of 1 pixel oversampled by a factor
    oversample (default 5).  The oversample value is forced to be odd.

    The approach is:
    
    #. Construct a 2D histogram count_diff in delta_x and delta_y from -maxdiff to maxdiff with a resolution 1/oversample pixels.  The position differences are computed for each source pair from the x1,y1 and x2,y2 lists.
    #. Smooth the histogram with a boxcar of size oversample, giving a bin resolution of 1 pixel
    #. Compute Hanning-filter smoothed version of the histogram.
    #. Find all significant peaks in the histogram.
    #. Filter out peaks near postarg offset (unless offset is near zero)
    #. Return the position of the best peak x, y and rotation angle. This returns zeros if the peak is judged not to be significant.
        
    :param x1: X coordinates to match (list 1) 
    :param y1: Y coordinates to match (list 1)
    :param x2: X coordinates to match (list 2) 
    :param y2: Y coordinates to match (list 2)
    :param maxdiff: Max x,y difference value that will be used to generate histogram
    :param postarg1: (dx,dy) for image 1
    :param postarg2: (dx,dy) for image 2
    :param verbose: Verbose output (True/False)? Default value = False
    :param oversample: oversampling factor (NOTE: should be odd value). Default value = 5
    :type x1: numpy.ndarray
    :type y1: numpy.ndarray
    :type x2: numpy.ndarray
    :type y2: numpy.ndarray
    :type maxdiff: Integer
    :type postarg1: (float, float)
    :type postarg2: (float, float)
    :type verbose: Boolean
    :type oversample: integer
    :returns: (numpy.ndarray, float) x, y position of the peak and rotation angle.
    """
    if x1.ndim != 1 or y1.ndim != 1 or x2.ndim != 1 or y2.ndim != 1:
        raise ValueError("x and y parameters must be 1-D arrays")
    n1 = len(x1)
    n2 = len(x2)
    if n1 != len(y1) or n2 != len(y2):
        raise ValueError("x and y arrays must be of equal length")

    # force oversample to be odd
    oversample = oversample + (1 - (oversample % 2))
    halfsample = (oversample-1)/2
    foversample = float(oversample)

    # add padding around edges so smoothing works better
    # midpoint of array
    xmid = maxdiff*oversample + halfsample
    # size of array
    nx = 2*xmid + 1
    ny = nx
    # trimmed size
    nxtrim = nx - 2*halfsample
    nytrim = ny - 2*halfsample

    # get all pairs that match in box
    p1, p2 = boxMatch(x1,y1,x2,y2,maxdiff)

    # make histograms for a range of rotations
    # rotation increment is determined by image size
    rotinc = 1.0/max(abs(x1).max(), abs(x2).max())
    # nrot = 51

    #XXX RLW, 2017 October 17
    #XXX this version has reverted to just doing rotation 0
    #XXX keeping the rest of the code for the future though
    nrot = 1
    #XXX RLW, 2017 October 17

    rotangle = (numpy.arange(nrot,dtype=float) - (nrot-1)/2) * rotinc
    crot = numpy.cos(rotangle)
    srot = numpy.sin(rotangle)
    count_diff = numpy.zeros((nrot,int(nytrim),int(nxtrim)), dtype=numpy.int64)
    for k in range(nrot):
        cx2 = x2*crot[k] - y2*srot[k]
        cy2 = y2*crot[k] + x2*srot[k]
        delta_x = x1[p1] - cx2[p2]
        delta_y = y1[p1] - cy2[p2]
        ii = (delta_x*oversample + xmid).astype(int)
        jj = (delta_y*oversample + xmid).astype(int)
        count_diff[k,:,:] = fast_boxsum(bincount2d(jj,ii,ny,nx,clipped=False), oversample)
    xmid = xmid - halfsample

    irot0 = (nrot-1)/2
    smax0 = count_diff[int(irot0)].max()
    smax = count_diff.max()
    # use zero-rotation histogram unless a rotated version is at least 1-sigma better
    sthresh = (smax+smax0)/2.0 + numpy.sqrt((smax+smax0)/2.0)
    if smax <= sthresh:
        if verbose and smax > smax0:
            log.info("Choosing zero-rot histogram peak {} over rotated max {} threshold={}".format(smax0,smax,sthresh))
        count_diff = count_diff[int(irot0)]
        rotfactor = 1
        rotangle = 0.0
    else:
        # pull out best slice of the histogram array
        irot, jmax, imax = numpy.unravel_index(numpy.argmax(count_diff), count_diff.shape)
        count_diff = count_diff[irot]
        rotfactor = 2*numpy.abs(irot-irot0) + 1
        rotangle = rotangle[irot]
        if verbose:
            log.info("Best histogram has rotation {} degrees".format(rotangle*180/numpy.pi))
            log.info("Rotated peak {} unrotated peak {}".format(smax,smax0))

    # find the best peak
    jmax, imax, npred = findbestpeak(count_diff, oversample, postarg1, postarg2, verbose=verbose)
    npred = npred*rotfactor
    maxcount = count_diff[jmax,imax]
    if verbose:
        log.info("{} {} {}".format(maxcount, imax, jmax))
        log.info('number of predicted peaks of this size = {}'.format(npred))
    # Is it significant?
    # Require probability > 1-e3 that one of the pixels iin the histogram has this max or bigger
    if npred >= 0.001:
        if verbose: log.info('Maximum not very significant: predicted {} peaks'.format(npred))
        return None
    x_base_offset = (imax-xmid)/foversample
    y_base_offset = (jmax-xmid)/foversample
    return (numpy.array((x_base_offset, y_base_offset)), rotangle)


def fast_boxsum(a, size):
    """
    Return box sum over region size**2 pixels while trimming off the (size-1)/2 pixels on each edge
    
    :param a: array of values to compute sums over
    :param size: size of sampling box used to compute sums
    :type a: numpy.ndarray
    :type size: integer
    :returns: numpy.ndarray (2-D) of sums
    """
    b = numpy.insert(a,0,0,axis=1).cumsum(axis=1)
    b = b[:,size:] - b[:,:-size]
    b = numpy.insert(b,0,0,axis=0).cumsum(axis=0)
    b = b[size:] - b[:-size]
    return b


def findbestpeak(h, size, postarg1, postarg2, minthresh=1.e-3, verbose=False):
    """
    Return the x,y location of the best histogram peak along with the number of predicted peaks of that amplitude
    
    :param h: 2-D histogram of counts
    :param size: number of samples per pixel (used in smoothing)
    :param postarg1: (dx,dy) for image 1
    :param postarg2: (dx,dy) for image 2
    :param minthresh: minimum probability threshold for peaks to consider. Default value = '1.0e-3'.
    :param verbose: if true, prints info
    :type h: numpy.ndarray
    :type size: integer
    :type postarg1: (float, float)
    :type postarg2: (float, float)
    :type minthresh: float
    :type verbose: boolean
    :returns: tuple (jmax, imax, npred)
    """

    # Find maximum value
    jmax, imax = numpy.unravel_index(numpy.argmax(h), h.shape)
    maxcount = h[jmax,imax]
    if maxcount == 0:
        if verbose:
            log.info("No matches in histogram")
        return (jmax, imax, 1.0)
    hmean = numpy.mean(h)

    # determine peak amplitude to get below minthresh random probability
    # The gammainc function gives the Poisson cumulative probability that N >= maxcount
    prob = scipy.special.gammainc(numpy.arange(maxcount)+1, hmean)
    # get slope to empirically improve probability estimate
    power = getprobcorr(prob, h, jmax, imax, size, verbose=verbose)
    if verbose:
        log.info("False-peak probability parameter {}".format(power))
    prob = prob**power
    w = numpy.where(prob <= minthresh)[0]
    if w.size == 0:
        # no significant peaks
        # just return the highest
        return (jmax, imax, prob.min()*h.size)
    thresh = w[0]+1
    # locate all peaks above threshold in Hanning-smoothed version of image
    sh = hanning_smooth(h, size)

    jj, ii = findpeaks(sh, thresh)

    # weed out peaks near postarg offset
    xmid = (h.shape[1]-1)/2
    ymid = (h.shape[0]-1)/2
    dx = size*(postarg2[0] - postarg1[0])
    dy = size*(postarg2[1] - postarg1[1])
    # do not apply postarg filtering if zero shift is included
    if dx**2+dy**2 > size**2:
        w = numpy.where((ii-dx-xmid)**2+(jj-dy-ymid)**2 > size**2)
        ss = ii.size
        ii = ii[w]
        jj = jj[w]
        if ss != ii.size and verbose:
            log.info("Filtered out {} peak near POSTARG position".format(ss-ii.size))
    else:
        if verbose:
            log.info("POSTARG near zero, no peaks excluded")

    if jj.size == 0:
        # no peaks found
        return (jmax, imax, 1.0)
    # sort by increasing distance
    distance = numpy.sqrt((jj-h.shape[0]/2)**2 + (ii-h.shape[1]/2)**2)
    index = numpy.argsort(distance)
    distance = distance[index]
    ii = ii[index]
    jj = jj[index]
    # compute local mean of h at each peak using annulus
    shalf = (size-1)/2
    r1 = 4*shalf
    r2 = 6*shalf
    ky, kx = numpy.ogrid[-r2:r2+1,-r2:r2+1]
    kr = kx*kx+ky*ky
    kernel = ((kr >= r1**2) & (kr <= r2**2)).astype(float)
    lmean = scipy.signal.convolve2d(h.astype(float), kernel, mode='same', boundary='fill', fillvalue=0)
    norm = scipy.signal.convolve2d(h*0+1.0, kernel, mode='same', boundary='fill', fillvalue=0)
    # use the global mean if it is higher than the local value
    numpy.clip(lmean/norm, hmean, None, out=lmean)
    hpeak = h[jj,ii]
    bmean = lmean[jj,ii]
    npred = numpy.pi*(distance+1)**2 * scipy.special.gammainc(hpeak, bmean)**power
    # sort by expected chance count, then by decreasing peak amplitude, then by increasing background
    # this handles the case where the predicted number is zero
    index = numpy.lexsort((bmean,-hpeak,npred))
    ii = ii[index]
    jj = jj[index]
    hpeak = hpeak[index]
    bmean = bmean[index]
    npred = npred[index]
    keep = (ii-ii[0])^2+(jj-jj[0])^2 > size^2
    keep[0] = True
    ii = ii[keep]
    jj = jj[keep]
    hpeak = hpeak[keep]
    bmean = bmean[keep]
    npred = npred[keep]
    if len(ii) == 1:
        if verbose:
            log.info("Only 1 peak left after weeding out nearby peaks")
        npbest = npred[0]
    else:
        if len(npred) > 1 and npred[1] < minthresh and hpeak[1] >= hpeak[0]-3*numpy.sqrt(hpeak[0]):
            # if second-highest peak is also significant, this peak is not reliable
            # return highest peak but with low significance
            npbest = 1.0
            if verbose:
                log.info("More than 1 significant peak, not reliable (top 2 npred = %e %e hpeak = %d %d)" % (npred[0], npred[1], hpeak[0], hpeak[1]))
        else:
            if len(npred) > 1 and npred[1] < minthresh and verbose:
                log.info("2nd peak judged not significant (top 2 npred = %e %e hpeak = %d %d)" % (npred[0], npred[1], hpeak[0], hpeak[1]))
            npbest = npred[0]
    return (jj[0], ii[0], npbest)


def findpeaks(image, thresh):
    """
    Return positions of all peaks in image above threshold thresh
    Based on `"detect_peaks" Stack Overflow discussion <https://stackoverflow.com/questions/3684484/peak-detection-in-a-2d-array/3689710#3689710>`_
    
    :param image: array of values to search
    :param thresh: threshold for peaks
    :type image: numpy.ndarray
    :type thresh: float
    :returns: index array (equivalent of where output)
    """
    # define an 8-connected neighborhood
    neighborhood = generate_binary_structure(2,2)
    # find local maximum for each pixel
    amax = maximum_filter(image, footprint=neighborhood)
    w = numpy.where((image == amax) & (image >= thresh))
    return w


def hanning_smooth(a, size):
    """
    Return Hanning-smoothed version of histogram
    
    :param a: array of values to smooth
    :param size: size of sampling box used to compute sums
    :type a: numpy.ndarray
    :type size: integer
    :returns: numpy.ndarray (2-D) with smoothed array
    """
    # create normalized 2-D Hanning filter
    hfilter = numpy.hanning(size+3)[1:-1]
    hfilter = numpy.outer(hfilter, hfilter)
    hfilter = hfilter/hfilter.sum()
    cim = scipy.signal.convolve2d(a, hfilter, mode='same', boundary='fill', fillvalue=0)
    # normalize for pixels off edge
    norm = scipy.signal.convolve2d(a*0+1.0, hfilter, mode='same', boundary='fill', fillvalue=0)
    return cim/norm

def getprobcorr(prob, h, jmax, imax, size, verbose=False):
    """Use measured distribution of histogram values to improve probability estimate

    Returns power <= 1 to raise original probability to get better value

    :param prob: array of probability estimates for counts from 1 to max(h)
    :param h: 2-D histogram of bin counts
    :param jmax: location of maximum value in h
    :param imax: location of maximum value in h
    :param size: number of samples per pixel
    :param verbose: if true, print info
    :type prob: array
    :type h: array
    :type jmax: int
    :type imax: int
    :type size: int
    :type verbose: boolean
    :returns: float power to improve probability
    """
    # remove points from histogram around peak
    keep = numpy.ones(h.shape, dtype=bool)
    shalf = int((size-1)/2)
    j1 = max(jmax-shalf,0)
    j2 = min(jmax+shalf, h.shape[0]-1)
    i1 = max(imax-shalf,0)
    i2 = min(imax+shalf, h.shape[1]-1)
    keep[j1:j2,i1:i2] = False
    hsub = h[numpy.where(keep)]
    # distribution of counts excluding the peak
    chist = numpy.bincount(hsub)
    # trim leading zeros (and always remove zero counts)
    hmin = max(numpy.argmax(chist != 0), 1)
    chist = chist[hmin:]
    # calculate reverse cumulative distribution
    csum = chist[::-1].cumsum()[::-1] / float(hsub.size)
    # note first element of prob is for 1 count, not zero
    xx = prob[hmin-1:hmin-1+csum.size]
    # avoid zero division (note these have zero weight too so they don't affect the weighted median)
    # also drop log(0) values
    w = numpy.where((xx != 1) & (xx != 0))[0]
    if w.size == 0:
        if verbose:
            log.info("No points to fit after removing zero weights")
        return 1.0
    xx = numpy.log(xx[w])
    csum = csum[w]
    # weighted median gives robust fit to power
    power = wtmedian(numpy.log(csum)/xx, -xx) 
    if power > 1:
        power = 1.0
    return power


def wtmedian(a, wt):
    """Computed weighted median of array a
    Use wt=1/sigma to get maximum likelihood estimator for
    Laplacian noise distribution with amplitude sigma.

    :param a: array with values
    :param wt: array with weights (non-negative with size as a)
    :type a: array_like
    :type wt: array_like
    :returns: a floating point weighted median value.
    """

    a = numpy.ravel(a)
    wt = numpy.ravel(wt).clip(min=0)
    if a.size != wt.size:
        raise ValueError("a and wt must be the same size")
    # ignore weights if all zero
    if wt.max() == 0:
        return numpy.median(a)

    index = numpy.argsort(a)
    wts = wt[index]
    wsum = wts.cumsum()

    # find first element in cumulative weight array >= sum/2
    wthresh = 0.5*wsum[-1]
    i = numpy.searchsorted(wsum, wthresh, side='left')
    if wsum[i] == wthresh:
        # Special case if sum is exactly half of total
        # Average with next non-zero weight (equivalent to median with
        # even number of points.)
        ihi = numpy.searchsorted(wsum, wthresh, side='right')
        return 0.5*(a[index[i]]+a[index[ihi]])
    else:
        return a[index[i]]


def bincount2d(jj, ii, ny, nx, clipped=False):

    """
    Fast 2-D histogram using numpy.bincount

    :param jj: 1-d array of y values
    :param ii: 1-d array of x values 
    :param ny: X size of output 2-d histogram
    :param nx: Y size of output 2-d histogram
    :param clipped: Clip output histogram (True/False)? Note: if clipped is True, ii/jj values are already clipped to range 0..nx-1,0..ny-1. Default value = False
    :type jj: numpy.ndarray
    :type ii: numpy.ndarray
    :type ny: integer 
    :type nx: integer
    :type clipped: Boolean
    :returns: (ny,nx) integer histogram with counts.
    """

    if clipped:
        index = ii + jj*nx
    else:
        ww = numpy.where((ii>=0) & (ii<nx) & (jj>=0) & (jj<ny))[0]
        index = ii[ww] + jj[ww]*nx
    return numpy.bincount(index.astype('int'), minlength=int(nx*ny)).reshape(int(ny),int(nx))


def catMatch(x1,y1, x2,y2, sep, offset=None, rotangle=None):

    """
    Routine to match two lists of objects by position using 2-D Cartesian distances

    Matches positions x1,y1 with positions in x2,y2, for matches within separation (sep).
    If more than one match is found, the nearest is returned.  Input catalogs need not be sorted.

    2017 June 15, Rick White
    Based loosely on Marcel Haas's translation of my IDL routine xymatch.pro
    
    :param x1: X components of catalog 1
    :param y1: Y components of catalog 1
    :param x2: X components of catalog 2
    :param y2: Y components of catalog 2
    :param sep: Maximum allowed separation (in pixels) between sources for said sources to be considered 'matching'
    :param offset: a 2-element array that (if specified) gets subtracted from the x1,y1 values
    :param rotangle: rotation in degrees applied to x1,y1 values
    :type x1: numpy.ndarray
    :type y1: numpy.ndarray
    :type x2: numpy.ndarray
    :type y2: numpy.ndarray
    :type sep: integer
    :type offset: numpy.ndarray
    :returns: an array of indices for 2nd list that correspond to the closest match (within sep) in cat2 for each object in cat1, so x1[i] matches x2[return_value[i]]. Note that objects in cat1 with no match within sep have indices -N-1 with N the length of cat2, so that IndexErrors will be raised if trying to assign these indices.
    """

    if x1.ndim != 1 or y1.ndim != 1 or x2.ndim != 1 or y2.ndim != 1:
        raise ValueError("x and y parameters must be 1-D arrays")
    if len(x1) != len(y1) or len(x2) != len(y2):
        raise ValueError("x and y arrays must be of equal length")

    # apply offset and rotation
    if rotangle is not None or offset is not None:
        if rotangle is None:
            rotangle = 0.0
        if offset is None:
            offset = (0.0, 0.0)
        elif len(offset) != 2:
            raise ValueError("offset must be a 2-element array")
        crot = numpy.cos(rotangle)
        srot = numpy.sin(rotangle)
        x1 = x1*crot + y1*srot - offset[0]
        y1 = y1*crot - x1*srot - offset[1]

    # Sort the arrays by increasing y-coordinate
    is1 = y1.argsort()
    x1 = x1[is1]
    y1 = y1[is1]
    is2 = y2.argsort()
    x2 = x2[is2]
    y2 = y2[is2]
    # find search limits in y2 for each object in y1
    # note this is designed to include the points that are exactly equal to maxdiff
    kvlo = y2.searchsorted(y1-sep,'left').clip(0, len(y2))
    kvhi = y2.searchsorted(y1+sep,'right').clip(kvlo, len(y2))

    nnomatch = 0
    n1 = len(x1)
    p2 = numpy.zeros(n1, dtype='int') - len(x2) - 1
    sepsq = sep**2
    for i in range(n1):
        y = y1[i]
        x = x1[i]
        klo = kvlo[i]
        khi = kvhi[i]
        dx = numpy.abs(x2[klo:khi] - x)
        w = (dx <= sep).nonzero()[0]
        if len(w) == 0:
            # Nothing matched
            nnomatch += 1
        else:
            distsq = (x - x2[klo+w])**2 + (y - y2[klo+w])**2
            if distsq.min() <= sepsq:
                p2[is1[i]] = is2[klo+w[distsq.argmin()]]
            else:
                nnomatch += 1
    return p2


def boxMatch(x1,y1, x2,y2, sep):

    """
    Routine to match two lists of objects by position using 2-D Cartesian distances

    Matches positions x1,y1 with positions in x2,y2, for matches within a box of size +- sep.
    *All* matches within the box are returned.  Input catalogs need not be sorted.

    2017 October 14, Rick White
    
    :param x1: X components of catalog 1
    :param y1: Y components of catalog 1
    :param x2: X components of catalog 2
    :param y2: Y components of catalog 2
    :param sep: Half-size of xy box in pixels between sources for said sources to be considered 'matching'
    :type x1: numpy.ndarray
    :type y1: numpy.ndarray
    :type x2: numpy.ndarray
    :type y2: numpy.ndarray
    :type sep: integer
    :returns: (p1, p2) index arrays into source lists with matching pairs
    """

    if x1.ndim != 1 or y1.ndim != 1 or x2.ndim != 1 or y2.ndim != 1:
        raise ValueError("x and y parameters must be 1-D arrays")
    if len(x1) != len(y1) or len(x2) != len(y2):
        raise ValueError("x and y arrays must be of equal length")

    # Sort the arrays by increasing y-coordinate
    is1 = y1.argsort()
    x1 = x1[is1]
    y1 = y1[is1]
    is2 = y2.argsort()
    x2 = x2[is2]
    y2 = y2[is2]
    # find search limits in y2 for each object in y1
    # note this is designed to include the points that are exactly equal to maxdiff
    kvlo = y2.searchsorted(y1-sep,'left').clip(0, len(y2))
    kvhi = y2.searchsorted(y1+sep,'right').clip(kvlo, len(y2))

    p1 = []
    p2 = []
    n1 = len(x1)
    for i in range(n1):
        y = y1[i]
        x = x1[i]
        klo = kvlo[i]
        khi = kvhi[i]
        dx = numpy.abs(x2[klo:khi] - x)
        w = (dx <= sep).nonzero()[0]
        if w.size > 0:
            p1.append(numpy.zeros(w.size,dtype=int)+is1[i])
            p2.append(is2[klo+w])
    if p1:
        return (numpy.hstack(p1), numpy.hstack(p2))
    else:
        return (numpy.array([],dtype=int), numpy.array([],dtype=int))


if __name__ == '__main__':
    import getopt

    def usage(msg=None):
        log.info("Usage: %s [-r xref,yref] coo1 coo2..." % sys.argv[0])
        if msg: log.info("{}".format(msg))
        sys.exit(1)

    try:
        opts, args = getopt.getopt(sys.argv[1:], "r:h")
    except getopt.error as e:
        usage(str(e))
    xref = 0.0
    yref = 0.0
    for opt, value in opts:
        if opt == "-r":
            try:
                xref, yref = list(map(float, value.split(',')))
            except ValueError as e:
                usage("Argument for -r must be comma-separated pair of floats")
        elif opt == "-h":
            usage()
        else:
            usage("Unknown option '%s'" % opt)
    msgunit = sys.stderr
    source_list_dict = {}
    for cooname in args:
        source_list_dict[cooname] = len(open(cooname).readlines())
    matched=run(source_list_dict, xref=xref, yref=yref)
    # print >> msgunit, 'matched keys'
    # print >> msgunit, matched.keys()
    # print "x1 y1 x2 y2"
    # print '\n'.join(matched[list2]['matchlist'])
