import string,os
import numpy as np
import ndimage

from pytools import asnutil,irafglob,parseinput
import pyfits

def parse_input(input,prodonly=False):    
    catlist = None
    
    if (isinstance(input, list) == False) and \
       ('_asn' in input or '_asc' in input) :
        # Input is an association table
        # Get the input files
        oldasndict = asnutil.readASNTable(input, prodonly=prodonly)
        filelist = [fileutil.buildRootname(fname) for fname in oldasndict['order']]
        
    elif (isinstance(input, list) == False) and \
       (input[0] == '@') :
        # input is an @ file
        f = open(input[1:])
        # Read the first line in order to determine whether
        # catalog files have been specified in a second column...
        line = f.readline()
        f.close()
        # Parse the @-file with irafglob to extract the input filename
        filelist = irafglob.irafglob(input, atfile=atfile_sci)
        print line
        # If there are additional columns for catalog files...
        if len(line.split()) > 1:
            # ...parse out the names of the catalog files as well 
            catlist = parse_atfile_cat(input)        
    else:
        #input is a string or a python list
        try:
            filelist, output = parseinput.parseinput(input)
            filelist.sort()
        except IOError: raise    

    return filelist,catlist
    
def atfile_sci(line):
    return line.split()[0]

def parse_atfile_cat(input):
    """ Return the list of catalog filenames specified as part of the input @-file
    """
    # input is an @ file
    f = open(input[1:])
    catlist = []
    for line in f.readlines():
        catlist.append(line.split()[1:])
    f.close()
    return catlist

#
# functions to help work with configobj input
#
def get_configobj_root(configobj):
    kwargs = {}
    for key in configobj:
        # Only copy in those entries which start with lower case letters
        # since sections are all upper-case for this task
        if key[0].islower(): kwargs[key] = configobj[key]
    return kwargs

# Object finding algorithm based on NDIMAGE routines
def ndfind(array,hmin,fwhm,sharplim=[0.2,1.0],roundlim=[-1,1],minpix=5):
    """ Source finding algorithm based on NDIMAGE routines
    """
    #cimg,c1 = idlgauss_convolve(array,sigma)
    #cimg = np.abs(ndimage.gaussian_laplace(array,fwhm))
    cimg = -1*ndimage.gaussian_laplace(array,fwhm)
    cimg = np.clip(cimg,0,cimg.max())
    
    cmask = cimg >= hmin
    # find and label sources
    ckern = ndimage.generate_binary_structure(2,1)
    clabeled,cnum = ndimage.label(cmask,structure=ckern)
    cobjs = ndimage.find_objects(clabeled)
    xpos = []
    ypos = []
    flux = []
    for s in cobjs:
        nmask = cmask[s].sum()
        if nmask >= minpix: # eliminate spurious detections
            yx = ndimage.center_of_mass(cimg[s]*cmask[s])
            # convert position to chip position in (0-based) X,Y
            xpos.append(yx[1]+s[1].start)
            ypos.append(yx[0]+s[0].start)
            flux.append((array[s]*cmask[s]).sum())
    # Still need to implement sharpness and roundness limits
    return np.array(xpos),np.array(ypos),np.array(flux)

def isfloat(value):
    """ Return True if all characters are part of a floating point value
    """
    tab = string.maketrans('','')
    # Check to see if everything but '-+.e' in string are numbers
    # as these characters would be used for floating-point/exponential numbers
    if (value.translate(tab.lower(),'-.+e').isdigit()):
        return True
    else:
        return False
        
def readcols(infile, cols=None):
    """ Function which reads specified columns from either FITS tables or 
        ASCII files
    """
    if infile.find('.fits') > 0:
        outarr = read_FITS_cols(infile,cols=cols)
    else:
        outarr = read_ASCII_cols(infile,cols=cols)
    return outarr

def read_FITS_cols(infile,cols=None):
    """ Read columns from FITS table
    """        
    ftab = pyfits.open(infile)
    extnum = 0
    extfound = False
    for extn in ftab:
        if extn.header.has_key('tfields'):
            extfound = True
            break
        extnum += 1
    if not extfound:
        print 'ERROR: No catalog table found in ',infile
        ftab.close()
        raise ValueError
    # Now, read columns from the table in this extension
    # if no column names were provided by user, simply read in all columns from table
    if cols[0] in [None,' ','','INDEF']:
        cols = ftab[extnum].data.names
    # Define the output 
    outarr = []
    for c in cols:
        outarr.append(ftab[extnum].data.field(c))
    
    ftab.close()
    return outarr

def read_ASCII_cols(infile,cols=[1,2,3,4]):
    """ Copied from 'reftools.wtraxyutils'.
        Input column numbers must be 1-based, or 'c'+1-based ('c1','c2',...')
    """
    # Convert input column names into column numbers 
    colnums = []
    for colname in cols:
        cnum = None
        if colname not in [None,""," ","INDEF"]:
            if isinstance(colname, str) and colname[0] == 'c':
                cname = colname[1:]
            else:
                cname = colname
            colnums.append(int(cname)-1)
    outarr = [] # initialize output result
    print 'colnums: ',colnums,cols
    # Open catalog file
    fin = open(infile,'r')
    for l in fin.readlines(): # interpret each line from catalog file
        l = l.strip()
        if len(l) == 0 or len(l.split()) < len(colnums) or (len(l) > 0 and l[0] == '#' or (l.find("INDEF") > -1)): continue
        for i in range(10):
            lnew = l.replace("  "," ")
            if lnew == l: break
            else: l = lnew
        lspl = lnew.split(" ")
        
        if len(outarr) == 0:
            if len(colnums) == 0: # No columns were specified, return them all
                colnums = range(len(lspl))
            for c in range(len(colnums)): outarr.append([])

        for c,n in zip(colnums,range(len(colnums))):
            if isfloat(lspl[c]):
                cval = float(lspl[c])
            else:
                cval = lspl[c]
                
            outarr[n].append(cval)
    fin.close()

    for n in range(len(colnums)):
        outarr[n] = np.array(outarr[n])

    return outarr            
 
def write_shiftfile(image_list,filename,outwcs='tweak_wcs.fits'):
    """ Write out a shiftfile for a given list of input Image class objects
    """       
    rows = ''
    nrows = 0
    for img in image_list:
        row = img.get_shiftfile_row()
        if row is not None:
            rows += row
            nrows += 1
    if nrows == 0: # If there are no fits to report, do not write out a file
        return
    
    # write out reference WCS now
    if os.path.exists(outwcs):
        os.remove(outwcs)
    p = pyfits.HDUList()
    p.append(pyfits.PrimaryHDU())
    p.append(createWcsHDU(image_list[0].refWCS))
    p.writeto(outwcs)
    
    # Write out shiftfile to go with reference WCS
    if os.path.exists(filename):
        os.remove(filename)
    f = open(filename,'w')
    f.write('# frame: output\n')
    f.write('# refimage: %s[wcs]\n'%outwcs)
    f.write('# form: delta\n')
    f.write('# units: pixels\n')
    f.write(rows)
    f.close()
    print 'Writing out shiftfile :',filename
    
def createWcsHDU(wcs):
    """ Generate a WCS header object that can be used to
        populate a reference WCS HDU.
        
        For most applications, 
        stwcs.wcsutil.HSTWCS.wcs2header() will work just as well.
    """

    hdu = pyfits.ImageHDU()
    hdu.header.update('EXTNAME','WCS')
    hdu.header.update('EXTVER',1)
    # Now, update original image size information
    hdu.header.update('WCSAXES',2,comment="number of World Coordinate System axes")
    hdu.header.update('NPIX1',wcs.naxis1,comment="Length of array axis 1")
    hdu.header.update('NPIX2',wcs.naxis2,comment="Length of array axis 2")
    hdu.header.update('PIXVALUE',0.0,comment="values of pixels in array")

    # Write out values to header...
    hdu.header.update('CD1_1',wcs.wcs.cd[0,0],comment="partial of first axis coordinate w.r.t. x")
    hdu.header.update('CD1_2',wcs.wcs.cd[0,1],comment="partial of first axis coordinate w.r.t. y")
    hdu.header.update('CD2_1',wcs.wcs.cd[1,0],comment="partial of second axis coordinate w.r.t. x")
    hdu.header.update('CD2_2',wcs.wcs.cd[1,1],comment="partial of second axis coordinate w.r.t. y")
    hdu.header.update('ORIENTAT',wcs.orientat,comment="position angle of image y axis (deg. e of n)")
    hdu.header.update('CRPIX1',wcs.wcs.crpix[0],comment="x-coordinate of reference pixel")
    hdu.header.update('CRPIX2',wcs.wcs.crpix[1],comment="y-coordinate of reference pixel")
    hdu.header.update('CRVAL1',wcs.wcs.crval[0],comment="first axis value at reference pixel")
    hdu.header.update('CRVAL2',wcs.wcs.crval[1],comment="second axis value at reference pixel")
    hdu.header.update('CTYPE1',wcs.wcs.ctype[0],comment="the coordinate type for the first axis")
    hdu.header.update('CTYPE2',wcs.wcs.ctype[1],comment="the coordinate type for the second axis")

    return hdu
    
#
# Code used for testing source finding algorithms
#
def idlgauss_convolve(image,fwhm):
    sigmatofwhm = 2*np.sqrt(2*np.log(2))
    radius = 1.5 * fwhm / sigmatofwhm # Radius is 1.5 sigma
    if radius < 1.0:
        radius = 1.0
        fwhm = sigmatofwhm/1.5
        print( "WARNING!!! Radius of convolution box smaller than one." )
        print( "Setting the 'fwhm' to minimum value, %f." %fwhm )
    sigsq = (fwhm/sigmatofwhm)**2 # sigma squared
    nhalf = int(radius) # Center of the kernel
    nbox = 2*nhalf+1 # Number of pixels inside of convolution box
    middle = nhalf # Index of central pixel
    
    kern_y, kern_x = np.ix_(np.arange(nbox),np.arange(nbox)) # x,y coordinates of the kernel
    g = (kern_x-nhalf)**2+(kern_y-nhalf)**2 # Compute the square of the distance to the center
    mask = g <= radius**2 # We make a mask to select the inner circle of radius "radius"
    nmask = mask.sum() # The number of pixels in the mask within the inner circle.
    g = np.exp(-0.5*g/sigsq) # We make the 2D gaussian profile
    
    ###
    # Convolving the image with a kernel representing a gaussian (which is assumed to be the psf)
    ###
    c = g*mask # For the kernel, values further than "radius" are equal to zero
    c[mask] = (c[mask] - c[mask].mean())/(c[mask].var() * nmask) # We normalize the gaussian kernel
    
    c1 = g[nhalf] # c1 will be used to the test the roundness
    sumc1 = c1.mean()
    sumc1sq = (c1**2).sum() - sumc1
    c1 = (c1-c1.mean())/((c1**2).sum() - c1.mean())

    h = ndimage.convolve(image,c,mode='constant',cval=0.0) # Convolve image with kernel "c"
    h[:nhalf,:] = 0 # Set the sides to zero in order to avoid border effects
    h[-nhalf:,:] = 0
    h[:,:nhalf] = 0
    h[:,-nhalf:] = 0
    
    return h,c1

def gauss_array(nx,ny=None,sigma_x=1.0,sigma_y=None,zero_norm=False):
    """ Computes the 2D Gaussian with size n1*n2.  
        Sigma_x and sigma_y are the stddev of the Gaussian functions.
        The kernel will be normalized to a sum of 1.
    """

    if ny == None: ny = nx
    if sigma_y == None: sigma_y = sigma_x
    
    xradius = nx//2
    yradius = ny//2
        
    # Create grids of distance from center in X and Y
    xarr = np.abs(np.arange(-xradius,xradius+1))
    yarr = np.abs(np.arange(-yradius,yradius+1))
    hnx = gauss(xarr,sigma_x)
    hny = gauss(yarr,sigma_y)
    hny = hny.reshape((ny,1))
    h = hnx*hny
    
    # Normalize gaussian kernel to a sum of 1
    h = h / np.abs(h).sum()
    if zero_norm:
        h -= h.mean()

    return h
  
def gauss(x,sigma):
    """ Compute 1-D value of gaussian at position x relative to center."""
    return np.exp(-np.power(x,2)/(2*np.power(sigma,2))) / (sigma*np.sqrt(2*np.pi))
