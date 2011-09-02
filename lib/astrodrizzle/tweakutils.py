import string,os
import numpy as np
import stsci.ndimage as ndimage

from stsci.tools import asnutil, irafglob, parseinput, fileutil
import pyfits
import astrolib.coords as coords

import pylab as pl
import stsci.imagestats as imagestats 

def parse_input(input, prodonly=False):    
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
            if input.find('*') > -1: # if wild-cards are given, sort for uniform usage
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
    
        This function provides a simple replacement for the DAOFIND task.
        
        Parameters
        ----------
        array : arr
            Input image as numpy array
        hmin  : float
            Limit for source detection in pixel values
        fwhm  : float
            Full-width half-maximum of the PSF in the image
        minpix : int
            Minimum number of pixels for any valid source
        sharplim : tuple
            [Not used at this time]
        roundlim : tuple
            [Not used at this time]
        
        Returns
        -------
        x  : arr
            Array of detected source X positions (in array coordinates, 0-based)
        y  : arr
            Array of detected source Y positions (in array coordinates, 0-based)
        flux : arr
            Array of detected source fluxes in pixel values
        id  : arr
            Array of detected source ID numbers
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
    return np.array(xpos),np.array(ypos),np.array(flux),np.arange(len(xpos))

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

def radec_hmstodd(ra,dec):
    """ Function to convert HMS values into decimal degrees.
        
        This function relies on the astrolib.coords package to perform the
        conversion to decimal degrees.  
        
        Parameters
        ----------
        ra  : list or array
            List or array of input RA positions 
        dec : list or array
            List or array of input Dec positions
            
        Returns
        -------
        pos : arr
            Array of RA,Dec positions in decimal degrees
            
        Notes
        -----
        Formats of ra and dec inputs supported::
        
            ["nn","nn","nn.nn"]
            "nn nn nn.nnn"
            "nn:nn:nn.nn"
            "nnH nnM nn.nnS" or "nnD nnM nn.nnS"
        
        See Also
        --------
        astrolib.coords
        
    """
    hmstrans = string.maketrans(string.letters,' '*len(string.letters))

    if isinstance(ra,list):
        rastr = ':'.join(ra)
    elif ra.find(':') < 0:
        # convert any non-numeric characters to spaces (we already know the units)
        rastr = string.translate(ra,hmstrans).strip()
        rastr = rastr.replace('  ',' ')
        # convert 'nn nn nn.nn' to final 'nn:nn:nn.nn' string
        rastr = rastr.replace(' ',':')
    else:
        rastr = ra

    if isinstance(dec,list):
        decstr = ':'.join(dec)
    elif dec.find(':') < 0:
        decstr = string.translate(dec,hmstrans).strip()
        decstr = decstr.replace('  ',' ')
        decstr = decstr.replace(' ',':')
    else:
        decstr = dec

    pos = coords.Position(rastr+' '+decstr,units='hours')
    return pos.dd()

def parse_colname(colname):
    """ Common function to interpret input column names provided by the user.
    
        This function translates column specification provided by the user 
        into a column number.

        Notes
        -----
        This function will understand the following inputs::
        
            '1,2,3' or   'c1,c2,c3' or ['c1','c2','c3']
            '1-3'   or   'c1-c3'
            '1:3'   or   'c1:c3'
            '1 2 3' or   'c1 c2 c3'
            '1'     or   'c1'
            1

        Parameters
        ----------
        colname :
            Column name or names to be interpreted
        
        Returns
        -------
        cols : list
            The return value will be a list of strings.

    """
    if isinstance(colname,list):
        cname = ''
        for c in colname:
            cname += str(c)+','
        cname = cname.rstrip(',')
    elif isinstance(colname,int) or colname.isdigit():
        cname = str(colname)
    else:   
        cname = colname
         
    if 'c' in cname[0]: cname = cname.replace('c','')
    
    ctok = None
    cols = None
    if '-' in cname:
        ctok = '-'
    if ':' in cname:
        ctok = ':'
    if ctok is not None:
        cnums = cname.split(ctok)
        c = range(int(cnums[0]),int(cnums[1])+1)
        cols = []
        for i in c:
            cols.append(str(i))
            
    if cols is None:
        ctok = ' '
        if ',' in cname:
            ctok = ','
        cols = cname.split(ctok)

    return cols

def readcols(infile, cols=None):
    """ Function which reads specified columns from either FITS tables or
        ASCII files

        This function reads in the columns specified by the user into numpy arrays
        regardless of the format of the input table (ASCII or FITS table).
        
        Parameters
        ----------
        infile : string
            Filename of the input file 
        cols   : string or list of strings
            Columns to be read into arrays
        
        Returns
        -------
        outarr : array
            Numpy array or arrays of columns from the table
    
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

def read_ASCII_cols(infile,cols=[1,2,3]):
    """ Interpret input ASCII file to return arrays for specified columns.
    
        Notes
        -----
        The specification of the columns should be expected to have lists for
        each 'column', with all columns in each list combined into a single entry.
        For example::
        
            cols = ['1,2,3','4,5,6',7]

        where '1,2,3' represent the X/RA values, '4,5,6' represent the Y/Dec values
        and 7 represents the flux value for a total of 3 requested columns of data
        to be returned.
        
        Returns
        -------
        outarr : list of arrays
            The return value will be a list of numpy arrays, one for each 'column'.
    """
    # build dictionary representing format of each row
    # Format of dictionary: {'colname':col_number,...}
    # This provides the mapping between column name and column number
    coldict = {}
    fin = open(infile,'r')
    flines = fin.readlines()
    fin.close()
    hcount = 0
    for l in flines: 
        if l[0] == '#': hcount += 1
    
    for l in flines: # interpret each line from catalog file
        if hcount == 1:
            # Parse colnames directly from column headers
            colnames = l.split()
            for name,i in zip(colnames,range(len(colnames))):
                name = name.replace('#','')
                coldict[name] = i
        elif l[0] == '#':
            continue
        else:
            # convert first row of data into column definitions using indices
            numcols = len(l.split())
            colnames = range(1,numcols+1)
            for name in colnames:
                coldict[str(name)] = name-1
            break
    numcols = len(cols)
    outarr = []
    for col in range(numcols):
        outarr.append([])
    convert_radec = False

    # Now, map specified columns to columns in file and populate output arrays
    # Open catalog file
    fin = open(infile,'r')
    for l in fin.readlines(): # interpret each line from catalog file
        if l[0] == '#':
            continue
        l = l.strip()
        # skip blank lines, comment lines, or lines with 
        # fewer columns than requested by user
        if len(l) == 0 or len(l.split()) < numcols or (
            len(l) > 0 and (l[0] == '#' or "INDEF" in l)
            ): continue
        lspl = l.split()
        nsplit = len(lspl)
        
        # For each 'column' requested by user, pull data from row
        for c,i in zip(cols,range(numcols)):
            cnames = parse_colname(c)
            if len(cnames) > 1:
                # interpret multi-column specification as one value
                outval = ''
                for cn in cnames:
                    cnum = coldict[cn]
                    cval = lspl[cnum]
                    outval += cval+' '
                outarr[i].append(outval)
                convert_radec = True
            else:                
                #pull single value from row for this column
                cnum = coldict[cnames[0]]
                if isfloat(lspl[cnum]):
                    cval = float(lspl[cnum])
                else:
                    cval = lspl[cnum]
                    # Check for multi-column values given as "nn:nn:nn.s"
                    if ':' in cval:
                        cval = cval.replace(':',' ')
                        convert_radec = True
                outarr[i].append(cval)
                
    fin.close()
    # convert multi-column RA/Dec specifications
    if convert_radec:
        outra = []
        outdec = []
        for ra,dec in zip(outarr[0],outarr[1]):
            radd,decdd = radec_hmstodd(ra,dec)
            outra.append(radd)
            outdec.append(decdd)
        outarr[0] = outra
        outarr[1] = outdec
        
    # convert all lists to numpy arrays
    for c in range(len(outarr)):
        outarr[c] = np.array(outarr[c])
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
    """ Generate a WCS header object that can be used to populate a reference WCS HDU.

        For most applications, stwcs.wcsutil.HSTWCS.wcs2header() 
        will work just as well.
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

def gauss_array(nx,ny=None,fwhm=1.0,sigma_x=None,sigma_y=None,zero_norm=False):
    """ Computes the 2D Gaussian with size nx*ny.
    
        Parameters
        ----------
        nx : int
        ny : int [Default: None]
            Size of output array for the generated Gaussian. If ny == None,
            output will be an array nx X nx pixels. 

        fwhm : float [Default: 1.0]
            Full-width, half-maximum of the Gaussian to be generated 

        sigma_x : float [Default:  None]
        sigma_y : float [Default:  None]
            Sigma_x and sigma_y are the stddev of the Gaussian functions.

        zero_norm : bool [Default:  False]
            The kernel will be normalized to a sum of 1 when True.
        
        Returns
        -------
        gauss_arr : array
            A numpy array with the generated gaussian function

    """

    if ny == None: ny = nx

    if sigma_x is None:
        if fwhm is None:
            print 'A value for either "fwhm" or "sigma_x" needs to be specified!'
            raise ValueError 
        else:
            # Convert input FWHM into sigma
            sigma_x = fwhm/(2*np.sqrt(2*np.log(2)))
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


#### Plotting Utilities for astrodrizzle
def make_vector_plot(coordfile,columns=[1,2,3,4],data=None,title=None, axes=None, every=1,
                    limit=None, xlower=None, ylower=None, output=None, headl=4,headw=3,
                    xsh=0.0,ysh=0.0,fit=None,scale=1.0,vector=True,textscale=5,append=False,linfit=False,rms=True):
    """ Convert a XYXYMATCH file into a vector plot or set of residuals plots.
    
        This function provides a single interface for generating either a vector
        plot of residuals or a set of 4 plots showing residuals.  The data being
        plotted can also be adjusted for a linear fit on-the-fly.  
        
        Parameters
        ----------
        coordfile : string
            Name of file with matched sets of coordinates. This input file can
            be a file compatible for use with IRAF's geomap. 
        columns : list [Default: [0,1,2,3]]
            Column numbers for the X,Y positions from each image
        data : list of arrays
            If specified, this can be used to input matched data directly
        title : string
            Title to be used for the generated plot
        axes : list
            List of X and Y min/max values to customize the plot axes
        every : int [Default: 1]
            Slice value for the data to be plotted
        limit : float
            Radial offset limit for selecting which sources are included in the plot
        xlower : float
        ylower : float
            Limit in X and/or Y offset for selecting which sources are included in the plot
        output : string
            Filename of output file for generated plot
        headl : int [Default: 4]
            Length of arrow head to be used in vector plot
        headw : int [Default: 3]
            Width of arrow head to be used in vector plot
        xsh : float
        ysh : float
            Shift in X and Y from linear fit to be applied to source positions 
            from the first image
        scale : float
            Scale from linear fit to be applied to source positions from the 
            first image
        fit : array
            Array of linear coefficients for rotation (and scale?) in X and Y from 
            a linear fit to be applied to source positions from the first image
        vector : bool [Default: True]
            Specifies whether or not to generate a vector plot. If False, task
            will generate a set of 4 residuals plots instead
        textscale : int [Default: 5]
            Scale factor for text used for labelling the generated plot
        append : bool [Default: False]
            If True, will overplot new plot on any pre-existing plot
        linfit : bool [Default: False]
            If True, a linear fit to the residuals will be generated and 
            added to the generated residuals plots
        rms : bool [Default: True]
            Specifies whether or not to report the RMS of the residuals as a
            label on the generated plot(s).
        
    
    """
    if data is None:
        data = readcols(coordfile,cols=columns)

    xy1x = data[0]
    xy1y = data[1]
    xy2x = data[2]
    xy2y = data[3]

    numpts = xy1x.shape[0]
    if fit is not None:
        xy1x,xy1y = apply_db_fit(data,fit,xsh=xsh,ysh=ysh)
        fitstr = '-Fit applied'
        dx = xy2x - xy1x
        dy = xy2y - xy1y
    else:
        dx = xy2x - xy1x - xsh
        dy = xy2y - xy1y - ysh
    # apply scaling factor to deltas
    dx *= scale
    dy *= scale

    print 'Total # points: ',len(dx)
    if limit is not None:
        indx = (np.sqrt(dx**2 + dy**2) <= limit)
        dx = dx[indx].copy()
        dy = dy[indx].copy()
        xy1x = xy1x[indx].copy()
        xy1y = xy1y[indx].copy()
    if xlower is not None:
        xindx = (np.abs(dx) >= xlower)
        dx = dx[xindx].copy()
        dy = dy[xindx].copy()
        xy1x = xy1x[xindx].copy()
        xy1y = xy1y[xindx].copy()
    print '# of points after clipping: ',len(dx)
    
    if output is not None:
        write_xy_file(output,[xy1x,xy1y,dx,dy])
        
    if not append:
        pl.clf()
    if vector: 
        dxs = imagestats.ImageStats(dx.astype(np.float32))
        dys = imagestats.ImageStats(dy.astype(np.float32))
        minx = xy1x.min()
        maxx = xy1x.max()
        miny = xy1y.min()
        maxy = xy1y.max()
        xrange = maxx - minx
        yrange = maxy - miny

        pl.quiver(xy1x[::every],xy1y[::every],dx[::every],dy[::every],\
                  units='y',headwidth=headw,headlength=headl)
        pl.text(minx+xrange*0.01, miny-yrange*(0.005*textscale),'DX: %f to %f +/- %f'%(dxs.min,dxs.max,dxs.stddev))
        pl.text(minx+xrange*0.01, miny-yrange*(0.01*textscale),'DY: %f to %f +/- %f'%(dys.min,dys.max,dys.stddev))
        pl.title(r"$Vector\ plot\ of\ %d/%d\ residuals:\ %s$"%(xy1x.shape[0],numpts,title))
    else:
        plot_defs = [[xy1x,dx,"X (pixels)","DX (pixels)"],\
                    [xy1y,dx,"Y (pixels)","DX (pixels)"],\
                    [xy1x,dy,"X (pixels)","DY (pixels)"],\
                    [xy1y,dy,"Y (pixels)","DY (pixels)"]]
        if axes is None:
            # Compute a global set of axis limits for all plots
            minx = xy1x.min()
            maxx = xy1x.max()
            miny = dx.min()
            maxy = dx.max()
            
            if xy1y.min() < minx: minx = xy1y.min()
            if xy1y.max() > maxx: maxx = xy1y.max()
            if dy.min() < miny: miny = dy.min()
            if dy.max() > maxy: maxy = dy.max()
        else:
            minx = axes[0][0]
            maxx = axes[0][1]
            miny = axes[1][0]
            maxy = axes[1][1]
        xrange = maxx - minx
        yrange = maxy - miny 
        
        for pnum,plot in zip(range(1,5),plot_defs):
            ax = pl.subplot(2,2,pnum)
            ax.plot(plot[0],plot[1],'.')
            if title is None:
                ax.set_title("Residuals [%d/%d]: No FIT applied"%(xy1x.shape[0],numpts))
            else:
                # This definition of the title supports math symbols in the title
                ax.set_title(r"$"+title+"$")
            pl.xlabel(plot[2])
            pl.ylabel(plot[3])
            lx=[ int((plot[0].min()-500)/500) * 500,int((plot[0].max()+500)/500) * 500]
            pl.plot([lx[0],lx[1]],[0.0,0.0],'k')
            pl.axis([minx,maxx,miny,maxy])
            if rms:
                pl.text(minx+xrange*0.01, maxy-yrange*(0.01*textscale),'RMS(X) = %f, RMS(Y) = %f'%(dx.std(),dy.std()))
            if linfit:
                lxr = int((lx[-1] - lx[0])/100)
                lyr = int((plot[1].max() - plot[1].min())/100)
                A = np.vstack([plot[0],np.ones(len(plot[0]))]).T
                m,c = np.linalg.lstsq(A,plot[1])[0]
                yr = [m*lx[0]+c,lx[-1]*m+c]
                pl.plot([lx[0],lx[-1]],yr,'r')
                pl.text(lx[0]+lxr,plot[1].max()+lyr,"%0.5g*x + %0.5g [%0.5g,%0.5g]"%(m,c,yr[0],yr[1]),color='r')
    
def apply_db_fit(data,fit,xsh=0.0,ysh=0.0):
    xy1x = data[0]
    xy1y = data[1]
    numpts = xy1x.shape[0]
    if fit is not None:
        xy1 = np.zeros((xy1x.shape[0],2),np.float64)
        xy1[:,0] = xy1x 
        xy1[:,1] = xy1y 
        xy1 = np.dot(xy1,fit)
        xy1x = xy1[:,0] + xsh
        xy1y = xy1[:,1] + ysh
    return xy1x,xy1y

def write_xy_file(outname,xydata,append=False,format=["%20.6f"]):
    if not isinstance(xydata,list):
        xydata = list(xydata)
    if not append:
        if os.path.exists(outname): os.remove(outname)
    fout1 = open(outname,'a+')
    for row in range(len(xydata[0][0])):
        outstr = ""
        for cols,fmts in zip(xydata,format):
            for col in range(len(cols)):
                outstr += fmts%(cols[col][row])
        fout1.write(outstr+"\n")
    fout1.close()
    print 'wrote XY data to: ',outname  
 
def find_xy_peak(img,center=None):
    """ Find the center of the peak of offsets
    """
    # find level of noise in histogram
    istats = imagestats.ImageStats(img.astype(np.float32),nclip=3)
    imgsum = img.sum()
    
    # clip out all values below mean+3*sigma from histogram
    imgc =img[:,:].copy()
    imgc[imgc < istats.mean+istats.stddev*3] = 0.0
    # identify position of peak
    yp0,xp0 = np.where(imgc == imgc.max())

    # take sum of 13x13 pixel box around peak
    xp_slice = (slice(int(yp0[0])-13,int(yp0[0])+15),
                slice(int(xp0[0])-13,int(xp0[0])+15))
    yp,xp = ndimage.center_of_mass(img[xp_slice])
    xp += xp_slice[1].start
    yp += xp_slice[0].start

    # compute S/N criteria for this peak: flux/sqrt(mean of rest of array)
    flux = imgc[xp_slice].sum()
    delta_size = img.size - imgc[xp_slice].size
    delta_flux = imgsum - flux
    if delta_size != 0.0 and delta_flux != 0.0:
        zpqual = flux/np.sqrt((imgsum - flux)/(delta_size))
    else:    
        zpqual = -1.0

    if center is not None:
        xp -= center[0]
        yp -= center[1]
    
    del imgc
    return xp,yp,flux,zpqual

def build_xy_zeropoint(imgxy,refxy,searchrad=3.0,histplot=False):
    """ Create a matrix which contains the delta between each XY position and 
        each UV position. 
    """
    print 'Computing initial guess for X and Y shifts...'
    xyshape = int(searchrad*2)+1
    zpmat = np.zeros([xyshape,xyshape],dtype=np.int32)

    for xy in imgxy:        
        deltax = xy[0] - refxy[:,0]
        deltay = xy[1] - refxy[:,1]
        xyind = np.bitwise_and(np.abs(deltax) < searchrad,np.abs(deltay) < searchrad)
        
        deltax = (deltax+searchrad)[xyind].astype(np.int32)
        deltay = (deltay+searchrad)[xyind].astype(np.int32)
        zpmat[deltay,deltax] += 1
    
    zpsum = zpmat.sum()
    xp,yp,flux,zpqual = find_xy_peak(zpmat,center=(searchrad,searchrad))
    print 'Found initial X and Y shifts of ',xp,yp,'\n     with significance of ',zpqual
    if histplot:
        zpstd = (((zpsum-flux)/zpmat.size)+zpmat.std())
        pl.clf()
        a=pl.imshow(zpmat,vmin=1,vmax=zpstd,interpolation='nearest')
        pl.gray()
        pl.colorbar()
        pl.title("Histogram of offsets: Peak S/N=%0.2f at (%0.4g, %0.4g)"%(zpqual,xp,yp))
        pl.plot(xp+searchrad,yp+searchrad,color='red',marker='+',markersize=24)
        pl.plot(searchrad,searchrad,color='yellow',marker='+',markersize=120)
        pl.text(searchrad,searchrad,"Offset=0,0",
                verticalalignment='bottom',color='yellow')
        a = raw_input("Press 'Enter' to continue...")

    return xp,yp,flux,zpqual
    
