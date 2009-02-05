#!/usr/bin/env python
"""
A library of utility functions

"""
import numpy as np
from pytools import asnutil,fileutil
import buildmask


def findrootname(filename):
    """
    return the rootname of the given file
    """

    puncloc = [filename.find(char) for char in string.punctuation]
    val = sys.maxint
    for num in puncloc:
        if num !=-1 and num < val:
            val = num
    return filename[0:val]


"""
These two functions are for reading in an "at" which contains
two columns of filenames, the first column is assumed to
be the science image and the second column is assumed
to be the IVM file that is associated with it
"""

def atfile_sci(filename):
    """
    return the filename of the science image 
    which is assumed to be the first word
    in the atfile the user gave
    """
    return filename.split()[0]

    
def atfile_ivm(filename):
    """
    return the filename of the IVM file
    which is assumed to be the second word
    in the atfile the user gave
    """
    return filename.split()[1]    
    
    
def printParams(paramDictionary):
    """ Print nicely the parameters from the dictionary
    """

    if (len(paramDictionary) == 0):
        print "\nNo parameters were supplied\n"
    else:
        keys=paramDictionary.keys()
        keys.sort()
        for key in keys:
            print key,":\t",paramDictionary[key]


def isASNTable(inputFilelist):
    """return TRUE if inputFilelist is a fits ASN file"""
    if ("_asn"  or "_asc") in inputFilelist:
        return True
    return False

def isCommaList(inputFilelist):
    """return True if the input is a comma separated list of names"""
    if "," in inputFilelist:
        return True
    return False        
  
def loadFileList(inputFilelist):
    """open up the @ file and read in the science and possible
      ivm filenames from the first two columns
    """
    f = open(inputFilelist[1:])
    # check the first line in order to determine whether
    # IVM files have been specified in a second column...
    lines = f.readline()
    f.close()
    
    # If there is a second column...
    if len(line.split()) == 2:
    	# ...parse out the names of the IVM files as well 
    	ivmlist = irafglob.irafglob(input, atfile=atfile_ivm) 
    
    # Parse the @-file with irafglob to extract the input filename
    filelist = irafglob.irafglob(input, atfile=atfile_sci)
    return filelist


def readCommaList(fileList):
    """ return a list of the files with the commas removed """
    names=fileList.split(',')
    fileList=[]
    for item in names:
        fileList.append(item)
    return fileList
    


def getInputAsList(input, output=None, ivmlist=None, prodonly=False):
    if (isinstance(input, list) == False) and \
       ('_asn' in input or '_asc' in input) :
       
        # Input is an association table
        # Get the input files, and run makewcs on them
        oldasndict = asnutil.readASNTable(input, prodonly=prodonly)

        if not output:
            output = oldasndict['output']

        filelist = [fileutil.buildRootname(fname) for fname in oldasndict['order']]
        
    elif (isinstance(input, list) == False) and \
       (input[0] == '@') :
        # input is an @ file
        f = open(input[1:])
        # Read the first line in order to determine whether
        # IVM files have been specified in a second column...
        line = f.readline()
        f.close()
        # Parse the @-file with irafglob to extract the input filename
        filelist = irafglob.irafglob(input, atfile=atfile_sci)
        # If there is a second column...
        if len(line.split()) == 2:
            # ...parse out the names of the IVM files as well 
            ivmlist = irafglob.irafglob(input, atfile=atfile_ivm)        
    else:
        #input is a string or a python list
        try:
            filelist, output = parseinput.parseinput(input, outputname=output)
            #filelist.sort()
        except IOError: raise
        
    return filelist, output

def runmakewcs(input):
    """
    Runs make wcs and recomputes the WCS keywords
    input: a list of files
    output: returns a list of names of the modified files
            (For GEIS files returns the translated names.)
    
    """
    newNames = updatewcs.updatewcs(input, checkfiles=False)
    
    return newNames


def update_input(filelist, ivmlist=None, removed_files=None):
    """
    Removes files flagged to be removed from the input filelist.
    Removes the corresponding ivm files if present.
    """
    newfilelist = []

    if removed_files == []:
        return filelist, ivmlist
    else:
        sci_ivm = zip(filelist, ivmlist)
        for f in removed_files:
            result=[sci_ivm.remove(t) for t in sci_ivm if t[0] == f ]
        ivmlist = [el[1] for el in sci_ivm] 
        newfilelist = [el[0] for el in sci_ivm] 
        return newfilelist, ivmlist 



def setOutputNames(filename):
    """
    Define the default output filenames for drizzle products,
    these are based on the original rootname of the image 
    
    filename should be just 1 filename, so call this in a loop
    for chip names contained inside a file
    
    """

    # Define FITS output filenames for intermediate products
    outFinal = filename+'_drz.fits'
    outSci = filename+'_sci.fits'
    outWeight = filename+'_weight.fits'
    outContext = filename+'_context.fits'
    outSky = filename + '_sky.fits'
    blotImage = filename + '_blot.fits'
    crImage = filename + '_cr.fits'
    outSingle = filename+'_single.fits'
    outSWeight = filename+'_wht.fits'
    outSContext = None
        
    fnames={'outFinal':outFinal,
            'outSci':outSci,
            'outWeight':outWeight,
            'outContext':outContext,
            'outSingle':outSingle,
            'outSWeight':outSWeight,
            'outSContext':outSContext,
            'blotImage':blotImage,
            'crImage':crImage,
            'outSingle':outSingle,
            'outSWeight':outSWeight,
            'outSContext':outSContext,
            'outSky':outSky}

    return fnames

####
#
# The following functions were required for use with the drizzling code
# and were copied in from pydrizzle_tng.py.
#
####

def build_mask(dqname,mask_name,detnum,dq_extn,bits,extver,instrument):
    '''Build masks for a given chip. This can be called for either single_mask
       or final_mask creation.    
    '''
    # The WFPC2 specific logic will need to be pulled in as
    # an InputImage class method.
    if bits != None:
        print 'building mask for ',dq_extn,',',extver,' for output ',mask_name
        if instrument != 'WFPC2':
            out_mask = buildmask.buildMaskImage(dqname,bits,mask_name,extname=dq_extn,extver=extver)
        else:
            out_mask = buildmask.buildShadowMaskImage(dqname,detnum,extver,mask_name, bitvalue=bits, binned=binned)

def get_detnum(hstwcs,filename,extnum):
    detnum = None
    binned = None
    if hstwcs.filename == filename and hstwcs.extver == extnum:
        detnum = hstwcs.chip
        binned = hstwcs.binned

    return detnum,binned

def get_exptime(header,primary_hdr):

    if primary_hdr.has_key('exptime'):
        exphdr = primary_hdr
    else:
        exphdr = header
        
    exptime = float(exphdr['EXPTIME'])
    if exptime == 0.: exptime = 1.0

    if exphdr.has_key('EXPSTART'):
        expstart = float(exphdr['EXPSTART'])
        expend = float(exphdr['EXPEND'])
    else:
        expstart = 0.
        expend = exptime

    return (exptime,expstart,expend)

def compute_texptime(imageObjectList):
    """
    Add up the exposure time for all the members in
    the pattern, since 'drizzle' doesn't have the necessary
    information to correctly set this itself.
    """
    expnames = []
    exptimes = []
    start = []
    end = []
    for img in imageObjectList:
        expnames += img.getKeywordList('_expname')
        exptimes += img.getKeywordList('_exptime')
        start += img.getKeywordList('_expstart')
        end += img.getKeywordList('_expend')
    
    exptime = 0.
    expstart = min(start)
    expend = max(end)
    exposure = None
    for n in range(len(expnames)):
        if expnames[n] != exposure:
            exposure = expnames[n]
            exptime += exptimes[n]

    return (exptime,expstart,expend)

def computeRange(corners):
    """ Determine the range spanned by an array of pixel positions. """
    _xrange = (np.minimum.reduce(corners[:,0]),np.maximum.reduce(corners[:,0]))
    _yrange = (np.minimum.reduce(corners[:,1]),np.maximum.reduce(corners[:,1]))
    return _xrange,_yrange

def getRotatedSize(corners,angle):
    """ Determine the size of a rotated (meta)image."""
    # If there is no rotation, simply return original values
    if angle == 0.:
        _corners = corners
    else:
        # Find center
        #_xr,_yr = computeRange(corners)
        #_cen = ( ((_xr[1] - _xr[0])/2.)+_xr[0],((_yr[1]-_yr[0])/2.)+_yr[0])
        _rotm = fileutil.buildRotMatrix(angle)
        # Rotate about the center
        #_corners = N.dot(corners - _cen,_rotm)
        _corners = np.dot(corners,_rotm)

    return computeRange(_corners)
