#!/usr/bin/env python
"""
A library of utility functions

"""
import numpy as np
import pyfits
from pytools import asnutil,fileutil
from pytools import cfgepar,cfgpars


__version__ = "0.1.0tng1"
__pyfits_version__ = pyfits.__version__
__numpy_version__ = np.__version__

def _ptime():
    import time

    # Format time values for keywords IRAF-TLM, and DATE
    _ltime = time.localtime(time.time())
    tlm_str = time.strftime('%H:%M:%S (%d/%m/%Y)',_ltime)
    #date_str = time.strftime('%Y-%m-%dT%H:%M:%S',_ltime)
    return tlm_str

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


def getDefaultConfigObj(taskname,configObj,input_dict={},loadOnly=True):
    """ Return default configObj instance for task updated 
        with user-specified values from input_dict.
        If configObj already defined, it will simply 
        return configObj unchanged. 
    """
    if configObj is None:
        # Uncomment the 'loadOnly' parameter to turn off GUI in this function.
        configObj = cfgepar.epar(taskname, loadOnly=loadOnly)
        # Merge in user-input into configObj
        #configObj.update(input_dict)
        if input_dict is not None and configObj is not None:
            mergeConfigObj(configObj,input_dict)
    return configObj

def mergeConfigObj(configObj,input_dict):
    for key in configObj:
        if isinstance(configObj[key],dict):
            keys = configObj[key].keys()
            for k in keys:
                for i in input_dict.keys():
                    if k == i:
                        configObj[key][k] = input_dict[i]
                        break
        else:
            for i in input_dict.keys():
                if i == key:
                    configObj[key] = input_dict[i]

def getSectionName(configObj,stepnum):
    """ Return section label based on step number.
    """
    for key in configObj.keys():
        if key.find('STEP '+str(stepnum)) >= 0:
            return key

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


####
#
# The following functions were required for use with the drizzling code
# and were copied in from pydrizzle_tng.py.
#
####

def countImages(imageObjectList):
    expnames = []
    for img in imageObjectList:
       expnames += img.getKeywordList('_expname')
    imgnames = []

    nimages = 0
    for e in expnames:
        if e not in imgnames:
            imgnames.append(e)
            nimages += 1
    return nimages

def get_detnum(hstwcs,filename,extnum):
    detnum = None
    binned = None
    if hstwcs.filename == filename and hstwcs.extver == extnum:
        detnum = hstwcs.chip
        binned = hstwcs.binned

    return detnum,binned

def get_exptime(header,primary_hdr):
    """shouldn't this just be defined in the instrument subclass of imageobject?"""

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

def createFile(dataArray=None, outfile=None, header=None):
    """Create a simple fits file for the given data array and header"""

    try:    
        assert(outfile != None), "Please supply an output filename for createFile"
        assert(dataArray != None), "Please supply a data array for createFiles"
    except AssertionError:
        raise AssertionError
        
    try:
        # Create the output file
        fitsobj = pyfits.HDUList()
        if (header != None):
            del(header['NAXIS1'])
            del(header['NAXIS2'])
            if header.has_key('XTENSION'):
                del(header['XTENSION'])
            if header.has_key('EXTNAME'):
                del(header['EXTNAME'])
            if header.has_key('EXTVER'):
                del(header['EXTVER'])

            if header.has_key('NEXTEND'):
                header['NEXTEND'] = 0

            hdu = pyfits.PrimaryHDU(data=dataArray,header=header)
            del hdu.header['PCOUNT']
            del hdu.header['GCOUNT']

        else:
            hdu = pyfits.PrimaryHDU(data=dataArray)

        fitsobj.append(hdu)
        fitsobj.writeto(outfile)

    finally:
        # CLOSE THE IMAGE FILES
        fitsobj.close()
        del fitsobj
