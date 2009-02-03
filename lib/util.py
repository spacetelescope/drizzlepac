#!/usr/bin/env python
"""
A library of utility functions

"""
from pytools import asnutil

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
