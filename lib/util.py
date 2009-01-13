#!/usr/bin/env python
"""
These are utility functions that can be used by any of the
classes

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
"""
    return filename.split()[0]

    
def atfile_ivm(filename):
"""
return the filename of the IVM file
which is assumed to be the second word
"""
    return filename.split()[1]    


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

""" get some better logic here"""
def isSingleFile(inputFilelist):
	if !isCommaList(inputFilelist):
    	if !isASNTable(inputFilelist):
        	if !isFilelist(inputFilelist):
            	return True
    return False

def isFilelist(inputFilelist):
"""return true if input is an @ style list"""
	if '@' in inputFilelist:
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
    

def readASNtable(filename):
"""
read the given association file (fits table)
and return a list of the science files
the asnutil function is used and
a list of filenames is returned

"""
	asndict=asnutil.readASNTable(filename, output=None, prodonly=False):	
	fileList=[]
    for member in asndict:
    	fileList.append(asndict[member]['MEMNAME'])       
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



