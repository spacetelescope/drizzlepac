#
#   Authors:    Christopher Hanley, Warren Hack, Megan Sosey
#   Program:    static_mask.py
#   Purpose:    Class that manages the creation of a global static
#               mask which is used to mask pixels that are some
#               sigma BELOW the mode computed for the image.

import numpy as np
from imagestats import ImageStats



def createStaticMask(imageSet,paramDict={},saveFile=True):
    """
    Manage the creation of the global static mask which
    masks pixels which are negative in the SCI array.
    A static mask numarray object gets created for each global
    mask needed, one for each chip from each instrument/detector.
    Each static mask array has type Int16, and resides in memory.

    The parameter 'badval' defaults to 64 and represents the
    DQ value used to mark pixels flagged as bad by this mask.

    The parameter 'goodval' defaults to 1.0 and represents the
    pixel value of the good pixels.

    """



def myCreateStaticMask(imageList=[],configObj={}):
    """user callable routine that accepts varied input"""

    #imageList here is assumed to be a python list of filenames
    #though I think this might be part of configObj too
    if len(imageList) == 0:
        print "Empty input image list given to static mask routine, checking configObj"
        if(len(configObj["imageList"]) == 0):
            print "no image list in configObject either"
            return ValueError
        else:
            imageList = configObj["imageList"]
            
    #make up a dictionary of the task parameter values
    paramDict=_setDefaults()
    
    #if configobj has been provided, then use those parameters instead
    if (len(configObj) != 0):
        #now replace the defaults with the user specs
        #this assumes only matching keys, though with this
        #syntax configObj could contain more than the standard keys
        #and they will get appended onto the paramDict
        for key in configObj:
            paramDict[key]=configObj[key]

    # Print out the parameters provided by the user
    print "\nUSER INPUT PARAMETERS for Static Mask:"
    util.printParams(paramDict)        

    #create image object    
    for image in imageList:
        imageSet=imageObject(image)
        print "\nWorking on: ",imageSet._filename
        imageSet.info()
        #call the real static mask routine
        createStaticMask(imageSet,paramDict,saveFile)
        imageSet.close()
    
    
def __init__ (self, badval=64, goodval=1.0, staticsig= 3.0):

    # For now, we don't use badval. It is supposed to
    # be used to flag back the DQ array of the input
    # images. This may lead to confusion when running
    # the task repeatedly over a set of images, whenever
    # additional images are included in the set each
    # time.

    self.static_sig = staticsig
    self.static_badval = badval
    self.static_goodval = goodval

    self.masklist = {}

def addMember(image,sci_arr,signature):
    """
    Combines the input image with the static mask that
    has the same signature.  The signature parameter
    consists of the tuple:
    (instrument/detector, (nx,ny), chip_id)
    """
    # If this is a new signature, create a new Static Mask file
    if not self.masklist.has_key(signature):
        self.masklist[signature] = self._buildMaskArray(signature)

    # Operate on input image DQ array to flag 'bad' pixels in the
    # global static mask
    stats = ImageStats(sci_arr,nclip=3,fields='mode')
    mode = stats.mode
    rms  = stats.stddev
    del stats

    print('  mode = %9f;   rms = %7f')  %  (mode,rms)
    #
    # The scale value (3.0) could potentially become a useful 
    # user settable parameter.  Consider for future revisions.
    # 29-April-2004 WJH/CJH/AMK
    #
    sky_rms_diff = mode - (self.static_sig*rms)

    np.bitwise_and(self.masklist[signature],np.logical_not(np.less( sci_arr, sky_rms_diff)),self.masklist[signature])

def _buildMaskArray(self,nx,ny):
    """ Creates numarray array for static mask array signature. """
    return np.ones([nx,ny],dtype=np.int16)

def getMask(image):
    """ Returns the appropriate StaticMask array for the image. """
    return image



#set up the default values for the task, this is still too clunky
#but perhaps could eventually be used to create the default configObj file?
def _setDefaults():

    badval=64
    goodval=1
    sigma=4.0

    paramDict={"goodval":goodval,"badval":badval,"sigma":sigma}

    return paramDict
