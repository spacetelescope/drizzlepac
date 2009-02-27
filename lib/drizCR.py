# DRIZ_CR  -- mask blemishes in dithered data by comparison of an image
#             with a model image and the derivative of the model image.
#
#
# Import external packages
import numpy as np
import convolve as NC
import pyfits
import os
import quickDeriv

__version__ = '1.1' #we should go through and update all these

__taskname__= "BigBlackBox.drizCR" #looks in BigBlackBox for sky.cfg
_step_num_ = 3  #this relates directly to the syntax in the cfg file

def getHelpAsString():
    """ I'm thinking we could just make a file called sky.help
    then use this function to read it into an array or list and return that?
    """
    helpString="Help string for drizCR will be here"

    return helpString

#this is the user access function
def drizCR(imageList=None,configObj=None, editpars=False, **inputDict):
    """
        create a median image from the seperately drizzled images   
    """
    inputDict["input"]=imageList        
    configObj = util.getDefaultConfigObj(__taskname__,configObj,inputDict,loadOnly=loadOnly(not editpars))
    if configObj is None:
        return

    run(configObj)
     

#this is the function that will be called from TEAL
def run(configObj):
 
    imgObjList,outwcs = processInput.setCommonInput(configObj,createOutwcs=False) #outwcs is not neaded here
    rundrizCR(imgObjList,configObj,saveFile=configObj["clean"])
    
    
#the final function that calls the workhorse  
def rundrizCR(imgObjList,configObj,saveFile=True):

    for image in imgObjList:
        for chip in range(1,image._numchips,1):
            _drizCr(image,chip,configObj,saveFile)



#the workhorse function
def _drizCr(sciImage=None,chip=None,configObj={},saveFiles=True):
    """mask blemishes in dithered data by comparison of an image
    with a model image and the derivative of the model image.

    sciImage is an imageObject which contains the science data
    blotImage is inferred from the sciImage object here which knows the name of its blotted image :)
    chip should be the science chip that corresponds to the blotted image that was sent
    configObj contains the user parameters
    dgMask is inferred from the sciImage object, the name of the mask file to combine with the generated Cosmic ray mask
    saveFiles saves intermediate files to disk
    
    here are the options you can override in configObj

    gain     = 7               # Detector gain, e-/ADU
    grow     = 1               # Radius around CR pixel to mask [default=1 for 3x3 for non-NICMOS]   
    ctegrow  = 0               # Length of CTE correction to be applied
    rn       = 5               # Read noise in electrons
    snr      = "4.0 3.0"       # Signal-to-noise ratio
    scale    = "0.5 0.4"       # scaling factor applied to the derivative
    backg    = 0              # Background value
    expkey   = "exptime"        # exposure time keyword
    
    blot images are saved out to simple fits files with 1 chip in them
    so for example in ACS, there will be 1 image file with 2 chips that is
    the original image and 2 blotted image files, each with 1 chip
    
    so I'm imagining calling this function twice, once for each chip,
    but both times with the same original science image file, output files
    and some input (output from previous steps) are referenced in the imageobject
    itself
    
    """    
    
    paramDict=setDefaults(configObj)
    grow=paramDict["grow"]
    ctegrow=paramDict["ctegrow"]
    
    scienceChip=sciImage[sciImage.scienceExt,chip]
    
    try:
        assert(chip != None), 'Please specify a chip to process for blotting'
        assert(sciImage !=None), 'Please specify a science image object for blotting'

    except AssertionError:
        raise AssertionError
        
    blotImageName=scienceChip.outputNames["blotImage"]
    blotDerivName=scienceChip.outputNames["blotDeriv"]
    crCorImage=scienceChip.outputNames["crcorImage"]
    crMaskImage=scienceChip.outputNames["crmaskImage"]
    ctedir=scienceChip.cte_dir
    
    #check that sciImage and blotImage are the same size?
    try:
        fileutil.checkFileExists(blotImageName)
    except IOError:
        print "Could not find the Blotted image on disk:",blotImageName
        raise IOError
   
  
    #grab the actual images from disk
    __inputImage=sciImage.getData(chip)
    
    try:
        __blotImage=fileutil.openImage(blotImageName,mode='readonly',writefits=False,memmap=0)
        __blotData=__blotImage.data
    except IOError:
        print "Problem opening blot images"
        return IOError
    
    #make the derivative blot image
    __blotDeriv = quickDeriv.qderiv(__blotData)
     
    
    #this grabs the original dq mask from the science image
    __dq = sciImage.maskExt + ',' + str(chip)
    __dqMask=sciImage.getData(exten=__dq)
    
    #parse out the SNR information
    __SNRList=(paramDict["snr"]).split()
    __snr1=float(__SNRList[0])
    __snr2=float(__SNRList[1])
    
    #parse out the scaling information 
    __scaleList = (paramDict["scale"]).split()
    __mult1 = float(__scaleList[0])
    __mult2 = float(__scaleList[1])
    
    __gain=scienceChip._gain
    __rn=scienceChip._rdnoise

    # Define output cosmic ray mask to populate
    __crMask = np.zeros(__inputImage.shape,dtype=np.uint8)

    # Determine a scaling factor depending on the units of the input image, "counts" or "cps"
    if (scienceChip._bunit== "counts"):
        __expmult = 1.
    elif(scienceChip._bunit=="cps"):
        __expmult = scienceChip._exptime
    else:
        print "drizCR found Unrecognized value for BUNIT:", scienceChip._bunit
        raise ValueError

##################   COMPUTATION PART I    ###################
    # Create a temporary array mask
    __t1 = np.absolute(__inputImage - __blotData)
    __ta = np.sqrt(__gain * np.absolute(__blotData * __expmult + paramDict["backg"] * __expmult) + __rn * __rn)
    __tb = ( __mult1 * __blotDeriv + __snr1 * __ta / __gain )
    del __ta
    __t2 = __tb / __expmult
    del __tb
    __tmp1 = np.logical_not(np.greater(__t1, __t2))
    del __t1
    del __t2

    # Create a convolution kernel that is 3 x 3 of 1's
    __kernel = np.ones((3,3),dtype=np.uint8)
    # Create an output tmp file the same size as the input temp mask array
    __tmp2 = np.zeros(__tmp1.shape,dtype=np.int16)
    # Convolve the mask with the kernel
    NC.convolve2d(__tmp1,__kernel,output=__tmp2,fft=0,mode='nearest',cval=0)
    del __kernel
    del __tmp1

##################   COMPUTATION PART II    ###################
    # Create the CR Mask
    __xt1 = np.absolute(__inputImage - __blotData)
    __xta = np.sqrt(__gain * np.absolute(__blotData * __expmult + paramDict["backg"] * __expmult) + __rn * __rn)
    __xtb = ( __mult2 *__blotDeriv + __snr2 * __xta / __gain )
    del __xta
    __xt2 = __xtb / __expmult
    del __xtb
    # It is necessary to use a bitwise 'and' to create the mask with numarray objects.
    __crMask = np.logical_not(np.greater(__xt1, __xt2) & np.less(__tmp2,9) )
    del __xt1
    del __xt2
    del __tmp2


##################   COMPUTATION PART III    ###################        
#flag additional cte 'radial' and 'tail' pixels surrounding CR pixels as CRs

    # In both the 'radial' and 'length' kernels below, 0->good and 1->bad, so that upon
    # convolving the kernels with __crMask, the convolution output will have low->bad and high->good 
    # from which 2 new arrays are created having 0->bad and 1->good. These 2 new arrays are then 'anded'
    # to create a new __crMask.

    # recast __crMask to int for manipulations below; will recast to Bool at end
    __crMask_orig_bool= __crMask.copy() 
    __crMask= __crMask_orig_bool.astype( np.int8 )

    # make radial convolution kernel and convolve it with original __crMask 
    cr_grow_kernel = np.ones((grow, grow))     # kernel for radial masking of CR pixel
    cr_grow_kernel_conv = __crMask.copy()   # for output of convolution
    NC.convolve2d( __crMask, cr_grow_kernel, output = cr_grow_kernel_conv)

    # make tail convolution kernel and convolve it with original __crMask
    cr_ctegrow_kernel = np.zeros((2*ctegrow+1,2*ctegrow+1))  # kernel for tail masking of CR pixel
    cr_ctegrow_kernel_conv = __crMask.copy()  # for output convolution 

    # which pixels are masked by tail kernel depends on sign of ctedir (i.e.,readout direction):
    if ( ctedir == 1 ):  # HRC: amp C or D ; WFC: chip = sci,1 ; WFPC2
        cr_ctegrow_kernel[ 0:ctegrow, ctegrow ]=1    #  'positive' direction
    if ( ctedir == -1 ): # HRC: amp A or B ; WFC: chip = sci,2
        cr_ctegrow_kernel[ ctegrow+1:2*ctegrow+1, ctegrow ]=1    #'negative' direction
    if ( ctedir == 0 ):  # NICMOS: no cte tail correction
        pass

    # do the convolution
    NC.convolve2d( __crMask, cr_ctegrow_kernel, output = cr_ctegrow_kernel_conv)    

    # select high pixels from both convolution outputs; then 'and' them to create new __crMask
    where_cr_grow_kernel_conv    = np.where( cr_grow_kernel_conv < grow*grow,0,1 )        # radial
    where_cr_ctegrow_kernel_conv = np.where( cr_ctegrow_kernel_conv < ctegrow, 0, 1 )     # length
    __crMask = np.logical_and( where_cr_ctegrow_kernel_conv, where_cr_grow_kernel_conv) # combine masks

    __crMask = __crMask.astype(np.uint8) # cast back to Bool

    del __crMask_orig_bool
    del cr_grow_kernel 
    del cr_grow_kernel_conv 
    del cr_ctegrow_kernel 
    del cr_ctegrow_kernel_conv
    del where_cr_grow_kernel_conv  
    del where_cr_ctegrow_kernel_conv 


    # Apply CR mask to the DQ array in place
    np.bitwise_and(__dqMask,__crMask,__dqMask)

    ####### Create the corr file
    __corrFile = np.zeros(__inputImage.shape,dtype=__inputImage.dtype)
    __corrFile = np.where(np.equal(__dqMask,0),__blotData,__inputImage)

    if(saveFiles):
        # Remove the existing cor file if it exists
        if(os.access(crcorimage, os.F_OK)):
            os.remove(crcorimage)
            print "Removing old corr file:",corrName 

        createFile(__corrFile,outfile=crcorimage,header=None)
    
    ######## Save the cosmic ray mask file to disk
    _cr_file = np.zeros(__inputImage.shape,np.uint8)
    _cr_file = np.where(__crMask,1,0).astype(np.uint8)

    
    if(saveFiles):
        # Remove the existing mask file if it exists
        if(os.access(crmaskimage, os.F_OK)):
            os.remove(crmaskimage)
            print "Removed old cosmic ray mask file:",crName 

        createFile(_cr_file, outfile=crmaskimage, header = None)
   
  
def createFile(dataArray=None, outfile=None, header=None):
    """Create a simple fits file for the give data array and header"""

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
            hdu = pyfits.PrimaryHDU(data=__corrFile)

        fitsobj.append(hdu)
        fitsobj.writeto(outfile)

    finally:
        # CLOSE THE IMAGE FILES
        fitsobj.close()
        del fitsobj,__corrFile
            
   
def setDefaults(configObj={}):
    """return a dictionary of the default parameters
        which also been updated with the user overrides
    """
    gain     = 7               # Detector gain, e-/ADU
    grow     = 1               # Radius around CR pixel to mask [default=1 for 3x3 for non-NICMOS]   
    ctegrow  = 0               # Length of CTE correction to be applied
    rn       = 5               # Read noise in electrons
    snr      = "4.0 3.0"       # Signal-to-noise ratio
    scale    = "0.5 0.4"       # scaling factor applied to the derivative
    backg    = 0              # Background value
    expkey   = "exptime"        # exposure time keyword
    
    paramDict={"gain":gain,
                "grow": grow,
                "ctegrow":ctegrow,
                "rn":rn,
                "snr":snr,
                "scale":scale,
                "backg":backg,
                "expkey":expkey}
    
    if (len(configObj) != 0):
        for key in configObj:
            paramDict[key]=configObj[key]
            
            
    return paramDict


