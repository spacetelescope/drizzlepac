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
import util
from pytools import fileutil


__version__ = '1.1' #we should go through and update all these

__taskname__= "BigBlackBox.drizCR" #looks in BigBlackBox for sky.cfg
_step_num_ = 6  #this relates directly to the syntax in the cfg file

def getHelpAsString():
    """ 
    return useful help from a file in the script directory called module.help
    """
    #get the local library directory where the code is stored
    localDir=os.path.split(__file__)
    helpfile=__taskname__.split(".")
    helpfile=localDir[0]+"/"+helpfile[1]+".help"
    
    if os.access(helpfile,os.R_OK):
        fh=open(helpfile,'r')
        ss=fh.readlines()
        fh.close()
        helpString=""
        for line in ss:
            helpString+=line
    else:    
        helpString=__doc__

    return helpString

#this is the user access function
def drizCR(input=None,configObj=None, editpars=False, **inputDict):
    """
        look for cosmic rays
    """
    print inputDict    
    inputDict["input"]=input    
    configObj = util.getDefaultConfigObj(__taskname__,configObj,inputDict,loadOnly=(not editpars))

    if editpars == False:
        run(configObj)
     

#this is the function that will be called from TEAL
def run(configObj):
 
    imgObjList,outwcs = processInput.setCommonInput(configObj,createOutwcs=False) #outwcs is not neaded here
    rundrizCR(imgObjList,configObj,saveFile=configObj["clean"])
    
    
#the final function that calls the workhorse  
def rundrizCR(imgObjList,configObj,saveFile=True):

    step_name = util.getSectionName(configObj,_step_num_)
    if not configObj[step_name]['driz_cr']:
        print 'Cosmic-ray identification (driz_cr) step not performed.'
        return
        
    for image in imgObjList:
#        for chip in range(1,image._numchips,1):
 #           print image,chip
#            _drizCr(image,chip,configObj,saveFile)
        _drizCr(image,configObj,saveFile)


#the workhorse function
def _drizCr(sciImage=None,configObj={},saveFile=True):
    """mask blemishes in dithered data by comparison of an image
    with a model image and the derivative of the model image.

    sciImage is an imageObject which contains the science data
    blotImage is inferred from the sciImage object here which knows the name of its blotted image :)
    chip should be the science chip that corresponds to the blotted image that was sent
    configObj contains the user parameters
    dgMask is inferred from the sciImage object, the name of the mask file to combine with the generated Cosmic ray mask
    saveFile saves intermediate files to disk
    
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
    
    step_name = util.getSectionName(configObj,_step_num_)

    grow=configObj[step_name]["driz_cr_grow"]
    ctegrow=configObj[step_name]["driz_cr_ctegrow"]
    
            
#    try:
#        assert(chip != None), 'Please specify a chip to process for blotting'
#        assert(sciImage != None), 'Please specify a science image object for blotting'

#    except AssertionError:
#        print "Problem with value of chip or sciImage to drizCR"
#        print sciImage
#        raise AssertionError
 
    for chip in range(1,sciImage._numchips+1,1):  
        exten=sciImage.scienceExt + ',' +str(chip)    
        scienceChip=sciImage[exten]

        if scienceChip.group_member:
            blotImageName=scienceChip.outputNames["blotImage"]
            crCorImage=scienceChip.outputNames["crcorImage"]
            crMaskImage=scienceChip.outputNames["crmaskImage"]
            ctedir=scienceChip.cte_dir

            #check that sciImage and blotImage are the same size?
            try:
                os.access(blotImageName,os.F_OK)
            except IOError:
                print "Could not find the Blotted image on disk:",blotImageName
                raise IOError


            #grab the actual image from disk
            __inputImage=sciImage.getData(exten)

            try:
                __blotImage=pyfits.open(blotImageName,mode="readonly")
            except IOError:
                print "Problem opening blot images"
                return IOError

            #make the derivative blot image
            __blotData=__blotImage[0].data #simple fits
            __blotDeriv = quickDeriv.qderiv(__blotData)
            __blotImage.close()

            #this grabs the original dq mask from the science image
            # This mask needs to take into account any crbits values 
            # specified by the user to be ignored. A call to the 
            # buildMask() method may work better here...
            #__dq = sciImage.maskExt + ',' + str(chip)
            #__dqMask=sciImage.getData(__dq)
            __dqMask = sciImage.buildMask(chip,configObj['crbit'])

            #parse out the SNR information
            __SNRList=(configObj[step_name]["driz_cr_snr"]).split()
            __snr1=float(__SNRList[0])
            __snr2=float(__SNRList[1])

            #parse out the scaling information 
            __scaleList = (configObj[step_name]["driz_cr_scale"]).split()
            __mult1 = float(__scaleList[0])
            __mult2 = float(__scaleList[1])

            __gain=scienceChip._effGain
            __rn=scienceChip._rdnoise
            __backg = scienceChip.subtractedSky

            # Define output cosmic ray mask to populate
            __crMask = np.zeros(__inputImage.shape,dtype=np.uint8)

            # Determine a scaling factor depending on the units of the input image, "counts" or "cps"
            if (scienceChip.in_units== "counts"):
                __expmult = 1.
            elif(scienceChip.in_units=="cps"):
                __expmult = scienceChip._exptime
            else:
                print "drizCR found Unrecognized value in input image for BUNIT:", scienceChip.in_units
                raise ValueError

        ##################   COMPUTATION PART I    ###################
            # Create a temporary array mask
            __t1 = np.absolute(__inputImage - __blotData)
            __ta = np.sqrt(__gain * np.absolute(__blotData * __expmult + __backg * __expmult) + __rn * __rn)
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
            __xta = np.sqrt(__gain * np.absolute(__blotData * __expmult + __backg * __expmult) + __rn * __rn)
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

            if(saveFile):
                # Remove the existing cor file if it exists
                if(os.access(crCorImage, os.F_OK)):
                    os.remove(crCorImage)
                    print "Removing old corr file:",crCorImage 

                util.createFile(__corrFile,outfile=crCorImage,header=None)

            ######## Save the cosmic ray mask file to disk
            _cr_file = np.zeros(__inputImage.shape,np.uint8)
            _cr_file = np.where(__crMask,1,0).astype(np.uint8)


            #if(saveFile):
            # Always write out crmaskimage, as it is required input for 
            # the final drizzle step. The final drizzle step combines this 
            # image with the DQ information on-the-fly.
            #
            # Remove the existing mask file if it exists
            if(os.access(crMaskImage, os.F_OK)):
                os.remove(crMaskImage)
                print "Removed old cosmic ray mask file:",crMaskImage 
        
            util.createFile(_cr_file, outfile=crMaskImage, header = None)

  
            
   
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


