#   PROGRAM: numcombine.py
#   AUTOHOR: Christopher Hanley (based upon original code of Anton Koekemoer)
#   DATE:    February 11, 2004
#   PURPOSE: Create an image combination algroithm that chooses between a minimum
#               or median image based up bad pixel identification.
#
#   HISTORY:
#      Version 0.1.0: Initial Development -- CJH -- 02/11/04
#      Version 0.1.1: Cast boxsize as type int.  This should always be the case
#       anyways but this protects against errors in the MDRIZTAB -- CJH/WJH -- 07/08/04
#      Version 0.1.2: Improve error message handing in the case where the boxcar
#       convolution step fails.  --CJH -- 10/13/04
#      Version 0.2.0: The creation of the median image will now more closesly replicate
#       the IRAF IMCOMBINE behavior of nkeep = 1 and nhigh = 1. -- CJH -- 03/29/05
from __future__ import division # confidence medium

import numpy as np
import stsci.convolve as NC

from stsci.image.numcombine import numCombine

from .version import *

class minmed:
    """ Create a median array, rejecting the highest pixel and computing the lowest valid pixel after mask application"""

    """
        # In this case we want to calculate two things:
        #   1) the median array, rejecting the highest pixel (thus running
        #      imcombine with nlow=0, nhigh=1, nkeep=1, using the masks)
        #   2) the lowest valid pixel after applying the masks (thus running
        #      imcombine with nlow=0, nhigh=3, nkeep=1, using the masks)
        #
        # We also calculate the sum of the weight files (to produce the total
        # effective exposure time for each pixel).
        #
        # The total effective background in the final image is calculated as follows:
        #   - convert background for each input image to counts/s (divide by exptime)
        #   - multiply this value by the weight image, to obtain the effective background
        #     counts (in DN) for each pixel, for each image
        #   - Add these images together, to obtain the total effective background
        #     for the combined image.
        #
        # Once we've made these two files, then calculate the SNR based on the
        # median-pixel image, and compare with the minimum.
    """

    """ In this version of the mimmed algorithm we assume that the units of all input data is electons."""
    def __init__(self,
            imageList,              # list of input data to be combined.
            weightImageList,        # list of input data weight images to be combined.
            readnoiseList,          # list of readnoise values to use for the input images.
            exposureTimeList,       # list of exposure times to use for the input images.
            backgroundValueList,    # list of image background values to use for the input images
            weightMaskList= None,   # list of imput data weight masks to use for pixel rejection.
            combine_grow = 1,       # Radius (pixels) for neighbor rejection
            combine_nsigma1 = 4,    # Significance for accepting minimum instead of median
            combine_nsigma2 = 3,     # Significance for accepting minimum instead of median
            fillval = False         # Turn on use of imedian/imean

            ):

        # Define input variables
        self.__imageList = imageList
        self.__weightImageList = weightImageList
        self.__weightMaskList = weightMaskList
        self.__exposureTimeList = exposureTimeList
        self.__readnoiseList = readnoiseList
        self.__backgroundValueList = backgroundValueList
        self.__numberOfImages = len( self.__imageList)
        self.__combine_grow = combine_grow
        self.__combine_nsigma1 = combine_nsigma1
        self.__combine_nsigma2 = combine_nsigma2

        if fillval:
            combtype_mean = 'imean'
            combtype_median = 'imedian'
        else:
            combtype_mean = 'mean'
            combtype_median = 'median'


        # Create a different median image based upon the number of images in the input list.
        __median_file = np.zeros(self.__imageList[0].shape,dtype=self.__imageList[0].dtype)
        if (self.__numberOfImages == 2):
            __tmp = numCombine(self.__imageList,numarrayMaskList=self.__weightMaskList,
                                 combinationType=combtype_mean,nlow=0,nhigh=0,
                                 nkeep=1,upper=None,lower=None)
            __median_file = __tmp.combArrObj
        else:
            # The value of NHIGH=1 will cause problems when there is only 1 valid
            # unmasked input image for that pixel due to a difference in behavior
            # between 'numcombine' and 'iraf.imcombine'.
            # This value may need to be adjusted on the fly based on the number of
            # inputs and the number of masked values/pixel.
            #
            __tmp = numCombine(self.__imageList,numarrayMaskList=self.__weightMaskList,
                                 combinationType=combtype_median,nlow=0,nhigh=1,
                                 nkeep=1,upper=None,lower=None)
            __median_file = __tmp.combArrObj

            if self.__weightMaskList in [None,[]]:
                self.__weightMaskList = [np.zeros(self.__imageList[0].shape,dtype=self.__imageList[0].dtype)]*len(self.__imageList)
            # The following section of code will address the problem caused by having
            # a value of nhigh = 1.  This will behave in a way similar to the way the
            # IRAF task IMCOMBINE behaves.  In order to accomplish this, the following
            # procedure will be followed:
            # 1) The input masks will be summed.
            # 2) The science data will be summed.
            # 3) In the locations of the summed mask where the sum is 1 less than the
            #    the total number of images, the value of that location in the summed
            #    sceince image will be used to replace the existing value in the
            #    existing __median_file.
            #
            # This procuedure is being used to prevent too much data from being thrown
            # out of the image.  Take for example the case of 3 input images.  In two
            # of the images the pixel locations have been masked out.  Now, if nhigh
            # is applied there will be no value to use for that position.  However,
            # if this new procedure is used that value in the resulting images will
            # be the value that was rejected by the nhigh rejection step.
            #

            # We need to make certain that "bad" pixels in the sci data are set to 0.  That way,
            # when the sci images are summed, the value of the sum will only come from the "good"
            # pixels.
            tmpList = []
            for image in xrange(len(self.__imageList)):
                tmp =np.where(self.__weightMaskList[image] == 1, 0, self.__imageList[image])
                tmpList.append(tmp)

            # Sum the mask files
            maskSum = self.__sumImages(self.__weightMaskList)
            # Sum the science images
            sciSum = self.__sumImages(tmpList)
            del(tmpList)
            # Use the summed sci image values in locations where the maskSum indicates
            # that there is only 1 good pixel to use.  The value will be used in the
            # __median_file image
            __median_file = np.where(maskSum == self.__numberOfImages-1,sciSum,__median_file)

        if self.__weightMaskList in [None,[]]:
            self.__weightMaskList = [np.zeros(self.__imageList[0].shape,dtype=self.__imageList[0].dtype)]*len(self.__imageList)
        # Sum the weightMaskList elements
        __maskSum = self.__sumImages(self.__weightMaskList)

        # Create the minimum image from the stack of input images.
        # Find the maximum pixel value for the image stack.
        _maxValue = -1e+9
        for image in self.__imageList:
            _newMax = image.max()
            if (_newMax > _maxValue):
                _maxValue = _newMax

        # For each image, set pixels masked as "bad" to the "super-maximum" value.
        for image in xrange(len(self.__imageList)):
            self.__imageList[image] = np.where(self.__weightMaskList[image] == 1,_maxValue+1,self.__imageList[image])

        # Call numcombine throwing out the highest N - 1 pixels.
        __tmp = numCombine(self.__imageList,numarrayMaskList=None,
                                 combinationType=combtype_median,nlow=0,nhigh=self.__numberOfImages-1,
                                 nkeep=1,upper=None,lower=None)
        __minimum_file = __tmp.combArrObj
        # Reset any pixl at _maxValue + 1 to 0.
        __minimum_file = np.where(__maskSum == self.__numberOfImages, 0, __minimum_file)

        # Scale the weight images by the background values and add them to the bk
        __backgroundFileList = []
        for image in xrange(len(self.__weightImageList)):
            __tmp = self.__weightImageList[image] * (self.__backgroundValueList[image]/(self.__exposureTimeList[image]))
            __backgroundFileList.append(__tmp)

        # Create an image of the total effective background (in DN) per pixel:
        # (which is the sum of all the background-scaled weight files)
        #
        __bkgd_file = self.__sumImages(__backgroundFileList)
        del(__backgroundFileList)

        #
        # Scale the weight mask images by the square of the readnoise values
        #
        __readnoiseFileList = []
        for image in xrange(len(self.__weightMaskList)):
            __tmp = np.logical_not(self.__weightMaskList[image]) * (self.__readnoiseList[image] * self.__readnoiseList[image])
            __readnoiseFileList.append(__tmp)

        # Create an image of the total readnoise**2 per pixel:
        # (which is the sum of all the input readnoise values)
        #
        __readnoise_file = self.__sumImages(__readnoiseFileList)
        del(__readnoiseFileList)

        # Create an image of the total effective exposure time per pixel:
        # (which is simply the sum of all the drizzle output weight files)
        #
        __weight_file = self.__sumImages(self.__weightImageList)


        # Scale up both the median and minimum arrays by the total effective exposure time
        # per pixel.
        #
        __minimum_file_weighted = __minimum_file * __weight_file
        __median_file_weighted = __median_file * __weight_file
        del(__weight_file)

        # Calculate the 1-sigma r.m.s.:
        #   variance = median_electrons + bkgd_electrons + readnoise**2
        #   rms = sqrt(variance)
        #   This image has units of electrons.
        #
        # make this the abs value so that negative numbers dont throw an exception?
        __rms_file = np.sqrt(__median_file_weighted + __bkgd_file + __readnoise_file)

        del __bkgd_file, __readnoise_file
        # For the median array, calculate the n-sigma lower threshold to the array
        # and incorporate that into the pixel values.
        #
        __median_rms_file = __median_file_weighted - (__rms_file * self.__combine_nsigma1)

        if self.__combine_grow != 0:
            #
            # Do a more sophisticated rejection: For all cases where the minimum pixel will
            # be accepted instead of the median, set a lower threshold for that pixel and the
            # ones around it (ie become less conservative in rejecting the median). This is
            # because in cases of triple-incidence cosmic rays, quite often the low-lying
            # outliers of the CRs can influence the median for the initial relatively high
            # value of sigma, so a lower threshold must be used to mnake sure that the minimum
            # is selected.
            #
            # This is done as follows:
            # 1) make an image which is zero everywhere except where the minimum will be accepted
            # 2) box-car smooth this image, to make these regions grow.
            # 3) In the file "median_rms_file_electrons", replace these pixels
            #     by median - combine_nsigma2 * rms
            #
            # Then use this image in the final replacement, in the same way as for the
            # case where this option is not selected.

            __minimum_flag_file = np.where(np.less(__minimum_file_weighted,__median_rms_file),1,0)

            # The box size value must be an integer.  This is not a problem since __combine_grow should always
            # be an integer type.  The combine_grow column in the MDRIZTAB should also be an integer type.
            __boxsize = int(2 * self.__combine_grow + 1)
            __boxshape = (__boxsize,__boxsize)
            __minimum_grow_file = np.zeros(self.__imageList[0].shape,dtype=self.__imageList[0].dtype)


            # If the boxcar convolution has failed it is potentially for two reasons:
            #   1) The kernel size for the boxcar is bigger than the actual image.
            #   2) The grow parameter was specified with a value < 0.  This would result
            #      in an illegal boxshape kernel.  The dimensions of the kernel box *MUST*
            #      be integer and greater than zero.
            #
            #   If the boxcar convolution has failed, try to give a meaningfull explanation
            #   as to why based upon the conditionals described above.

            if (__boxsize <= 0):
                errormsg1 =  "############################################################\n"
                errormsg1 += "# The boxcar convolution in minmed has failed.  The 'grow' #\n"
                errormsg1 += "# parameter must be greater than or equal to zero. You     #\n"
                errormsg1 += "# specified an input value for the 'grow' parameter of:    #\n"
                errormsg1 += "        combine_grow: " + str(self.__combine_grow)+'\n'
                errormsg1 += "############################################################\n"
                raise ValueError,errormsg1
            if (__boxsize > self.__imageList[0].shape[0]):
                errormsg2 =  "############################################################\n"
                errormsg2 += "# The boxcar convolution in minmed has failed.  The 'grow' #\n"
                errormsg2 += "# parameter specified has resulted in a boxcar kernel that #\n"
                errormsg2 += "# has dimensions larger than the actual image.  You        #\n"
                errormsg2 += "# specified an input value for the 'grow' parameter of:    #\n"
                errormsg2 += "        combine_grow: " +str(self.__combine_grow)+'\n'
                errormsg2 += "############################################################\n"
                print self.__imageList[0].shape
                raise ValueError,errormsg2

            # Attempt the boxcar convolution using the boxshape based upon the user input value of "grow"
            NC.boxcar(__minimum_flag_file,__boxshape,output=__minimum_grow_file,mode='constant',cval=0)

            del(__minimum_flag_file)

            __temp1 = (__median_file_weighted - (__rms_file * self.__combine_nsigma1))
            __temp2 = (__median_file_weighted - (__rms_file * self.__combine_nsigma2))
            __median_rms2_file = np.where(np.equal(__minimum_grow_file,0),__temp1,__temp2)
            del(__temp1)
            del(__temp2)
            del(__rms_file)
            del(__minimum_grow_file)

            # Finally decide whether to use the minimim or the median (in counts/s),
            # based on whether the median is more than 3 sigma above the minimum.
            #
            self.combArrObj = __tmp
            self.combArrObj = np.where(np.less(__minimum_file_weighted,__median_rms2_file),
                                        __minimum_file,
                                        __median_file)

        else:
            # Finally decide whether to use the minimim or the median (in counts/s),
            # based on whether the median is more than 3 sigma above the minimum.
            #
            self.combArrObj = __tmp
            self.combArrObj = np.where(np.less(__minimum_file_weighted,__median_rms_file),
                                        __minimum_file,
                                        __median_file)

        # Set fill regions to a pixel value of 0.
        self.combArrObj = np.where(__maskSum == self.__numberOfImages, 0, self.combArrObj)

#        self.out_file1 = __median_rms2_file.copy()
#        self.out_file2 = __minimum_file_weighted.copy()



    def __sumImages(self,numarrayObjectList):
        """ Sum a list of numarray objects. """
        if numarrayObjectList in [None,[]]:
            return None
        __sum = np.zeros(numarrayObjectList[0].shape,dtype=numarrayObjectList[0].dtype)
        for image in numarrayObjectList:
            __sum += image
        return __sum
