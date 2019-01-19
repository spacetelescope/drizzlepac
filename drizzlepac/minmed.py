"""
Create an image combination algroithm that chooses between a minimum
or median image based up bad pixel identification.

:Authors: Christopher Hanley

:License: :doc:`LICENSE`

"""
#   PROGRAM: numcombine.py
#   AUTHOR: Christopher Hanley (based upon original code of Anton Koekemoer)
#   DATE:    February 11, 2004
#   PURPOSE: Create an image combination algroithm that chooses between
#            a minimum or median image based up bad pixel identification.
#   HISTORY:
#      Version 0.1.0: Initial Development -- CJH -- 02/11/04
#      Version 0.1.1: Cast boxsize as type int.  This should always be the case
#        anyways but this protects against errors in the
#        MDRIZTAB -- CJH/WJH -- 07/08/04
#      Version 0.1.2: Improve error message handing in the case where
#        the boxcar convolution step fails.  --CJH -- 10/13/04
#      Version 0.2.0: The creation of the median image will now more closesly
#        replicate the IRAF IMCOMBINE behavior of nkeep = 1 and nhigh = 1.
#        -- CJH -- 03/29/05
#      Version 0.3.0: Rewritten to optimize the code and to bring
#        code up to modern standards.-- Mihai Cara -- 02/19/2018
import warnings
import numpy as np
from scipy import signal
from stsci.image.numcombine import numCombine, num_combine
from .version import *

class minmed:
    """ **DEPRECATED** Create a median array, rejecting the highest pixel and
    computing the lowest valid pixel after mask application

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

    In this version of the mimmed algorithm we assume that the units of all
    input data is electons.
    """
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

        warnings.warn("The 'minmed' class is deprecated and may be removed"
                      " in a future version. Use 'min_med()' instead.",
                      DeprecationWarning)

        # Define input variables
        self._imageList = imageList
        self._weightImageList = weightImageList
        self._weightMaskList = weightMaskList
        self._exposureTimeList = exposureTimeList
        self._readnoiseList = readnoiseList
        self._backgroundValueList = backgroundValueList
        self._numberOfImages = len( self._imageList)
        self._combine_grow = combine_grow
        self._combine_nsigma1 = combine_nsigma1
        self._combine_nsigma2 = combine_nsigma2

        if fillval:
            combtype_mean = 'imean'
            combtype_median = 'imedian'
        else:
            combtype_mean = 'mean'
            combtype_median = 'median'


        # Create a different median image based upon the number of images in the input list.
        median_file = np.zeros(self._imageList[0].shape,dtype=self._imageList[0].dtype)
        if (self._numberOfImages == 2):
            tmp = numCombine(self._imageList,numarrayMaskList=self._weightMaskList,
                                 combinationType=combtype_mean,nlow=0,nhigh=0,
                                 nkeep=1,upper=None,lower=None)
            median_file = tmp.combArrObj
        else:
            # The value of NHIGH=1 will cause problems when there is only 1 valid
            # unmasked input image for that pixel due to a difference in behavior
            # between 'numcombine' and 'iraf.imcombine'.
            # This value may need to be adjusted on the fly based on the number of
            # inputs and the number of masked values/pixel.
            #
            tmp = numCombine(self._imageList,numarrayMaskList=self._weightMaskList,
                                 combinationType=combtype_median,nlow=0,nhigh=1,
                                 nkeep=1,upper=None,lower=None)
            median_file = tmp.combArrObj

            if self._weightMaskList in [None,[]]:
                self._weightMaskList = [np.zeros(self._imageList[0].shape,dtype=self._imageList[0].dtype)]*len(self._imageList)
            # The following section of code will address the problem caused by having
            # a value of nhigh = 1.  This will behave in a way similar to the way the
            # IRAF task IMCOMBINE behaves.  In order to accomplish this, the following
            # procedure will be followed:
            # 1) The input masks will be summed.
            # 2) The science data will be summed.
            # 3) In the locations of the summed mask where the sum is 1 less than the
            #    the total number of images, the value of that location in the summed
            #    sceince image will be used to replace the existing value in the
            #    existing median_file.
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
            for image in range(len(self._imageList)):
                tmp =np.where(self._weightMaskList[image] == 1, 0, self._imageList[image])
                tmpList.append(tmp)

            # Sum the mask files
            maskSum = self._sumImages(self._weightMaskList)
            # Sum the science images
            sciSum = self._sumImages(tmpList)
            del(tmpList)
            # Use the summed sci image values in locations where the maskSum indicates
            # that there is only 1 good pixel to use.  The value will be used in the
            # median_file image
            median_file = np.where(maskSum == self._numberOfImages-1,sciSum,median_file)

        if self._weightMaskList in [None,[]]:
            self._weightMaskList = [np.zeros(self._imageList[0].shape,dtype=self._imageList[0].dtype)]*len(self._imageList)
        # Sum the weightMaskList elements
        maskSum = self._sumImages(self._weightMaskList)

        # Create the minimum image from the stack of input images.
        # Find the maximum pixel value for the image stack.
        maxValue = -1e+9
        for image in self._imageList:
            newMax = image.max()
            if (newMax > maxValue):
                maxValue = newMax

        # For each image, set pixels masked as "bad" to the "super-maximum" value.
        for image in range(len(self._imageList)):
            self._imageList[image] = np.where(self._weightMaskList[image] == 1,maxValue+1,self._imageList[image])

        # Call numcombine throwing out the highest N - 1 pixels.
        tmp = numCombine(self._imageList,numarrayMaskList=None,
                                 combinationType=combtype_median,nlow=0,nhigh=self._numberOfImages-1,
                                 nkeep=1,upper=None,lower=None)
        minimum_file = tmp.combArrObj
        # Reset any pixl at maxValue + 1 to 0.
        minimum_file = np.where(maskSum == self._numberOfImages, 0, minimum_file)

        # Scale the weight images by the background values and add them to the bk
        backgroundFileList = []
        for image in range(len(self._weightImageList)):
            tmp = self._weightImageList[image] * (self._backgroundValueList[image]/(self._exposureTimeList[image]))
            backgroundFileList.append(tmp)

        # Create an image of the total effective background (in DN) per pixel:
        # (which is the sum of all the background-scaled weight files)
        #
        bkgd_file = self._sumImages(backgroundFileList)
        del(backgroundFileList)

        #
        # Scale the weight mask images by the square of the readnoise values
        #
        readnoiseFileList = []
        for image in range(len(self._weightMaskList)):
            tmp = (np.logical_not(self._weightMaskList[image]) *
                   (self._readnoiseList[image] * self._readnoiseList[image]))
            readnoiseFileList.append(tmp)

        # Create an image of the total readnoise**2 per pixel:
        # (which is the sum of all the input readnoise values)
        #
        readnoise_file = self._sumImages(readnoiseFileList)
        del(readnoiseFileList)

        # Create an image of the total effective exposure time per pixel:
        # (which is simply the sum of all the drizzle output weight files)
        #
        weight_file = self._sumImages(self._weightImageList)


        # Scale up both the median and minimum arrays by the total effective exposure time
        # per pixel.
        #
        minimum_file_weighted = minimum_file * weight_file
        median_file_weighted = median_file * weight_file
        del(weight_file)

        # Calculate the 1-sigma r.m.s.:
        #   variance = median_electrons + bkgd_electrons + readnoise**2
        #   rms = sqrt(variance)
        #   This image has units of electrons.
        #
        # make this the abs value so that negative numbers dont throw an exception?
        rms_file2 = np.fmax(
            median_file_weighted + bkgd_file + readnoise_file,
            np.zeros_like(median_file_weighted)
        )
        rms_file = np.sqrt(rms_file2)

        del bkgd_file
        del readnoise_file
        # For the median array, calculate the n-sigma lower threshold to the array
        # and incorporate that into the pixel values.
        #
        median_rms_file = median_file_weighted - (rms_file * self._combine_nsigma1)

        if self._combine_grow != 0:
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

            minimum_flag_file = np.where(np.less(minimum_file_weighted,median_rms_file), 1, 0)

            # The box size value must be an integer.  This is not a problem since __combine_grow should always
            # be an integer type.  The combine_grow column in the MDRIZTAB should also be an integer type.
            boxsize = int(2 * self._combine_grow + 1)
            boxshape = (boxsize, boxsize)
            minimum_grow_file = np.zeros(self._imageList[0].shape,dtype=self._imageList[0].dtype)


            # If the boxcar convolution has failed it is potentially for two reasons:
            #   1) The kernel size for the boxcar is bigger than the actual image.
            #   2) The grow parameter was specified with a value < 0.  This would result
            #      in an illegal boxshape kernel.  The dimensions of the kernel box *MUST*
            #      be integer and greater than zero.
            #
            #   If the boxcar convolution has failed, try to give a meaningfull explanation
            #   as to why based upon the conditionals described above.

            if (boxsize <= 0):
                errormsg1 =  "############################################################\n"
                errormsg1 += "# The boxcar convolution in minmed has failed.  The 'grow' #\n"
                errormsg1 += "# parameter must be greater than or equal to zero. You     #\n"
                errormsg1 += "# specified an input value for the 'grow' parameter of:    #\n"
                errormsg1 += "        combine_grow: " + str(self._combine_grow)+'\n'
                errormsg1 += "############################################################\n"
                raise ValueError(errormsg1)
            if (boxsize > self._imageList[0].shape[0]):
                errormsg2 =  "############################################################\n"
                errormsg2 += "# The boxcar convolution in minmed has failed.  The 'grow' #\n"
                errormsg2 += "# parameter specified has resulted in a boxcar kernel that #\n"
                errormsg2 += "# has dimensions larger than the actual image.  You        #\n"
                errormsg2 += "# specified an input value for the 'grow' parameter of:    #\n"
                errormsg2 += "        combine_grow: " +str(self._combine_grow)+'\n'
                errormsg2 += "############################################################\n"
                print(self._imageList[0].shape)
                raise ValueError(errormsg2)

            # Attempt the boxcar convolution using the boxshape based upon the user input value of "grow"
            ker = np.ones(boxshape) / float(boxsize**2)
            minimum_grow_file = signal.convolve2d(minimum_flag_file, ker,
                                                  boundary='fill', mode='same')

            del(minimum_flag_file)

            temp1 = (median_file_weighted - (rms_file * self._combine_nsigma1))
            temp2 = (median_file_weighted - (rms_file * self._combine_nsigma2))
            median_rms2_file = np.where(np.equal(minimum_grow_file, 0), temp1, temp2)
            del(temp1)
            del(temp2)
            del(rms_file)
            del(minimum_grow_file)

            # Finally decide whether to use the minimim or the median (in counts/s),
            # based on whether the median is more than 3 sigma above the minimum.
            #
            self.combArrObj = np.where(
                np.less(minimum_file_weighted, median_rms2_file),
                minimum_file,
                median_file
            )

        else:
            # Finally decide whether to use the minimim or the median (in counts/s),
            # based on whether the median is more than 3 sigma above the minimum.
            #
            self.combArrObj = np.where(
                np.less(minimum_file_weighted, median_rms_file),
                minimum_file,
                median_file
            )

        # Set fill regions to a pixel value of 0.
        self.combArrObj = np.where(maskSum == self._numberOfImages, 0, self.combArrObj)

#        self.out_file1 = median_rms2_file.copy()
#        self.out_file2 = minimum_file_weighted.copy()


    def _sumImages(self,numarrayObjectList):
        """ Sum a list of numarray objects. """
        if numarrayObjectList in [None, []]:
            return None

        tsum = np.zeros(numarrayObjectList[0].shape, dtype=numarrayObjectList[0].dtype)

        for image in numarrayObjectList:
            tsum += image

        return tsum


def min_med(images, weight_images, readnoise_list, exptime_list,
            background_values, weight_masks=None, combine_grow=1,
            combine_nsigma1=4, combine_nsigma2=3, fillval=False):
    """ Create a median array, rejecting the highest pixel and
    computing the lowest valid pixel after mask application.

    .. note::
        In this version of the mimmed algorithm we assume that the units of
        all input data is electons.

    Parameters
    ----------
    images : list of numpy.ndarray
        List of input data to be combined.

    weight_images : list of numpy.ndarray
        List of input data weight images to be combined.

    readnoise_list : list
        List of readnoise values to use for the input images.

    exptime_list : list
        List of exposure times to use for the input images.

    background_values : list
        List of image background values to use for the input images.

    weight_masks : list of numpy.ndarray, None
        List of imput data weight masks to use for pixel rejection.
        (Default: `None`)

    combine_grow : int
        Radius (pixels) for neighbor rejection. (Default: 1)

    combine_nsigma1 : float
        Significance for accepting minimum instead of median. (Default: 4)

    combine_nsigma2 : float
        Significance for accepting minimum instead of median. (Default: 3)

    fillval : bool
        Turn on use of imedian/imean. (Default: `False`)

    Returns
    -------
    combined_array : numpy.ndarray
        Combined array.

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
    # The total effective background in the final image is calculated as
    # follows:
    #   - convert background for each input image to counts/s
    #     (divide by exptime)
    #   - multiply this value by the weight image, to obtain the effective
    #     background counts (in DN) for each pixel, for each image
    #   - Add these images together, to obtain the total effective background
    #     for the combined image.
    #
    # Once we've made these two files, then calculate the SNR based on the
    # median-pixel image, and compare with the minimum.

    nimages = len(images)
    combtype_median = 'imedian' if fillval else 'median'
    images = np.asarray(images)
    weight_images = np.asarray(weight_images)

    if weight_masks == [] or weight_masks is None:
        weight_masks = None
        mask_sum = np.zeros(images.shape[1:], dtype=np.int16)
        all_bad_idx = np.array([], dtype=np.int)
        all_bad_idy = np.array([], dtype=np.int)
    else:
        weight_masks = np.asarray(weight_masks, dtype=np.bool)
        mask_sum = np.sum(weight_masks, axis=0, dtype=np.int16)
        all_bad_idx, all_bad_idy = np.where(mask_sum == nimages)

    # Create a different median image based upon the number of images in the
    # input list.
    if nimages == 2:
        median_file = num_combine(
            images,
            masks=weight_masks,
            combination_type='imean' if fillval else 'mean',
            nlow=0, nhigh=0, lower=None, upper=None
        )

    else:
        # The value of NHIGH=1 will cause problems when there is only 1 valid
        # unmasked input image for that pixel due to a difference in behavior
        # between 'num_combine' and 'iraf.imcombine'.
        # This value may need to be adjusted on the fly based on the number of
        # inputs and the number of masked values/pixel.
        #
        median_file = num_combine(
            images,
            masks=weight_masks,
            combination_type=combtype_median,
            nlow=0, nhigh=1, lower=None, upper=None
        )

        # The following section of code will address the problem caused by
        # having a value of nhigh = 1.  This will behave in a way similar to
        # the way the IRAF task IMCOMBINE behaves.  In order to accomplish
        # this, the following procedure will be followed:
        # 1) The input masks will be summed.
        # 2) The science data will be summed.
        # 3) In the locations of the summed mask where the sum is 1 less than
        #    the total number of images, the value of that location in the
        #    summed science image will be used to replace the existing value
        #    in the existing median_file.
        #
        # This procedure is being used to prevent too much data from being
        # thrown out of the image. Take for example the case of 3 input images.
        # In two of the images the pixel locations have been masked out.
        # Now, if nhigh is applied there will be no value to use for that
        # position.  However, if this new procedure is used that value in
        # the resulting images will be the value that was rejected by the
        # nhigh rejection step.

        # We need to make certain that "bad" pixels in the sci data are set to
        # 0. That way, when the sci images are summed, the value of the sum
        # will only come from the "good" pixels.
        if weight_masks is None:
            sci_sum = np.sum(images, axis=0)
            if nimages == 1:
                median_file = sci_sum

        else:
            sci_sum = np.sum(images * np.logical_not(weight_masks), axis=0)
            # Use the summed sci image values in locations where the mask_sum
            # indicates that there is only 1 good pixel to use. The value will
            # be used in the median_file image
            idx = np.where(mask_sum == (nimages - 1))
            median_file[idx] = sci_sum[idx]

    # Create the minimum image from the stack of input images.
    if weight_masks is not None:
        # make a copy of images to avoid side-effect of modifying input
        # argument:
        images = images.copy()
        images[weight_masks] = np.nan
        images[:, all_bad_idx, all_bad_idy] = 0
        minimum_file = np.nanmin(images, axis=0)
    else:
        minimum_file = np.amin(images, axis=0)

    # Scale the weight images by the background values and add them to the bk
    # Create an image of the total effective background (in DN) per pixel:
    # (which is the sum of all the background-scaled weight files)
    s = np.asarray([bv / et for bv, et in
                    zip(background_values, exptime_list)])
    bkgd_file = np.sum(weight_images * s[:, None, None], axis=0)

    # Scale the weight mask images by the square of the readnoise values.
    # Create an image of the total readnoise**2 per pixel
    # (which is the sum of all the input readnoise values).
    if weight_masks is None:
        rdn2 = sum((r**2 for r in readnoise_list))
        readnoise_file = rdn2 * np.ones_like(images[0])

    else:
        readnoise_file = np.sum(
            np.logical_not(weight_masks) *
            (np.asarray(readnoise_list)**2)[:, None, None],
            axis=0
        )

    # Create an image of the total effective exposure time per pixel:
    # (which is simply the sum of all the drizzle output weight files)
    weight_file = np.sum(weight_images, axis=0)

    # Scale up both the median and minimum arrays by the total effective
    # exposure time per pixel.
    minimum_file_weighted = minimum_file * weight_file
    median_file_weighted = median_file * weight_file
    del weight_file

    # Calculate the 1-sigma r.m.s.:
    #   variance = median_electrons + bkgd_electrons + readnoise**2
    #   rms = sqrt(variance)
    #   This image has units of electrons.
    #
    # make this the abs value so that negative numbers dont throw an exception?
    rms_file2 = np.fmax(
        median_file_weighted + bkgd_file + readnoise_file,
        np.zeros_like(median_file_weighted)
    )
    rms_file = np.sqrt(rms_file2)
    del bkgd_file, readnoise_file

    # For the median array, calculate the n-sigma lower threshold to the array
    # and incorporate that into the pixel values.
    median_rms_file = median_file_weighted - rms_file * combine_nsigma1

    if combine_grow != 0:
        # Do a more sophisticated rejection: For all cases where the minimum
        # pixel will be accepted instead of the median, set a lower threshold
        # for that pixel and the ones around it (ie become less conservative
        # in rejecting the median). This is because in cases of
        # triple-incidence cosmic rays, quite often the low-lying outliers
        # of the CRs can influence the median for the initial relatively high
        # value of sigma, so a lower threshold must be used to mnake sure that
        # the minimum is selected.
        #
        # This is done as follows:
        # 1) make an image which is zero everywhere except where the minimum
        #    will be accepted
        # 2) box-car smooth this image, to make these regions grow.
        # 3) In the file "median_rms_file_electrons", replace these pixels
        #     by median - combine_nsigma2 * rms
        #
        # Then use this image in the final replacement, in the same way as for
        # the case where this option is not selected.
        minimum_flag_file = np.less(minimum_file_weighted,
                                    median_rms_file).astype(np.float64)

        # The box size value must be an integer. This is not a problem since
        # __combine_grow should always be an integer type. The combine_grow
        # column in the MDRIZTAB should also be an integer type.
        boxsize = int(2 * combine_grow + 1)
        boxshape = (boxsize, boxsize)
        minimum_grow_file = np.zeros_like(images[0])

        # If the boxcar convolution has failed it is potentially for
        # two reasons:
        #   1) The kernel size for the boxcar is bigger than the actual image.
        #   2) The grow parameter was specified with a value < 0.  This would
        #      result in an illegal boxshape kernel. The dimensions of the
        #      kernel box *MUST* be integer and greater than zero.
        #
        #   If the boxcar convolution has failed, try to give a meaningfull
        #   explanation as to why based upon the conditionals described above.
        if boxsize <= 0:
            errormsg1 = "############################################################\n"
            errormsg1 += "# The boxcar convolution in minmed has failed.  The 'grow' #\n"
            errormsg1 += "# parameter must be greater than or equal to zero. You     #\n"
            errormsg1 += "# specified an input value for the 'grow' parameter of:    #\n"
            errormsg1 += "        combine_grow: " + str(combine_grow)+'\n'
            errormsg1 += "############################################################\n"
            raise ValueError(errormsg1)

        if boxsize > images.shape[1]:
            errormsg2 = "############################################################\n"
            errormsg2 += "# The boxcar convolution in minmed has failed.  The 'grow' #\n"
            errormsg2 += "# parameter specified has resulted in a boxcar kernel that #\n"
            errormsg2 += "# has dimensions larger than the actual image.  You        #\n"
            errormsg2 += "# specified an input value for the 'grow' parameter of:    #\n"
            errormsg2 += "        combine_grow: " + str(combine_grow) + '\n'
            errormsg2 += "############################################################\n"
            print(images.shape[1:])
            raise ValueError(errormsg2)

        # Attempt the boxcar convolution using the boxshape based upon the user
        # input value of "grow"
        ker = np.ones((boxsize, boxsize)) / float(boxsize**2)
        minimum_grow_file = signal.convolve2d(minimum_flag_file, ker,
                                              boundary='fill', mode='same')

        median_rms_file = np.where(
            np.equal(minimum_grow_file, 0),
            median_file_weighted - rms_file * combine_nsigma1,
            median_file_weighted - rms_file * combine_nsigma2
        )
        del rms_file, minimum_grow_file

    # Finally decide whether to use the minimim or the median (in counts/s),
    # based on whether the median is more than 3 sigma above the minimum.
    combined_array = np.where(
        np.less(minimum_file_weighted, median_rms_file),
        minimum_file,
        median_file
    )
    # Set fill regions to a pixel value of 0.
    combined_array[all_bad_idx, all_bad_idy] = 0

    return combined_array
