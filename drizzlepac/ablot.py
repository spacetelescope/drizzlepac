"""
In this step the median image now gets blotted back to create median-cleaned
images which can be compared directly with each input image to identify
cosmic-rays.

:Authors: Warren Hack

:License: :doc:`/LICENSE`

"""

import os
import logging
import numpy as np
from stsci.tools import fileutil
from . import outputimage
from . import wcs_functions
from . import util
import stwcs
from stwcs import distortion

log = logging.getLogger(__name__)

try:
    from . import cdriz
except ImportError:
    cdriz = None
    log.error('\n Coordinate transformation and image resampling library NOT found!')
    log.error('\n Please check the installation of this package to insure C code was built successfully.')
    raise ImportError

from . import __version__

__all__ = ["blot", "runBlot"]

__taskname__ = "ablot"
STEP_NUM = 5
PROCSTEPS_NAME = "Blot"


def blot(
    data,
    reference,
    outdata,
    configObj=None,
    wcsmap=wcs_functions.WCSMap,
    editpars=False,
    **input_dict,
):
    """
    The median image is the combination of the WCS aligned input images
    that have already had the distortion model applied. Taking the median
    of the aligned images allows for a statistical rejection of bad pixels
    from the image stack. The resulting median image can then be input for
    the blot task with the goal of creating 'cleaned' versions of the input
    images at each of their respective dither locations. These "blotted" images
    can then be directly compared to the original distorted input images for
    detection of image artifacts (i.e. bad-pixels, hot pixels, and cosmic-rays)
    whose locations will be saved to the output badpixel masks.

    Aside from the input parameters, this step only requires
    opening the single median image created from all the input images.
    A distorted version of the median image corresponding
    to each input 'chip' (extension) is written as output from this
    step as separate simple FITS images.

    For more information on the science applications of the blot task,
    see the `DrizzlePac Handbook <http://drizzlepac.stsci.edu>`_.

    Parameters
    ----------
    data : str
        Input distortion-corrected (median or drizzled) image to be used as the
        source for creating blotted images.

    reference : str
        Filename of image to read to define the blotted WCS; image with distortion
        to be matched by output blotted image.

    outdata : str
        Filename for output blotted image.

    configObj : object, optional
        Contains all the parameters needed to control the blot operation.


    Notes
    -----

    The following parameters are contained in the ``configObj`` object.

    coeffs : bool (Default Value = True)
        This parameters specifies whether or not to use the header-based distortion
        coefficients when creating the blotted, distorted image.  If False, no
        distortion will be applied at all, effectively working as a cut-out operation.

    interp : str{'nearest', 'linear', 'poly3', 'poly5', 'sinc'} (Default = 'poly5')
        This parameter defines the method of interpolation to be used when
        blotting drizzled images back to their original WCS solution.
        Valid options include:

        - **nearest**: Nearest neighbor
        - **linear**: Bilinear interpolation in x and y
        - **poly3**: Third order interior polynomial in x and y
        - **poly5**: Fifth order interior polynomial in x and y
        - **sinc**: Sinc interpolation (accurate but slow)

        The 'poly5' interpolation method has been chosen as the default because
        it is relatively fast and accurate.

        If 'sinc' interpolation is selected, then the value of the parameter
        for ``blot_sinscl`` will be used to specify the size of the sinc
        interpolation kernel.

    sinscl : float (Default Value = 1.0)
        Size of the sinc interpolation kernel in pixels.

    stepsize : int (Default Value = 10)
        Number of pixels for WCS interpolation.  The distortion model will be sampled
        exactly and completely every ``stepsize`` pixel with bi-linear interpolation
        being used to compute the distortion for intermediate pixels. This optimization
        speeds up the computation significantly when ``stepsize`` >> 1 at the expense
        of interpolation errors for intermediate pixels.

    addsky : bool (Default Value = Yes)
        Add back a sky value using the ``MDRIZSKY`` value from the header.
        If 'Yes' (``True``), the ``blot_skyval`` parameter is ignored.

    skyval : float (Default Value = 0.0)
        This is a user-specified custom sky value to be added to the blot image.
        This is only used if ``blot_addsky`` is 'No' (``False``).

    in_units : str{'cps', 'counts'} (Default Value= 'cps')
        Units of input (drizzled) image.
        Valid options are **'cps'** and **'counts'**.

    out_units : str{'cps', 'counts'} (Default Value = 'counts')
        Units of the ouput (blotted) image.
        Valid options are **'cps'** and **'counts'**.

    expkey : str (Default Value = 'exptime)
        Name of keyword to use to extract exposure time value, which will be used to
        scale the blotted image to the final output flux values when ``out_units`` is
        set to **counts**.

    expout : str or float (Default Value = 'input')
        Value of exposure time to use in scaling the output blotted image when
        ``out_units`` is set to **counts**. If set to **'input'**, the value will be
        read in from the input image header keyword specified by ``expkey``.

        .. note:: The following parameters, when set, will override any value determined
            from ``refimage`` if a reference image was specified.

    outscale : float, optional
        Absolute size of output pixels in arcsec/pixel.

    orient : float
        Orientation of output (PA of Y axis, N through E)

    raref : float
        RA of reference point on output image (CRVAL1, degrees).

    decref : float
        Dec of reference point on output image (CRVAL2, degrees)

    xrefpix : float
        Reference pixel X position on output (CRPIX1).

    yrefpix : float
        Reference pixel Y position on output (CRPIX2)

    outnx : float
        Size of output image's X-axis (pixels).

    outny : float
        Size of output image's Y-axis (pixels).


    These tasks are designed to work together seemlessly when run in the full
    ``AstroDrizzle`` interface. More advanced users may wish to create specialized
    scripts for their own datasets, making use of only a subset of the
    predefined ``AstroDrizzle`` tasks, or add additional processing, which may
    be usefull for their particular data. In these cases, individual access to
    the tasks is important.

    Something to keep in mind is that the full ``AstroDrizzle`` interface will
    make backup copies of your original files and place them in the ``OrIg/``
    directory of your current working directory. If you are working with
    the stand alone interfaces, it is assumed that the user has already
    taken care of backing up their original datafiles as the input file
    with be directly altered.


    Examples
    --------
    1. Basic example of how to call :py:func:`blot` yourself from a Python
       command line, using the default parameter settings::

        >>> from drizzlepac import ablot
        >>> ablot.blot()

    2. Re-create the ``SCI,1`` chip from ``j8c0d1bwq_flc.fits`` using the
       ``AstroDrizzle`` median image ``adriz_aligned_wcs_f814w_med.fits``::

        >>> from drizzlepac import ablot
        >>> from stsci.tools import teal
        >>> blotobj = teal.load('ablot')  # get default values
        >>> ablot.blot('adriz_aligned_wcs_f814w_med.fits',
        ...            'j8c0d1bwq_flc.fits[sci,1]',
        ...            'aligned_f814w_sci1_blot.fits',
        ...            configObj=blotobj)
    """

    if input_dict is None:
        input_dict = {}
    input_dict["data"] = data
    input_dict["reference"] = reference
    input_dict["outdata"] = outdata

    # gets configObj defaults using EPAR/TEAL.
    #
    # Also insure that the input_dict (user-specified values) are folded in
    # with a fully populated configObj instance.
    configObj = util.getDefaultConfigObj(__taskname__, configObj,
                                         input_dict, loadOnly=(not editpars))
    if configObj is None:
        return

    if not editpars:
        run(configObj, wcsmap=wcsmap)


def run(configObj, wcsmap=None):
    """
    Run the blot task based on parameters provided interactively by the user.

    """

    # Insure all output filenames specified have .fits extensions
    if configObj["outdata"][-5:] != ".fits":
        configObj["outdata"] += ".fits"

    scale_pars = configObj["Data Scaling Parameters"]
    user_wcs_pars = configObj["User WCS Parameters"]

    # PyFITS can be used here as it will always operate on
    # output from PyDrizzle (which will always be a FITS file)
    # Open the input (drizzled?) image
    _fname, _sciextn = fileutil.parseFilename(configObj["data"])
    _inimg = fileutil.openImage(_fname, memmap=False)
    _expin = fileutil.getKeyword(configObj["data"], scale_pars["expkey"], handle=_inimg)

    # Return the PyFITS HDU corresponding to the named extension
    _scihdu = fileutil.getExtn(_inimg, _sciextn)
    _insci = _scihdu.data.copy()

    _inexptime = 1.0
    if scale_pars["in_units"] == "counts":
        if scale_pars["expkey"] in _inimg["PRIMARY"].header:
            _inexptime = _inimg["PRIMARY"].header[scale_pars["expkey"]]
        elif "DRIZEXPT" in _inimg["PRIMARY"].header:
            # Try keyword written out by new 'drizzle' if no valid 'expkey' was given
            _inexptime = _inimg["PRIMARY"].header["DRIZEXPT"]
        else:
            raise ValueError('No valid exposure time keyword could be found '
                             'for input %s' % configObj['data'])
    # always convert input to 'cps' for blot() algorithm
    if _inexptime != 0.0 or _inexptime != 1.0:
        np.divide(_insci, _inexptime, _insci)

    _inimg.close()
    del _inimg

    # read in WCS from source (drizzled) image
    source_wcs = stwcs.wcsutil.HSTWCS(configObj["data"])
    if source_wcs.wcs.is_unity():
        log.warning("No valid WCS found for input drizzled image: {}!".format(configObj['data']))

    # define blot_wcs
    blot_wcs = None
    _refname, _refextn = fileutil.parseFilename(configObj["reference"])
    if os.path.exists(_refname):
        # read in WCS from pre-existing output image
        blot_wcs = stwcs.wcsutil.HSTWCS(configObj["reference"])
        if blot_wcs.wcs.is_unity():
            log.warning("No valid WCS found for output image: {} !".format(configObj['reference']))

    # define blot WCS based on input images or specified reference WCS values
    if user_wcs_pars["user_wcs"]:
        blot_wcs = wcs_functions.build_hstwcs(
            user_wcs_pars["raref"],
            user_wcs_pars["decref"],
            user_wcs_pars["xrefpix"],
            user_wcs_pars["yrefpix"],
            int(user_wcs_pars["outnx"]),
            int(user_wcs_pars["outny"]),
            user_wcs_pars["outscale"],
            user_wcs_pars["orient"],
        )
        configObj["coeffs"] = None

    # If blot_wcs is still not defined at this point, we have a problem...
    if blot_wcs is None:
        blot_wcs = stwcs.distortion.utils.output_wcs([source_wcs], undistort=False)

    out_wcs = blot_wcs.copy()
    # perform blotting operation now
    _outsci = do_blot(_insci, source_wcs, out_wcs, _expin, coeffs=configObj['coeffs'],
                    interp=configObj['interpol'], sinscl=configObj['sinscl'],
            stepsize=configObj['stepsize'], wcsmap=wcsmap)
    # create output with proper units and exptime-scaling
    if scale_pars["out_units"] == "counts":
        if scale_pars["expout"] == "input":
            _outscale = fileutil.getKeyword(
                configObj["reference"], scale_pars["expkey"]
            )
            # _outscale = _expin
        else:
            _outscale = float(scale_pars['expout'])
        log.debug("Output blotted images scaled by exptime of {}".format(_outscale))
        np.multiply(_outsci, _outscale, _outsci)

    # Add sky back in to the blotted image, as specified by the user
    if configObj["addsky"]:
        skyval = _scihdu.header["MDRIZSKY"]
    else:
        skyval = configObj['skyval']
    log.debug("Added {} counts back in to blotted image as sky.".format(skyval))
    _outsci += skyval

    del _scihdu

    # Write output Numpy objects to a PyFITS file
    # Blotting only occurs from a drizzled SCI extension
    # to a blotted SCI extension...
    outputimage.writeSingleFITS(
        _outsci, blot_wcs, configObj["outdata"], configObj["reference"]
    )


#
#### Top-level interface from inside AstroDrizzle
#
def runBlot(imageObjectList, output_wcs, configObj={},
            wcsmap=wcs_functions.WCSMap, procSteps=None):
    """
    Top-level interface for running the blot step from within AstroDrizzle.

    This function performs the blot operation on a list of input images to create
    distorted versions that can be compared with the original input images for
    cosmic ray detection. It takes the median (drizzled) image and "blots" it
    back to match the original distorted geometry of each input image.

    Parameters
    ----------
    imageObjectList : list
        List of imageObject instances to be processed through the blot step.
        Each imageObject contains the input image data and associated metadata.

    output_wcs : WCS object
        World Coordinate System object that defines the coordinate transformation
        for the output (median/drizzled) image that will be blotted back to the
        input image geometries.

    configObj : dict, optional
        Configuration object containing all the parameters needed to control
        the blot operation, including interpolation method, sky handling, and
        other blot-specific settings. Default is an empty dictionary.

    wcsmap : WCSMap, optional
        Custom mapping class to use for coordinate transformations between
        the drizzled and blotted image coordinate systems.

    procSteps : ProcessingSteps, optional
        Object used to track and log the progress of processing steps within
        the AstroDrizzle pipeline. If provided, the blot step will be logged.
        Default is None.

    Notes
    -----
    This function serves as the high-level interface called by AstroDrizzle to
    perform the blot operation. It checks whether the blot step should be
    performed based on the configuration settings, and if so, calls the
    lower-level ``run_blot`` function to perform the actual blotting operation.
    """
    if procSteps is not None:
        procSteps.addStep(PROCSTEPS_NAME)
        if not imageObjectList:
            procSteps.endStep(PROCSTEPS_NAME, reason="skipped", delay_msg=True)

    blot_name = util.getSectionName(configObj, STEP_NUM)

    # This can be called directly from MultiDrizle, so only execute if
    # switch has been turned on (no guarantee MD will check before calling).
    if configObj[blot_name]["blot"]:
        paramDict = buildBlotParamDict(configObj)

        log.debug(f"USER INPUT PARAMETERS for {PROCSTEPS_NAME} Step:")
        util.printParams(paramDict, log=log)

        run_blot(imageObjectList, output_wcs.single_wcs, paramDict, wcsmap=wcsmap)
    else:
        log.debug('Blot step not performed.')
        return

    if procSteps is not None:
        procSteps.endStep(PROCSTEPS_NAME)


# Run 'drizzle' here...
#
def buildBlotParamDict(configObj):
    blot_name = util.getSectionName(configObj, STEP_NUM)

    paramDict = {'blot_interp':configObj[blot_name]['blot_interp'],
                'blot_sinscl':configObj[blot_name]['blot_sinscl'],
                'blot_addsky':configObj[blot_name]['blot_addsky'],
                'blot_skyval':configObj[blot_name]['blot_skyval'],
                'coeffs':configObj['coeffs']}
    return paramDict


def _setDefaults(configObj={}):
    """set up the default parameters to run drizzle
    build,single,units,wt_scl,pixfrac,kernel,fillval,
    rot,scale,xsh,ysh,blotnx,blotny,outnx,outny,data
    """

    paramDict = {
        "build": True,
        "single": True,
        "in_units": "cps",
        "wt_scl": 1.0,
        "pixfrac": 1.0,
        "kernel": "square",
        "fillval": 999.0,
        "rot": 0.0,
        "scale": 1.0,
        "xsh": 0.0,
        "ysh": 0.0,
        "blotnx": 2048,
        "blotny": 2048,
        "outnx": 4096,
        "outny": 4096,
        "data": None,
        "driz_separate": True,
        "driz_combine": False,
    }

    if len(configObj) != 0:
        for key in configObj.keys():
            paramDict[key] = configObj[key]

    return paramDict


def run_blot(imageObjectList, output_wcs, paramDict, wcsmap=wcs_functions.WCSMap):
    """
    run_blot(imageObjectList, output_wcs, paramDict, wcsmap=wcs_functions.WCSMap)

    Perform the blot operation on the list of images.
    """
    # Insure that input imageObject is a list
    if not isinstance(imageObjectList, list):
        imageObjectList = [imageObjectList]
    #
    # Setup the versions info dictionary for output to PRIMARY header
    # The keys will be used as the name reported in the header, as-is
    #
    _versions = {
        "AstroDrizzle": __version__,
        "PyFITS": util.__fits_version__,
        "Numpy": util.__numpy_version__,
    }

    _hdrlist = []

    for img in imageObjectList:

        for chip in img.returnAllChips(extname=img.scienceExt):
            log.debug(f'Blot: creating blotted image: {chip.outputNames["data"]}')

            #### Check to see what names need to be included here for use in _hdrlist
            chip.outputNames["driz_version"] = _versions["AstroDrizzle"]
            outputvals = chip.outputNames.copy()
            outputvals.update(img.outputValues)
            outputvals["blotnx"] = chip.wcs.naxis1
            outputvals["blotny"] = chip.wcs.naxis2
            _hdrlist.append(outputvals)

            plist = outputvals.copy()
            plist.update(paramDict)

            # PyFITS can be used here as it will always operate on
            # output from PyDrizzle (which will always be a FITS file)
            # Open the input science file
            medianPar = "outMedian"
            outMedianObj = img.getOutputName(medianPar)
            if img.inmemory:
                outMedian = img.outputNames[medianPar]
                _fname, _sciextn = fileutil.parseFilename(outMedian)
                _inimg = outMedianObj
            else:
                outMedian = outMedianObj
                _fname, _sciextn = fileutil.parseFilename(outMedian)
                _inimg = fileutil.openImage(_fname, memmap=False)

            # Return the PyFITS HDU corresponding to the named extension
            _scihdu = fileutil.getExtn(_inimg, _sciextn)
            _insci = _scihdu.data.copy()
            _inimg.close()
            del _inimg, _scihdu

            _outsci = do_blot(
                _insci,
                output_wcs,
                chip.wcs,
                chip._exptime,
                coeffs=paramDict["coeffs"],
                interp=paramDict["blot_interp"],
                sinscl=paramDict["blot_sinscl"],
                wcsmap=wcsmap,
            )
            # Apply sky subtraction and unit conversion to blotted array to
            # match un-modified input array
            if paramDict["blot_addsky"]:
                skyval = chip.computedSky
            else:
                skyval = paramDict["blot_skyval"]
            _outsci /= chip._conversionFactor
            if skyval is not None:
                _outsci += skyval
                log.debug(f"Applying sky value of {skyval:0.6f} to blotted "
                          f"image {chip.outputNames['data']}")

            # Write output Numpy objects to a PyFITS file
            # Blotting only occurs from a drizzled SCI extension
            # to a blotted SCI extension...

            _outimg = outputimage.OutputImage(
                _hdrlist, paramDict, build=False, wcs=chip.wcs, blot=True
            )
            _outimg.outweight = None
            _outimg.outcontext = None
            outimgs = _outimg.writeFITS(
                plist["data"],
                _outsci,
                None,
                versions=_versions,
                blend=False,
                virtual=img.inmemory,
            )

            img.saveVirtualOutputs(outimgs)
            # _buildOutputFits(_outsci,None,plist['outblot'])
            _hdrlist = []

            del _outsci

        del _outimg


def do_blot(
    source,
    source_wcs,
    blot_wcs,
    exptime,
    coeffs=True,
    interp="poly5",
    sinscl=1.0,
    stepsize=10,
    wcsmap=None,
):
    """Core functionality of performing the 'blot' operation to create a single
    blotted image from a single source image.
    All distortion information is assumed to be included in the WCS specification
    of the 'output' blotted image given in 'blot_wcs'.

    This is the simplest interface that can be called for stand-alone
    use of the blotting function.

    Parameters
    ----------
    source
        Input numpy array of undistorted source image in units of 'cps'.
    source_wcs
        HSTWCS object representing source image distortion-corrected WCS.
    blot_wcs
        (py)wcs.WCS object representing the blotted image WCS.
    exptime
        exptime to use for scaling output blot image. A value of 1 will
        result in output blot image in units of 'cps'.
    coeffs
        Flag to specify whether or not to use distortion coefficients
        associated with blot_wcs. If False, do not apply any distortion
        model.
    interp
        Form of interpolation to use when blotting pixels. Valid options::

            "nearest","linear","poly3", "poly5"(default), "spline3", "sinc"
    sinscl
        Scale for sinc interpolation kernel (in output, blotted pixels)
    stepsize
        Number of pixels for WCS interpolation
    wcsmap
        Custom mapping class to use to provide transformation from
        drizzled to blotted WCS.  Default will be to use
        `~drizzlepac.wcs_functions.WCSMap`.

    """
    _outsci = np.zeros(blot_wcs.array_shape, dtype=np.float32)

    # Now pass numpy objects to callable version of Blot...
    build = False
    misval = 0.0
    kscale = 1.0

    xmin = 1
    ymin = 1
    xmax, ymax = source_wcs.pixel_shape

    # compute the undistorted 'natural' plate scale for this chip
    if coeffs:
        wcslin = distortion.utils.make_orthogonal_cd(blot_wcs)
    else:
        wcslin = blot_wcs
        blot_wcs.sip = None
        blot_wcs.cpdis1 = None
        blot_wcs.cpdis2 = None
        blot_wcs.det2im = None

    if wcsmap is None and cdriz is not None:
        """
        Use default C mapping function.
        """
        log.debug('Using default C-based coordinate transformation...')
        mapping = cdriz.DefaultWCSMapping(
            blot_wcs,
            source_wcs,
            blot_wcs.pixel_shape[0],
            blot_wcs.pixel_shape[1],
            stepsize,
        )
        pix_ratio = source_wcs.pscale / wcslin.pscale
    else:
        #
        ##Using the Python class for the WCS-based transformation
        #
        # Use user provided mapping function
        log.debug('Using coordinate transformation defined by user...')
        if wcsmap is None:
            wcsmap = wcs_functions.WCSMap
        wmap = wcsmap(blot_wcs, source_wcs)
        mapping = wmap.forward
        pix_ratio = source_wcs.pscale / wcslin.pscale

    t = cdriz.tblot(
        source,
        _outsci,
        xmin,
        xmax,
        ymin,
        ymax,
        pix_ratio,
        kscale,
        1.0,
        1.0,
        "center",
        interp,
        exptime,
        misval,
        sinscl,
        1,
        mapping,
    )
    del mapping

    return _outsci
