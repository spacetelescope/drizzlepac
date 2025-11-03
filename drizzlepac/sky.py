#!/usr/bin/env python
"""
This step measures, subtracts and/or equalizes the sky from each
input image while recording the subtracted value in the image header.

:Authors: Christopher Hanley, Megan Sosey, Mihai Cara

:License: :doc:`/LICENSE`

"""
import os, sys

import numpy as np

from stsci.tools import fileutil, logutil
from stsci.tools.bitmask import interpret_bit_flags
import stsci.imagestats as imagestats

from stsci.skypac.skymatch import skymatch
from stsci.skypac.utils import MultiFileLog, ext2str, \
     file_name_components, in_memory_mask, temp_mask_file, openImageEx
from stsci.skypac.parseat import FileExtMaskInfo, parse_at_file

from . import processInput

from . import util
from . import __version__


__taskname__= "drizzlepac.sky" #looks in drizzlepac for sky.cfg
__all__ = ['sky']
STEP_NUM = 2  #this relates directly to the syntax in the cfg file
PROCSTEPS_NAME = "Subtract Sky"


log = logutil.create_logger(__name__, level=logutil.logging.NOTSET)


# this is the user access function
def sky(input=None,outExt=None,configObj=None, group=None, editpars=False, **inputDict):
    """
    Function for computing and subtracting (or equalizing/matching) the backgroud
    in input images. The algorithm for sky subtraction can be selected through
    the ``skymethod`` parameter. This function will update the ``MDRIZSKY`` keyword
    in the headers of the input files.

    Sky subtraction is generally recommended for optimal flagging and removal of
    CR's when the sky background is more than a few electrons.  However, some
    science applications may require the sky to not be removed, allowing for the
    final drizzle step to be performed with no sky subtraction. If you turn off
    sky subtraction, you should also set drizzle.pixfrac to 1, otherwise
    variations in sky between images will add noise to your data.

    In addition to the "pure" sky computation, this task can be used for sky
    "equalization", that is, it can match sky values in the images
    that are part of a mosaic.

    For cameras with multiple detectors (such as ACS/WFC, WFPC2, or WFC3),
    the sky values in each exposure are first measured separately for
    the different detectors. These different values are then compared,
    and the lowest measured sky value is used as the estimate for all of the
    detectors for that exposure. This is based on the premise that for large
    extended or bright targets, the pixel intensity distribution in one or more
    of the detectors may be significantly skewed toward the bright end by the
    target itself, thereby overestimating the sky on that detector. If the other
    detector is less affected by such a target, then its sky value  will be lower,
    and can therefore also be substituted as the sky value for the detector
    with the bright source. The input file's primary headers is updated with the
    computed sky value.

    For more information on the science applications of the sky task,
    see the `DrizzlePac Handbook: <http://drizzlepac.stsci.edu>`_.


    Parameters
    ----------

    input : str or list of str (Default = None)
        A Python list of image filenames, or just a single filename.

    outExt : str (Default = None)
        The extension of the output image. If the output already exists
        then the input image is overwritten.

    configObj : configObject (Default = None)
        An instance of ``configObject``

    group : int (Default = None)
        The group of the input image.

    editpars : bool (Default = False)
        A parameter that allows user to edit input parameters by hand in the GUI.

    inputDict : dict, optional
        An optional list of parameters specified by the user.

        .. note::
            These are parameters that ``configObj`` should contain by default. These
            parameters can be altered on the fly using the ``inputDict``. If ``configObj``
            is set to None and there is no ``inputDict`` information, then the values
            for the parameters will be pulled from the default configuration files
            for the task.

        Table of optional parameters that should be in ``configobj`` and can also be
        specified in ``inputDict``.

        ===============   ===================================================================
        Name              Definition
        ===============   ===================================================================
        skyuser           KEYWORD in header which indicates a sky subtraction value to use.
        skymethod         Sky computation method
        skysub            Perform sky subtraction?
        skywidth          Bin width of histogram for sampling sky statistics (in sigma)
        skystat           Sky correction statistics parameter
        skylower          Lower limit of usable data for sky (always in electrons)
        skyupper          Upper limit of usable data for sky (always in electrons)
        skyclip           Number of clipping iterations
        skylsigma         Lower side clipping factor (in sigma)
        skyusigma         Upper side clipping factor (in sigma)
        skymask_cat       Catalog file listing image masks
        use_static        Use static mask for skymatch computations?
        sky_bits          Bit flags for identifying bad pixels in DQ array
        skyuser           KEYWORD indicating a sky subtraction value if done by user
        skyfile           Name of file with user-computed sky values
        in_memory         Optimize for speed or for memory use
        ===============   ===================================================================

        These optional parameters are described in more detail below in the
        "Other Parameters" section.

    Notes
    -----

    skysub : bool (Default = Yes)
        Turn on or off sky subtraction on the input data. When ``skysub`` is set
        to ``no``, then ``skyuser`` field will be enabled and if user specifies a header
        keyword showing the sky value in the image, then that value will be used for
        CR-rejection but it will not be subtracted from the (drizzled) image data.
        If user sets ``skysub`` to ``yes`` then ``skyuser`` field will be disabled
        (and if it is not empty - it will be ignored) and user can use one of the
        methods available through the ``skymethod`` parameter to compute the sky
        or provide a file (see ``skyfile`` parameter) with values that should be
        subtracted from (single) drizzled images.

    skymethod : {'localmin', 'globalmin+match', 'globalmin', 'match'}, optional (Default = 'localmin')

        Select the algorithm for sky computation:

        * *localmin* : compute a common sky for all members of *an exposure*
          (see NOTES below). For a typical use, it will compute
          sky values for each chip/image extension (marked for sky
          subtraction in the :py:obj:`input` parameter) in an input image,
          and it will subtract the previously found minimum sky value
          from all chips (marked for sky subtraction) in that image.
          This process is repeated for each input image.

          .. note::
              This setting is recommended when regions of overlap between images
              are dominated by "pure" sky (as opposite to extended, diffuse
              sources). This is similar to the "skysub" algorithm used in previous
              versions of astrodrizzle.

        * *globalmin* : compute a common sky value for all members of
          *all exposures* (see NOTES below). It will compute
          sky values for each chip/image extension (marked for sky
          subtraction in the ``input`` parameter) in all input
          images, find the minimum sky value, and then it will
          subtract the same minimum sky value from all chips
          (marked for sky subtraction) in all images. This method *may*
          useful when input images already have matched background values.

        * *match* : compute differences in sky values between images
          in common (pair-wise) sky regions. In this case computed sky values
          will be relative (delta) to the sky computed in one of the
          input images whose sky value will be set to (reported to be) 0.
          This setting will "equalize" sky values between the images in
          large mosaics. However, this method is not recommended when used
          in conjunction with `AstroDrizzle <http://stsdas.stsci.edu/stsci_python_sphinxdocs_2.13/drizzlepac/astrodrizzle.html>`_
          because it computes relative sky values while ``AstroDrizzle`` needs
          "measured" sky values for median image generation and CR rejection.

        * *globalmin+match* : first find a minimum "global" sky value
          in all input images and then use 'match' method to
          equalize sky values between images.

          .. note::
              This is the *recommended* setting for images
              containing diffuse sources (e.g., galaxies, nebulae)
              covering significant parts of the image.

    skywidth : float, optional (Default Value = 0.1)
        Bin width, in sigma, used to sample the distribution of pixel flux values 
        in order to compute the sky background statistics.


    skystat : {'median', 'mode', 'mean'}, optional (Default Value = 'median')
        Statistical method for determining the sky value from the image pixel values.


    skylower : float, optional (Default Value = INDEF)
        Lower limit of usable pixel values for computing the sky. This value 
        should be specified in the units of the input image.


    skyupper : float, optional (Default Value = INDEF)
        Upper limit of usable pixel values for computing the sky. This value 
        should be specified in the units of the input image.


    skyclip : int, optional (Default Value = 5)
        Number of clipping iterations to use when computing the sky value.


    skylsigma : float, optional (Default Value = 4.0)
        Lower clipping limit, in sigma, used when computing the sky value.


    skyusigma : float, optional (Default Value = 4.0)
        Upper clipping limit, in sigma, used when computing the sky value.


    skymask_cat : str, optional (Default Value = '')
        File name of a catalog file listing user masks to be used with images.


    use_static : bool, optional (Default Value = True)
        Specifies whether or not to use static mask to exclude masked image 
        pixels from sky computations.


    sky_bits : int, None, optional (Default = 0)
        Integer sum of all the DQ bit values from the input image's DQ array 
        that should be considered "good" when building masks for sky computations. 
        For example, if pixels in the DQ array can be combinations of 1, 2, 4, 
        and 8 flags and one wants to consider DQ "defects" having flags 2 and 4 
        as being acceptable for sky computations, then ``sky_bits`` should be 
        set to 2+4=6. Then a DQ pixel having values 2,4, or 6 will be considered 
        a good pixel, while a DQ pixel with a value, e.g., 1+2=3, 4+8=12, etc. 
        will be flagged as a "bad" pixel.

        Default value (0) will make *all* non-zero pixels in the DQ mask to be
        considered "bad" pixels, and the corresponding image pixels will not be 
        used for sky computations.

        Set ``sky_bits`` to ``None`` to turn off the use of image's DQ array 
        for sky computations.

        .. note::
            DQ masks (if used), *will be* combined with user masks specified 
            in the input @-file.


    skyfile : str, optional (Default Value = '')
        Name of file containing user-computed sky values to be used with each input
        image. This ASCII file should only contain 2 columns: image filename in
        column 1 and sky value in column 2. The sky value should be provided in
        units that match the units of the input image and for multi-chip images,
        the same value will be applied to all chips.


    skyuser : str (Default = '')
        Name of header keyword which records the sky value already subtracted
        from the image by the user. The ``skyuser`` parameter is ignored when
        ``skysub`` is set to ``yes``.

        .. note::
            When ``skysub``=``no`` and ``skyuser`` field is empty, then ``AstroDrizzle``
            will assume that sky background is 0.0 for the purpose of cosmic-ray
            rejection.


    in_memory : bool, optional (Default Value = False)
        Specifies whether to optimize execution for speed (maximum memory usage) or
        use a balanced approach in which a minimal amount of image data is kept in
        memory and retrieved from disk as needed. The default setting is
        recommended for most systems.


    **Further Notes:**

    :py:func:`sky` provides new algorithms for sky value computations
    and enhances previously available algorithms used by, e.g.,
    `Astrodrizzle <http://stsdas.stsci.edu/stsci_python_sphinxdocs_2.13/drizzlepac/astrodrizzle.html>`_.

    First, the standard sky computation algorithm
    (see skymethod = 'localmin') was upgraded to be able to use
    DQ flags and user supplied masks to remove "bad" pixels from being
    used for sky statistics computations.

    Second, two new methods have been introduced: ``'globalmin'`` and
    ``'match'``, as well as a combination of the two -- ``'globalmin+match'``.

    - The ``'globalmin'`` method computes the minimum sky value across *all*
      chips in *all* input images. That sky value is then considered to be
      the background in all input images.

    - The ``'match'`` algorithm is somewhat similar to the traditional sky
      subtraction method (skymethod ='localmin') in the sense that it
      measures the sky independently in input images (or detector chips). The
      major differences are that, unlike the traditional method,

        * ``'match'`` algorithm computes *relative* sky values with regard
          to the sky in a reference image chosen from the input list
          of images; *and*

        * Sky statistics is computed only in the part of the image
          that intersects other images.

    This makes ``'match'`` sky computation algorithm particularly useful
    for "equalizing" sky values in large mosaics in which one may have
    only (at least) pair-wise intersection of images without having
    a common intersection region (on the sky) in all images.

    The ``'match'`` method works in the following way: for each pair
    of intersecting images, an equation is written that
    requires that average surface brightness in the overlapping part of
    the sky be equal in both images. The final system of equations is then
    solved for unknown background levels.

    .. warning::

        Current algorithm is not capable of detecting cases when some groups of
        intersecting images (from the input list of images) do not intersect
        at all other groups of intersecting images (except for the simple
        case when *single* images do not intersect any other images). In these
        cases the algorithm will find equalizing sky values for each group.
        However since these groups of images do not intersect each other,
        sky will be matched only within each group and the "inter-group"
        sky mismatch could be significant.

        Users are responsible for detecting such cases and adjusting processing
        accordingly.

    .. warning::

        Because this method computes *relative sky values* compared to a
        reference image (which will have its sky value set to 0), the sky
        values computed with this method usually are smaller than the
        "absolute" sky values computed, e.g., with the ``'localmin'``
        algorithm. Since `AstroDrizzle <http://stsdas.stsci.edu/stsci_python_sphinxdocs_2.13/drizzlepac/astrodrizzle.html>`_ expects
        "true" (as opposite to *relative*) sky values in order to
        correctly compute the median image or to perform cosmic-ray
        detection, this algorithm in not recommended to be used *alone*
        for sky computations to be used with ``AstroDrizzle``.

        For the same reason, IVM weighting in ``AstroDrizzle`` should **not**
        be used with ``'match'`` method: sky values reported in ``MDRIZSKY``
        header keyword will be relative sky values (sky offsets) and derived
        weights will be incorrect.

    The ``'globalmin+match'`` algorithm combines ``'match'`` and
    ``'globalmin'`` methods in order to overcome the limitation of the
    ``'match'`` method described in the note above: it uses ``'globalmin'``
    algorithm to find a baseline sky value common to all input images
    and the ``'match'`` algorithm to "equalize" sky values in the mosaic.
    Thus, the sky value of the "reference" image will be equal to the
    baseline sky value (instead of 0 in ``'match'`` algorithm alone)
    making this method acceptable for use in conjunction with
    ``AstroDrizzle``.

    **Glossary:**

    *Exposure* -- a *subset* of FITS image extensions in an input image
    that correspond to different chips in the detector used to acquire
    the image. The subset of image extensions that form an exposure
    is defined by specifying extensions to be used with input images
    (see parameter ``input``).

    See help for :py:func:`~stsci.skypac.parseat.parse_at_line` for details
    on how to specify image extensions.

    **Footprint** -- the outline (edge) of the projection of a chip or
    of an exposure on the celestial sphere.

    .. note::

        * Footprints are managed by the
          `spherical_geometry.polygon.SphericalPolygon
          <https://spherical-geometry.readthedocs.io/en/latest/api/spherical_geometry.polygon.SphericalPolygon.html>`_
          class.

        * Both footprints *and* associated exposures (image data, WCS
          information, and other header information) are managed by the
          :py:class:`~stsci.skypac.skyline.SkyLine` class.

        * Each :py:class:`~stsci.skypac.skyline.SkyLine` object contains one or more
          :py:class:`~stsci.skypac.skyline.SkyLineMember` objects that manage
          both footprints *and* associated *chip* data that form an exposure.

    **Remarks:**
    
    * :py:func:`sky` works directly on *geometrically distorted*
      flat-fielded images thus avoiding the need to perform an additional
      drizzle step to perform distortion correction of input images.    
      Initially, the footprint of a chip in an image is approximated by a
      2D planar rectangle representing the borders of chip's distorted
      image. After applying distortion model to this rectangle and
      projecting it onto the celestial sphere, it is approximated by
      spherical polygons. Footprints of exposures and mosaics are
      computed as unions of such spherical polygons while overlaps
      of image pairs are found by intersecting these spherical polygons.

    **Limitations and Discussions:**
    Primary reason for introducing "sky match" algorithm was to try to
    equalize the sky in large mosaics in which computation of the
    "absolute" sky is difficult due to the presence of large diffuse
    sources in the image. As discussed above, :py:func:`sky`
    accomplishes this by comparing "sky values" in a pair of images in the
    overlap region (that is common to both images). Quite obviously the
    quality of sky "matching" will depend on how well these "sky values"
    can be estimated. We use quotation marks around *sky values* because
    for some image "true" background may not be present at all and the
    measured sky may be the surface brightness of large galaxy, nebula, etc.

    Here is a brief list of possible limitations/factors that can affect
    the outcome of the matching (sky subtraction in general) algorithm:

    * Since sky subtraction is performed on *flat-fielded* but
        *not distortion corrected* images, it is important to keep in mind
        that flat-fielding is performed to obtain uniform surface brightness
        and not flux. This distinction is important for images that have
        not been distortion corrected. As a consequence, it is advisable that
        point-like sources be masked through the user-supplied mask files.
        Alternatively, one can use ``upper`` parameter to limit the use of
        bright objects in sky computations.

    * Normally, distorted flat-fielded images contain cosmic rays. This
        algorithm does not perform CR cleaning. A possible way of minimizing
        the effect of the cosmic rays on sky computations is to use
        clipping (``nclip`` > 0) and/or set ``upper`` parameter to a value
        larger than most of the sky background (or extended source) but
        lower than the values of most CR pixels.

    * In general, clipping is a good way of eliminating "bad" pixels:
        pixels affected by CR, hot/dead pixels, etc. However, for
        images with complicated backgrounds (extended galaxies, nebulae,
        etc.), affected by CR and noise, clipping process may mask different
        pixels in different images. If variations in the background are
        too strong, clipping may converge to different sky values in
        different images even when factoring in the "true" difference
        in the sky background between the two images.

    * In general images can have different "true" background values
        (we could measure it if images were not affected by large diffuse
        sources). However, arguments such as ``lower`` and ``upper`` will
        apply to all images regardless of the intrinsic differences
        in sky levels.

    **How to use the tasks stand alone interface in your own scripts:**
    These tasks are designed to work together seemlessly when run in the
    full ``AstroDrizzle`` interface. More advanced users may wish to create
    specialized scripts for their own datasets, making use of only a subset
    of the predefined ``AstroDrizzle`` tasks, or add additional processing,
    which may be usefull for their particular data. In these cases,
    individual access to the tasks is important.

    Something to keep in mind is that the full ``AstroDrizzle`` interface will
    make backup copies of your original files and place them in the ``OrIg/``
    directory of your current working directory. If you are working with
    the stand alone interfaces, it is assumed that the user has already taken
    care of backing up their original datafiles as the input file with be
    directly altered.

    Examples
    --------
    Basic example of how to call sky yourself from a Python command line,
    this example will use the default parameter settings and subtract a sky
    value from each ``*flt.fits`` image in the current directory,
    saving the output file with the extension of "mysky":

    >>> from drizzlepac import sky
    >>> sky.sky('*flt.fits',outExt='mysky')
    """

    if input is not None:
        inputDict['input']=input
        inputDict['output']=None
        inputDict['updatewcs']=False
        inputDict['group']=group
    else:
        print("Please supply an input image", file=sys.stderr)
        raise ValueError

    configObj = util.getDefaultConfigObj(__taskname__,configObj,inputDict,loadOnly=(not editpars))
    if configObj is None:
        return

    if not editpars:
        run(configObj,outExt=outExt)


# this is the function that will be called from TEAL
def run(configObj,outExt=None):

    #now we really just need the imageObject list created for the dataset
    filelist,output,ivmlist,oldasndict=processInput.processFilenames(configObj['input'],None)

    imageObjList=processInput.createImageObjectList(filelist,instrpars={},group=configObj['group'])

    #set up the output names, if no extension given the default will be used
    #otherwise, the user extension is used and if the file already exists it's overwritten
    saveFile = False
    if(outExt not in [None,'','None']):
        saveFile = True
        for image in imageObjList:
            outsky = image.outputNames['outSky']
            if outExt not in outsky:
                outsky = outsky.replace("sky",outExt)
                image.outputNames['outSky']=outsky
                log.info(outsky)

    subtractSky(imageObjList,configObj,saveFile=saveFile)


# this is the workhorse looping function
def subtractSky(imageObjList,configObj,saveFile=False,procSteps=None):
    # if neither 'skyfile' nor 'skyuser' are specified, subtractSky will
    # call _skymatch to perform "sky background matching". When 'skyuser'
    # is specified, subtractSky will call the old _skysub.

    if procSteps is not None:
        procSteps.addStep(PROCSTEPS_NAME)

    # General values to use
    step_name = util.getSectionName(configObj, STEP_NUM)
    paramDict = configObj[step_name]

    if not util.getConfigObjPar(configObj, 'skysub'):
        log.info('Sky Subtraction step not performed.')
        _addDefaultSkyKW(imageObjList)
        if 'skyuser' in paramDict and not util.is_blank(paramDict['skyuser']):
            kwd = paramDict['skyuser'].lstrip()
            if kwd[0] == '@':
                # user's sky values are in a file:
                log.info("Retrieving user computed sky values from file '{}'"
                         .format(kwd[1:]))
                _skyUserFromFile(imageObjList, kwd[1:],apply_sky=False)
            else:
                # user's sky values are stored in a header keyword:
                log.info("Retrieving user computed sky values from image "
                         "headers ")
                log.info("recorded in the '{:s}' header keywords."
                         .format(paramDict['skyuser']))
                for image in imageObjList:
                    log.info('Working on sky for: %s' % image._filename)
                    _skyUserFromHeaderKwd(image, paramDict)
        else:
            # reset "computedSky" chip's attribute:
            for image in imageObjList:
                numchips    = image._numchips
                extname     = image.scienceExt
                for extver in range(1, numchips + 1, 1):
                    chip = image[extname, extver]
                    if not chip.group_member:
                        continue
                    chip.computedSky = None

        if procSteps is not None:
            procSteps.endStep(PROCSTEPS_NAME)
        return

    #get the sub-dictionary of values for this step alone and print them out
    log.info('USER INPUT PARAMETERS for Sky Subtraction Step:')
    util.printParams(paramDict, log=log)
    if 'skyfile' in paramDict and not util.is_blank(paramDict['skyfile']):
        _skyUserFromFile(imageObjList,paramDict['skyfile'])
    else:
        # in_memory:
        if 'in_memory' in configObj:
            inmemory = configObj['in_memory']
        elif len(imageObjList) > 0 and imageObjList[0].inmemory is not None:
            inmemory = imageObjList[0].inmemory
        else:
            inmemory = False
        # clean:
        if 'STATE OF INPUT FILES' in configObj and \
           'clean' in configObj['STATE OF INPUT FILES']:
            clean = configObj['STATE OF INPUT FILES']['clean']
        else:
            clean = True

        _skymatch(imageObjList, paramDict, inmemory, clean, log)

    if procSteps is not None:
        procSteps.endStep(PROCSTEPS_NAME)


def _skymatch(imageList, paramDict, in_memory, clean, logfile):
    # '_skymatch' converts input imageList and other parameters to
    # data structures accepted by the "skymatch" package.
    # It also creates a temporary mask by combining 'static' mask,
    # DQ image, and user-supplied mask. The combined mask is then
    # passed to 'skymatch' to be used for excluding "bad" pixels.

    #header keyword that contains the sky that's been subtracted
    skyKW = "MDRIZSKY"

    nimg = len(imageList)
    if nimg == 0:
        log.info("Skymatch needs at least one image to perform{0} \
                    sky matching. Nothing to be done.",os.linesep)
        return

    # create a list of input file names as provided by the user:
    user_fnames   = []
    loaded_fnames = []
    filemaskinfos = nimg * [ None ]
    for img in imageList:
        user_fnames.append(img._original_file_name)
        loaded_fnames.append(img._filename)

    # parse sky mask catalog file (if any):
    catfile = paramDict['skymask_cat']
    if catfile:
        #extname = imageList[0].scienceExt
        #assert(extname is not None and extname != '')
        catfile = catfile.strip()
        mfindx = parse_at_file(fname = catfile,
                               default_ext = ('SCI','*'),
                               default_mask_ext = 0,
                               clobber      = False,
                               fnamesOnly   = True,
                               doNotOpenDQ  = True,
                               match2Images = user_fnames,
                               im_fmode     = 'update',
                               dq_fmode     = 'readonly',
                               msk_fmode    = 'readonly',
                               logfile      = MultiFileLog(console=True),
                               verbose      = True)
        for p in mfindx:
            filemaskinfos[p[1]] = p[0]

    # step through the list of input images and create
    # combined (static + DQ + user supplied, if any) mask, and
    # create a list of FileExtMaskInfo objects to be passed
    # to 'skymatch' function.
    #
    # This needs to be done in several steps, mostly due to the fact that
    # the mask catalogs use "original" (e.g., GEIS, WAIVER FITS) file names
    # while ultimately we want to open the version converted to MEF. Second
    # reason is that we want to combine user supplied masks with DQ+static
    # masks provided by astrodrizzle.
    new_fi = []
    sky_bits = interpret_bit_flags(paramDict['sky_bits'])
    for i in range(nimg):
        # extract extension information:
        extname = imageList[i].scienceExt
        extver  = imageList[i].group
        if extver is None:
            extver = imageList[i].getExtensions()
        assert(extname is not None and extname != '')
        assert(extver)

        # create a new FileExtMaskInfo object
        fi = FileExtMaskInfo(default_ext=(extname,'*'),
                             default_mask_ext=0,
                             clobber=False,
                             doNotOpenDQ=True,
                             fnamesOnly=False,
                             im_fmode='update',
                             dq_fmode='readonly',
                             msk_fmode='readonly')

        # set image file and extensions:
        fi.image = loaded_fnames[i]
        extlist  = [ (extname,ev) for ev in extver ]
        fi.append_ext(extlist)

        # set user masks if any (this will open the files for a later use):
        fi0 = filemaskinfos[i]
        if fi0 is not None:
            nmask = len(fi0.mask_images)
            for m in range(nmask):
                mask = fi0.mask_images[m]
                ext  = fi0.maskext[m]
                fi.append_mask(mask, ext)
        fi.finalize()

        # combine user masks with static masks:
        assert(len(extlist) == fi.count)

        masklist = []
        mextlist = []

        for k in range(fi.count):
            if fi.mask_images[k].closed:
                umask = None
            else:
                umask = fi.mask_images[k].hdu[fi.maskext[k]].data
            (mask, mext) = _buildStaticDQUserMask(imageList[i], extlist[k],
                               sky_bits, paramDict['use_static'],
                               fi.mask_images[k], fi.maskext[k], in_memory)

            masklist.append(mask)
            mextlist.append(mext)

        # replace the original user-supplied masks with the
        # newly computed combined static+DQ+user masks:
        fi.clear_masks()
        for k in range(fi.count):
            if in_memory and mask is not None:
                # os.stat() on the "original_fname" of the mask will fail
                # since this is a "virtual" mask. Therefore we need to compute
                # mask_stat ourselves. We will simply use id(data) for this:
                mstat = os.stat_result((0,id(mask.hdu)) + 8*(0,))
                fi.append_mask(masklist[k], mextlist[k], mask_stat=mstat)
            else:
                fi.append_mask(masklist[k], mextlist[k])
            if masklist[k]:
                masklist[k].release()
        fi.finalize()

        new_fi.append(fi)

    try:
        # Run skymatch algorithm:
        skymatch(new_fi,
                 skymethod   = paramDict['skymethod'],
                 skystat     = paramDict['skystat'],
                 lower       = paramDict['skylower'],
                 upper       = paramDict['skyupper'],
                 nclip       = paramDict['skyclip'],
                 lsigma      = paramDict['skylsigma'],
                 usigma      = paramDict['skyusigma'],
                 binwidth    = paramDict['skywidth'],
                 skyuser_kwd = skyKW,
                 units_kwd   = 'BUNIT',
                 readonly    = not paramDict['skysub'],
                 dq_bits     = None,
                 optimize    = 'inmemory' if in_memory else 'balanced',
                 clobber     = True,
                 clean       = clean,
                 verbose     = True,
                 flog        = MultiFileLog(console = True, enableBold = False),
                 _taskname4history = 'AstroDrizzle')
    except Exception:
        if 'match' in paramDict['skymethod']:  # This catches 'match' and 'globalmin+match'
            new_method = 'globalmin' if 'globalmin' in paramDict['skymethod'] else 'localmin'

            # revert to simpler sky computation algorithm
            log.warning('Reverting sky computation to "localmin" from "{}'.format(paramDict['skymethod']))
            skymatch(new_fi,
                     skymethod=new_method,
                     skystat=paramDict['skystat'],
                     lower=paramDict['skylower'],
                     upper=paramDict['skyupper'],
                     nclip=paramDict['skyclip'],
                     lsigma=paramDict['skylsigma'],
                     usigma=paramDict['skyusigma'],
                     binwidth=paramDict['skywidth'],
                     skyuser_kwd=skyKW,
                     units_kwd='BUNIT',
                     readonly=not paramDict['skysub'],
                     dq_bits=None,
                     optimize='inmemory' if in_memory else 'balanced',
                     clobber=True,
                     clean=clean,
                     verbose=True,
                     flog=MultiFileLog(console=True, enableBold=False),
                     _taskname4history='AstroDrizzle')
        else:
            raise

    # Populate 'subtractedSky' and 'computedSky' of input image objects:
    for i in range(nimg):
        assert(not new_fi[i].fnamesOnly and not new_fi[i].image.closed)
        image = imageList[i]
        skysubimage = new_fi[i].image.hdu
        numchips = image._numchips
        extname = image.scienceExt
        assert(os.path.samefile(image._filename, skysubimage.filename()))

        for extver in range(1, numchips + 1, 1):
            chip = image[extname, extver]
            if not chip.group_member:
                continue
            subtracted_sky = skysubimage[extname, extver].header.get(skyKW, 0.)
            chip.subtractedSky = subtracted_sky
            chip.computedSky = subtracted_sky

    # clean-up:
    for fi in new_fi:
        fi.release_all_images()

def _buildStaticDQUserMask(img, ext, sky_bits, use_static, umask,
                           umaskext, in_memory):
    # creates a temporary mask by combining 'static' mask,
    # DQ image, and user-supplied mask.

    def merge_masks(m1, m2):
        if m1 is None: return m2
        if m2 is None: return m1
        return np.logical_and(m1, m2).astype(np.uint8)

    mask = None

    # build DQ mask
    if sky_bits is not None:
        mask = img.buildMask(img[ext]._chip,bits=sky_bits)

    # get correct static mask mask filenames/objects
    staticMaskName = img[ext].outputNames['staticMask']
    smask = None
    if use_static:
        if img.inmemory:
            if staticMaskName in img.virtualOutputs:
                smask = img.virtualOutputs[staticMaskName].data
        else:
            if staticMaskName is not None and os.path.isfile(staticMaskName):
                sm, dq = openImageEx(
                    staticMaskName,
                    mode='readonly',
                    memmap=False,
                    saveAsMEF=False,
                    clobber=False,
                    imageOnly=True,
                    openImageHDU=True,
                    openDQHDU=False,
                    preferMEF=False,
                    verbose=False
                )
                if sm.hdu is not None:
                    smask = sm.hdu[0].data
                    sm.release()
            else:
                log.warning("Static mask for file \'{}\', ext={} NOT FOUND." \
                            .format(img._filename, ext))
        # combine DQ and static masks:
        mask = merge_masks(mask, smask)

    # combine user mask with the previously computed mask:
    if umask is not None and not umask.closed:
        if mask is None:
            # return user-supplied mask:
            umask.hold()
            return (umask, umaskext)
        else:
            # combine user mask with the previously computed mask:
            dtm  = umask.hdu[umaskext].data
            mask = merge_masks(mask, dtm)

    if mask is None:
        return (None, None)
    elif mask.sum() == 0:
        log.warning("All pixels masked out when applying DQ, " \
                    "static, and user masks!")

    # save mask to a temporary file:
    (root,suffix,fext) = file_name_components(img._filename)
    if in_memory:
        tmpmask = in_memory_mask(mask)
        strext = ext2str(ext, compact=True, default_extver=None)
        tmpmask.original_fname = "{1:s}{0:s}{2:s}{0:s}{3:s}" \
            .format('_', root, suffix, 'in-memory_skymatch_mask')
    else:
        (tmpfname, tmpmask) = temp_mask_file(mask, root,
            prefix='', suffix='skymatch_mask', ext=ext,
            randomize_prefix=False)
        img[ext].outputNames['skyMatchMask'] = tmpfname

    return (tmpmask, 0)

# this function applies user supplied sky values from an input file
def _skyUserFromFile(imageObjList, skyFile, apply_sky=None):
    """
    Apply sky value as read in from a user-supplied input file.

    """
    skyKW="MDRIZSKY" #header keyword that contains the sky that's been subtracted

    # create dict of fname=sky pairs
    skyvals = {}
    if apply_sky is None:
        skyapplied = False # flag whether sky has already been applied to images
    else:
        skyapplied = apply_sky

    for line in open(skyFile):
        if apply_sky is None and line[0] == '#' and 'applied' in line:
            if '=' in line: linesep = '='
            if ':' in line: linesep = ':'
            appliedstr = line.split(linesep)[1].strip()
            if appliedstr.lower() in ['yes','true','y','t']:
                skyapplied = True
                print('...Sky values already applied by user...')

        if not util.is_blank(line) and line[0] != '#':
            lspl = line.split()
            svals = []
            for lvals in lspl[1:]:
                svals.append(float(lvals))
            skyvals[lspl[0]] = svals

    # Apply user values to appropriate input images
    for imageSet in imageObjList:
        fname = imageSet._filename
        numchips=imageSet._numchips
        sciExt=imageSet.scienceExt
        if fname in skyvals:
            print("    ...updating MDRIZSKY with user-supplied value.")
            for chip in range(1,numchips+1,1):
                if len(skyvals[fname]) == 1:
                    _skyValue = skyvals[fname][0]
                else:
                    _skyValue = skyvals[fname][chip-1]

                chipext = '%s,%d'%(sciExt,chip)
                _updateKW(imageSet[chipext],fname,(sciExt,chip),skyKW,_skyValue)

                # Update internal record with subtracted sky value
                #
                # .computedSky:   value to be applied by the
                #                 adrizzle/ablot steps.
                # .subtractedSky: value already (or will be by adrizzle/ablot)
                #                 subtracted from the image
                if skyapplied:
                    imageSet[chipext].computedSky = None # used by adrizzle/ablot
                else:
                    imageSet[chipext].computedSky = _skyValue
                imageSet[chipext].subtractedSky = _skyValue
                print("Setting ",skyKW,"=",_skyValue)
        else:
            print("*"*40)
            print("*")
            print("WARNING:")
            print("    .... NO user-supplied sky value found for ",fname)
            print("    .... Setting sky to a value of 0.0! ")
            print("*")
            print("*"*40)

def _skyUserFromHeaderKwd(imageSet,paramDict):
    """
    subtract the sky from all the chips in the imagefile that imageSet represents

    imageSet is a single imageObject reference
    paramDict should be the subset from an actual config object

    """
    _skyValue=0.0    #this will be the sky value computed for the exposure
    skyKW="MDRIZSKY" #header keyword that contains the sky that's been subtracted

    #just making sure, tricky users and all, these are things that will be used
    #by the sky function so we want them defined at least
    try:
        assert imageSet._numchips > 0, "invalid value for number of chips"
        assert imageSet._filename != '', "image object filename is empty!, doh!"
        assert imageSet._rootname != '', "image rootname is empty!, doh!"
        assert imageSet.scienceExt !='', "image object science extension is empty!"

    except AssertionError:
        raise AssertionError

    numchips=imageSet._numchips
    sciExt=imageSet.scienceExt

    # User Subtraction Case, User has done own sky subtraction,
    # so use the image header value for subtractedsky value
    skyuser=paramDict["skyuser"]

    if skyuser != '':
        print("User has computed their own sky values...")

        if skyuser != skyKW:
            print("    ...updating MDRIZSKY with supplied value.")
            for chip in range(1,numchips+1,1):
                chipext = '%s,%d'%(sciExt,chip)
                if not imageSet[chipext].group_member:
                    # skip extensions/chips that will not be processed
                    continue
                try:
                    _skyValue = imageSet[chipext].header[skyuser]
                except:
                    print("**************************************************************")
                    print("*")
                    print("*  Cannot find keyword ",skyuser," in ",imageSet._filename)
                    print("*")
                    print("**************************************************************\n\n\n")
                    raise KeyError

                _updateKW(imageSet[sciExt+','+str(chip)],
                          imageSet._filename,(sciExt,chip),skyKW,_skyValue)

                # Update internal record with subtracted sky value
                imageSet[chipext].subtractedSky = _skyValue
                imageSet[chipext].computedSky = None
                print("Setting ",skyKW,"=",_skyValue)

# this is the main function that does all the real work in computing the
# statistical sky value for each image (set of chips)
# mcara: '_skySub' is obsolete now:
#        was replaced with '_skyUserFromHeaderKwd' and '_skymatch'
def _skySub(imageSet,paramDict,saveFile=False):
    """
    subtract the sky from all the chips in the imagefile that imageSet represents

    imageSet is a single imageObject reference
    paramDict should be the subset from an actual config object
    if saveFile=True, then images that have been sky subtracted are saved to a predetermined output name
    else, overwrite the input images with the sky-subtracted results

    the output from sky subtraction is a copy of the original input file
    where all the science data extensions have been sky subtracted

    """

    _skyValue=0.0    #this will be the sky value computed for the exposure
    skyKW="MDRIZSKY" #header keyword that contains the sky that's been subtracted

    #just making sure, tricky users and all, these are things that will be used
    #by the sky function so we want them defined at least
    try:
        assert imageSet._numchips > 0, "invalid value for number of chips"
        assert imageSet._filename != '', "image object filename is empty!, doh!"
        assert imageSet._rootname != '', "image rootname is empty!, doh!"
        assert imageSet.scienceExt !='', "image object science extension is empty!"

    except AssertionError:
        raise AssertionError

    numchips=imageSet._numchips
    sciExt=imageSet.scienceExt

    # User Subtraction Case, User has done own sky subtraction,
    # so use the image header value for subtractedsky value
    skyuser=paramDict["skyuser"]

    if skyuser != '':
        print("User has computed their own sky values...")

        if skyuser != skyKW:
            print("    ...updating MDRIZSKY with supplied value.")
            for chip in range(1,numchips+1,1):
                try:
                    chipext = '%s,%d'%(sciExt,chip)
                    _skyValue = imageSet[chipext].header[skyuser]

                except:
                    print("**************************************************************")
                    print("*")
                    print("*  Cannot find keyword ",skyuser," in ",imageSet._filename)
                    print("*")
                    print("**************************************************************\n\n\n")
                    raise KeyError

                _updateKW(imageSet[sciExt+','+str(chip)],imageSet._filename,(sciExt,chip),skyKW,_skyValue)

                # Update internal record with subtracted sky value
                imageSet[chipext].subtractedSky = _skyValue
                imageSet[chipext].computedSky = None
                print("Setting ",skyKW,"=",_skyValue)

    else:
        # Compute our own sky values and record the values for use later.
        # The minimum sky value from all the  science chips in the exposure
        # is used as the reference sky for each chip

        log.info("Computing minimum sky ...")
        minSky=[] #store the sky for each chip
        minpscale = []

        for chip in range(1,numchips+1,1):
            myext=sciExt+","+str(chip)

            #add the data back into the chip, leave it there til the end of this function
            imageSet[myext].data=imageSet.getData(myext)

            image=imageSet[myext]
            _skyValue= _computeSky(image, paramDict, memmap=False)
            #scale the sky value by the area on sky
            # account for the case where no IDCSCALE has been set, due to a
            # lack of IDCTAB or to 'coeffs=False'.
            pscale=imageSet[myext].wcs.idcscale
            if pscale is None:
                log.warning("No Distortion coefficients available...using "
                            "default plate scale.")
                pscale = imageSet[myext].wcs.pscale
            _scaledSky=_skyValue / (pscale**2)
            #_skyValue=_scaledSky
            minSky.append(_scaledSky)
            minpscale.append(pscale)

        _skyValue = min(minSky)

        _reportedSky = _skyValue*(minpscale[minSky.index(_skyValue)]**2)
        log.info("Minimum sky value for all chips %s" % _reportedSky)

        #now subtract that value from all the chips in the exposure
        #and update the chips header keyword with the sub
        for chip in range(1,numchips+1,1):
            image=imageSet[sciExt,chip]
            myext = sciExt+","+str(chip)
            # account for the case where no IDCSCALE has been set, due to a
            # lack of IDCTAB or to 'coeffs=False'.
            idcscale = image.wcs.idcscale
            if idcscale is None: idcscale = image.wcs.pscale
            _scaledSky=_skyValue * (idcscale**2)
            image.subtractedSky = _scaledSky
            image.computedSky = _scaledSky
            log.info("Using sky from chip %d: %f\n" % (chip,_scaledSky))
            ###_subtractSky(image,(_scaledSky))
            # Update the header so that the keyword in the image is
            #the sky value which should be subtracted from the image
            _updateKW(image,imageSet._filename,(sciExt,chip),skyKW,_scaledSky)


###############################
##  Helper functions follow  ##
###############################

def _computeSky(image, skypars, memmap=False):

    """
    Compute the sky value for the data array passed to the function
    image is a fits object which contains the data and the header
    for one image extension

    skypars is passed in as paramDict

    """
    #this object contains the returned values from the image stats routine
    _tmp = imagestats.ImageStats(image.data,
            fields      = skypars['skystat'],
            lower       = skypars['skylower'],
            upper       = skypars['skyupper'],
            nclip       = skypars['skyclip'],
            lsig        = skypars['skylsigma'],
            usig        = skypars['skyusigma'],
            binwidth    = skypars['skywidth']
            )

    _skyValue = _extractSkyValue(_tmp,skypars['skystat'].lower())
    log.info("    Computed sky value/pixel for %s: %s "%
             (image.rootname, _skyValue))

    del _tmp

    return _skyValue


def _extractSkyValue(imstatObject,skystat):
    if (skystat =="mode"):
        return imstatObject.mode
    elif (skystat == "mean"):
        return imstatObject.mean
    else:
        return imstatObject.median


def _subtractSky(image,skyValue,memmap=False):
    """
    subtract the given sky value from each the data array
    that has been passed. image is a fits object that
    contains the data and header for one image extension
    """
    try:
        np.subtract(image.data,skyValue,image.data)

    except IOError:
        print("Unable to perform sky subtraction on data array")
        raise IOError


def _updateKW(image, filename, exten, skyKW, Value):
    """update the header with the kw,value"""
    # Update the value in memory
    image.header[skyKW] = Value

    # Now update the value on disk
    if isinstance(exten,tuple):
        strexten = '[%s,%s]'%(exten[0],str(exten[1]))
    else:
        strexten = '[%s]'%(exten)
    log.info('Updating keyword %s in %s' % (skyKW, filename + strexten))
    fobj = fileutil.openImage(filename, mode='update', memmap=False)
    fobj[exten].header[skyKW] = (Value, 'Sky value computed by AstroDrizzle')
    fobj.close()

def _addDefaultSkyKW(imageObjList):
    """Add MDRIZSKY keyword to "commanded" SCI headers of all input images,
        if that keyword does not already exist.
    """
    skyKW = "MDRIZSKY"
    Value = 0.0
    for imageSet in imageObjList:
        fname = imageSet._filename
        numchips=imageSet._numchips
        sciExt=imageSet.scienceExt
        fobj = fileutil.openImage(fname, mode='update', memmap=False)
        for chip in range(1,numchips+1,1):
            ext = (sciExt,chip)
            if not imageSet[ext].group_member:
                # skip over extensions not used in processing
                continue
            if skyKW not in fobj[ext].header:
                fobj[ext].header[skyKW] = (Value, 'Sky value computed by AstroDrizzle')
                log.info("MDRIZSKY keyword not found in the %s[%s,%d] header."%(
                            fname,sciExt,chip))
                log.info("    Adding MDRIZSKY to header with default value of 0.")
        fobj.close()

# this is really related to each individual chip
# so pass in the image for that chip, image contains header and data
def getreferencesky(image,keyval):

    _subtractedSky=image.header[keyval]
    _refplatescale=image.header["REFPLTSCL"]
    _platescale=image.header["PLATESCL"]

    return (_subtractedSky * (_refplatescale / _platescale)**2 )
