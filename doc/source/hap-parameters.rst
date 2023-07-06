HAP Parameters
-----------------------------------

.. _hap-parameters:

Shown below are the parameters that are used by the Hubble Advanced Products. We include the parameter, the default value for WFC3 processing, and a description of that parameter.

.. run_hap_processing
.. identified in json files. 


Initialize HAP
^^^^^^^^^^^^^^

\_mdriztab_btn\_: Update From MDRIZTAB
    text asdlkfjasldkj falskdjf laskdjf laskd flaksjdfl kas dlfkjal sdfka sdlfkj alsdkjf laksjd laksjdlfka sldk

runfile: "astrodrizzle.log"
    sdflaskj dflaksj dflkasjd flkasjd lfkjasldk flaskjd f

wcskey: ""

proc_unit: "native"

coeffs: true

context: true

group: ""

build: true

crbit: 4096

stepsize: 10

resetbits: "4096"

num_cores: 1

in_memory: false,

Instrument Parameters
^^^^^^^^^^^^^^^^^^^^^

gnkeyword: "ATODGNA,ATODGNB,ATODGNC,ATODGND"

rnkeyword: "READNSEA,READNSEB,READNSEC,READNSED"

expkeyword: "EXPTIME"

gain: ""

rdnoise: ""

exptime: ""

State of input files
^^^^^^^^^^^^^^^^^^^^

restore: false

preserve: false

overwrite: false

clean: true

Step 1: Static mask
^^^^^^^^^^^^^^^^^^^

static: true

static_sig: 4.0

Step 2: Sky Subtraction
^^^^^^^^^^^^^^^^^^^^^^^

skysub: false

skymethod: "match"

skystat: "median"

skywidth: 0.1

skylower: -100.0

sky_bits: "16"

skyupper: null

skyclip: 5

skylsigma: 4.0

skyusigma: 4.0

skymask_cat: ""

use_static: true

skyfile: ""

skyuser: ""

Step 3: Drizzle Separate images
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

driz_separate : bool, default=False
    This parameter specifies whether or not to drizzle each input image onto separate output images. The separate output images will all have the same WCS as the final combined output frame. These images are used to create the median image, needed for cosmic ray rejection.

driz_sep_bits: int, Default = "16"
    Integer sum of all the DQ bit values from the input image's DQ array that should be considered 'good' when building the weighting mask. This can also be used to reset pixels to good if they had been flagged as cosmic rays during a previous run of ``AstroDrizzle``, by adding the value 4096 for ACS and WFPC2 data. Please see the section on Selecting the ``Bits`` Parameter for a more detailed discussion.

driz_sep_kernel : str, default=square
    Used for the initial separate drizzling operation only, this parameter specifies the form of the kernel function used to distribute flux onto the separate output images. The current options are: 'square', 'point', 'gaussian', 'turbo', 'tophat', and 'lanczos3'. See adrizzle.help for more details. 

driz_sep_wt_scl: float, default=image header exposure time 
    This parameter specifies the weighting factor for input image. If ``driz_sep_wt_scl``\ =\ ``exptime``, then the scaling value will be set equal to the exposure time found in the image header. The use of the default value is recommended for producing optimal behavior for most scenarious. It is possible to set ``wt_scl``\ =\ 'expsq' for weighting by the square of the exposure time, which is optimal for read-noise dominated images.

driz_sep_pixfrac: float, default=1.0
    Fraction by which input pixels are "shrunk" before being drizzled onto the output image grid, given as a real number between 0 and 1. This specifies the size of the footprint, or "dropsize", of a pixel in units of the input pixel size. If ``pixfrac`` is set to less than 0.001, the kernel parameter will be reset to 'point' for more efficient processing. In the step of drizzling each input image onto a separate output image, the default value of 1.0 is best in order to ensure that each output drizzled image is fully populated with pixels from the input image. For more information, see the help for the ``drizzle`` task.

.. null?!?!?!?!?

driz_sep_fillval: int or INDEF, default = null
    Value to be assigned to output pixels that have zero weight, or that receive flux from any input pixels during drizzling. This parameter corresponds to the ``fillval`` parameter of the ``drizzle`` task. If the default of ``INDEF`` is used, and if the weight in both the input and output images for a given pixel are zero, then the output pixel will be set to the value it would have had if the input had a non-zero weight. Otherwise, if a numerical value is provided (e.g. 0), then these pixels will be set to that value.

driz_sep_compress: bool, default=False
    Whether to use compression when writing out product.

Step 3a: Custom WCS for Separate Outputs
""""""""""""""""""""""""""""""""""""""""

driz_sep_wcs: bool, default=False
    Define custom WCS for seperate output images?

driz_sep_refimage: str, default=""
    Reference image from which a WCS solution can be obtained.

driz_sep_rot : float, default=null
    Position Angle of output image's Y-axis relative to North.  A value of 0.0 would orient the final output image to be North up.  The default of ``INDEF`` specifies that the images will not be rotated,  but will instead be drizzled in the default orientation for the camera  with the x and y axes of the drizzled image corresponding approximately  to the detector axes. This conserves disk space, as these single  drizzled images are only used in the intermediate step of creating  a median image.

driz_sep_scale : float, default=null
    Linear size of the output pixels in arcseconds/pixel for each separate  drizzled image (used in creating the median for cosmic ray rejection).  The default value of ``INDEF`` specifies that the undistorted pixel  scale for the first input image will be used as the pixel scale for  all the output images.

driz_sep_outnx : int, default=null
    Size, in pixels, of the X axis in the output images that each input  will be drizzled onto. If no value is specified, the smallest size that  can accommodate the full dithered field will be used.

driz_sep_outny : int, default=null
    Size, in pixels, of the Y axis in the output images that each input  will be drizzled onto. If no value is specified, the smallest size  that can accommodate the full dithered field will be used.

driz_sep_ra : float, default=null
    Right ascension (in decimal degrees) specifying the center of the output  image. If this value is not designated, the center will automatically  be calculated based on the distribution of image dither positions.

driz_sep_dec : float, default=null
    Declination (in decimal degrees) specifying the center of the output  image. If this value is not designated, the center will automatically  be calculated based on the distribution of image dither positions.

driz_sep_crpix1: float, or null, default=null
    Reference pixel X position on output (CRPIX1).

driz_sep_crpix2: float, or null, default=null
    Reference pixel Y position on output (CRPIX2).

Step 4: Create Median Image
^^^^^^^^^^^^^^^^^^^^^^^^^^^

median: false
    Create a median image?

median_newmasks: true
    Create new masks when doing the median?

combine_type: "minmed"
    Type of combine operation. "minmed","iminmed","median","mean","imedian","imean","sum".

combine_nlow: 0
    Minmxa, number of low pixels to reject.

combine_nhigh: 1
    Minmxa, number of high pixels to reject.

combine_maskpt: 0.3
    Percentage of weight image value below which it is flagged as a bad pixel.

combine_nsigma: "4 3"
    Significance for accepting minimum instead of median.


combine_lthresh: null
    Lower threshold for clipping input pixel values.

combine_hthresh: null
    Upper threshold for clipping input pixel values.

combine_grow: 1
    Radius (pixels) for neighbor rejection.

combine_bufsize: null
    Size of buffer(in Mb) for each input image.


Step 5: Blot back the median image
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

blot: false

blot_interp: "poly5"

blot_sinscl: 1.0

blot_addsky: true

blot_skyval: 0.0


Step 6: Remove cosmic rays with deriv, driz_cr
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

driz_cr: false

driz_cr_snr: "5.0 4.0"

driz_cr_grow: 1

driz_cr_ctegrow: 0

driz_cr_scale: "3.0 2.4"

driz_cr_corr: false

Step 7: Drizzle final combined image
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

driz_combine: true

final_pixfrac: 1.0

final_fillval: null

final_bits: "65535"

final_maskval: null

final_wht_type: "EXP"

final_kernel: "square"

final_wt_scl: "exptime"

final_units: "cps"

Step 7a: Custom WCS for final output
""""""""""""""""""""""""""""""""""""

final_wcs: true

final_rot: 0.0

final_refimage: ""

final_scale: null

final_outnx: null

final_outny: null

final_ra: null

final_dec: null

final_crpix1: null

final_crpix2: null
