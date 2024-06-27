HAP Parameters
-----------------------------------

.. _hap-parameters:

Shown below are the parameters that are used by the Hubble Advanced Products. We include the parameter, the default value for WFC3 processing, and a description of that parameter.

.. run_hap_processing
.. identified in json files. 


Initialize HAP
^^^^^^^^^^^^^^

\_mdriztab_btn\_: str (*default=Update From MDRIZTAB*)
    Immediately read and use the values from the MDRIZTAB.

runfile: str (*default="astrodrizzle.log"*)
    File for logging the processing.

wcskey: str (*default=""*)
    WCS version to use in processing

proc_unit: str (*default="native"*)
    Units used during processing: "native", "electrons"

coeffs: bool (*default=True*)
    Use header-based distortion coefficients?

context: bool (*default=True*)
    Create context image during final drizzle?

group: str (*default=""*)
    Single extension or group to be combined/cleaned.

build: bool (*default=True*)
    Create multi-extension output file for final drizzle?

crbit: int (*default=4096*)
    Bit value for CR ident. in DQ array.

stepsize: int (*default=10*)
    Step size for drizzle coordinate computation.

resetbits: int (*default="4096"*)
    Bit values to reset in all input DQ arrays.

num_cores: int (*default=1*)
    Max CPU cores to use (n<2 disables, None = auto-decide).

in_memory: bool (*default=False*)
    Process everything in memory to minimize disk I/O?

Instrument Parameters
^^^^^^^^^^^^^^^^^^^^^

.. or float?

gnkeyword: str (*default="ATODGNA,ATODGNB,ATODGNC,ATODGND"*)
    .. the default readnoise/gain value? what are the options?

rnkeyword: str (*default="READNSEA,READNSEB,READNSEC,READNSED"*)
    Keyword used to specify a value, which is used to override the instrument specific default readnoise values. The value is assumed to be in units of electrons. This parameter should not be populated if the ``rdnoise`` parameter is in use.

expkeyword: str (*default="EXPTIME"*)
    Keyword used to specify a value, which is used to override the default exposure time image header values. The value is assumed to be in units of seconds. This parameter should not be populated if the ``exptime`` parameter is in use.

gain: str (*default=""*)
    Value used to override instrument specific default gain values. The value is assumed to be in units of electrons/count. This parameter should not be populated if the ``gainkeyword`` parameter is in use.

rdnoise: str (*default=""*)
    Value used to override instrument specific default readnoise values. The value is assumed to be in units of electrons. This parameter should not be populated if the ``rnkeyword`` parameter is in use.

exptime: str (*default=""*)
    .. ?

State of input files
^^^^^^^^^^^^^^^^^^^^

restore: bool (*default=False*)
    Copy input files FROM archive directory for processing?

preserve: bool (*default=False*)
    Copy input files to archive directory, if not already archived?

overwrite: bool (*default=False*)
    Copy input files into archive, overwriting if required?

clean: bool (*default=True*)
    Delete temporary files after completion?

Step 1: Static mask
^^^^^^^^^^^^^^^^^^^

static: bool (*default=True*)
    Create static bad-pixel mask from the data?

static_sig: float (*default=4.0*)
    Sigma*rms below mode to clip for static mask

Step 2: Sky Subtraction
^^^^^^^^^^^^^^^^^^^^^^^

skysub: bool (*default=False*)
    Turn on or off sky subtraction on the input data. When ``skysub`` is set  to ``no``, then ``skyuser`` field will be enabled and if user specifies a  header keyword showing the sky value in the image, then that value will  be used for CR-rejection but it will not be subtracted from the (drizzled)  image data. If user sets ``skysub`` to ``yes`` then ``skyuser`` field will be  disabled (and if it is not empty - it will be ignored) and user can use  one of the methods available through the ``skymethod`` parameter to  compute the sky or provide a file (see ``skyfile`` parameter) with values  that should be subtracted from (single) drizzled images.

skymethod: str (*default="match"*)
    Sky computation method: "globalmin+match","localmin", "globalmin", "match". See astrodrizzle.help for more details.

skystat: str (*default="median"*)
    Statistical method for determining the sky value from the image pixel values: "median","mode","mean".

skywidth: float (*default=0.1*)
    Bin width of histogram for sampling sky statistics (in sigma)

skylower: float (*default=-100.0*)
    Lower limit of usable data for sky (always in electrons)

sky_bits: str (*default="16"*)
    Integer mask bit values considered good pixels in DQ array

skyupper: int or null (*default=null*)
    Upper limit of usable data for sky (always in electrons)

skyclip: int (*default=5*)
    Number of clipping iterations

skylsigma: float (*default=4.0*)
    Lower side clipping factor (in sigma)

skyusigma: float (*default=4.0*)
    Upper side clipping factor (in sigma)

skymask_cat: str (*default=""*)
    Catalog file listing image masks

use_static: bool (*default=True*)
    Use static mask for skymatch computations?

skyfile: str (*default""*)
    Name of file with user-computed sky values to be subtracted

skyuser: str (*default""*)
    KEYWORD indicating a sky subtraction value if done by user

Step 3: Drizzle Separate images
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

driz_separate : bool (*default=False*)
    This parameter specifies whether or not to drizzle each input image onto separate output images. The separate output images will all have the same WCS as the final combined output frame. These images are used to create the median image, needed for cosmic ray rejection.

driz_sep_bits: int (*default="16"*)
    Integer sum of all the DQ bit values from the input image's DQ array that should be considered "good" when building the weighting mask. This can also be used to reset pixels to good if they had been flagged as cosmic rays during a previous run of ``AstroDrizzle``, by adding the value 4096 for ACS and WFPC2 data. Please see the section on Selecting the ``Bits`` Parameter for a more detailed discussion.

driz_sep_kernel: str (*default="turbo"*)
    Used for the initial separate drizzling operation only, this parameter specifies the form of the kernel function used to distribute flux onto the separate output images. The current options are: "square", "point", "turbo", "gaussian", and "lanczos3". The latter two options ("gaussian" and "lanczos3") are not guaranteed to conserve flux, but may produce reasonable results; understand the effects of these kernels before using them. A former option "tophat" has been removed as it was found to produce poor results. See adrizzle.help for more details. 

driz_sep_wt_scl: float (*default=exposure time (from image header)*)
    This parameter specifies the weighting factor for input image. If ``driz_sep_wt_scl``\ =\ ``exptime``, then the scaling value will be set equal to the exposure time found in the image header. The use of the default value is recommended for producing optimal behavior for most scenarious. It is possible to set ``wt_scl``\ =\ "expsq" for weighting by the square of the exposure time, which is optimal for read-noise dominated images.

driz_sep_pixfrac: float (*default=1.0*)
    Fraction by which input pixels are "shrunk" before being drizzled onto the output image grid, given as a real number between 0 and 1. This specifies the size of the footprint, or "dropsize", of a pixel in units of the input pixel size. If ``pixfrac`` is set to less than 0.001, the kernel parameter will be reset to "point" for more efficient processing. In the step of drizzling each input image onto a separate output image, the default value of 1.0 is best in order to ensure that each output drizzled image is fully populated with pixels from the input image. For more information, see the help for the ``drizzle`` task.

.. null?!?!?!?!?

driz_sep_fillval: int or INDEF (*default=null*)
    Value to be assigned to output pixels that have zero weight, or that receive flux from any input pixels during drizzling. This parameter corresponds to the ``fillval`` parameter of the ``drizzle`` task. If the default of ``INDEF`` is used, and if the weight in both the input and output images for a given pixel are zero, then the output pixel will be set to the value it would have had if the input had a non-zero weight. Otherwise, if a numerical value is provided (e.g. 0), then these pixels will be set to that value.

driz_sep_compress: bool (*default=False*)
    Whether to use compression when writing out product.

Step 3a: Custom WCS for Separate Outputs
""""""""""""""""""""""""""""""""""""""""

driz_sep_wcs: bool (*default=False*)
    Define custom WCS for separate output images?

driz_sep_refimage: str (*default=""*)
    Reference image from which a WCS solution can be obtained.

driz_sep_rot : float or null (*default=null*)
    Position Angle of output image's Y-axis relative to North. A value of 0.0 would orient the final output image to be North up. The default of ``INDEF`` specifies that the images will not be rotated, but will instead be drizzled in the default orientation for the camera with the x and y axes of the drizzled image corresponding approximately to the detector axes. This conserves disk space, as these single drizzled images are only used in the intermediate step of creating a median image.

driz_sep_scale : float or null (*default=null*)
    Linear size of the output pixels in arcseconds/pixel for each separate drizzled image (used in creating the median for cosmic ray rejection). The default value of ``INDEF`` specifies that the undistorted pixel scale for the first input image will be used as the pixel scale for all the output images.

driz_sep_outnx : int or null (*default=null*)
    Size, in pixels, of the X axis in the output images that each input will be drizzled onto. If no value is specified, the smallest size that can accommodate the full dithered field will be used.

driz_sep_outny : int or null (*default=null*)
    Size, in pixels, of the Y axis in the output images that each input will be drizzled onto. If no value is specified, the smallest size that can accommodate the full dithered field will be used.

driz_sep_ra : float or null (*default=null*)
    Right ascension (in decimal degrees) specifying the center of the output image. If this value is not designated, the center will automatically be calculated based on the distribution of image dither positions.

driz_sep_dec : float or null (*default=null*)
    Declination (in decimal degrees) specifying the center of the output image. If this value is not designated, the center will automatically be calculated based on the distribution of image dither positions.

driz_sep_crpix1: float or null (*default=null*)
    Reference pixel X position on output (CRPIX1).

driz_sep_crpix2: float or null (*default=null*)
    Reference pixel Y position on output (CRPIX2).

Step 4: Create Median Image
^^^^^^^^^^^^^^^^^^^^^^^^^^^

median: bool (*default=False*)
    Create a median image?

median_newmasks: bool (*default=True*)
    Create new masks when doing the median?

combine_type: str (*default="minmed"*)
    Type of combine operation. "minmed","iminmed","median","mean","imedian","imean","sum".

combine_nlow: int (*default=0*)
    Minmxa, number of low pixels to reject.

combine_nhigh: int (*default=1*)
    Minmxa, number of high pixels to reject.

combine_maskpt: float (*default=0.3*)
    Percentage of weight image value below which it is flagged as a bad pixel.

combine_nsigma: str (*default="4 3"*)
    Significance for accepting minimum instead of median.


combine_lthresh: ??? (*default=null*)
    Lower threshold for clipping input pixel values.

combine_hthresh: ??? (*default=null*)
    Upper threshold for clipping input pixel values.

combine_grow: int (*default=1*)
    Radius (pixels) for neighbor rejection.

combine_bufsize: ??? (*default=null*)
    Size of buffer(in Mb) for each input image.


Step 5: Blot back the median image
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

blot: bool (*default=False*)
    Blot the median back to the input frame?

blot_interp: str (*default="poly5"*)
    Interpolant (nearest,linear,poly3,poly5,sinc)

blot_sinscl: float (*default=1.0*)
    Scale for sinc interpolation kernel

blot_addsky: bool (*default=True*)
    Add sky using MDRIZSKY value from header?

blot_skyval: float (*default=0.0*)
    Custom sky value to be added to blot image

Step 6: Remove cosmic rays with deriv, driz_cr
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

driz_cr: bool (*default=False*)
    Perform CR rejection with deriv and driz_cr?

driz_cr_snr: str (*default="5.0 4.0"*)
    Driz_cr.SNR parameter*
    
driz_cr_grow: int (*default=1*)
    Driz_cr_grow parameter

driz_cr_ctegrow: int (*default=0*)
    Driz_cr_ctegrow parameter

driz_cr_scale: str (*default="3.0 2.4"*)
    Driz_cr.scale parameter

driz_cr_corr: bool (*default=False*)
    Create CR cleaned _crclean file and a _crmask file?

Step 7: Drizzle final combined image
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

driz_combine: bool (*default=True*)
    Perform final drizzle image combination?

final_pixfrac: float (*default=1.0*)
    Linear size of drop in input pixels

final_fillval: int (*default=null*)
    Value to be assigned to undefined output points

final_bits: str (*default="65535"*)
    Integer mask bit values considered good

final_maskval: ??? (*default=null*)
    Value to be assigned to regions outside SCI image

final_wht_type: str (*default="EXP"*)
    Type of weighting for final drizzle

final_kernel: str (*default="square"*)
    Shape of kernel function

final_wt_scl: str (*default="exptime"*)
    Weighting factor for input data image

final_units: str (*default="cps"*)
    Units for final drizzle image (counts or cps)

Step 7a: Custom WCS for final output
""""""""""""""""""""""""""""""""""""

final_wcs: bool (*default=True*)
    Define custom WCS for final output image?

final_rot: float (*default=0.0*)
    Position Angle of drizzled image's Y-axis w.r.t. North (degrees)

final_refimage: str (*default=""*)
    Reference image from which to obtain a WCS

final_scale: int (*default=null*)
    Absolute size of output pixels in arcsec/pixel

final_outnx: int (*default=null*)
    Size of FINAL output frame X-axis (pixels)

final_outny: int (*default=null*)
    Size of FINAL output frame Y-axis (pixels)

final_ra: float (*default=null*)
    right ascension output frame center in decimal degrees

final_dec: float (*default=null*)
    declination output frame center in decimal degrees

final_crpix1: ??? (*default=null*)
    Reference pixel X position on output (CRPIX1)

final_crpix2: ??? (*default=null*)
    Reference pixel Y position on output (CRPIX2)


Alignment
^^^^^^^^^

General
"""""""

MIN_FIT_MATCHES: int (*default=4*)
    Minimum number of acceptable cross-matches for a good fit.

MAX_FIT_RMS: int (*default=20*)
    Not currently in use.

MAX_SOURCES_PER_CHIP: int (*default=250*)
    Not currently in use. 


run_align (*primarily in align.py*)
"""""""""""""""""""""""""""""""""""

update_hdr_wcs: bool (*default=True*)
    Not currently in use. Write newly computed WCS information to image headers?

catalog_list: list of strings (*default=["GAIAeDR3", "GSC242", "2MASS"]*)
    Not currently in use. Set of astrometric catalogs which should be used as references for fitting the input images. A separate fit will be performed for each catalog specified. The catalog name will also be used as part of the output ``WCSNAME`` value for the fit determined from that catalog.

fit_algorithm_list_ngt1: list of strings (*default=["match_relative_fit", "match_2dhist_fit", "match_default_fit"]*)
    Not currently in use.

fit_algorithm_list_ng1: list of strings (*default=["match_2dhist_fit", "match_default_fit"]*)
    Not currently in use.

MIN_CATALOG_THRESHOLD: int (*default=3*)
    Not currently in use.

MIN_OBSERVABLE_THRESHOLD: int (*default=10*)
    Not currently in use. The minimum number of observed sources required to continue fitting. If below this threshold, the code will return a status=1 and try with another catalog. 

MAX_FIT_LIMIT: int (*default=150*)
    Not currently in use.

mosaic_catalog_list: list of strings (*default=["GAIAeDR3", "GSC242", "2MASS"]*)
    List of available catalogs for aligning for both pipeline and SVM products. The code will go through each catalog in this order. 

mosaic_fit_list: list of strings (*default=["match_relative_fit", "match_2dhist_fit", "match_default_fit"]*)
    List of available fit algorithms for aligning for both pipeline and SVM products; match_default_fit relative alignment without using 2dhist and different throusholds (see json configuration files).

mosaic_fitgeom_list: dict (*default={"rshift": 10, "rscale": 10, "general": 6}*)
    The different fit geometries tried in alignment as well as their minobj value which specifies the number of matched sources required for a successful fit. For pipeline products, the fitgeometry value is ignored and defaults to a fit geometry of ``rscale``. The fitgeom for the pipeline products is specified as a default in *align_utils.perform_fit*. The value for minobj specified here, however, is used for the pipeline products. 

fit quality categories
""""""""""""""""""""""

    * 1 = valid solution with rms < 10 mas

    * 2 = Valid but compromised solution with rms < 10 mas

    * 3 = Valid solution with RMS >= 10 mas

    * 4 = Valid but compromised solution with RMS >= 10 mas

    * 5 = Not valid solution


generate_source_catalogs (*primarily in align_utils.py*)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""

box_size: int (*default=13*)
    The size of the box used for calculating the 2D Background of the catalog "white light" image along each axis in pixels. In nothing is specified, a default of BKG_BOX_SIZE=27 is used. 

win_size: int (*default=3*)
    The size of the 2D filter to apply to the background of the catalog "white light" image. If nothing is specified, a default of BKG_FILTER_SIZE=3 is used.

nsigma: float (*default=3.0*)
    The number of signma used for setting the level of the average background rms. If nothing is specified, a value of 5.0 is used.

centering_mode: str (*default="starfind"*)
    Algorithm to use when computing the positions of the detected sources. Options are "segmentaton" or ``starfind``. Centering will only take place after ``threshold`` has been determined, and sources are identified using segmentation. Centering using ``segmentation`` will rely on ``photutils.segmentation.SourceCatalog`` to generate the properties for the source catalog. Centering using ``starfind`` will use ``photutils.detection.IRAFStarFinder`` to characterize each source in the catalog.

bkg_estimator: str (*default="MedianBackground"*)
    Name of the algorithm to use for background calculation in *photutils.Background2D*.

rms_estimator: str (*default="StdBackgroundRMS"*)
    Name of the rms estimation algorithm used by *photutils.Background2D*.

num_sources: int (*default=250*)
    Not currently in use. Maximum number of brightest/faintest sources to return in catalog.

deblend: bool (*default=false*)
    Specify whether or not to apply photutils deblending algorithm when evaluating each of the identified segments (sources) from the chip.

fwhmpsf: float (*default=0.13*)
    The full width at half maximum of the PSF in arcseconds used for the starfind algorithm.

classify: bool (*default=false*)
    Not currently in use as cosmic rays are being removed before segmentation. Specifies whether or not to apply photutils classification algorithm when evaluating each of the identified segments (sources) from the chip.

threshold: float (*default=-1.1*)
    Value from the image which serves as the limit for determining sources. If None, compute a default value of (background+5*rms(background)). If threshold < 0.0, use absolute value as scaling factor for default value. If nothing is specified a default value of None is used to compute the background. 



generate_astrometric_catalog
""""""""""""""""""""""""""""

gaia_only: bool (*default=false*)
    Not currently in use. Specify whether or not to only use sources from GAIA in output catalog.

existing_wcs: Class Instance (*default=null*)
    Existing WCS object specified by the user


perform_fit (*primarily external in tweakwcs.matchutils.XYXYMatch*)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
For match_relative_fit, match_default_fit, and match_2dhist_fit, the following parameters are used:

fitgeom": "rscale",
    As used above, this is ignored for pipeline products. 

searchrad: float (*default=125*)
    The search radius for a match (in units of the tangent plane).

separation: float (*default=4.0*)
    The minimum separation in the tangent plane (in units of the tangent plane) for sources in the image and reference catalogs in order to be considered to be distinct sources. Objects closer together than ``separation`` distance are removed from the image and reference coordinate catalogs prior to matching.

tolerance: float (*default=2.0*)
    The matching tolerance (in units of the tangent plane) after applying an initial solution derived from the ``triangles`` algorithm.

use2dhist: bool(*default=true*)
    Use 2D histogram to find initial offset?


determine_fit_quality
"""""""""""""""""""""

MIN_CATALOG_THRESHOLD": int (*default=3*)
    The minimum number of catalog sources required to continue fitting. If below this threshold, the code will return a fit_quality=5 and try with another catalog. 

MIN_OBSERVABLE_THRESHOLD": int (*default=4*)
    If the number of observed sources is below this threshold, the code ends alignment and defers to an *a priori* solution.

MIN_CROSS_MATCHES": int (*default=3*)
    Not currently in use. 

MIN_FIT_MATCHES": int (*default=4*)
    Not currently in use.

MAX_FIT_RMS": float (*default=20*)
    Not currently in use. Maximum RMS value for a fit to be considered good. Currently a warning is printed but nothing is done with this parameter. 

MAX_FIT_LIMIT": int (*default=150*)
    The maximum allowable RMS value for a fit to be considered good. If not, the fit is considered compromised.

MAX_SOURCES_PER_CHIP: int (*default=250*)
    Maximum number of brightest sources per chip which will be used for cross-matching and fitting.

MAS_TO_ARCSEC: float (*default=1000*)
    Conversion factor from milliarcseconds to arcseconds.

GOOD_FIT_QUALITY_VALUES: int (*default=[1, 2, 3, 4]*)
    The fit_quality (see above) flag values that are allowable for a successful fit. 