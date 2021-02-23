.. _running-astrodrizzle:

********************
Running Astrodrizzle
********************

``runastrodriz`` is a module to control operation of astrodrizzle which removes distortion and combines HST images in the pipeline.


Typical Usage
=============

    >>> runastrodriz.py [-fhibn] inputFilename [newpath]


Alternative Usage
=================

    >>> python
    >>> from wfc3tools import runastrodriz
    >>> runastrodriz.process(inputFilename,force=False,newpath=None,inmemory=False)


GUI Usage under Python
======================

    >>> python
    >>> from stsci.tools import teal
    >>> import wfc3tools
    >>> cfg = teal.teal('runastrodriz')

PyRAF Usage
===========

    >>> epar runastrodriz



Options
=======

If the '-i' option gets specified, no intermediate products will be written out
to disk. These products, instead, will be kept in memory. This includes all
single drizzle products (*single_sci* and *single_wht*), median image,
blot images, and crmask images.  The use of this option will therefore require
significantly more memory than usual to process the data.

If a value has been provided for the newpath parameter, all processing will be
performed in that directory/ramdisk.  The steps involved are:

    * create a temporary directory under that directory named after the input file
    * copy all files related to the input to that new directory
    * change to that new directory and run astrodrizzle
    * change back to original directory
    * move (not copy) ALL files from temp directory to original directory
    * delete temp sub-directory

The '-b' option will run this task in BASIC mode without creating headerlets
for each input image.

The '-n' option allows the user to specify the number of cores to be used in
running AstroDrizzle.

.. note:: This value will be forced to a value of '1' (one) on Windows systems due to
  exceptions caused by threaded logging under Windows.  Future versions will
  lift this enforced restriction on Windows systems once issues with logging are
  resolved.


.. _runastrodriz-description:

Pipeline Astrometric Calibration Description
=============================================
A lot of effort goes into trying to determine the WCS solution which provides
closest agreement with the GAIA astrometry frame.  Combining the
images successfully requires the most accurate alignment possible between all
the input exposures taking into account the best available distortion model
for the input exposures.  Being able to compare the combined image to other exposures
taken of the same field at other times, or even with other telescopes, requires
that the WCS be defined to come as close to the GAIA frame as possible.  The
processing done by `runastrodriz` attempts to not only apply the most current
distortion model, but also the best available pre-computed GAIA-based WCS
solutions while proceeding to determine it's own solution to the available GAIA
reference stars for the field-of-view.  It also performs numerous checks to see
which of these solutions results in the most accurately aligned images to each
other.  The processing continued to attempt an alignment to GAIA using sources
it identifies from the observations and checks to see if any fit was successful.
Finally, it selects the WCS solution most closely aligned to GAIA that as the
basis for creating the final, distorted-corrected, combined
drizzle products for the set of exposures being processed.

.. note:: The API for the code used for these operations are described in more detail
          :ref:`advanced_products_api`.

Overview
--------
The overall logic implemented by `runastrodriz` to generate the final set of
drizzle products involves creating multiple sets of drizzle products and ultimately
selecting one as the 'best'.  The basic steps are outlined here, with following
sections providing additional detail on each step.

    #. Initialization

      a) Get any defined environment variables to control processing
      b) Interpret input file
      c) Make sure all input files are in a local directory
      d) Check for DRIZCORR calibration switch in input files
      e) Create name for output log files
      f) Define lists of individual input files to be processed

    #. Run updatewcs without astrometry database update on all input exposures (FLCs? and FLTs)

    #. Generate initial default products and perform verification

        a) perform cosmic-ray identification and generate drizzle products using astrodrizzle for all sets of inputs
        b) verify relative alignment with focus index after masking out CRs
        c) if alignment fails, update trailer file with failure information

    #. If alignment is verified,

        a) copy inputs to separate sub-directory for processing
        b) run updatewcs to get a priori updates

          * apply 'best' apriori (not aposteriori) solution

        c) generate drizzle products for all sets of inputs (FLC and/or FLT)
        d) verify alignment using focus index on FLC or, if no FLC, FLT products
        e) if alignment fails, update trailer file with info on failure
        f) if product alignment verified:

            * copy all drizzle products to parent directory
            * copy updated input exposures to parent directory

    #. If a posteriori correction enabled,

        a) copy all inputs to separate sub-directory for processing
        b) run align to align the images
        c) generate drizzle products for all sets of inputs (FLC and/or FLT) without CR identification
        d) verify alignment using focus index on FLC or, if no FLC, FLT products
        e) determine similarity index relative to pipeline default product
        f) if either focus or similarity indicates a problem, update trailer file with info on failure
        g) if product alignment verified:

           * copy all drizzle products to parent directory
           * copy updated input exposures to parent directory

    #. Remove all processing sub-directories


Initialization
--------------

Environment Variables
^^^^^^^^^^^^^^^^^^^^^^
The pipeline processing code starts out by looking to see whether the user has defined any processing behavior through the use of these environment variables:

  * **'ASTROMETRY_COMPUTE_APOSTERIORI'**: This environment variable specifies whether or not to attempt an *a posteriori* alignment where the code looks for sources in each of the images and uses those positions to perform relative alignment between the images and then fit those images to the GAIA frame.
  * **'ASTROMETRY_APPLY_APRIORI'**: This environment variable turns on/off application of any pre-defined(*a priori*) WCS solution found in the astrometry database.
  * **'ASTROMETRY_STEP_CONTROL' [DEPRECATED, do not use]**: Old variable replaced by 'ASTROMETRY_APPLY_APRIORI'.

Values that can be provided for setting these variables are:

  * 'on', 'yes', 'true': Any of these values will turn **on** the processing controlled by the variable
  * 'off', 'no', 'false': Any of these values will turn **off** the processing controleed by the variable

By default, all the processing steps are turned **on** during pipeline processing in order to maximize the chances of aligning the data as closely as possible to the absolute astrometry standard coordinate system defined through the use of the GAIA catalogs.  However, these controls are provided to support those observations which would not be suitable for such alignment, including observations of single sources.

Input Data
^^^^^^^^^^^
The processing code needs to be told what data to process, and for `runastrodriz`, a single input filename is all that **can** be provided.  This single input will be either:

  * the name of an association table for a whole set of input exposures with a filename that looks like **'<rootname>_asn.fits'**, where <rootname> is the designation for the association, such as *'ie6d07030_asn.fits'*.
  * the name of a single (uncalibrated) exposure with a filename that looks like **'<rootname>_raw.fits'**.

This one input filename, though, will simply provide the code with the information it needs to find all the calibrated input exposures which need to have their distortion-models updated and applied.  The whole set of input files required for processing includes:

  * ASN (``*_asn.fits``) files: These small FITS tables provide the relationship between the input exposures and the output products with the output filenames defined in the table.  There will NOT be an ASN table for exposures which were taken by themselves (called 'singletons').
  * RAW (``*_raw.fits``) files: Not processed directly, but required in order to get the intended value of the `DRIZCORR` calibration switch.  The ASN files also only give the rootname, and with the possibility of multiple suffixes (_flt, _flc,...) for calibrated products, the code starts with the _raw files to insure that what is specified in the ASN table is actually present and has been calibrated before processing.
  * FLT/FLC (``*_flt.fits`` or ``*_flc.fits``) files: These are the non-CTE-corrected (_flt) and CTE-corrected (_flc) calibrated exposures to be processed.

The FLT/FLC files will be the ones that actually get processed and updated with the new distortion models and WCSs, while the others allow the code to know what FLT/FLC files should be included in the processing.  This allows for multiple associations of data to live in the same directory and not interfere with each other as they are re-processed.  That can be useful when interested in combining data from multiple visits, for example.

.. warning::  Should any of these files not be available (found in the local directory), the code will raise an Exception when trying to run 'drizzlepac.astrodrizzle.AstroDrizzle' on the data.  The message will indicate what file was missing with something like: **"Exception: File ie6d07ujq_flt not found."**

Calibration Switches
^^^^^^^^^^^^^^^^^^^^
This processing serves as an official calibration step defined for HST data through the use of the **DRIZCORR** header keyword.  This keyword can be found along with all the other calibration switches in the PRIMARY header (extension 0) of the exposures FITS file. A quick way to view this (or any keyword) value would be with:

.. code:: python

    from astropy.io import fits
    val = fits.getval('ie6d07ujq_flt.fits', 'drizcorr')


This switch must be set to 'PERFORM' in order to allow the processing to be done. Processing will be completely skipped should the value of this switch in the '_raw.fits' file be set to 'OMIT'.

Log Files
^^^^^^^^^^
A number of log files, or 'trailer' files, get generated during processing, and their filenames get defined as early as possible in the processing.  The primary file will be a file with a '.tra' extension and should have the same '<rootname>' as the input file used to start the code.  For example, if you were to reprocess 'ie6d07030_asn.fits', you would end up with a trailer file with the name 'ie6d07030.tra'.

This log file contains the messages generated from performing all the updates to the distortion model, updates from the astrometry database (if any), and all the image combinations performed by 'AstroDrizzle()' to create the final set of calibrated, drizzled exposures.  Should any problems arise when during the processing, the log can provide the error messages and tracebacks to determine what went wrong.


Data to be Processed
^^^^^^^^^^^^^^^^^^^^^
Once the code has performed all the initialization, it prepares the processing by defining what files need to be combined together from the input files it can find.  This includes looking for CTE-corrected versions of the calibrated exposures (FLC files) as well as all the non-CTE-corrected files (FLT files) and creating a separate list of each type.  Many types of data do not get CTE-corrected by the instruments calibration software, such as calacs.e or calwf3.e, and so no list of FLC files will be made.  This will tell the code that it only needs to process the FLT files by themselves.  If FLC files are found, all updates to the astrometry and WCS will be performed on those files and the results then get copied into the FLT file headers upon completion of the processing.


Update the WCS
----------------
The first operation on the calibrated input files focuses on applying the calibrations
for the distortion model to the WCS.  This operation gets performed using the
`updatewcs` task using the syntax:

.. code:: python

    from stwcs.updatewcs import updatewcs
    updatewcs(calfiles_flc, use_db=False)

where `calfiles_flc` is the list of CTE-corrected FLC files or in the case there are
no CTE-corrected files, the list of calibrated FLT files.  Crucially, the use
of `use_db=False` forces `updatewcs` to only apply the distortion model to the
default WCS to create what is referred to as the **pipeline-default WCS**.  This
WCS has a `WCSNAME` associated with it that has the format ``IDC_<rootname>`` where
``<rootname>`` is the rootname of the `IDCTAB` reference files applied to the WCS.

This default WCS serves as the basis for all subsequent processing as the code
tries to determine the WCS which is aligned most closely to the GAIA astrometric
coordinate system.



Generate the initial default products
--------------------------------------
The instrument teams have calibrated the distortion models extremely well for nearly
all imaging modes with the latest calibration model being applied to the WCS keywords
when the observations were updated in the previous step.  The observations at this
point represent what the best calibration of the pointings as observed by the
telescope.  The accuracy of the guiding allows for sub-pixel alignment of the
observations for most of the data and this step applies the distortion model to
generate the 'pipeline-default' drizzle products.

The default products get generated using the ``astrodrizzle`` task.  This initial
run relies on a couple of default settings to generate the default drizzle products;
namely,

  * reads and applies default parameter settings from MDRIZTAB specified in observation header
  * uses ``resetbits=4096``
  * runs with ``crbit=4096`` to define cosmic-rays/bad-pixels with DQ flag of 4096

Identify Cosmic-Rays
^^^^^^^^^^^^^^^^^^^^
Generating these drizzle products serves as the initial attempt to identify and to flag
bad-pixels or cosmic-rays in each of the observations.  Assuming the relative
alignment of the initial pointing by the telescope is good (aligned to <0.1 pixels),
most of the cosmic-rays will be successfully identified at this point by flagging those
pixels with a value of 4096 in the DQ array for each chip.  This will
make it easier to find sources and confirm alignment without having to weed through
so many false sources.  However, there are times when the default alignment by
the telescope was not maintained which can result in all sources (real and cosmic-rays
alike) to be flagged, so subsequent steps can reset the DQ bits from 4096 to 0
while processing the data again with `astrodrizzle` using different WCS solutions.

These initial products will only be generated for the CTE-corrected versions of
the observations (``*_flc.fits`` or FLC files) if they are present, and the standard
calibrated versions of the observations (``*flt.fits`` or FLT files) otherwise.

Verifying Alignment
^^^^^^^^^^^^^^^^^^^
The relative alignment of these pipeline-default products relies entirely on the
guiding accuracy of the telescope.  Unfortunately, there are times when guiding
problems impact the observations. These guiding errors can occur due to any of
several reasons, including but not limited to:

  * re-acquisition of a different guide star from one orbit to another, usually as a result of using a close binary that was not previously identified in the guide star catalog
  * high slew rate due to only guiding on gyros due to problems with acquiring guide stars
  * spurious guiding problems due to the aging telescope and guiding systems


**Computing the Focus Index**

Verifying whether or not we can identify any problems with the relative alignment
for these products starts by measuring the focus index for the drizzled products.
The focus index was based on using the properties of the Laplacian of Gaussian (LoG)
operator as an edge detector.  See http://alumni.media.mit.edu/~maov/classes/vision09/lect/09_Image_Filtering_Edge_Detection_09.pdf for background on the Laplacian of Gaussian
operator and its use in image filtering.  The index that has been implemented is based
on the maximum value of the LoG operation on each drizzled product.

The process for computing this index is:

  * use the drizzled product, with as many cosmic-rays removed as possible, as the input
  * mask out all the saturated sources as well as possible
  * apply the LoG operator to the image
  * pick out the pixel with the maximum value to serve as the value of the focus index

This measurement process gets applied to the total drizzle product for an association,
as well as the drizzle product for each input exposure as well,
known as 'single drizzled' products.  The single drizzled products represent the
optimal focus since there is only a single exposure with only telescope focus
changes affecting the image focus value.  The range of values from the single drizzled
products establishes the distribution of 'good' focus values that gets used to
evaluate whether the total drizzle product passes focus verification.  This range
of values comes as a result of the changing focus of
the telescope from one exposure to another and to a lesser extent the effect of noise
in low-S/N observations.

A Z-score then gets computed for the focus index value of each single drizzle
product.  In simplest terms, the Z-score is a measure of how many sigma above or
below the population mean a measured valued is.  The actual
computation is:

.. code:: python

    from scipy.stats as st

    p = st.norm.cdf(x=val, loc=mean, scale=sigma)
    z_score = st.norm.ppf(p)

A Z-score then gets computed for the focus index value derived from the total
drizzle product.  If this score falls within the range of values defined by the
single drizzle focus index Z-score values, this WCS solution is considered to
have passed the 'focus verification' check.

**Computing the Similarity Index**

In addition to the focus index, a similarity index can also be computed between
the single drizzle products (again treated as 'truth') and the total drizzle
product.  The function used to compute this is the ``max_overlap_diff`` function
in ``astrometric_utils``.  The similarity index gets computed only for the region of maximum
overlap of all the input exposures.  This region of overlap gets determined
using the ``SkyFootprint`` class from the ``cell_utils`` module.  Should an input
exposure not overlap the regions where most of the exposures overlap, then the region
which overlaps at least 1 other exposure will be used for computing the index.

Point sources are detected in the selected region of overlap with a mask being
generated for each source containing a value of 1 for the point source and 0 for
the background.  The sources are identified in the single drizzle image overlap
region and the total drizzle product overlap region.  These single drizzle mask
then gets subtracted from the total drizzle mask, then scaled by the number of
non-zero pixels in the single drizzle mask resulting in a Hamming distance between
the two images.  The Hamming distance, simply put, provides the percentage of
differences pixel-by-pixel between two arrays as described in the
`scipy package spatial.distance <https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.hamming.html>`_.
This distance then gets scaled by the relative exposure time of the
single drizzle image to account for uncertainties introduced by readout noise, low
S/N detection of sources and other variances due to exposure time.

We then compute a variant of the Mean Squared Error (MSE) algorithm used in the
AmphiIndex image comparison code used for comparing images taken of amphibians.
One description of how the MSE measures the similarity between images can be found at
`https://www.pyimagesearch.com/2014/09/15/python-compare-two-images/
<https://www.pyimagesearch.com/2014/09/15/python-compare-two-images/>`_.
This similarity index is sensitive to small offsets between exposures, as well
as differences in noise, overall S/N, and even presence of cosmic-rays.
In contrast, the Hamming-distance is not as sensitive to noise.  Therefore, we
compare the MSE similarity with the Hamming
distance and take the minimum of the two values as a more robust measure of the
similarity of the images.  Both values share one key characteristic: values > 1.0
indicate more pixels are different than similar.  The code takes the maximum value
of the similarity indices computed for the total drizzle product compared to
all the single drizzle products as the final measure of the similarity.  If this
value is less than 1.0, then this WCS is considered to have passed the similarity
check.


Updating the Trailer File
^^^^^^^^^^^^^^^^^^^^^^^^^^
Associations where there are problems with the alignment will cause this verification
to fail since the sources will not be 'as sharp' based on the LoG operator.  As a
result, it can flag situations where even sub-pixel offsets down less than 0.5 pixels
are identified. For the default pipeline alignment, failure at this point is only
noted in the log with the hope that later alignment efforts will resolve the
problem affecting the original input data as noted in this check.


Applying A Priori WCS Solutions
-------------------------------
A priori WCS solutions defined for use with HST data refer to improvements to the
WCS solutions that were pre-computed.  As of 2020, there were 2 primary sources
of `a priori` WCS solutions:

    * **GSC240**:  correcting the previous guide star coordinates to the GAIA frame
    * **HSC30**: corrections derived using the Hubble Source Catalog(HSC) coordinates cross-matched to the GAIA catalog

The updated ``a priori`` solutions are stored as ``headerlets`` in an astrometry database.
The headerlet format allows them to be applied directly to the exposure using the
STWCS package while requiring very little storage space (typically, < 120Kb per
headerlet). More details on the ``headerlet`` can be found at `https://stwcs.readthedocs.io/en/latest/headerlet.html <https://stwcs.readthedocs.io/en/latest/headerlet.html>`_.

These solutions get applied through the use of the `updatewcs <https://stwcs.readthedocs.io/en/latest/updatewcs.html>`_ task in STWCS.  This task not only recomputes the PRIMARY WCS (one used by DS9 for
coordinates), but also queries the astrometry database to append all additional updated WCS solutions
as headerlet extensions based on the IDCTAB specified in the image header. The astrometry database may
also have solutions based on additional IDCTAB solutions, but those will only be applied if ``updatewcs``
gets run manually with a non-default value for the ``all_wcs`` parameter.

In the process of modifying the file, ``updatewcs`` also insures that there
are no duplicate solutions based on the ``HDRNAME`` keyword unless otherwise specified by the user.
Duplicate solutions can come from any source, even inadvertantly by the user when performing image
alignment on their own, so removing duplicates insures that the file does not get cluttered with
unnecessary extensions. This also highlights the need to insure that all new WCS solutions get
provided with unique ``HDRNAME``, and preferably ``WCSNAME`` also, keyword values.  This will insure
that the headerlet module does not thrown an Exception when trying to work with these alternate WCS
headerlet extensions.

Astrometry Database
^^^^^^^^^^^^^^^^^^^^
A publicly accessible database has been established to serve as a repository of
``a priori`` WCS solutions (full descriptions of which are found in the following
sections) as well as pipeline-generated ``a posteriori`` WCS solutions.
This database can be accessed through functions provided by the `STWCS updatewcs.astrometry_utils
module <https://stwcs.readthedocs.io/en/latest/astrometry_utils.html>`_.  This can result
in several WCS solutions being available for each exposure, with one set of solutions for
each distortion model that has been in use for these instruments since we initialized the
database in early 2019.  The functions in the ``stwcs.updatewcs.astrometry_utils`` module will
allow someone to determine the full list of WCSs available for a given exposure and have them applied
as desired to a given exposure.


Supporting New IDCTABs
^^^^^^^^^^^^^^^^^^^^^^^
Calibrations of the distortion model for each instrument evolves over time due to changes in the telescope as
well as improvements in the modeling of the distortion, including better understanding of the time-dependent
aspects of the distortion model.  These new models get provided as new versions of the
``IDCTAB`` reference file, along with the ``D2IMFILE`` and ``NPOLFILE``.  The astrometry database
contains ``a priori`` WCSs which represent the WCS for each
exposure based on the coordinates of the guide stars used for the exposure after updating
their coordinates to ones determined from the GAIA catalogs. However, they were originally computed
based on the ``IDCTAB`` reference file in use when the database was first established.

If the ``IDCTAB`` specified in
the image is not found in any of the WCSs in the database, the ``a priori`` WCS based on that ``IDCTAB`` get
determined by the ``updatewcs`` task. It starts by querying the guide-star web
interface to retrieve the corrections from the original guide star coordinates using the `stwcs.updatewcs.astrometry_utils.find_gsc_offset function <https://stwcs.readthedocs.io/en/latest/astrometry_utils.html>`_.  These offsets can also evolve as new GAIA catalogs are released to provide more accurate coordinates for the guide stars.  These offsets are then used to correct the reference point of the pipeline-default WCS based on the new ``IDCTAB`` using the ``apply_new_apriori`` method in the `AstrometryDB class <https://stwcs.readthedocs.io/en/latest/astrometry_utils.html#stwcs.updatewcs.astrometry_utils.AstrometryDB.apply_new_apriori>`_.
This method uses the same lines of code used to populate the astrometry database with the original set
of ``a priori`` WCS solutions. This not only insures that there is
always a GAIA-based WCS available for all exposures, but it does it using the best available information.
This new ``a priori`` WCS not only gets added to the image as an alternate (and perhaps PRIMARY) WCS, but it
also gets written out as a headerlet as an archive of the new WCS.


GSC240: GAIA and the HST Guide Stars
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Observations taken prior to October 2017 used guide star coordinates which were
based on guide star coordinates derived primarily from ground-based observations.
This resulted in an uncertainty of 1 arcsecond in the absolute pointing of the
telescope for any given observation.  The development and availability of the
space-based GAIA astrometric catalog finally allowed for the guide star coordinates
to be known to better than 10 milli-arcseconds in 2015 with proper motion uncertainties
increasing by 5 milli-arcseconds per year on average.  The GAIA astrometric catalog
was then cross-matched to the HST guide star catalog used for pointing the telescope,
and corrections were determined. These corrections were then applied to every HST
observation taken before Oct 2017 as if the telescope used the GAIA coordinates originally to
generate updated WCS solutions to describe the GAIA-based pointing.  These updated
WCS solutions were labelled with 'GSC240' in the WCSNAME and stored in an
astrometry database to be applied on-demand to all observations taken before Oct 2017.

These solutions will not result in perfect alignment to the GAIA catalog, due to
temporal uncertainties in the calibration of the instrument's field of view relative
to the FGS's used to point and to guide the telescope during the observations.  This
uncertainty can be up to 0.5 arcseconds, but it still represents a significant improvement
in the absolute astrometry from the 1-sigma of 1 arcsecond for previous WCS solutions.

All observations
taken after Oct 2017 already used guide-star coordinates based on GAIA, so no new
WCS was needed as it would simply be the same as the pipeline default WCS.  However, if
``updatewcs`` computes the new ``a priori`` WCS on-the-fly for a new IDCTAB for observations
taken after Oct 2017, it will be given the 'GSC240' (or newer) label in the `WCSNAME` to indicate
the type of WCS being applied to the image.

HSC30: Hubble Source Catalog WCSs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The Hubble Source Catalog(HSC) (https://archive.stsci.edu/hst/hsc/) developed a comprehensive
catalog of a majority of the sources observed in Hubble data.  This catalog was
then cross-matched to the GAIA catalog to determine improved positions for those
sources.  By using the updated positions from Version 3.0 of the HSC and comparing them to the original
positions based on the pipeline default WCS solutions, updates were derived for
all observations with sources from the HSC.  The updates were then used to recompute
the WCS solutions for those observations which were labelled as 'HSC30' in the WCSNAME and
stored in the astrometry database.

Separate Directories
^^^^^^^^^^^^^^^^^^^^
One mechanism used to enable comparisons of various WCS solutions is to keep
copies of the observations with different types of WCS solutions in separate
directories.  Up until this point in the processing, the data has been processed
in the directory where the processing was started.  In order to keep the ``a priori``
solutions separate, a sub-directory gets created with name based on the association
table rootname or the rootname of the single exposure being processed using the
convention:  `<rootname>_apriori`.  All the FLC (or FLT, if no FLC files are present),
and ASN file (if processing an association) are copied from the main directory into
the new sub-directory and the process moves to the sub-directory to continue its
processing.

Applying A Priori Solutions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Application of ``a priori`` WCS solutions computed in previous STScI automated
calibration (pipeline) processing also occurs when running the ``updatewcs``
task with ``use_db=True`` (the default setting).  This queries the astrometry
database and retrieves the headerlets for all the ``a priori`` solutions.  Only
those WCSs based on the currently specified IDCTAB will be retained unless the user
requests that all solutions be kept.

The database reports what solution is flagged as the ``best`` solution, which will
typically result in the closest alignment to GAIA and will be the previously
computed ``a posteriori`` solution if available.  All the retrieved headerlets get appended
as new extensions to the observations FITS file, then the database WCS solution flagged as ``best``
gets applied to replace the active or primary WCS in the observation after saving
a copy of the original primary WCS.  However, this solution only gets used
to replace the current PRIMARY WCS in the SCI header if it was based on the
same IDCTAB as currently specified in the image primary header.
The other solutions returned by the database are retained, but not applied at this point
to enable the user to switch between them later as appropriate for their work.

When performing the standard processing with ``runastrodriz``, we are only interested in seeing
whether there are any issues in applying the pre-defined ``a priori`` corrections, while
also setting the standard for the relative alignment for comparison with any new ``a posteriori``
fit that may be determined later in the processing.

Generating A Priori Products
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The FLC images updated with the ``a priori`` WCS solutions now get combined using
``astrodrizzle``.  If the pipeline default focus verification succeeded, then
``resetbits`` will be set to 0 so that the previous DQ flags can be used.  If the
verification failed, though, ``resetbits`` gets set to 4096 so that the cosmic-rays
can be identified and flagged fresh based on the alignment provided by the ``a priori``
WCS solutions.

This processing will result in a total combined drizzle product based on the
``a priori`` solution.

Evaluating Alignment
^^^^^^^^^^^^^^^^^^^^^
Confirming that the relative alignment between the images in the association was
maintained with the ``a priori`` WCS now can be done.  Although the ``a priori``
WCS solutions are vetted for accuracy, HST has taken a few hundred thousand
different exposures in dozens of configurations and not all of those exposures were
taken exactly as planned.  Therefore, considerable effort goes into trying to verify
that the alignment between the images has been maintained.

This verification starts by computing the focus index and similarity values for the total
drizzle product and the single drizzle products using the same code used to verify
the pipeline default WCS drizzle product.  It then extends to include computing
the similarity index between the ``a priori`` drizzle products and the pipeline
default drizzle products.  This will attempt to measure whether or not the ``a priori``
alignment is significantly different than the presumably good pipeline default
alignment.  Once again, if the similarity index is less than 1, the ``a priori``
alignment is considered to be successful.

Keeping the A Priori Alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Should all the verification steps indicate a successful alignment, the ``a priori``
WCS solution should be retained as an improved WCS solution over the pipeline
default WCS.  This gets done by simply copying the calibrated images which have been
updated with the WCS solution (both the FLC and FLT images) from the ``<rootname>_apriori``
sub-directory to the main processing directory.  This will replace the FLC and FLT
files with the pipeline default solutions so that should no other WCS prove to be
better, the ``a priori`` WCS solution will end up being used to generate the final
drizzle products which get archived and provided to the end-user.


Performing An A Posteriori Alignment
-------------------------------------
The ultimate goal of this processing would be to have the input observations
aligned as closely to an astrometric standard coordinate system as much as
possible.  The highest quality, highest precision astrometric catalog available
would be the GAIA astrometric catalog and this processing seeks to align HST
observations as closely to that catalog's coordinate system.

The ``a priori`` solutions provide an update to the astrometry based
on either the guide stars used (the ``GSC240`` and related solutions) or manually
verified alignment of sources from the observations field-of-view performed using
the Hubble Source Catalog (the ``HSC30`` solution).  Unfortunately, both of
these types of solutions fail to account for sources of astrometric error which
can still affect the observations and result in offsets from the GAIA system due to
updates in the distortion calibration for the instruments or uncertainties in
the position of the detectors field-of-view relative to the Fine Guidance Sensors
(FGS) and the guide stars used for taking the observaitons.

The only way to correct for those effects remains to identify sources from the
observations and perform a fit to the GAIA catalog directly.  This is called an
``a posteriori`` solution when it can be done successfully.  However, this can
only be performed for observations which contain enough detectable sources,
specifically sources found in the GAIA catalog.  Not all observations meet this
criteria either due to exposure time (too long or too short), wavelength of
observation, filter bandpass (narrowband vs wide-band) and even number of sources
in the field.  This processing code makes no assumptions about the possibility of
success and tries to perform this ``a posteriori`` fit on all observations.

Copying the Observations
^^^^^^^^^^^^^^^^^^^^^^^^^
Copies of the observaions are made in a sub-directory named after the input
file used to start the processing with the convention:

   <rootname>_aposteriori

For example, if the association **icw402010_asn.fits** was being processed, this
directory would be named **icw402010_aposteriori**.

All the calibrated FLC and/or FLT images along with the ASN file are copied into
this sub-directory.  These files, at this point, have the best available WCS at
this time which is most likely an ``a priori`` solution.  This improves the
chance that the ``a posteriori`` fit will work by minimizing the offset from GAIA
which needs to be searched to find a cross-match with the GAIA sources in the
field-of-view.

Aligning the Observations
^^^^^^^^^^^^^^^^^^^^^^^^^^
The alignment process gets performed using the ``perform_align()`` function from
the ``drizzlepac/align`` module. This function performs the following steps in
an attempt to perform an ``a posteriori`` fit to GAIA:

    * Evaluates all the input observations to identify any which can not be
      aligned, such as GRISM or SCAN mode observations.  For a full description
      of all the type of observations that can be filtered out, see
      :ref:`analyze_api`.
    * Compute a 2D background for all the observations using ``photutils``
    * Determine a PSF kernel from the detectable sources in the image, if possible.
    * Segments the image after applying the 2D background to identify as many
      sources as possible above a threshold using ``photutils.segmentation``
    * Performs source centering using ``photutils.DAOStarFinder``
    * Keeps the position of the single brightest source nearest the center of
      the segment as the catalog position for each segment's object.
    * Checks whether there are enough sources to potentially get a viable linear
      fit.

        * If not, the attempt at an ``a posteriori`` fit quits without updating
          the WCS of the input files.

    * Queries the GAIA DR2 catalog through the STScI web service to obtain a catalog
      of GAIA sources that overlap the field-of-view of the combined set of
      observations. This catalog will serve as the **reference catalog** for the
      fitting process.

        * If there are not enough GAIA sources overlapping these observations,
          then the fit attempt quits without updating the WCS of the input
          files.

    * Provide the source catalogs for each input image, each input images's WCS,
      and the GAIA reference catalog to function ``align_wcs()`` in the ``tweakwcs``
      package.

        * This function cross-matches the source catalog from each image with
          the GAIA catalog and performs an **rscale** linear fit (as defined by
          ``runastrodriz``), then updates the input WCS with the results of the
          fit upon success.  See the `tweakwcs readthedocs pages
          <https://tweakwcs.readthedocs.io/en/latest/imalign.html>`_ for more
          details.
        * The function ``align_wcs`` is first called without using the GAIA
          reference catalog in order to perform a relative alignment between the observations.
        * The function ``align_wcs`` is then called with the GAIA catalog as
          the reference in order to finally perform a single fit to the GAIA catalog
          for all the observations at the same time.

    * Evaluate the success/failure state of the fit and the quality of any
      successful fit.
    * Repeat the fit with ``tweakwcs.align_wcs`` with other GAIA catalogs;
      including GAIA DR1 or any others specified for use in ``runastrodriz`` itself.
    * Select the fit to the GAIA catalog which results in the lowest RMS.

        * Some fields are dominated by external galaxies with no proper motion for
          which GAIA DR1 without proper motions provides the best fit (lowest RMS).
        * Other fields are dominated by local galactic stars with appreciable
          proper motions best accounted for (still with some error) by the
          GAIA DR2 catalog with its proper motions.

    * Keep the WCS's updated with the **best** solution and update the **WCSNAME**
      keyword for those WCSs to reflect the type of fit that was successful and
      the catalog that was used.

        * The naming convention is more fully described on the
          `Drizzlepac Astrometry description
          <https://drizzlepac.readthedocs.io/en/latest/astrometry.html>`_.

The result of this lengthy process is a set of WCS objects which have been
updated with a fit to a GAIA catalog representing an ``a posteriori`` solution.


Generate the Aligned Drizzle Products
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Successful alignment of the WCSs to a GAIA catalog means that these ``a posteriori``
updated exposures can be combined to create a drizzled product using ``AstroDrizzle``.


Verify the A Posteriori Alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
These newly updated drizzle products still need to be evaluated to insure that the
fit performed to GAIA maintained relative alignment between the images as well.
Mis-alignment of the images to each other can result from too few sources being
used for the fit imprinting the errors in those source positions on the relative
alignment.  The verification used is the same focus and similarity checks that were
performed on the ``a priori`` updated drizzle products and even the pipeline
default drizzle products.

A Posteriori Failure
^^^^^^^^^^^^^^^^^^^^
At any number of points throughout this computation and verification, it could
end up quitting and flagging this attempt as a failure. If this happens, no
updated WCS solutions get created or saved and processing returns to the parent
directory while deleting the entire ``<dataset>_aposteriori`` directory along
with all the mis-aligned or un-alignable files.  This allows the processing to
revert to using the previously verified WCS solutions as the ``best`` WCS solution
available for these observations.

A Posteriori Success
^^^^^^^^^^^^^^^^^^^^
Successfully fitting to GAIA can only be declared after the verification process
returned values indicating good alignment in the drizzle product.  The processing
would then copy these ``a posteriori``-updated input exposures from the sub-directory
these computations were being performed in based up to the parent directory to
replace the previously updated versions of the input files.  This entire sub-directory
then gets deleted, unless the processing was being run in debug mode.

Creation of Final Aligned Products
----------------------------------
The starting directory now contains updated input FLC/FLT files based on WCSs which
have been verified to have maintained relative alignment and with alignment as close
to the GAIA astrometric coordinate system as possible.  These exposures get
processed by ``AstroDrizzle`` to create the final, combined drizzle products for
the user and for archiving at STScI in the Mikulski Archive for Space Telescopes (MAST).
These products include the calibrated drizzle(DRZ) products as well as any
CTE-corrected drizzle(DRC) products depending on what input exposures are
available.
