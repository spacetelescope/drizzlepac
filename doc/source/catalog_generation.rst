.. _catalog_generation:

==================
Catalog Generation
==================

Point (Aperture) Photometric Catalog Generation
================================================

1: Source Detection
-------------------

1.1: Important Clarifications
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
As previously discussed in :ref:`singlevisit`, AstroDrizzle creates a single multi-filter, detector-level 
drizzle-combined image for source identification and one or more detector/filter-level drizzle-combined images 
(depending on
which filters were used in the dataset) for photometry. The same set of sources identified in the
multi-filter detection image is used to measure photometry for each filter. We use this method to maximize the
signal across all available wavelengths at the source detection stage, thus providing photometry with the
best quality source list across all available input filters.

It should also be stressed here that the point and segment photometry source list generation algorithms
identify source catalogs independently of each other and DO NOT use a shared common source catalog for
photometry.

1.2: Preliminaries
^^^^^^^^^^^^^^^^^^^^

1.2.1: Generation of Pixel Masks
""""""""""""""""""""""""""""""""""
Every multi-filter, detector-level drizzle-combined image is associated with a boolean footprint mask which 
defines the illuminated (True) and non-illuminated (False) portions of the image based upon its constituent 
exposures and the corresponding WCS solution.  The boundary of the illuminated portion
is iteratively eroded or contracted to minimize the impact of regions where signal
quality is known to be degraded, and thereby, could affect source identification and subsequent 
photometric measurements.  The erosion depth is approximately ten pixels and is done by using the 
`ndimage.binary_erosion <https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.binary_erosion.html>`_ scipy tool.
For computational convenience, an inverse footprint mask is also available for functions
which utilize masks to indicate pixels which should be ignored during processing of the 
input data.


1.2.2: Detection Image Background Determination
""""""""""""""""""""""""""""""""""""""""""""""""
For consistency, the same background and background RMS images are used by both the point and
segment algorithms.
To ensure optimal source detection, the multi-filter detection image must be background-subtracted. 
In order to accommodate the different types of detectors, disparate signal levels, and highly varying 
astronomical image content, three background computations are used as applicable.  The first category 
of background definition is a special case situation, and if it is found to be applicable to the detection 
image, the background and background RMS images are defined and no further background evaluation is done. 
It has been observed that some background regions of ACS SBC drizzle-combined 
detection images, *though the evaluation is done for all instrument detection images*,
are measured to have values identically
equal to zero.  If the number of identically zero pixels in the footprint portion of the detection image
exceeds a configurable percentage threshold value (default is 25%), then a two-dimensional background image 
is constructed and set to the constant of zero. Its companion constructed RMS image set to the RMS
value computed for the non-zero pixels which reside within the footprint portion of the image.

If the special background determination category above is not applicable, then sigma-clipped statistics are
computed for the detection image using the `astropy.stats.sigma_clipped_stats <https://docs.astropy.org/en/stable/api/astropy.stats.sigma_clipped_stats.html>`_ 
Astropy tool. This algorithm uses the detection image and its inverse footprint mask, as well
as the specification for the number of standard deviations and the maximum number of iterations
to compute the mean, median, and rms of the
sigma-clipped data.  The specification for the number of standard deviations and the maximum number
of iterations are configurable values which are set to 3.0 and 3 by default, respectively.

At this point the Pearson's second coefficient of skewness is computed.

.. math::
    skewness = 3.0 * (mean - median) / rms 

The skewness compares a sample distribution with a normal distribution where the
larger the absolute value of the skewess, the more the sample distribution differs from
a normal distribution. The skewness is computed in this context to aid in determining 
whether it is worth computing the background by other means.

In order to ensure the sigma-clipped statistics are reasonable, a negative mean or median value is
reset to zero, and a minimum RMS value is computed for comparison based upon FITS keyword values in 
the detection image header.  For the CCD detectors only, a-to-d gain (ATODGN), read noise
(READNSE), number of drizzled images (NDRIZIM), and total exposure time (TEXPTIME) are employed
to compute a minimum RMS.  Once viable background and background RMS values are determined, 
two-dimensional images matching the dimensions of the detection image are constructed.
Through configuration settings, a user can specify the sigma-clipped statistics algorithm be
used to compute the background and RMS images, though the special case of identically zero 
background data will be evaluated and will supersede the user request when applicable.

The final background determination algorithm which is the
`photutils.background.Background2d <https://photutils.readthedocs.io/en/stable/api/photutils.background.Background2D.html>`_ Astropy tool
is only invoked if the special case identically zero algorithm has not been applied,
the user has not requested that only the sigma-clipped statistics algorithm be computed, and if the 
skewness value derived using the sigma-clipped statistics is less than a pre-defined and configurable
threshold (default value of 0.5).

The `photutils.background.Background2d <https://photutils.readthedocs.io/en/stable/api/photutils.background.Background2D.html>`_ 
algorithm uses
sigma-clipped statistics to determine background and RMS values across the image. An initial low-resolution
estimate of the background is performed by computing sigma-clipped median values in 27x27 pixel boxes across
the image. This low-resolution background image is then median-filtered using a 3x3 pixel sample window to
correct for local small-scale overestimates and/or underestimates. 

Once a background and RMS image are determined using this final technique, a preliminary 
background-subtracted image is computed so it can be evaluated for the percentage of negative
values in the illuminated portion of the image. If the percentage of negative values exceeds a
configurable and defined threshold, the computation of the background and RMS image from this
algorithm are discarded.  The background and RMS images computed using the sigma-clipped statistics in
technique two, with its associated updates, are ultimately chosen as the images to use.

It should be noted these are configurable values.
Our catalogs use these values deeming them to be the best for the general situation, but users can tune these values to
optimize for their own data. To this end, users can adjust parameter values "bkg_box_size" and/or
"bkg_filter_size" in the <instrument>_<detector>_catalog_generation_all.json files in the following path:
/drizzlepac/pars/hap_pars/default_parameters/<instrument>/<detector>/.

1.3: Source Identification with DAOStarFinder
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
We use the `photutils.detection.DAOStarFinder <https://photutils.readthedocs.io/en/stable/api/photutils.detection.DAOStarFinder.html>`_ Astropy tool to identify sources in the background-subtracted
multi-filter detection image. Regions flagged in the previously created bad pixel mask are ignored by
DAOStarFinder. This algorithm works by identifying local brightness maxima with roughly gaussian
distributions whose peak values are above a predefined minimum threshold. Full details of the process are
described in `Stetson 1987; PASP 99, 191 <http://adsabs.harvard.edu/abs/1987PASP...99..191S>`_.
The exact set of input parameters fed into DAOStarFinder is detector-dependent. The parameters can be found in
the <instrument>_<detector>_catalog_generation_all.json files mentioned in the previous section.

2: Aperture Photometry Measurement
------------------------------------

2.1: Flux determination
^^^^^^^^^^^^^^^^^^^^^^^^
Aperture photometry is then preformed on the previously identified sources using a pair of concentric
photometric apertures. The sizes of these apertures depend on the specific detector being used, and are
listed below in table 1:

.. table:: Table 1: Aperture photometry aperture sizes

    +-------------------+----------------------------+----------------------------+
    |Instrument/Detector|Inner aperture size (arcsec)|Outer aperture size (arcsec)|
    +===================+============================+============================+
    |ACS/HRC            |0.03                        |0.125                       |
    +-------------------+----------------------------+----------------------------+
    |ACS/SBC            |0.07                        |0.125                       |
    +-------------------+----------------------------+----------------------------+
    |ACS/WFC	        |0.05                        |0.15                        |
    +-------------------+----------------------------+----------------------------+
    |WFC3/IR	        |0.15                        |0.45                        |
    +-------------------+----------------------------+----------------------------+
    |WFC3/UVIS          |0.05                        |0.15                        |
    +-------------------+----------------------------+----------------------------+

Raw (non-background-subtracted) flux values are computed by summing up the enclosed flux within the two specified
apertures using the `photutils.aperture.aperture_photometry <https://photutils.readthedocs.io/en/stable/api/photutils.aperture.aperture_photometry.html>`_
tool. Input values are detector-dependent, and can be found in the \*_catalog_generation_all.json files described above
in section 1.3.

Local background values are computed based on the 3-sigma-clipped mode of pixel values present in a circular annulus
with an inner radius of 0.25 arcseconds and an outer radius of 0.50 arcseconds surrounding each identified source. This
local background value is then subtracted from the raw inner and outer aperture flux values to compute the
background-subtracted inner and outer aperture flux values found in the output .ecsv catalog file by the formula

.. math::
    f_{bgs} = f_{raw} - f_{bg} \cdot a

where
    * :math:`f_{bgs}` is the background-subtracted flux, in electrons per second
    * :math:`f_{raw}` is the raw, non-background-subtracted flux, in electrons per second
    * :math:`f_{bg}` is the per-pixel background flux, in electrons per second per pixel
    * :math:`a` is the area of the photometric aperture, in pixels

The overall standard deviation and mode values of pixels in the background annulus are also reported for each
identified source in the output .ecsv catalog file in the “STDEV” and “MSKY” columns respectively (see Section 3 for
more details).

2.2: Calculation of photometric errors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
2.2.1: Calculation of flux uncertainties
"""""""""""""""""""""""""""""""""""""""""
For every identified source, the `photutils.aperture_photometry() <https://photutils.readthedocs.io/en/stable/api/photutils.aperture.aperture_photometry.html>`_
tool calculates standard deviation values for each aperture based on a 2-dimensional RMS array computed using the
`photutils.background.Background2d()  <https://photutils.readthedocs.io/en/stable/api/photutils.background.Background2D.html>`_
tool that we previously utilized to compute the 2-dimensional background array in order to background-subtract the
detection image for source identification. We then compute the final flux errors as seen in the output .ecsv catalog
file using the following formula:

.. math::
    \Delta f = \sqrt{\frac{\sigma^2 }{g}+(a\cdot\sigma_{bg}^{2})\cdot (1+\frac{a}{n_{sky}})}

where
    * :math:`{\Delta} f`  is the flux uncertainty, in electrons per second
    * :math:`{\sigma}` is the standard deviation of photometric aperture signal, in counts per second
    * :math:`{g}` is effective gain in electrons per count
    * :math:`{a}` is the photometric aperture area, in pixels
    * :math:`{\sigma_{bg}}` is standard deviation of the background
    * :math:`{n_{sky}}` is the sky annulus area, in pixels

2.2.2: Calculation of ABmag uncertainties
"""""""""""""""""""""""""""""""""""""""""""
Magnitude error calculation comes from computing :math:`{\frac{d(ABMAG)}{d(flux)}}`. We use the following formula:

.. math::
    \Delta mag_{AB} = 1.0857 \cdot  \frac{\Delta f}{f}

where
    * :math:`{\Delta mag_{AB}}` is the uncertainty in AB magnitude
    * :math:`{\Delta f}` is the flux uncertainty, in electrons per second
    * :math:`{f}` is the flux, in electrons per second

2.3: Calculation of concentration index (CI) values and flag values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
2.3.1: Calculation of concentration index (CI) values
""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The Concentration index is a measure of the "sharpness" of a given source’s PSF, and computed with the following
formula:

.. math::
    CI = m_{inner} - m_{outer}

where
    * :math:`{CI}` is the concentration index, in AB magnitude
    * :math:`{m_{inner}}` is the inner aperture AB magnitude
    * :math:`{m_{outer}}` is the outer aperture AB magnitude

We use the concentration index to automatically classify each identified photometric source as either a point source
(i.e. stars), an extended source (i.e. galaxies, nebulosity, etc.), or as an “anomalous” source (i.e. saturation,
hot pixels, cosmic ray hits, etc.). This designation is described by the value in the "flags" column

2.3.2: Determination of flag values
"""""""""""""""""""""""""""""""""""""
The flag value associated with each source provides users with a means to distinguish between legitimate point sources,
legitimate extended sources, and scientifically dubious sources (those likely impacted by low signal to noise, detector
artifacts, saturation, cosmic rays, etc.). The values in the “flags” column of the catalog are a sum of a one or more of
these values. Specific flag values are defined below in table 2:

.. table:: Table 2: Flag definitions

    +------------+-----------------------------------------------------------+
    | Flag value | Meaning                                                   |
    +============+===========================================================+
    | 0          | Point source :math:`{(CI_{lower} < CI < CI_{upper})}`     |
    +------------+-----------------------------------------------------------+
    | 1          | Extended source :math:`{(CI > CI_{upper})}`               |
    +------------+-----------------------------------------------------------+
    | 2          | Bit value 2 not used in ACS or WFC3 sourcelists           |
    +------------+-----------------------------------------------------------+
    | 4          | Saturated Source                                          |
    +------------+-----------------------------------------------------------+
    | 8          | Faint Detection Limit                                     |
    +------------+-----------------------------------------------------------+
    | 16         | Hot pixels :math:`{(CI < CI_{lower})}`                    |
    +------------+-----------------------------------------------------------+
    | 32         | False Detection: Swarm Around Saturated Source            |
    +------------+-----------------------------------------------------------+
    | 64         | False detection due proximity of source to image edge     |
    |            | or other region with a low number of input images         |
    +------------+-----------------------------------------------------------+

2.3.2.1: Assignment of flag values 0 (point source), 1 (extended source), and 16 (hot pixels)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Assignment of flag values 0 (point source), 1 (extended source), and 16 (hot pixels) are determined purely based on the
concentration index (CI) value. The majority of commonly used filters for all ACS and WFC3 detectors have
filter-specific CI threshold values that are automatically set at run-time. However, if filter-specific CI threshold
values cannot be found, default instrument/detector-specific CI limits are used instead.  Instrument/detector/filter
combinations that do not have filter-specific CI threshold values are listed below in table 3 and  the default CI
values are listed below in table 4.

.. table:: Table 3: Instrument/detector/filter combinations that **do not** have filter-specific CI threshold values

    +------------------------+---------------------------------------------------+
    | Instrument/Detector    | Filters without specifically defined CI limits    |
    +========================+===================================================+
    | ACS/HRC                | F344N                                             |
    +------------------------+---------------------------------------------------+
    | ACS/SBC                | All ACS/SBC filters                               |
    +------------------------+---------------------------------------------------+
    | ACS/WFC                | F892N                                             |
    +------------------------+---------------------------------------------------+
    | WFC3/IR                | None                                              |
    +------------------------+---------------------------------------------------+
    | WFC3/UVIS              | None                                              |
    +------------------------+---------------------------------------------------+

.. note:: As photometry is not performed on observations that utilized grisms, prisms, polarizers, ramp filters, or quad filters, these elements were omitted from the above list.

.. table:: Table 4: Default concentration index threshold values

    +---------------------+----------------------+----------------------+
    | Instrument/Detector | :math:`{CI_{lower}}` | :math:`{CI_{upper}}` |
    +=====================+======================+======================+
    | ACS/HRC             | 0.9                  | 1.6                  |
    +---------------------+----------------------+----------------------+
    | ACS/SBC             | 0.15                 | 0.45                 |
    +---------------------+----------------------+----------------------+
    | ACS/WFC             | 0.9                  | 1.23                 |
    +---------------------+----------------------+----------------------+
    | WFC3/IR             | 0.25                 | 0.55                 |
    +---------------------+----------------------+----------------------+
    | WFC3/UVIS           | 0.75                 | 1.0                  |
    +---------------------+----------------------+----------------------+

2.3.2.2: Assignment of flag value 4 (Saturated Source)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""
A flag value of 4 is assigned to sources that are saturated. The process of identifying saturated sources starts by
first transforming the input image XY coordinates of all pixels flagged as saturated in the data quality arrays of each
input flc/flt.fits images (the images drizzled together to produce the drizzle-combined filter image being used to
measure photometry) from non-rectified, non-distortion-corrected coordinates to the rectified, distortion-corrected
frame of reference of the filter-combined image. We then identify impacted sources by cross-matching this list of
saturated pixel coordinates against the positions of sources in the newly created source catalog and assign flag values
where necessary.

2.3.2.3: Assignment of flag value 8 (faint detection limit)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
A flag value of 8 is assigned to sources whose signal to noise ratio is below a predefined value. We define sources as
being above the faint object limit if the following is true:

.. math::
    \Delta ABmag_{outer} \leq  \frac{2.5}{snr \cdot log(10))}

Where
    * :math:`{\Delta ABmag_{outer}}` is the outer aperture AB magnitude uncertainty
    * :math:`{snr}` is the signal to noise ratio, which is 1.5 for ACS/WFC and 5.0 for all other detectors.

2.3.2.4: Assignment of flag value 32 (false detection: swarm around saturated source)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The source identification routine has been shown to identify false sources in regions near bright or saturated
sources, and in image artifacts associated with bright or saturated sources, such as diffraction spikes, and in the
pixels surrounding saturated PSF where the brightness level “plateaus” at saturation. We identify impacted sources by
locating all sources within a predefined radius of a given source and checking if the brightness of each of these
surrounding sources is less than a radially-dependent minimum brightness value defined by a pre-defined stepped
encircled energy curve. The parameters used to determine assignment of this flag are instrument-dependent, can be found
in the “swarm filter” section of the \*_quality_control_all.json files in the path described above in section 1.3.


2.3.2.5: Assignment of flag value 64 (False detection due proximity of source to image edge or other region with a low number of input images)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Sources flagged with a value of 64 are flagged as “bad” because they are inside of or in close proximity to regions
characterized by low or null input image contribution. These are areas where for some reason or another, very few or no
input images contributed to the pixel value(s) in the drizzle-combined image.
We identify sources impacted with this effect by creating a two-dimensional weight image that maps the number of
contributing exposures for every pixel. We then check each source against this map to ensure that all sources and flag
appropriately.

3: The output catalog file
---------------------------
3.1: Filename format
^^^^^^^^^^^^^^^^^^^^^^
Source positions and photometric information are written to a .ecsv (Enhanced Character Separated Values) file. The
naming of this file is fully automatic and follows the following format:
<TELESCOPE>_<PROPOSAL ID>_<OBSERVATION SET ID>_<INSTRUMENT>_<DETECTOR>_
<FILTER>_<DATASET NAME>_<CATALOG TYPE>.ecsv

So, for example if we have the following information:
    * Telescope = HST
    * Proposal ID = 98765
    * Observation set ID = 43
    * Instrument = acs
    * Detector = wfc
    * Filter name = f606w
    * Dataset name = j65c43
    * Catalog type = point_cat

The resulting auto-generated catalog filename will be:
    * hst_98765_43_acs_wfc_f606w_j65c43_point-cat.ecsv

3.2: File format
^^^^^^^^^^^^^^^^^
The .ecsv file format is quite flexible and allows for the storage of not only character-separated datasets, but also
metadata. The first section (lines 4-17) contains a mapping that defines the datatype, units, and formatting
information for each data table column. The second section (lines 19-27) contains information explaining STScI’s use
policy for HAP data in refereed publications. The third section (lines 28-48) contains relevant image metadata. This
includes the following items:

    * WCS (world coordinate system) name
    * WCS (world coordinate system) type
    * Proposal ID
    * Image filename
    * Target name
    * Observation date
    * Observation time
    * Instrument
    * Detector
    * Target right ascension
    * Target declination
    * Orientation
    * Aperture right ascension
    * Aperture declination
    * Aperture position angle
    * Exposure start (MJD)
    * Total exposure duration in seconds
    * CCD Gain
    * Filter name
    * Total Number of sources in catalog

The next section (lines 50-66) contains important notes regarding the coordinate systems used, magnitude system used,
apertures used, concentration index definition and flag value definitions:

    * X, Y coordinates listed below use are zero-indexed (origin = 0,0)
    * RA and Dec values in this table are in sky coordinates (i.e. coordinates at the epoch of observation and fit to GAIADR1 (2015.0) or GAIADR2 (2015.5)).
    * Magnitude values in this table are in the ABMAG system.
    * Inner aperture radius in pixels and arcseconds (based on detector platescale)
    * Outer aperture radius in pixels and arcseconds (based on detector platescale)
    * Concentration index (CI) formulaic definition
    * Flag value definitions

Finally, the last section contains the catalog of source locations and photometry values. It should be noted that the
specific columns and their ordering were deliberately chosen to facilitate a 1:1 exact mapping to the_daophot.txt
catalogs produced by Hubble Legacy Archive. As this code was designed to be the HLA's replacement, we sought to
minimize any issues caused by the transition. The column names are as follows (Note that this is the same left-to-right
ordering in the .ecsv file as well):

    * X-Center: 0-indexed X-coordinate position
    * Y-Center: 0-indexed Y-coordinate position
    * RA: Right ascension (sky coordinates), in degrees
    * DEC: Declination (sky coordinates), in degrees
    * ID: Object catalog index number
    * MagAp1: Inner aperture brightness, in AB magnitude
    * MagErrAp1: Inner aperture brightness uncertainty, in AB magnitude
    * MagAp2: Outer aperture brightness, in AB magnitude
    * MagErrAp2: Outer aperture brightness uncertainty, in AB magnitude
    * MSkyAp2: Outer aperture background brightness, in AB magnitude
    * StdevAp2: Standard deviation of the outer aperture background brightness, in AB magnitude
    * FluxAp2: Outer aperture flux, in electrons/sec
    * CI: Concentration index (MagAp1 – MagAp2), in AB magnitude
    * Flags: See Section 2.3.2 for flag value definitions

3.3 Rejection of Cosmic-Ray Dominated Catalogs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Not all sets of observations contain multiple overlapping exposures in the same filter. This makes it impossible
to ignore all cosmic-rays that have impacted those single exposures.  The contributions of cosmic-rays often
overwhelm any catalog generated from those single exposures making recognizing astronomical sources almost
impossible amongst the noise of all the cosmic-rays.  As a result, those catalogs can not be trusted.  In an
effort to only publish catalogs which provide the highest science value, criteria developed by the Hubble Legacy
Archive (HLA) has been implemented to recognize those catalogs dominated by cosmic-rays and not provided as an
output product.

.. note ::
  This rejection criteria is NOT applied to WFC3/IR or ACS/SBC data since they are not affected by cosmic-rays
  in the same way as the other detectors.

3.3.1 Single-image CR Rejection Algorithm
"""""""""""""""""""""""""""""""""""""""""""
An algorithm has been implemented to identify and ignore cosmic-rays in single exposures.  This algorithm has
been used for ignoring cosmic-rays during the image alignment code used to determine the *a posteriori*
alignment to GAIA.

This algorithm starts by evaluating the central moments of all sources from the segment catalog.
Any source where the maximum central moment (as determined by
`photutils.segmentation.SourceProperties <https://photutils.readthedocs.io/en/stable/api/photutils.segmentation.SourceProperties.html#photutils.segmentation.SourceProperties>`_
is 0 for both X and Y moments gets identified as cosmic-rays.  This indicates that the source has a
concentration of flux greater than a point-source and most probably represents a 'head-on cosmic-ray'.

In addition to these 'head-on cosmic-rays', 'glancing cosmic-rays' produce streaks across the detector.
Those are identified by identifying sources with a minimum width (semiminor_axis) less than the FWHM of a point source
and an elongation > 2.  The width and elongation are also properties defined by
`photutils.segmentation.SourceProperties <https://photutils.readthedocs.io/en/stable/api/photutils.segmentation.SourceProperties.html#photutils.segmentation.SourceProperties>`_.
The combination of these criteria allows for the identification of a vast majority of cosmic-rays.  The DQ array
of the single exposure then gets updated to flag those pixels identified as cosmic-rays based on these criteria.
These DQ flags are then ONLY applied when creating the TotalProduct to limit the contribution of cosmic-rays
from the total detection image.  These flags are NOT used to generate any other product in order to avoid
affecting the photometry or astrometry of any source from the total detection image any more than necessary.

3.3.2 Rejection Criteria
"""""""""""""""""""""""""
The rejection criteria has been defined so that if either the point source catalog or the segment catalog fails,
then both catalogs are rejected and deleted.

In its simplest form the criteria for rejection is:
        n_cat < thresh
where:
        thresh = crfactor * (n1_residual * n1_exposure_time)**2 / texptime
and:
        n_cat    : Number of good point and extended sources in the catalog (flag < 2)
        crfactor : Number of expected cosmic-rays per second across the entire detector
        n1_exposure_time : amount of exposure time for all single filter exposures
        texptime : Total exposure time of the combined drizzle product
        n1_residual : Remaining fraction of cosmic-rays after applying single-image CR removal

The value of `crfactor` should be adjusted for sub-arrays to account for the smaller area being read out, but
that logic has not yet been implemented.  The values used in the processing of single-visit mosaics are:

    segment-catalog crfactor : 300
    point-catalog crfactor   : 150

These numbers are deliberately set high to be conservative about which catalogs to keep.  The CR rate varies
with position in the orbit, and these are set high enough that it is rare for approved catalogs to be dominated
by CRs (even though they can obviously have some CRs included.)

Finally, the `n1_residual` term gets set as a configuration parameter with a default value of 5% (0.05).  This
indicates that the single-image cosmic-ray identification process was expected to leave 5% of the cosmic-rays
unflagged. This process can be affected by numerous factors, and having this as a user settable parameter allows
the user to account for these effects when reprocessing the data manually.  Pipeline processing, though, may
still be subject to situations where this process does not do as well which can result in a catalog with a
higher than expected contribution of cosmic-rays.  Should this number of sources trigger the rejection criteria,
these catalogs will be rejected and not written out.

Also note that we reject both the point and segment catalogs if either one fails this test.  The reasoning
behind that is that since the catalogs are based on the same image, it is unlikely that one catalog will be
good and the other contaminated.

Should the catalogs fail this test, neither type of catalogs will be written out to disk for this visit.


Segment Photometric Catalog Generation
=======================================
Michele's documentation goes here!
