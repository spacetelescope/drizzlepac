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
As previously discussed in :ref:`singlevisit`, AstroDrizzle creates a single multi-filter detector-level drizzle-combined
image for source identification and one or more detector/filter-level drizzle-combined images (depending on
which filters were used in the dataset) for photometry. The same set of sources identified in the
multi-filter detection image is used to measure photometry for each filter. We use this method to maximize the
signal across all available wavelengths at the source detection stage, thus providing photometry with the
best quality source list across all available input filters.

It should also be stressed here that the point and segment photometry source list generation algorithms
identify source catalogs independently of each other and DO NOT use a shared common source catalog for
photometry.

1.2: Preliminaries
^^^^^^^^^^^^^^^^^^^^

1.2.1: Generation of the Bad Pixel Mask
""""""""""""""""""""""""""""""""""""""""""""""""
Before any source identification takes place, a bad pixel mask is created to identify regions of the
detection image where signal quality is known to be degraded. These are areas near the edge of the image,
areas with little to no input image contribution, and areas that contain saturated pixels. To minimize the
impact of these regions on source identification and subsequent photometric measurements, the regions flagged
in this bad pixel mask are iteratively “grown” for 10 steps using the `ndimage.binary_dilation <https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.binary_dilation.html>`_ scipy tool.
Pixels in the immediate vicinity of a given masked region may also be impacted to some degree. As we cannot
be fully certain that these pixels are or are not impacted, or to the degree of the impact, they are all
flagged.

1.2.2: Detection Image Background Subtraction
""""""""""""""""""""""""""""""""""""""""""""""""
To ensure optimal source detection, the multi-filter detection image is background-subtracted. We computed a
2-dimensional background image using the `photutils.background.Background2d <https://photutils.readthedocs.io/en/stable/api/photutils.background.Background2D.html>`_ Astropy tool. This algorithm uses
sigma-clipped statistics to determine background and RMS values across the image. An initial low-resolution
estimate of the background is performed by computing sigma-clipped median values in 27x27 pixel boxes across
the image. This low-resolution background image is then median-filtered using a 3x3 pixel sample window to
correct for local small-scale overestimates and/or underestimates. It should be noted these are configurable values.
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

Segment Photometric Catalog Generation
=======================================
Michele's documentation goes here!