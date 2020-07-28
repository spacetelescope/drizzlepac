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
As described in the previous step, AstroDrizzle creates a single multi-filter detector-level drizzle-combined
image for source identification and one or more detector/filter-level drizzle-combined images (depending on
which filters were used in the dataset) for photometry. The same set of sources identified in the
multi-filter detection image is used to measure photometry for each filter. We use method to maximize the
signal across all available wavelengths at the source detection stage, thus providing photometry with the
best quality source list across all available input filters.

It should also be stressed here that the point and segment photometry source list generation algorithms
identify source catalogs independently of each other and DO NOT use a shared common source catalog for
photometry.

1.2: Preliminaries
^^^^^^^^^^^^^^^^^^^^

1.2.1: Generation of the Bad Pixel Mask
""""""""""""""""""""""""""""""""""""""""""""""""
Before any source identification takes place, we created a bad pixel mask to identify regions of the
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
correct for local small-scale overestimates and/or underestimates.

1.3: Source Identification with DAOStarFinder
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
We use the `photutils.detection.DAOStarFinder <https://photutils.readthedocs.io/en/stable/api/photutils.background.Background2D.html>`_ Astropy tool to identify sources in the background-subtracted
multi-filter detection image. Regions flagged in the previously created bad pixel mask are ignored by
DAOStarFinder. This algorithm works by identifying local brightness maxima with a roughly gaussian
distributions whose peak values are above a predefined minimum threshold. Full details of the process are
described in `Stetson 1987; PASP 99, 191 <http://adsabs.harvard.edu/abs/1987PASP...99..191S>`_.
The exact set of input parameters fed into DAOStarFinder is detector-dependent. The parameters can be found in
the instrument>_<detector>_catalog_generation_all.json files in the following path:
/drizzlepac/pars/hap_pars/default_parameters/<instrument>/<detector>/.


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
tool. Input values are detector-dependent, and can be found in the \*_catalog_generation_all.json files described above in section 1.3.

Local background values are computed based on the 3-sigma-clipped mode of pixel values present in a circular annulus
with an inner radius of 0.25 arcseconds and an outer radius of 0.50 arcseconds surrounding each identified source. This
local background value is then subtracted from the raw inner and outer aperture flux values to compute the
background-subtracted inner and outer aperture flux values found in the output .ecsv catalog file by the formula

.. math::

    f_{bgs}= f_{raw} - f_{bg} \cdot a,

where *f*\ :sub:`bgs`\  is the background-subtracted flux, *f*\ :sub:`raw`\  is the raw, unbackground-subtracted flux,
and *a* is the area of the photometric aperture. The overall standard deviation and mode values of pixels in the
background annulus are also reported for each identified source in the output .ecsv catalog file in the “STDEV” and
“MSKY” columns respectively (see Section 3 for more details).

Segment Photometric Catalog Generation
=======================================
Michele's documentation goes here!