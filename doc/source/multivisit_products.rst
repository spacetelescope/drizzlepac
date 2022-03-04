.. _multivisit_products:

=============================
Multi-visit Mosaic Products
=============================

Multi-Visit Mosaic (MVM)
    A Multi-Visit Mosaic (MVM) is a single image (product) made by combining all observations taken of the same part of the sky.

Observations taken of the same part of the sky over all the years that HST has been operational can enable unique science
due to the high-resolution of the HST cameras.  Generating useful mosaics from these observations, though, requires
solving a number of key problems; namely,

  * aligning all the images to the same coordinate system
  * define the size on the sky of each mosaic
  * defining what exposures should go into each mosaic

MVM processing implemented as part of the Hubble Advanced Products (HAP) pipeline generates new products based on
solutions implemented for these critical issues.  These new products, SkyCell layers, are unlike other HST images
due to the fact that
they consist of many exposures taken at different times, sometimes years apart.  The format of these products and
the new aspects of these products are described in the following sections.


SkyCell Layers
===============
The most basic MVM product would be the SkyCell layer as described in `Defining SkyCell Layers`_ section.  These layers
represent all the exposures taken in a given detector/filter combination in that SkyCell (position on the sky).  As a
example, the exposures for sky cell **p1889x07y19** define 9 separate layers; namely,

  * **wfc3_uvis_f475w (0.04"/pixel)** :  hst_skycell-p1889x07y19_wfc3_uvis_f475w_all_drz.fits
  * **wfc3_ir_f105w_coarse  (0.12"/pixel)** : hst_skycell-p1889x07y19_wfc3_ir_f105w_coarse-all_drz.fits
  * **wfc3_ir_f105w  (0.04"/pixel)** : hst_skycell-p1889x07y19_wfc3_ir_f105w_all_drz.fits
  * **wfc3_ir_f125w_coarse  (0.12"/pixel)** : hst_skycell-p1889x07y19_wfc3_ir_f125w_coarse-all_drz.fits
  * **wfc3_ir_f125w  (0.04"/pixel)** : hst_skycell-p1889x07y19_wfc3_ir_f125w_all_drz.fits
  * **wfc3_ir_f160w_coarse  (0.12"/pixel)** : hst_skycell-p1889x07y19_wfc3_ir_f160w_coarse-all_drz.fits
  * **wfc3_ir_f160w  (0.04"/pixel)** : hst_skycell-p1889x07y19_wfc3_ir_f160w_all_drz.fits
  * **acs_wfc_f850lp  (0.04"/pixel)** : hst_skycell-p1889x07y19_acs_wfc_f850lp_all_drc.fits
  * **acs_wfc_f775w  (0.04"/pixel)** : hst_skycell-p1889x07y19_acs_wfc_f775w_all_drc.fits

A full SkyCell would cover an area on the sky of approximately 0.2\deg x 0.2\deg with a WCS defined as:

.. code-block::

    Number of WCS axes: 2
    CTYPE : 'RA---TAN'  'DEC--TAN'
    CRVAL : 180.0  26.0
    CRPIX : 96492.0  -160812.0
    CD1_1 CD1_2  : -1.1111111111111112e-05  0.0
    CD2_1 CD2_2  : 0.0  1.1111111111111112e-05
    NAXIS : 21954  21954


SkyCell Subarray Specification
-------------------------------
SkyCell layers would normally result in arrays that take up 1.8Gb each, or nearly 4Gb for the entire FITS file.  In
addition, most of a typical SkyCell layer will be empty due to the small size of most detectors, especially the WFC3/IR,
ACS/HRC and ACS/SBC detectors which are only 1024x1024 arrays.  Each SkyCell then gets evaluated to only define a common
subarray size that covers ALL HST data for the SkyCell regardless of the layer, then uses that to define the smallest
WCS specification possible for the SkyCell.  The data from **p1889x07y19** actually only covers about 1/4 of the entire
SkyCell resulting in a WCS defined as:

.. code:: python

    >>> wcs1889 = HSTWCS('hst_skycell-p1889x07y19_acs_wfc_f775w_all_drc.fits', ext=1)
    >>> print(wcs1889)

    WCS Keywords
    Number of WCS axes: 2
    CTYPE : 'RA---TAN'  'DEC--TAN'
    CRVAL : 180.0  26.0
    CRPIX : 94360.0  -171093.0
    CD1_1 CD1_2  : -1.1111111111111e-05  0.0
    CD2_1 CD2_2  : 0.0  1.11111111111111e-05
    NAXIS : 10840  11672

We can see how the exposures land in the SkyCell as defined with this subarray WCS.

.. list-table::

  * - .. figure:: images/skycell-p1889x07y19_f775w_full.jpg
         :figwidth: 95%
         :alt: SkyCell p1889x07y19 WFC3/UVIS F775W layer.

         All the WFC3/UVIS F775W exposures that overlap SkyCell **p1889x07y19**.

    -  .. figure:: images/skycell-p1889x07y19_f105w_full.jpg
          :figwidth: 95%
          :alt: SkyCell p1889x07y19 WFC3/IR F105W layer.

          All the WFC3/IR F105W exposures that overlap SkyCell **p1889x07y19**.


These figures demnstrate how the SkyCell subarray has been defined to only cover the exposures that overlap this
SkyCell while minimizing the amount of empty space in these layers making these products as small as possible without
resorting to compression.

.. note::
  These layers are defined as subarrays of the entire SkyCell WCS, which are just subarrays of the ProjectionCell.
  As a result, any position in these layers refers to the exact same position on the sky
  as defined in the SkyCell or ProjectionCell.  A source at (x,y)=(7936, 11133) will have the same sky coordinates
  regardless of what SkyCell layer it was measured in, including an array defined for the entire SkyCell.

