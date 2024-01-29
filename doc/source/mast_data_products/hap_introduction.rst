Introduction
------------
.. _hap-introduction:

The ``drizzlepac`` package can be used for many purposes, all related to aligning and combining images to create products which can provide the deepest available views of the data.  Combining the data with ``drizzlepac`` relies on the WCS solution specified in the input image headers.  These WCS solutions are expected to align the images to each other (relative astrometry) as well as align the image to the correct position on the sky (absolute astrometry).  The telemetry from HST allows the relative astrometry to be known extremely accurately (sub-milli-arcsecond level) when all images use the same guide stars and when the images were taken in the same visit.  However, data taken at different times using different guide stars have historically had errors in the alignment with a sigma of 1 arc-second (or more).  As a result, corrections to the alignment need to be made in order to successfully combine the images.

The code being used in the automated HST calibration pipeline to generate the products served by the HST Archive through the MAST Portal relies on multiple methods to do the best job possible in aligning and combining data regardless of whether they came from the same visit (or instrument) or not.

On December 17, 2020, MAST began production of new ACS and WFC3 products in the HST data calibration pipeline: Hubble Legacy Archive (HLA)-style mosaics comprising the data from a single HST visit which are aligned to a common astrometric reference frame. These are the first of two types of **Hubble Advanced Products (HAP)** that will be produced in the HST data pipeline and made available through the `MAST Discovery Portal <https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html>`_.

Three levels of products are available as part of this release:

  * Exposure level products contain data from a single HST exposure.
  * Filter level products are produced from all exposures in a visit with a common filter.
  * Total level products combine all exposures from a visit for a specific detector and are intended as a detection image for producing catalogs.

The **HAP Single-Visit Mosaics (SVMs)** differ from the standard HST drizzled data products, which are aligned filter-by-filter to Gaia.  SVM data products, on the other hand, are all drizzled onto the same north-up pixel grid and may have improved relative alignment across filters within a given visit, enabling easy comparison of the images through multiple filters or for images to be combined to create color mosaics. When possible, sources in the images have been aligned directly to the Gaia source catalog to improve the WCS of the images. SVM data products with both relative alignment (by filter) and absolute alignment to Gaia will contain the string 'FIT_SVM_GAIA' in the 'WCSNAME' keyword in the science extension of the image header. More discussion on HAP alignment, may be found on the webpage `Improvements in HST Astrometry <https://outerspace.stsci.edu/pages/viewpage.action?spaceKey=HAdP&title=Improvements+in+HST+Astrometry>`_.

Combining data across visits to create mosaics for every observed location on the sky gets performed using
**HAP Multi-Visit Mosaic (MVM)** processing.
