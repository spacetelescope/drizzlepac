.. drizzlepac api documentation master file, created by
   sphinx-quickstart on Wed Sep 15 14:22:04 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. tocdepth: 2

Drizzlepac
==========
.. _drizzlepac:

The Drizzlepac package is a collection packages for aligning and combining astronomical images and coordinate transformations. Drizzlepac includes ``AstroDrizzle``, ``Tweakreg``, and the Hubble Archival Program (HAP) processing pipline, responsible for the creation of Hubble Single Visit Mosaics (SVMs) and Multi Visit Mosaics (MVMs).

.. toctree:: 
    :maxdepth: 2

    getting-started

Drizzling images
================
.. _drizzling-images:

* introduction
* Primary User Interface: AstroDrizzle()

Tweakreg & Tweakback
====================
.. _tweakreg-tweakback:

* image registration tasks
* Coordinate transformation tasks

Pipeline processing
===================
.. _pipeline-processing:

* introduction pipeline processing
* overview of pipeline
* Running astrodriz

Hubble Archival Products
========================
.. _hubble-archival-products:

* introduction
* Astrometry and Advanced Pipeline Products
* Improvements to HST astrometry
* Astrometry and headerlets
* HAP parameters
* Glossary of HAP terms

SVM
~~~
.. _SVM:

* processing
* catalog generation
* Classes to manage Catalogs and WCS’s

MVM
~~~
.. _MVM:

* processing
* interfaces
* utilities

Data products
=============
.. _data-products:

* Drizzle
* SVM
* MVM

API
===
.. _API:

* Drizzle (Image drizzling step)
* Pipeline processing
* Tweakreg, Tweakback
* SVM
* MVM
* Dependencies
* Licence








     
.. users can create custom pipelines 
.. ``AstroDrizzle`` from `~drizzlepac.astrodrizzle`


.. Contents:


..    :maxdepth: 1

..    astrodrizzle
..    imageobject
..    process
..    static
..    sky
..    adrizzle
..    median
..    ablot
..    drizcr
..    util
..    tweakback
..    LICENSE

.. dependencies
.. This package relies on the STWCS package in order to provide the support for the WCS-based distortion models and alignment of the input images.



.. DrizzlePac Release Notes
.. ------------------------
.. The code for this package gets released through a number of methods: namely,
.. the use of the package for pipeline and archive processing of ACS and WFC3
.. data, SSB's semi-annual public release of the stsci_python package, and a
.. weekly beta release of the development version.  The following notes
.. provide some details on what has been revised for each version.


.. Image Registration Tasks
.. ------------------------
.. A number of tasks have been developed to support the registration
.. of images.  This includes documentation for the replacement task
.. for IRAF's ``tweakshifts``, currently named ``TweakReg``, along
.. with tasks for updating the WCS in HST images and performing
.. photometry equalization for WFPC2 data.

.. In addition, several tasks have been developed to perform
.. coordinate transformations that take into account the full
.. distortion model specified in HST image headers.

.. The full description of all of these tasks have been added
.. to the :ref:`imagereg` page.

.. .. toctree::
..    :maxdepth: 1

..    image_registration


.. Reproducing Pipeline Processing
.. -------------------------------
.. The task 'runastrodriz' can be used to reproduce the same Drizzle processing that gets performed on HST data when retrieving data from the HST archive.

.. .. toctree::
..    :maxdepth: 1

..    Running Astrodriz <runastrodriz.rst>
..    API for Astrometry Code <astrometry_api.rst>


.. Astrometry and Advanced Pipeline Products
.. ------------------------------------------
.. The ``drizzlepac`` package can be used for many purposes, all related to aligning and combining images to create products which can provide the deepest available views of the data.  Combining the data with ``drizzlepac`` relies on the WCS solution specified in the input image headers.  These WCS solutions are expected to align the images to each other (relative astrometry) as well as align the image to the correct position on the sky (absolute astrometry).  The telemetry from HST allows the relative astrometry to be known extremely accurately (sub-milli-arcsecond level) when all images use the same guide stars and when the images were taken in the same visit.  However, data taken at different times using different guide stars have historically had errors in the alignment with a sigma of 1 arc-second (or more).  As a result, corrections to the alignment need to be made in order to successfully combine the images.

.. The code being used in the automated HST calibration pipeline to generate the products served by the HST Archive through the MAST Portal relies on multiple methods to do the best job possible in aligning and combining data regardless of whether they came from the same visit (or instrument) or not.

.. On December 17, 2020, MAST began production of new ACS and WFC3 products in the HST data calibration pipeline: Hubble Legacy Archive (HLA)-style mosaics comprising the data from a single HST visit which are aligned to a common astrometric reference frame. These are the first of two types of **Hubble Advanced Products (HAP)** that will be produced in the HST data pipeline and made available through the `MAST Discovery Portal <https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html>`_.

.. Three levels of products are available as part of this release:

..   * Exposure level products contain data from a single HST exposure.
..   * Filter level products are produced from all exposures in a visit with a common filter.
..   * Total level products combine all exposures from a visit for a specific detector and are intended as a detection image for producing catalogs.

.. The **HAP Single Visit Mosaics (SVMs)** differ from the standard HST drizzled data products, which are aligned filter-by-filter to Gaia.  SVM data products, on the other hand, are all drizzled onto the same north-up pixel grid and may have improved relative alignment across filters within a given visit, enabling easy comparison of the images through multiple filters or for images to be combined to create color mosaics. When possible, sources in the images have been aligned directly to the Gaia source catalog to improve the WCS of the images. SVM data products with both relative alignment (by filter) and absolute alignment to Gaia will contain the string 'FIT_SVM_GAIA' in the 'WCSNAME' keyword in the science extension of the image header. More discussion on HAP alignment, may be found on the webpage `Improvements in HST Astrometry <https://outerspace.stsci.edu/pages/viewpage.action?spaceKey=HAdP&title=Improvements+in+HST+Astrometry>`_.

.. .. toctree::
..   :maxdepth: 1

..   Astrometry and Headerlets <astrometry.rst>

..   singlevisit.rst
..   runsinglehap.rst

..   catalog_generation
..   catalogs


.. Combining data across visits to create mosaics for every observed location on the sky gets performed using
.. **HAP Multi-Visit Mosaic (MVM)** processing.

.. .. toctree::
..   :maxdepth: 1

..   Multi-Visit Mosaic Processing <multivisit.rst>
..   Multi-Visit Mosaic Products <multivisit_products.rst>
..   runmultihap.rst
..   multivisit_api.rst

..   makecustommosaic.rst
..   mvmutilities.rst
..   mvmutilities_api.rst


.. Indices and tables
.. ==================
.. .. toctree::
..   :maxdepth: 1

..   Glossary of HAP Terms <hap_glossary.rst>

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
