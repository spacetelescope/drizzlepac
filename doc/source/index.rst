.. drizzlepac api documentation master file, created by
   sphinx-quickstart on Wed Sep 15 14:22:04 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to drizzlepac's API documentation!
===========================================
This package supports the use of ``AstroDrizzle`` as an integrated set of
modules that can be run in an automated manner to combine images, along
with other tasks to support image alignment and coordinate transformations
with distortion included.  The version of ``DrizzlePac`` described here
implements a single task to run the entire ``AstroDrizzle`` processing
pipeline, while also providing the framework for users to create their own
custom pipeline based on the modules in this package merged with their own
custom code if desired. These pages document what functions and classes are
available for use under Python while providing the syntax for calling those
functions from Python tasks.

Full documentation of how to run the primary ``AstroDrizzle`` and ``TweakReg``
tasks, along with fully worked examples, can be found in the
`DrizzlePac Handbook <http://drizzlepac.stsci.edu>`_.

This package relies on the STWCS package in order to provide
the support for the WCS-based distortion models and alignment
of the input images.

Contents:

.. toctree::
   :maxdepth: 2

   astrodrizzle
   imageobject
   process
   static
   sky
   adrizzle
   median
   ablot
   drizcr
   util
   tweakback
   LICENSE

DrizzlePac Release Notes
------------------------
The code for this package gets released through a number of methods: namely,
the use of the package for pipeline and archive processing of ACS and WFC3
data, SSB's semi-annual public release of the stsci_python package, and a
weekly beta release of the development version.  The following notes
provide some details on what has been revised for each version.

.. toctree::
    :maxdepth: 1

    CHANGELOG

Image Registration Tasks
------------------------
Documentation for the replacement task for IRAF's ``tweakshifts``,
currently named ``TweakReg``, has been added to this package.
These new modules describe how to run the new ``TEAL``-enabled task,
as well as use the classes in the task to generate catalogs interactively
for any chip and work with that catalog. The current implementation of this
code relies on a very basic source finding algorithm loosely patterned
after the DAOFIND algorithm and does not provide all the same features
or outputs found in DAOFIND. The fitting algorithm also reproduces the
fitting performed by IRAF's ``geomap`` in a limited fashion; primarily,
it only performs fits equivalent to ``geomap``'s 'shift' and 'rscale'
solutions. These algorithms will be upgraded as soon as replacements
are available.

.. toctree::
   :maxdepth: 2

   tweakreg
   refimagefindpars
   imagefindpars
   image
   catalogs
   catalog_generation
   wcscorr
   tweakutils
   updatehdr
   mapreg
   photeq
   pixreplace

Coordinate Transformation Tasks
-------------------------------
These tasks support transformations of source positions to and from
distorted and drizzled images.

.. toctree::
   :maxdepth: 1

   pixtopix
   pixtosky
   skytopix


ACS Header Update Task
----------------------
A task, 'updatenpol', has been written to automate the updating of ACS image headers with the filename of the appropriate NPOLFILE based on the DGEOFILE specified in the image header.  This task should be used to update all ACS images prior to processing them with 'astrodrizzle'.

.. toctree::
   :maxdepth: 2

   updatenpol


Reproducing Pipeline Processing
-------------------------------
The task 'runastrodriz' can be used to reproduce the same Drizzle processing that gets performed on HST data when retrieving data from the HST archive.

.. toctree::
   :maxdepth: 2

   Running Astrodriz <runastrodriz.rst>
   API for Astrometry Code <astrometry_api.rst>


Astrometry and Advanced Pipeline Products
------------------------------------------
The `drizzlepac` package can be used for many purposes, all related to aligning and combining images to create products which can provide the deepest available views of the data.  Combining the data with `drizzlepac` relies on the WCS solution specified in the input image headers.  These WCS solutions are expected to align the images to each other (relative astrometry) as well as align the image to the correct position on the sky (absolute astrometry).  The telemetry from HST allows the relative astrometry to be known extremely accurately (sub-milli-arcsecond level) when all images use the same guide stars and when the images were taken in the same visit.  However, data taken at different times using different guide stars have historically had errors in the alignment with a sigma of 1 arc-second (or more).  As a result, corrections to the alignment need to be made in order to successfully combine the images.

The code being used in the automated HST calibration pipeline to generate the products served by the HST Archive through the MAST Portal relies on multiple methods to do the best job possible in aligning and combining data regardless of whether they came from the same visit (or instrument) or not.

On December 17, 2020, MAST began production of new ACS and WFC3 products in the HST data calibration pipeline: Hubble Legacy Archive (HLA)-style mosaics comprising the data from a single HST visit which are aligned to a common astrometric reference frame. These are the first of two types of **Hubble Advanced Products (HAP)** that will be produced in the HST data pipeline and made available through the `MAST Discovery Portal <https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html>`_.

Three levels of products are available as part of this release:

  * Exposure level products contain data from a single HST exposure.
  * Filter level products are produced from all exposures in a visit with a common filter.
  * Total level products combine all exposures from a visit for a specific detector and are intended as a detection image for producing catalogs.

The **HAP Single Visit Mosaics (SVMs)** differ from the standard HST drizzled data products, which are aligned filter-by-filter to Gaia.  SVM data products, on the other hand, are all drizzled onto the same north-up pixel grid and may have improved relative alignment across filters within a given visit, enabling easy comparison of the images through multiple filters or for images to be combined to create color mosaics. When possible, sources in the images have been aligned directly to the Gaia source catalog to improve the WCS of the images. SVM data products with both relative alignment (by filter) and absolute alignment to Gaia will contain the string 'FIT_SVM_GAIA' in the 'WCSNAME' keyword in the science extension of the image header. More discussion on HAP alignment, may be found on the webpage `Improvements in HST Astrometry <https://outerspace.stsci.edu/pages/viewpage.action?spaceKey=HAdP&title=Improvements+in+HST+Astrometry>`_.

.. toctree::
  :maxdepth: 2

  Astrometry and Headerlets <astrometry.rst>

  singlevisit.rst
  runsinglehap.rst

  Glossary of HAP Terms <hap_glossary.rst>


Indices and tables
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
