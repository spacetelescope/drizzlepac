.. drizzlepac api documentation master file, created by
   sphinx-quickstart on Wed Sep 15 14:22:04 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to drizzlepac's API documentation!
===========================================
This package supports the use of AstroDrizzle as an integrated set of modules that can be run in an automated manner to combine images, along with other tasks to support image alignment and coordinate transformations with distortion included.  The version of DrizzlePac described here implements a single task to run the entire AstroDrizzle processing pipeline, while also providing the framework for users to create their own custom pipeline based on the modules in this package merged with their own custom code if desired. These pages document what functions and classes are available for use under Python while providing the syntax for calling those functions from Python tasks.

Full documentation of how to run the primary AstroDrizzle and tweakreg tasks, along with fully worked examples, can be found in the `DrizzlePac Handbook available at the http://drizzlepac.stsci.edu <http://drizzlepac.stsci.edu>`_.
 
This package relies on the STWCS and PyWCS packages in order to provide the support for the WCS-based distortion models and alignment of the input images.


Contents:
 
.. toctree::
   :maxdepth: 2

   astrodrizzle
   imageobject
   process 
   static
   sky
   drizzle
   median
   blot
   drizcr
   util
   
DrizzlePac Release Notes
------------------------
The code for this package gets released through a number of methods: namely, the use of the package for pipeline and archive processing of ACS and WFC3 data, SSB's semi-annual public release of the stsci_python package, and a weekly beta release of the development version.  The following notes provide some details on what has been revised for each version. 

.. toctree::
    :maxdepth: 1
    
    release_notes

Image Registration Tasks
------------------------
Documentation for the replacement task for IRAF's `tweakshifts`, currently named `tweakreg`, has been added to this package. These new modules describe how to run the new TEAL-enabled task, as well as use the classes in the task to generate catalogs interactively for any chip and work with that catalog. The current implementation of this code relies on a very basic source finding algorithm loosely patterned after the DAOFIND algorithm and does not provide all the same features or outputs found in DAOFIND. The fitting algorithm also reproduces the fitting performed by IRAF's `geomap` in a limited fashion; primarily, it only performs fits equivalent to `geomap`'s 'shift' and 'rscale' solutions. These algorithms will be upgraded as soon as replacements are available. 

.. toctree::
   :maxdepth: 2
   
   tweakreg
   imagefindpars
   image
   catalogs
   wcscorr
   tweakutils
   updatehdr

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


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

 
 
