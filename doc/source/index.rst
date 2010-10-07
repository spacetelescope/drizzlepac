.. betadrizzle documentation master file, created by
   sphinx-quickstart on Wed Sep 15 14:22:04 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to betadrizzle's documentation!
=======================================
This package supports the use of MultiDrizzle as an integrated set of modules that can be run in an automated manner to combine images.  The version of MultiDrizzle described here implements a single task to run the entire MultiDrizzle processing pipeline, while also providing the framework for users to create their own custom pipeline based on the modules in this package merged with their own custom code if desired. 

This package relies on the STWCS and PyWCS packages in order to provide the support for the WCS-based distortion models and alignment of the input images.


Contents:
 
.. toctree::
   :maxdepth: 2

   betadrizzle
   imageobject
   process 
   static
   sky
   drizzle
   median
   blot
   drizcr
   util
   
Image Registration Tasks
------------------------
Documentation for the replacement task for IRAF's `tweakshifts`, currently named `tweakreg`, has been added to this package. These new modules describe how to run the new TEAL-enabled task, as well as use the classes in the task to generate catalogs interactively for any chip and work with that catalog. The current implementation of this code relies on a very basic source finding algorithm loosely patterned after the DAOFIND algorithm and does not provide all the same features or outputs found in DAOFIND. The fitting algorithm also reproduces the fitting performed by IRAF's `geomap` in a limited fashion; primarily, without iterations of outliers and only to perform fits equivalent to `geomap`'s 'shift' and 'rscale' solutions. These algorithms will be upgraded as soon as replacements are available. 

.. toctree::
   :maxdepth: 2
   
   userint 
   image
   catalogs
   wcscorr
   
Building a TEAL Interface for Tasks
-----------------------------------
.. toctree::

   teal_guide

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

 
 
