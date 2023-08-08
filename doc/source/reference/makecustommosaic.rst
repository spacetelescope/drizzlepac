.. _makecustommosaic_api:

==========================
API for make_custom_mosaic
==========================
.. toctree::
   :maxdepth: 2

The task ``make_custom_mosaic`` serves as the primary interface for processing data from a multiple visits that span
multiple skycells into a uniform set of images.

Frequently Asked Questions
--------------------------
- Why would a user want to run this code?
    - Standard multi-visit mosaic (MVM) processing only processes a single skycell per run. Datasets with images that
      span one or more skycell boundaries require multiple MVM processing runs and produce multiple overlapping product
      images that are segmented at skycell boundaries. This code simplifies this process. It creates seamless
      multi-visit mosaics for input datasets that span one or more skycell boundaries.

- What are the required and opitonal inputs? What do they do?
    - The only required to make_custom_mosaic.py is either the name of a text file containing a list of calibrated ACS
      or WFC _flc.fits or _flt.fits files process or a search search sting that will be used to  find input files.
      Optional input arguments control the level of verbosity in the log statements that are displayed to the screen
      and written to log files, the formatting of  output filenames, and the method used to align input images. These
      inputs are discussed in more detail in the **USAGE** section.

- When the code is run, what output products should be expected? Do they differ depending on the inputs?
    - As this script depends on ``hapmultisequencer`` to perform the actual image processing, users can expect the
      output products of make_custom_mosaics.py to be functionally identical to standard MVM products. As with standard
      MVM processing, the only variation in output images occurs when processing WFC3/IR data, which produces output
      imagery products at the standard platescale of 0.04 arcseconds per pixel, as well as at the native WFC3/IR
      detector platescale of 0.12 arcseconds per pixel.

.. automodule:: drizzlepac.make_custom_mosaic
    :members: