.. _running-astrodrizzle:

********************
Running Astrodrizzle
********************

``runastrodriz`` is a module to control operation of astrodrizzle which removes distortion and combines HST images in the pipeline.


Typical Usage
=============

    >>> runastrodriz.py [-fhibn] inputFilename [newpath]


Alternative Usage
=================

    >>> python
    >>> from wfc3tools import runastrodriz
    >>> runastrodriz.process(inputFilename,force=False,newpath=None,inmemory=False)


GUI Usage under Python
======================

    >>> python
    >>> from stsci.tools import teal
    >>> import wfc3tools
    >>> cfg = teal.teal('runastrodriz')

PyRAF Usage
===========

    >>> epar runastrodriz



Options
=======

If the '-i' option gets specified, no intermediate products will be written out
to disk. These products, instead, will be kept in memory. This includes all
single drizzle products (*single_sci* and *single_wht*), median image,
blot images, and crmask images.  The use of this option will therefore require
significantly more memory than usual to process the data.

If a value has been provided for the newpath parameter, all processing will be
performed in that directory/ramdisk.  The steps involved are:

    * create a temporary directory under that directory named after the input file
    * copy all files related to the input to that new directory
    * change to that new directory and run astrodrizzle
    * change back to original directory
    * move (not copy) ALL files from temp directory to original directory
    * delete temp sub-directory

The '-b' option will run this task in BASIC mode without creating headerlets
for each input image.

The '-n' option allows the user to specify the number of cores to be used in
running AstroDrizzle.


.. _runastrodriz-description:

Description
===========
A lot of effort goes into trying to determine the WCS solution which provides
closest agreement with the GAIA astrometry frame.  Combining the 
images successfully requires the most accurate alignment possible between all 
the input exposures taking into account the best available distortion model
for the input exposures.  Being able to compare the combined image to other exposures
taken of the same field at other times, or even with other telescopes, requires
that the WCS be defined to come as close to the GAIA frame as possible.  The 
processing done by `runastrodriz` attempts to not only apply the most current 
distortion model, but also the best available pre-computed GAIA-based WCS 
solutions while proceeding to determine it's own solution to the available GAIA
reference stars for the field-of-view.  It then performs numerous checks to see 
which of these solutions results in the most accurately aligned images to each 
other, then to GAIA and selects that WCS solution.  This selected WCS solution 
serves as the basis for creating the final, distorted-corrected, combined 
drizzle products for the set of exposures being processed. 

Overview
--------
The overall logic implemented by `runastrodriz` to generate the final set of 
drizzle products involves creating multiple sets of drizzle products and ultimately
selecting one as the 'best'.  The basic steps are outlined here, with following 
sections providing additional detail on each step.

    #. Initialization
    
      a) Get any defined environment variables to control processing
      b) Interpret input file
      c) Make sure all input files are in a local directory
      d) Check for DRIZCORR calibration switch in input files
      e) Create name for output log files
      f) Define lists of individual input files to be processed

    #. Run updatewcs without astrometry database update on all input exposures (FLCs? and FLTs)

    #. Generate initial default products and perform verification

        a) perform cosmic-ray identification and generate drizzle products using astrodrizzle for all sets of inputs
        b) verify relative alignment with focus index after masking out CRs
        c) copy all drizzle products to parent directory
        d) if alignment fails, update trailer file with failure information

    #. If alignment is verified,

        a) copy inputs to separate sub-directory for processing
        b) run updatewcs to get a priori updates
        
          * apply 'best' apriori (not aposteriori) solution

        c) generate drizzle products for all sets of inputs (FLC and/or FLT) without CR identification
        d) verify alignment using focus index on FLC or, if no FLC, FLT products
        e) if alignment fails, update trailer file with info on failure
        f) if product alignment verified:
        
            * copy all drizzle products to parent directory
            * copy updated input exposures to parent directory

    #. If a posteriori correction enabled,

        a) copy all inputs to separate sub-directory for processing
        b) run align to align the images
        c) generate drizzle products for all sets of inputs (FLC and/or FLT) without CR identification
        d) verify alignment using focus index on FLC or, if no FLC, FLT products
        e) determine similarity index relative to pipeline default product
        f) if either focus or similarity indicates a problem, update trailer file with info on failure
        g) if product alignment verified:
        
           * copy all drizzle products to parent directory
           * copy updated input exposures to parent directory

    #. Remove all processing sub-directories


Initialization
--------------

Environment Variables
^^^^^^^^^^^^^^^^^^^^^^
The pipeline processing code starts out by looking to see whether the user has defined any processing behavior through the use of these environment variables:

  * **'ASTROMETRY_COMPUTE_APOSTERIORI'**: This environment variable specifies whether or not to attempt an *a posteriori* alignment where the code looks for sources in each of the images and uses those positions to perform relative alignment between the images and then fit those images to the GAIA frame.  
  * **'ASTROMETRY_APPLY_APRIORI'**: This environment variable turns on/off application of any pre-defined(*a priori*) WCS solution found in the astrometry database.  
  * **'ASTROMETRY_STEP_CONTROL' [DEPRECATED, do not use]**: Old variable replaced by 'ASTROMETRY_APPLY_APRIORI'.   

Values that can be provided for setting these variables are:

  * 'on', 'yes', 'true': Any of these values will turn **on** the processing controlled by the variable
  * 'off', 'no', 'false': Any of these values will turn **off** the processing controleed by the variable

By default, all the processing steps are turned **on** during pipeline processing in order to maximize the chances of aligning the data as closely as possible to the absolute astrometry standard coordinate system defined through the use of the GAIA catalogs.  However, these controls are provided to support those observations which would not be suitable for such alignment, including observations of single sources.

Input Data
^^^^^^^^^^^
The processing code needs to be told what data to process, and for `runastrodriz`, a single input filename is all that **can** be provided.  This single input will be either:

  * the name of an association table for a whole set of input exposures with a filename that looks like **'<rootname>_asn.fits'**, where <rootname> is the designation for the association, such as *'ie6d07030_asn.fits'*.  
  * the name of a single (uncalibrated) exposure with a filename that looks like **'<rootname>_raw.fits'**.

This one input filename, though, will simply provide the code with the information it needs to find all the calibrated input exposures which need to have their distortion-models updated and applied.  The whole set of input files required for processing includes:

  * ASN (``*_asn.fits``) files: These small FITS tables provide the relationship between the input exposures and the output products with the output filenames defined in the table.  There will NOT be an ASN table for exposures which were taken by themselves (called 'singletons').  
  * RAW (``*_raw.fits``) files: Not processed directly, but required in order to get the intended value of the `DRIZCORR` calibration switch.  The ASN files also only give the rootname, and with the possibility of multiple suffixes (_flt, _flc,...) for calibrated products, the code starts with the _raw files to insure that what is specified in the ASN table is actually present and has been calibrated before processing.
  * FLT/FLC (``*_flt.fits`` or ``*_flc.fits``) files: These are the non-CTE-corrected (_flt) and CTE-corrected (_flc) calibrated exposures to be processed.
  
The FLT/FLC files will be the ones that actually get processed and updated with the new distortion models and WCSs, while the others allow the code to know what FLT/FLC files should be included in the processing.  This allows for multiple associations of data to live in the same directory and not interfere with each other as they are re-processed.  That can be useful when interested in combining data from multiple visits, for example.  

.. warning::  Should any of these files not be available (found in the local directory), the code will raise an Exception when trying to run 'drizzlepac.astrodrizzle.AstroDrizzle' on the data.  The message will indicate what file was missing with something like: **"Exception: File ie6d07ujq_flt not found."**

Calibration Switches
^^^^^^^^^^^^^^^^^^^^
This processing serves as an official calibration step defined for HST data through the use of the **DRIZCORR** header keyword.  This keyword can be found along with all the other calibration switches in the PRIMARY header (extension 0) of the exposures FITS file. A quick way to view this (or any keyword) value would be with:

.. code:: python

    from astropy.io import fits
    val = fits.getval('ie6d07ujq_flt.fits', 'drizcorr')
    

This switch must be set to 'PERFORM' in order to allow the processing to be done. Processing will be completely skipped should the value of this switch in the '_raw.fits' file be set to 'OMIT'.

Log Files
^^^^^^^^^^
A number of log files, or 'trailer' files, get generated during processing, and their filenames get defined as early as possible in the processing.  The primary file will be a file with a '.tra' extension and should have the same '<rootname>' as the input file used to start the code.  For example, if you were to reprocess 'ie6d07030_asn.fits', you would end up with a trailer file with the name 'ie6d07030.tra'.  

This log file contains the messages generated from performing all the updates to the distortion model, updates from the astrometry database (if any), and all the image combinations performed by 'AstroDrizzle()' to create the final set of calibrated, drizzled exposures.  Should any problems arise when during the processing, the log can provide the error messages and tracebacks to determine what went wrong.


Data to be Processed
^^^^^^^^^^^^^^^^^^^^^
Once the code has performed all the initialization, it prepares the processing by defining what files need to be combined together from the input files it can find.  This includes looking for CTE-corrected versions of the calibrated exposures (FLC files) as well as all the non-CTE-corrected files (FLT files) and creating a separate list of each type.  Many types of data do not get CTE-corrected by the instruments calibration software, such as calacs.e or calwf3.e, and so no list of FLC files will be made.  This will tell the code that it only needs to process the FLT files by themselves.  If FLC files are found, all updates to the astrometry and WCS will be performed on those files and the results then get copied into the FLT file headers upon completion of the processing.  


Updating the WCS
----------------
The first operation on the calibrated input files focuses on applying the calibrations
for the distortion model to the WCS.  This operation gets performed using the 
`updatewcs` task using the syntax:

.. code:: python

    from stwcs.updatewcs import updatewcs
    updatewcs(calfiles_flc, use_db=False)
    
where `calfiles_flc' is the list of CTE-corrected FLC files or in the case there are
no CTE-corrected files, the list of calibrated FLT files.  Crucially, the use
of `use_db=False` forces `updatewcs` to only apply the distortion model to the
default WCS to create what is referred to as the **pipeline-default WCS**.  This
WCS has a `WCSNAME` associated with it that has the format ``IDC_<rootname>`` where
``<rootname>`` is the rootname of the `IDCTAB` reference files applied to the WCS. 

This default WCS serves as the basis for all subsequent processing as the code
tries to determine the WCS which is aligned most closely to the GAIA astrometric
coordinate system.  

 













