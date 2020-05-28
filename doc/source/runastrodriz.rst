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
        e) if alignment verified, copy updated input exposures to parent directory

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


