.. _release_2.0.0_notes:

**************************************
DrizzlePac v2.1.0 Release Notes
**************************************
This version builds upon the major set of changes implemented in v2.0.0 by not
only fixing some bugs, but also cleaning up/changing/revising some APIs and 
docstrings.  The complete list of changes includes:

- The 'updatewcs' parameter was removed from both the 'astrodrizzle' and 'tweakreg' interactive TEAL interfaces.  The 'updatewcs' parameter can still be used with the Python interface for both the astrodrizzle.Astrodrizzle() and tweakreg.TweakReg() functions instead of calling the 'stwcs.updatewcs.updatewcs()' function separately before running 'astrodrizzle' or 'tweakreg'. 
- The stand-alone interface for the blot routine (ablot.blot()) has been revised to work seamlessly with astrodrizzle-generated products while being more obvious how to call it correctly. The help file for this task was also heavily revised to document all the input parameters and to provide an example of how to use the task.
- Both astrodrizzle and tweakreg now return an output CD matrix which has identical cross-terms indicating the same scale and orientation in each axis. This relies on a revision to the stwcs.distortion.utils.output_wcs() function.
- The user interfaces to all 3 coordinate transformation tasks now use 'coordfile' as the input file of coordinates to transform. The use of 'coords' has been deprecated, but still can be used if needed. However, use of 'coordfile' will always override any input provided simultaneously with 'coords' parameter.  Help files have been updated to document this as clearly as possible for users. 
- User-provided list of input catalogs no longer needs to be matched exactly with input files. As long as all input images are included in input catalog list in any order, tweakreg will apply the correct catalog to the correct file.
- Tweakreg has been updated to correctly and fully apply source selection criteria for both input source catalogs and reference source catalogs based on fluxmin,fluxmax and nbright for each.
- All use of keyword deletion has been updated in drizzlepac (and fitsblender) to avoid warnings from astropy.
- All 3 coordinate transformation tasks rely on the input of valid WCS information for the calculations. These tasks now warn the user when it could not find a valid WCS and instead defaulted to using a unity WCS, so that the user can understand what input needs to be checked/revised to get the correct results.
- Exclusion/inclusion region files that can be used with 'tweakreg' can only be specified in image coordinates, not sky coordinates and will only support files written out using
DS9-compatible format. 
- The reported fit as written out to a file has been slightly modified to report more appropriate numbers of significant digits for the results. 
- Use of astrolib.coords was removed from drizzlepac and replaced by use of astropy functions instead. This eliminated one more obsolete dependency in our software.
- Code was revised to rely entirely on astropy.wcs instead of stand-alone pywcs.
- Code was revised to rely entirely on astropy.io.fits instead of stand-alone pyfits.


**************************************
DrizzlePac v2.0.0 Release Notes
**************************************
This version encompasses a large number of updates and revisions to the DrizzlePac code, including several parameter name changes and additions.  The scope of these changes indicates the level of effort that went into improving the DrizzlePac code to make it easier and more productive for users. 

Summary of Revisions
=====================
The most significant updates to the DrizzlePac code include:

  - The Python code has been updated to work identically (without change) under both Python 2.7 and Python 3.x.
  - Implementing sky matching, a new algorithm for matching the sky across a set of images being combined by astrodrizzle 
  - Updating tweakreg to now align full mosaics where some images may not overlap others in the mosaic
  - Added the option to write out single drizzle step images as compressed images (to save disk space for large mosaics, and I/O time for single drizzle step)
  - Improved tweakreg residual plots visually while allowing them to be written out automatically when tweakreg gets run in non-interactive mode
  - Renamed parameters in tweakreg and imagefind to eliminate name clashes
  - Added option to select sources based on sharpness/roundness when tweakreg searches for sources
  - Added support for exclusion and inclusion regions arbitrary shape/size when tweakreg searches for sources
  - Added a full set of source detection parameters for reference image to support multi-instrument alignment in tweakreg
  - Added support for new (simpler, more robust) ACS calibration of time-dependent distortion
  - A full 6-parameter general linear fit can now be performed using tweakreg, in addition to shift and rscale
  - Cleaned up logic for sky-subtraction: user can now turn off sky-subtraction with skysub=no, and still specify a user-defined sky value as the skyuser keyword.  This will reduce(eliminate?) the need to manually set MDRIZSKY=0. 
  
In addition to these major updates/changes, numerous smaller bugs were fixed and other revisions were implemented which affected a small portion of the use cases, such as:
  - headerlet code now accepts lists of files to be updated
  - source sky positions (RA and Dec) now included in match file
  - DQ flags can now be taken into account when performing source finding in tweakreg
  - all intermediate files generated by astrodrizzle will now be removed when using 'clean'='yes'
  - a problem was fixed that caused createMedian to crash where there were no good pixels in one of the images (when they did not overlap)
  - interpretation of shiftfile now improved to handle arbitrarily-long filenames, rather than being limited to 24 character filenames
  - documentation has been updated, sometimes with a lot more extensive descriptions

This version of DrizzlePac also requires use of the latest release version of astropy primarily for WCS and FITS I/O support.


