.. _release_notes:

**************************************
DrizzlePac Release Notes
**************************************
The code for this package gets released through a number of methods: namely,
  - the use of the package for pipeline and archive processing of ACS and WFC3 data,
  - SSB's semi-annual `public release of the stsci_python package <http://www.stsci.edu/institute/software_hardware/pyraf/stsci_python/installation>`_, and
  - a weekly beta release of the development version as part of the `IRAFX download <http://stsdas.stsci.edu/irafx/>`_.

The following notes provide some details on what has been revised for each version in
reverse chronological order (most recent version at the top of the list).

DrizzlePac(astrodrizzle) v1.1.15(30-Dec-2013)
-------------------------------------------------
**Publicly Released through PyPI:** Jan 14, 2014

**available under SSBX/IRAFX starting:** Jan 6, 2014

- [Bug Fix] Files created or updated by drizzlepac, fitsblender, or STWCS tasks, e.g. tweakreg or apply_headerlet, will now ensure that the NEXTEND keyword value correctly reflects the number of extensions in the FITS file upon completion.


DrizzlePac(astrodrizzle) v1.1.14dev(21-Oct-2013) in IRAFX
---------------------------------------------------------
**Installed in OPUS:** Dec 11, 2013

**available starting:** Oct 28, 2013

- [Bug Fix] DQ arrays in input images now get updated with cosmic-ray masks computed by astrodrizzle when run with the parameter 'in_memory=True'. This restored the cosmic-ray masks detected during pipeline processing.


DrizzlePac(astrodrizzle) v1.1.13dev(11-Oct-2013) in IRAFX
---------------------------------------------------------
**available starting:** Oct 21, 2013

- Tweakreg can now be run in 'batch' mode.  This allows the user to generate plots and have them saved to disk automatically without stopping processing and requiring any user input.


DrizzlePac(astrodrizzle) v1.1.12dev(05-Sep-2013) in IRAFX
---------------------------------------------------------
**available starting:** Sept 9, 2013

This version fixed a couple of bugs in astrodrizzle; namely,
  - Logic was updated to support pixfrac = 0.0 without crashing. Ths code will now automatically reset the kernel to 'point' in that case.
  - Astrodrizzle now forcibly removes all OPUS WCS keywords from drizzle product headers.
  - Default rules for generating drizzle product headers (as used in the archive) were modified to add definitions for 'float_one', 'int_one', 'zero' that generate output values of 1.0, 1, and 0 (zero) respectively for use as keyword values. This allows the LTM* rules to replace 'first' with 'float_one' so that the physical and image coordinates for drizzle products are consistent.

Additionally, changes were made to STWCS for reprocessing use:
  - Problems with using apply_headerlet_as_primary() from the STWCS package on WFPC2 data have been corrected in this revision.


DrizzlePac(astrodrizzle) v1.1.11dev(05-Jul-2013) in IRAFX
---------------------------------------------------------
**available starting:** July 15, 2013

- AstroDrizzle now can process all STIS data without crashing.


DrizzlePac(astrodrizzle) v1.1.10dev(06-Feb-2013) in IRAFX
---------------------------------------------------------
**available starting:** May 6, 2013

- The output drizzle image header no longer contains references to D2IM arrays. This allows tweakreg to work with drizzled images as input where 2-D D2IM corrections were needed.
- Deprecated references to PyFITS .has_key() methods were also removed from the entire package, making it compatible with PyFITS 3.2.x and later.


DrizzlePac(astrodrizzle) v1.1.8dev(06-Feb-2013) in IRAFX
--------------------------------------------------------
**available starting:** Feb 11, 2013

- Fixed a bug in astrodrizzle which caused blot to raise an exception when using 'sinc' interpolation.
- Cleaned up the logic for writing out the results from the pixtopix, pixtosky, and skytopix tasks to avoid an Exception when a list of inputs are provided and no output file is specified.
- A new parameter was added to the tweakback task to allow a user to specify the value of WCSNAME when updating the FLT images with a new solution from a DRZ image header.
- Code in tweakback for updating the header with a new WCS will now automatically generate a unique WCSNAME if the there is a WCS solution in the FLT headers with the default or user-defined value of WCSNAME.


DrizzlePac(astrodrizzle) v1.1.7dev(18-Dec-2012) in IRAFX
--------------------------------------------------------
**available starting:** Feb 4, 2013

- Updated astrodrizzle to work with input images which do not have WCSNAME defined. This should make it easier to support non-HST input images in the future.
- cleared up confusion between flux parameters in imagefindpars and catalog inputs in tweakreg.
- turned of use of fluxes for trimming input source catalogs when no flux column can be found in input source catalogs


DrizzlePac(astrodrizzle) v1.1.7dev(18-Dec-2012) in IRAFX
--------------------------------------------------------
**available starting:** Dec 10, 2012

- Update tweakreg 2d histogram building mode to correctly find the peak when all the inputs match with the same offset (no spurious sources in either source catalog).
- Fixed a bug so that Ctrl-C does not cause an exception when used while tweakreg is running
- revised the source finding logic to ignore sources near the image edge, a change from how daofind works (daofind expands the image with blanks then fits anyway)
- created a new function to apply the nsigma separation criteria to (try to) eliminate duplicate entries for the same source from the source list. It turns out daofind does have problems with reporting some duplicate sources as well. This function does not work perfectly, but works to remove nearly all (if not all) duplicates in most cases.

DrizzlePac(astrodrizzle) v1.1.7dev(8-Jan-2012) in IRAFX
--------------------------------------------------------
**available starting:** Jan 14, 2013

- Bug fixed in updatehdr module to allow shiftfiles without RMS columns to work as inputs to manually apply shifts to headers of input images
- Revised astrodrizzle to update WCS of all input images BEFORE checking whether or not they are valid. This ensures that all files provided as input to astrodrizzle in the pipeline have the headers updated with the distortion model and new WCS.
- Images with NGOODPIX=0 now identified for WFC3 and WFPC2 inputs, so they can be ignored during astrodrizzle processing.
- Replaced 2d histogram building code originally written in Python with a C function that run about 4x faster.


DrizzlePac(astrodrizzle) v1.1.6dev(5-Dec-2012) in IRAFX
-------------------------------------------------------
**available starting:** Dec 10, 2012

- tweakreg v1.1.0 source finding algorithm now runs many times faster (no algorithmic changes). No changes have been made yet to speed up the 2d histogram source matching code.
- The 'pixtopix' task was updated to make the 'outimage' parameter optional by using the input image as the default. This required no API changes, but the help files were updated
- Very minor update to guard against MDRIZTAB being specified without any explicit path.
- Update astrodrizzle to correctly report the exposure time, exposure start, and exposure end for the single drizzle products, in addition to insuring the final drizzle values remain correct.
- astrodrizzle also includes initial changes to safeguard the C code from getting improperly cast values from the configObj(TEAL) input.

DrizzlePac(astrodrizzle) v1.1.5dev(23-Oct-2012) in IRAFX
--------------------------------------------------------
**available starting:** Oct 29, 2012

- Scaling of sky array for WFC3/IR IVM generation now correct
- template mask files for WFPC2 no longer generated so that WFPC2 data can now be processed using num_cores > 1 (parallel processing)
- interpretation of the 'group' parameter fixed to support a single integer, a comma-separated list of integers or a single 'sci,<n>' value. The values correspond to the FITS extension number of the extensions that should be combined. This fix may also speed up the initialization step as more direct use of pyfits was implemented for the interpretation of the 'group' parameter.

DrizzlePac(astrodrizzle) v1.1.1(31-Aug-2012) in HST Archive
-----------------------------------------------------------
**available starting:** Sept 26, 2012

The HST Archive and operational calibration pipeline started using this version on Sept 26, 2012.

DrizzlePac(astrodrizzle) v1.1.4dev(20-Sep-2012) in IRAFX
--------------------------------------------------------
**available starting:** Sept 24, 2012

- Bug fixed to allow use of final_wht_type=IVM for processing WFPC2 data
- Revised Initialization processing to speed it up by using more up-to-date, direct pyfits calls.

DrizzlePac(astrodrizzle) v1.1.3(7-Sep-2012) in IRAFX
-----------------------------------------------------
**available starting:** Sept 17, 2012

- Fixed the logic so that crclean images always get created regardless of the value of the 'clean' parameter.

DrizzlePac(astrodrizzle) v1.1.2(5-Sep-2012) in IRAFX
-----------------------------------------------------
**available starting:** Sept 10, 2012

- Remove the restriction of only being able to process images which have WCSNAME keyword as imposed by r15631. The removal of this restriction will now allow for processing of non-updated input files with updatewcs=False for cases where no distortion model exists for the data (as required by CADC).
- Added log statements reporting what sky value was actually used in the drizzle and blot steps

DrizzlePac(astrodrizzle) v1.1.1(30-Aug-2012) in IRAFX
-----------------------------------------------------
**available starting:** Sept 3, 2012

- Major revision to astrodrizzle allowing the option to process without writing out any intermediate products to disk. The intermediate products remain in memory requiring significantly more memory than usual. This improves the overall processing time by eliminating as much disk activity as possible as long as the OS does not start disk swapping due to lack of RAM.
- revised to turn off 'updatewcs' when coeffs=False(no) so that exposures with filter combinations not found in the IDCTAB will not cause an error

DrizzlePac(astrodrizzle) v1.0.7(21-Aug-2012) in IRAFX
-----------------------------------------------------
**available starting:** Aug 27, 2012

- Fixes problems with missing single_sci images.
- Static mask step revised to skip updates to static mask if all pixel data falls within a single histogram bin. This avoids problems with masking out entire images, which happens if low S/N SBC data is processed with static_mask=yes.


DrizzlePac(astrodrizzle) v1.0.6(14-Aug-2012) in IRAFX
-----------------------------------------------------
**available starting:** Aug 20, 2012

Use of IVM for final_wht now correct, as previous code used wrong inputs when IVM weighting was automatically generated by astrodrizzle.

DrizzlePac(astrodrizzle) v1.0.5(8-Aug-2012) in IRAFX
----------------------------------------------------
**available starting:** Aug 13, 2012

- Completely removed the use of the TIME arrays for weighting IR drizzle products so that the photometry for saturated sources in drizzled products now comes out correct.
- Corrected a problem with astrodrizzle which affected processing of WFPC2 data where CRPIX2 was not found when creating the output single sci image.

stsci_python v2.13 [Includes astrodrizzle v1.0.2(13-July-2012)]
---------------------------------------------------------------
**available starting:** Aug 3, 2012

The complete version of stsci_python can be downloaded from `our download page <http://www.stsci.edu/institute/software_hardware/pyraf/stsci_python/current/stsci-python-download>`_

- `stsci_python v2.13 Release Notes <http://www.stsci.edu/institute/software_hardware/pyraf/stsci_python/release-notes/releasenotes.2.13>`_

- `Old stsci_python release notes <http://www.stsci.edu/institute/software_hardware/pyraf/stsci_python/release-notes>`_


DrizzlePac(astrodrizzle) v1.0.1(20-June-2012)
---------------------------------------------
**Used in archive/pipeline starting:** July 10, 2012

Pipeline and archive started processing ACS data with this version.

DrizzlePac(astrodrizzle) v1.0.0(25-May-2012)
--------------------------------------------
**Used in archive/pipeline starting:** June 6, 2012

Pipeline and archive first started using astrodrizzle by processing WFC3 images.
