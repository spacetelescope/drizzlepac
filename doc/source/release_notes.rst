.. _release_notes:

**************************************
DrizzlePac Release Notes 
**************************************
The code for this package gets released through a number of methods: namely,
  - the use of the package for pipeline and archive processing of ACS and WFC3 data, 
  - SSB's semi-annual `public release of the stsci_python package <http://www.stsci.edu/institute/software_hardware/pyraf/stsci_python/current/stsci-python-download>`_, and 
  - a weekly beta release of the development version as part of the `IRAFX download <http://stsdas.stsci.edu/irafx/>`_.  
  
The following notes provide some details on what has been revised for each version in
reverse chronological order (most recent version at the top of the list).

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

