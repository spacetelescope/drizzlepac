.. _release_notes:

========================
DrizzlePac Release Notes
========================

The version of DrizzlePac can be identified using ::

  > python
  >>> import drizzlepac
  >>> drizzlepac.__version__

The following notes provide some details on what has been revised for each
version in reverse chronological order (most recent version at the top
of the list).  The number at the end of each item is the Github Pull Request (PR)
number of the code change for that issue.  These PRs can be viewed at:

    https://github.com/spacetelescope/drizzlepac/pulls

3.10.0 (14-Jul-2025)
====================

- Overwrote the MVM "alignment" configuration files with the SVM
  files to keep the information in-sync.  [#2028]

- Clarified the sigma value used to compute the threshold above which
  sources are detected for the segmentation catalog when using the
  Gaussian or RickerWavelet smoothing kernel.  The value has corrected
  in the output Segmentation catalogs and given greater visibility in
  the trailer log files.  [#2027]

- Updated the multiplicative values in the catalog configuration files
  which are used in conjunction with the computed image RMS to derive
  a threshold above which sources are detected. The Point catalog uses
  the variable "nsigma".  The Segmentation catalog uses the variables
  "segm_nsigma" and "rw2d_nsigma" when using the Gaussian or
  RickerWavelet smoothing kernels, respectively. [#2026]

- Updated the catalog configuration files to remove obsolete variables,
  "fwhm" and "TWEAK_THRESHOLD", from the sourcex and dao sections, respectively.
  The "TWEAK_FWHMPSF" variable/value now resides in the general section of the
  file as it applies to both catalogs.  Modified the catalog_utils.py module
  so the Segmentation catalog now reports the proper value for the Gaussian
  Filter FWHM which is used to smooth the total detection image. [#2024]

- Corrected the use of a string comparison to "asn" in the build_poller_table
  routine of the poller_utils.py module as these characters are
  a valid portion of the root of an ipppssoot filename (e.g., j6kasn01q).
  The comparison is now done against "_asn" when looking for association
  names in order to perform the proper actions. [#2019]

- Added parameter setting, sub_shape, to the IterativePSFPhotometry invocation
  to ensure a rectangular shape around the center of a star is defined when
  subtracting PSF models. [#2014]

- Resolved the issue of duplicate "ID"s in the rows of the Total Point catalog.
  For "point" source identification, looping is done over a list of weight masks,
  tp_masks, when invoking the "finder" algorithms. Tables are returned with a
  unique "id" number for each row/source in the table. When tp_masks > 1, the
  returned tables are stacked to generate a final table. The "id" number was not
  updated to reflect the stacking of the tables which created rows with the same
  "id".  This error did not cause the code to fail, but it did generate a garbled
  table. [#2007]

- Removed the extra column in the Point source identifcation table when using the
  DAOStarFinder or IRAFStarFinder utilities.  The extra column, daofind_mag, was
  added in Photutils 2.0. [#2006]

- Dropped support for Python v3.10 due to conflict with upgrading to
  Photutils v2.2.0. [#1987]

- Added documentation for the alignment logic and the selection of the SVM
  reference image. [#1967]

- Implemented an additional RMS determination for the background based
  upon the Median Absolute Deviation (MAD) algorithm. The MAD algorithm
  is now one of three ways the RMS is computed under the umbrella of
  computing the background of the input image.  The largest of the RMS
  values is ultimately used for further computation.  Removed the obsolete
  "bthresh" variable in instrument/detector "catalog" JSON files.  Updated
  the "bthresh" variable in only the ACS WFC "quality" JSON file to be "5.0",
  matching the other detector files. [#1978]

- Deprecated the TEAL GUI; TEAL is still used for loading configuration 
  files. [#1975]

- Fixed the crfactor designation for the WFPC2 detector (PC) which caused the
  the computation for rejecting catalog creation based on expected cosmic ray
  detections to fail ONLY for WFPC2.  Also, updated the WFPC2 cr_residual factor
  from 0.0 to 0.05 as it had never be set correctly.  Created a PyTest for
  WFPC2 SVM processing. [#1969]

- Updated the Pyproject.toml file to force use of Photutils v2.0.0 or greater.
  This update is in support of the change addressed by #1950. [#1966]

- Set non-positive catalog fluxes to nans to remove warnings for dividing by 
  zero and calculating the log of negative numbers. [#1959]

- Added a check to make sure that the pre-alignment WCS solutions from the astrometry 
  database are within a reasonable distance of the header target positions. [#1958]

- Removed deprecated parameter, edge_method, from the instantiation of a 
  Background2D.  The default for this value is now always equal to "pad"
  which was the setting in use in our code. [#1957]

- Removed python<3.13 restriction and remove some warnings. [#1936]

- build and test with latest supported version of Python [#1955]

- Added a third test using the size in pixels of the largest determined
  segment as a discriminant with regard to switching from a Gaussian to a
  Ricker-Wavelet kernel.  The kernel is used to convolve with the input image
  to increase the signal-to-noise ratio for the detection of faint sources. [#1953]

- Replaced deprecated class IntegratedGaussianPRF with CircularGaussianSigmaPRF.
  [#1950]

- Updated path for regression test results on artifactory. [#1933]

- Added new header keywords and match requirements for relative fitting. [#1860]

- Updated alignment parameters. [#1932]

- Implemented fixes to address uneven detection thresholds in the HAP catalogs
  due to bugs in the function, make_wht_masks, which intends to create weight
  masks covering the full drizzled output footprint. A somewhat related bug in
  the compute_threshold method associated with only the Segment catalog was also
  addressed.  The "scale factor" which causes the RMS computation to be too small
  was deleted.  The RMS computation for the Point and Segment catalogs is now the
  same. [#1939]


3.9.1 (30-Jan-2025)
===================

- Further updates done to address the deprecated Photutils functionality as the
  original changes did not produce results at least as good as the results
  generated by the previous Photutils functionality.  [#1934]
  

3.9.0 (16-Dec-2024)
===================

- **This version used by operations but does not generate HAP products (SVM/MVM).**

- Include a minimum RMS value for the SBC detector, as is done for the other
  detectors, as there seems to be a lot of noise in the source catalogs due to 
  a low detection threshold. [#1908]

- Force an exit with a return code, KEYWORD_UPDATE_PROBLEM, in try/exception block
  when invoking refine_product_headers in hapsequencer.py and hapmultisequencer.py.
  If the FITS header keywords are not properly updated, this can cause errors during
  CAOM ingest. [#1911]

- Introduce warnings for fits extensions with science data of all zeros, and ensure 
  data with zeros in all science extensions are not processed. [#998]

- Change to the algorithm which chooses which background determination algorithm to
  use for processing when working on the output source catalogs.  If the RMS from
  the Photutils Background2D is greater than the RMS from the astropy.stats
  sigma_clipped_stats, then the background algorithm is now set to use the
  sigma_clipped_stats. [#1904]

- Properly account for illuminated and non-illuminated pixels through the use
  of masks when determining the background algorithm to use when working on the
  output source catalogs. Use of the proper masks is particularly important for
  the ACS/SBC data where the pixel is in the illuminated zone, but its value may
  be zero. SBC is a MAMA photo-counting detector. (#1894)

- Modifications to support an upgrade to Photutils v1.13.0. Changes were made
  to accommodate new APIs, modified low-level functionality, and address columns
  of a table in get_cutouts() by name rather than position to ensure the correct
  data is acquired.  Support is now for versions of Photutils>=1.10.0.  [#1844]

- Added documentation describing regression tests. [#1881]

- Addressed additional issues related to numpy 2.0 scalar promotion. [#1875]

- Update to HDRTABLE for MVM products to include SVM rootname and SVM creation date. [#1846]

- Added python 3.12 to testing matrix for Jenkins and github actions. [#1843]

- ``manageInputCopies`` now copies successfully even if the original files were
  defined by full paths rather than being in the current working directory. [#1835]


3.8.0
=====

- Version not released; internal testing only. 

3.7.1.1 (1-Oct-2024)
====================

- Improved S_REGION using simplify-polygon, eorions, and dilation. [#1323] 


3.7.1 (12-Aug-2024)
===================
- Avoid applying the estimated cosmic ray vs real sources threshold for the
  ACS/SBC and WFC3/IR detectors. [#1858]

- Corrected the way the n1_exposure_time and tot_exposure_time values
  are computed as these values are used in the computation for rejecting
  catalog creation based on expected cosmic ray detections.  Generalized
  the crfactor dictionary for all detectors. Ensure if any catalog type
  is rejected, all the catalog types are rejected. [#1853]

- Modified the call to the hamming function in the deconvolve_utils.py module
  as SciPy deprecated the way window filtering functions can be invoked. These
  functions can no longer be imported from the scipy.signal namespace but need
  to be accessed via scipy.signal.windows. [#1848]

- Corrected the way that the number of constituent images are accumulated
  per pixel by ensuring each contributing pixel has a finite value and
  is not zero. [#1820]

- Within the HAP configuration files, increased the minimum number of matches
  for a successful "rscale" fit from 6 to 10, and removed "shift" as a fit geometry
  option. [#1823].

- Removed the use of a custom smoothing kernel based upon actual image
  data as a poorly determined kernel can ultimately cause poor source
  position determination.  The default kernel has been set to a
  Gaussian with default dimensions of 11 x 11 pixels. [#1805]

- Addressed bugs caught by SonarQube static code analysis.  Interface
  changes listed here: Added missing input data parameter to the create_output
  calls, Added missing log level to run function, Removed the deprecated
  parameter, dao_threshold, from astrometric_utils.py/extract_sources, removed
  "ivmlist" parameter from the interface of multiple functions in processInput.py
  as it is an output parameter (buildFileListOrig, buildFileList, checkMultipleFiles,
  and process_input), and addressed missing parameters in the calls to
  get_ci_info and get_ci_from_file.. [#1802]

- Exclude single filter images from the generation of the total detection
  image to minimize cosmic ray contamination, unless there are only single
  filter images in the visit. [#1797]

- Implemented a series of bug fixes for the segmentation catalog [#1793]
- Define the threshold image to be (nsigma * background_rms).
- Fixed bug in the generation of the threshold image - ensure the final
  threshold is built up properly by using the weight mask for the region
  in question.
- Pass the background image to detect_segments() so the convolved image can be
  background subtracted.
- For the detection of sources, background subtract the input image for both the
  Gaussian and RickerWavelet kernels.  Do not do any clipping on the background
  subtracted image.
- Update configuration files for the RickerWavelet2DKernel: source_box is now 6
  and rw2d_nsigma is now 3.
- Fixed a bug in the computation of the "biggest source".

- Created a new method, ricker_matched_kernel(), to generate the RickerWavelet2DKernel
  properly. Sigma is now provided, versus the FWHM, to the RickerWavelet2dKernel
  constructor, and the normalization is handled by the new method where the
  normalization causes the RickerWavelet core to match the Gaussian core.  [#1791]

- Added contributors guide to readthedocs. [#1787]

- Removed "tophat" as a kernel option, added warnings for "gaussian" and "lanczos3"
  that they may not be conserving flux. [#1786]

- Updated config json to exclude bad pixels in single WFC3/IR SVM processing. [#1783]

- Bug fix for mdriztab=True option in Astrodrizzle previously overwriting user inputs. [#1774]

- Reverted PR #1222 allowing pixels to be filled with available data where WHT=0. [#1767]

- Force the identified bad rows to be removed from the total (aka white light)
  source catalog before the corresponding bad segments are removed from the
  segmentation image. [#1771]

- Improved calculation of S_REGION using dialation and erosion. [#1762]

- Skycell added to flt(c) and drz(c) science headers for the pipeline and svm products. [#1729]


3.7.0 (02-Apr-2024)
===================

- Update project.toml file to specify numpy>=1.18,  <2.0 [#1743]

- Update project.toml file to specify python_requires>=3.10 [#1737]

- Github branch "master" renamed to main. [#1725]

- Clean up spacing in toml file to eliminate improper spacing to
  avoid decprecation warning [#1731]

- Clean up YAML diagram in of workflows area [#1728]

- Updated installation instructions and small text changes [#1727]

- Remove outdated references of Pyraf and change to Python [#1726]

- Fix to add "stregion" to the requirements-dev.txt file to fix the build
  error under Python 3.12. [#1714]

- Reorganized the readthedocs documentation with the help of various STScI
  staff. [#1717]

- Updates requirements-dev.txt to not install eggs that cause problems
  for the regression tests [#1721]

- Regression Testing: allow "dev" jobs to fail [#1718]

- Initial setup for Architectural Design Records used to keep track of top-level
  thinking behind the code. [#1697]


3.6.2 (27-Nov-2023)
===================

- At this time pin Astrocut to versions <=0.9 to avoid conflicts with urllib3
  package.  [#1689]

- Added functionality to allow the use of a two-column poller file. This is used
  to update the WFPC2 SVM aperture header keywords from the values in the poller
  file. [#1683]

- Removed the version restriction on matplotlib. [#1649]

- Forced a preferential order on the final selection of the WCS solution
  from the common pool of solutions among all input exposurea.  All input images
  need to have the same WCSNAME (same WCS solution) when performing pipeline
  alignment to avoid imprinting differences from one catalog to another on the
  final fit and destroying the relative alignment. [#1645, #1638]

- Redesigned the overall structure of the documentation, readthedocs, for the
  package. [#1620]

- Addressed a bug in the calculation of measurements for each detected source
  in the filter catalogs. The detection catalog, based upon the "total" image,
  is now used in the correct manner to define the source centroids and shape
  properties.  In addition, these properties are used to perform aperture
  photometry. [#1614]

- Updated the HAP drizzle parameters for WFPC2. The primary change includes
  changing skymethod='localmin' from the prior 'match' which did not work well
  for the overlapping chips. [#1617]

- Corrected reference catalog weights from being proportional to sigma to
  the proper 1/sigma**2. [#1616]

- Removed the use of the shadow mask as an initial step in addressing the WFPC2
  chip gaps [#1551]

- Fixed a bug in processing of the ``group`` argument due to which the code
  would crash when ``group`` would be an integer number or a list of numbers.
  Also, added support for specifying extensions as tuples of
  ``(extname, extver)``. [#1612]


3.6.1 (15-Jun-2023)
===================

- Fixed an incompatiblity in the ``minmed`` code for cosmic ray rejection
  with the ``numpy`` version ``>=1.25``. [#1573]

- Fixed projection cell identification in overlapping regions. [#1572]

- Force the version of matplotlib to be <= 3.6.3 as the newer versions of
  the library cause problems with the calcloud preview generation. [#1571]

3.6.0 (12-Jun-2023)
===================

- Modified the pyproject.toml file to ensure the tweakwcs version is greater
  than 0.8.2 as the issue of taking a very long time to compute the bounding
  polygon now defaults to an approximate method which is significantly faster.
  [#1565]

- Modified Projection Cell 0 declination coordinate of the center to be
  -89.999999999997 and the Projection Cell 2643 declination coordinate to
  be 89.999999999997 to shift the WCS CRVAL position slightly off the pole.
  [#1560]

- Modified the criteria for the rejection of catalogs based upon the cosmic
  ray criterion.  An empty catalog (n_sources=0) should not be rejected by the
  CR contamination.  Also, if a catalog is empty, it should not trigger the
  rejection of the other "type" of catalog (type=point vs segment). [#1559]

- For WFPC2 datasets which turn out to have no viable data to process and
  a manifest file has been requested, force an empty manifest file to be
  generated and issue the exit code NO_VIABLE_DATA (65). [#1550]

- Protect against writing the S_REGION keyword in intentionally empty DRZ/DRC
  files in ``processinput.process`` to avoid messy crash. [#1547]

- Fix a bug in ``processinput.buildFileListOrig`` due to which astrodrizzle
  might crash when ``updatewcs`` is set to ``True``. [#1549]

- Turn off use of ``verify_guiding()`` for WFPC2 images only as its use
  incorrectly recognizes diffraction spikes from saturated stars as evidence
  of loss of lock and flags those exposures as 'bad'. [#1511]

- Ensure processing of all IMAGETYP=EXT WFPC2 targets. [#1505]

- Properly identify neighbor Projection Cells which overlap input
  exposures. [#1503]

- Updates identify and remove any WFPC2 calibration exposures that
  cannot be processed during standard pipeline alignment and drizzling.
  The list of recognized calibration target names was updated to
  accommodate WFPC2 and to identify exposures to be skipped and deleted
  after converting the D0M images into FLT images. [#1514]

- Compute a default kernel for use with astrometric_utils.extract_sources()
  function when the kernel parameter is None.  The default kernel is based on
  the fwhm parameter of the same function. [#1519]

- Address many ReadTheDocs issues. [#1521 - #1529]

- Write the EXPNAME keyword to the ACS SVM and MVM headers to avoid errors
  and enforce consistency with WFC3. [#1530]

- Properly populate the S_REGION keyword with a closed polygon for the
  pipeline FLT/FLC images. [#1533]

- Compute the S_REGION values for pipeline drizzled products. [#1535]

- Ensure the DATE keyword is written to the primary header of all output
  drizzled products. The DATE represents the date the file was written.
  [#1537]

- Update to ensure the SVM FLT/FLC files all contain the S_REGION keyword
  and the value of the keyword is a closed polygon. [#1536]

3.5.1 (08-Feb-2023)
===================

- Turn on use of ``verify_guiding()`` to ignore exposures where guide star
  lock was lost and the stars are trailed. [#1443]

- Ensure when no sources are found and the variable thresh is zero, the
  ``verify_crthesh()`` properly indicates the catalog failed the CR threshold.
  [#1450]

- Added informational text when the catalog service fails (e.g., service cannot
  be reached or the request was somehow malformed) to make the default response
  more helpful. The request specification is also sent to the log, so the user
  can see what was actually requested. [#1451]

- Protect against there being no sources left to measure
  the properties after cleaning cosmic rays from the input
  in ``verify_guiding()``.
  [#1466]

- Check the SCI extension(s) of the output FLT/FLC and DRZ/DRC files.  If the active
  WCS solution is 'a priori', delete the following keywords if they are associated
  with the active WCS as they are residue from a previous 'a posteriori' solution:
  NMATCHES, RMS_RA/RMS_DEC, FITGEOM, and CRDER1/CRDER2. Ensure the WCSTYPE is based
  upon the active WCSNAME to clean up any confusion.
  [#1465]

- Protect against inability to find a FWHM due to a fitting problem. [#1467]

- Implement photometric equalization for standard pipeline processing
  (runastrodriz) of WFPC2 data. [#1471]

- Update required to the compute_2d_background() function of the astrometric_utils
  module to accommodate changes in the PhotUtils API. [#1480]

3.5.0 (10-Oct-2022)
====================

- Introduced a new ``apply_tweak()`` function as a replacement to the
  ``tweakback()``. ``apply_tweak()`` preserves the functionality of ``tweakback``
  with a re-designed API. Existing ``tweakback`` was deprecated. [#1372]

- Updated segmentation source catalog generation to use ICRS as input RADESYS
  when input images have an unsupported REFFRAME value (like OTHER or B1950). [#1423]

- Refactored code to work with changes in ``tweakwcs`` version 0.8.0. [#1430]

- Ignore non-CTE-corrected exposures when SVM or MVM products also include
  CTE-corrected exposures as inputs. [#1433]


3.4.3 (24-Aug-2022)
===================
This release includes includes updates for these features in addition to various bug fixes:
  - Initial support for aligning and creating SVM and MVM products for WFPC2 data
    based on unoptimized processing parameters
  - Python 3.10 support
  - Photutils 1.4.0 (and newer) support
  - Updated documentation on SVM processing and output mosaics

The list of specific changes for the significant issues includes:

- Fixed skycell size in pixels as quoted in the documentation. (#1387)
- Ensure Ramp filter data is not used for MVM processing (#1393)
- Added requested values and clarification text regarding photometry to the catalogs (#1390)
- Modified the docstring which defines the HAPLEVEL and its associated meaning (#1395)
- Modified the "exposure level" products to have a HAPLEVEL = 1 (#1398)
- Get full S_REGION outline (#1401)
- Update readthedocs for SVM catalog generation (#1400)
- Delete all reference catalogs during SVM processing (#1409)
- Update runastrodriz to work with WFPC2 data as singletons (#1412)
- Revert sky matching to use local sky minimization upon any error (#1411)
- Update SVM to support processing WFPC2 exposures (#1418)
- Add support for Python 3.10 (#1420)
- Add WFPC2 support for MVM processing (#1422)
- Support additional RADESYS options for input files (#1423)
- Ensure the gain variables are defined for all detectors (#1425)
- Essentially remove restriction on PhotUtils package version (#1426)


3.4.2 (27-May-2022)
===================
This release addresses a number of issues related to SVM and MVM processing.

- Reset tasknames to work with TEAL (#1285)
- Protect computations when photflam is equal to 0.0 (#1295)
- MVM: Define MVM-specific processing parameters for drizzling (#1277)
- Remove IPPPSSOO keyword from MVM product headers (again) (#1297)
- Fix problem with astropy 5.0 table interpretation (#1292)
- Statistics for SVM and MVM  (#1300)
- SVM: add/remove/update Astrodrizzle Parameter files (#1303)
- Explicitly update boolean column in ASN tables (#1307)
- Synchronize output WCS specifications for SVM processing (#1312)
- Smooth out determination of S_REGION vertices (#1315)
- Ensure units of catalog variables comply with Astropy (#1316)
- Apply default alignment fit parameters for zero exptime exposures (#1319)
- Fix bug caused by Astropy Tables being interpreted as QTables (#1320)
- Revise logic for when mask keywords are computed (#1323)
- Restrict version of Photutils to < 1.4.0. (#1326)
- Add MEANWHT and MEDWHT keywords to drizzle products (#1324, #1349)
- Add documentation describing mvm products and artifacts (#1322)
- Add release notes for 3.4.1final (#1328)
- Fix typo in ACS MVM header rules file (#1332)
- Update astropy min version to 5.0.4 (#1335)
- Avoid archiving duplicate WCS solutions in SVM processing (#1333)
- Update installation dependencies for fitsblender and skypac (#1354)
- Flag and ignore bad images based on detecting linear features (#1351)
- Improve algorithm for identifying and filtering large segments (#1357)
- Carry over IDCSCALE keyword when updating WCS to match Grism WCS (#1355)
- Ignore MVM layers with no overlapping exposures (#1360)
- Update crder units (#1362)
- This change addresses bugs associated with the big_segments attribute of the segmentation image (#1365)
- Update the WFC3 rules files (#1366)
- Only allow "verify_guiding" check for MVM processing (#1368)
- Fix the size of the HAPEXPNAME column in the HDRTAB of the MVM output DRZ/DRZ file (#1371)
- Pass along default WCSNAME (#1370)
- Re-design tweakback (#1372)
- Bugfix: point-cat-fxm files being left around (#1369)

3.4.1 (5-Apr-2022)
==================
This release addresses issues found in v3.4.0.  The most significant
issues were:

- Add documentation describing mvm products and artifacts (#1322)

- Revise logic for when mask keywords are computed (#1323)

- Restrict version of Photutils to < 1.4.0. (#1326)

- Add MEANWHT and MEDWHT keywords to drizzle products (#1324)

- Modify the units of the catalog variables so they are astropy-compatible (#1318)

- Smooth out determination of S_REGION vertices (#1315)

- Apply default alignment fit parameters for zero exptime exposures (#1319)

- fix for tasknames to once again work with TEAL (#1289)

- Revise code to properly support Astropy v5.0 (#1286 , #1290 , #1292, #1296, #1307)

- Protect computations in catalog generation when photflam is equal to 0.0 (#1295)

- Define MVM-specific and SVM-specific processing parameters for drizzling (#1277, #1303)

- Remove IPPPSSOO keyword from header of output SVM or MVM drizzle products (#1297)

- Insure correct statistics are reported in MVM headers (#1300)




3.4.0 (7-Mar-2022)
==================
This major release adds support for multi-visit mosaic (MVM) processing, in
addition to including numerous revisions to try to align more datasets
successfully to GAIA during pipeline and single-visit mosaic (SVM) processing.
Multi-visit mosaics (MVM) introduce the concept of SkyCells with new code added to define
them.  SkyCells are subarrays of pre-defined tangent planes spaced regularly
on the sky as standardized definitions of mosaics to be created
from all HST observations taken of each part of the sky.

New features added in this version include:

- Support for creating MVMs as generated
  by the 'drizzlepac/hapmultisequencer.py' module or using the
  new command-line task ``runmultihap``.

- Tools for generating cutouts of MVM products found in the
  ``drizzlepac/haputils/hapcut_utils.py`` module.

The most significant revisions and bug fixes that affect
output products of this version of the code include:

- Detect extension name from WFPC2 flat-field files. [#1193]

- Refactored the build system to be PEP-517 ad PEP-518 compliant. [#1244]

- Fixed a bug in the drizzle algorithm due to which input pixels with
  zero weights may still contribute to the output image. [#1222]

- Added Sphinx documentation describing tools used for working with
  MVM products. [#1144, #1150]

- Changed names of "ISO" columns in Segmentation catalog to be unique [#1155]

- Add WCS keyword values to catalog metadata [#1160]

- Enforced a minimum number of cross-matches for alignment to be 4 sources [#1187, #1218]

- Revised 2D background determination for smaller detectors to improve source
  detection during alignment. [#1187]

- Create empty catalogs when exposures are effectively blank. [#1199]

- Cut processing time from days to minutes for exposures of crowded fields of
  faint sources or fields dominated by a single large extended source.  [#1198]

- Report correct value of NMATCHES keyword for number of sources actually
  used in alignment fit to GAIA. [#1217]

- Prevent older distortion models from overriding new distortion models
  when performing a posteriori alignment to GAIA. [#1220]

- Add explicit dependency on spherical-geometry package. [#1232]

- Update how make_poller_files.py generates visit numbers. [#1221]

- Insure both FLT and FLC headers have same a posteriori fit keywords. [#1238]

- MVM: Make tool to quantify quality of GAIA alignment generic for general use. [#1241]

- Fix logic to not align grism data in standard pipeline. [#1243]

- Remove nictools as a dependency for this package. [#1245]

- RickerWavelet Kernel for SBC to separate crowded PSFS needs to have
  dimensions which are odd [#1246]

- Refine headers for filter and total products to allow keywords like IPPPSSOO and ASN_ID
  which only apply to single exposures
  (or data from the same ASN) to be removed from SVM filter and total drizzle products and
  from MVM layers drizzle products  [#1249]

- Remove logic from align that related to checking for alignment results in align.py
  when it was not necessary so that more data can successfully align to GAIA. [#1250]

- Add support for using astropy 5.0. [#1280]


3.3.1 (19-Nov-2021)
===================
This version provides bug fixes primarily
for the single-visit mosaic (SVM) processing.

- Insure a compatible version of photutils gets installed. [#1151]

- Improve handling of segmentation catalog generation for
  mostly or completely blank images. [#1152]

- Changed default floating point value in catalogs
  from -9999.9 to -9999.0.  [#1165]

- Avoid creating an empty manifest file when no images
  get drizzled by SVM processing, unless the visit was
  comprised solely of Grism/Prism data. [#1174, #1181]

- Update total catalog to only remove sources which were
  not measured successfully in any filter. [#1175]

- Fix the units of a few variables in the output Point and
  Segmentation catalogs [#1178]


3.3.0 (28-Sep-2021)
===================

This version includes all the functionality needed to generate
source catalogs, both point source and extended (segment) source
catalogs, during single-visit mosaic (SVM) processing.  In fact,

- Updated code to work with Python >= 3.7
- **GAIAeDR3** catalog now the initial catalog of choice for a posteriori alignment
  during standard pipeline processing, as well as for SVM/MVM processing.
- SVM/MVM processing will loop over catalogs, fit methods and fit geometries in
  looking for a successful fit, using the first successful fit it computes.

  - CATALOGS used: **GAIAeDR3**, **GSC242**, **2MASS** (in this order)
  - methods: relative, image-by-image
  - geometries: **rscale**, **rshift**, **shift** (each with different minimum cross-matches)

- SVM processing will always generate both point source and extended source catalogs, even
  if the catalogs contain no rows of sources and measurements.

  - point source catalog will be generated using TinyTim PSF-based detection
  - extended source (segment) catalog will only have sources larger
    than the PSF kernel deblended.
  - catalog columns will closely resemble the Hubble Legacy Archive (HLA) catalogs columns

- Grism/Prism exposures do not get aligned, but instead get the WCS correction from direct images
- Added logic to handle visits where there are only Grism/Prism exposures with no direct images
- ``S_REGION`` keyword:

  - added to FLT/FLC file headers
  - revised region computation to match closely the actual exposure footprint within mosaic

- Always runs ``updatewcs`` on input files to insure pipeline-default WCSs are always present

  - Add ``WCSNAME=OPUS`` if no ``IDCTAB`` WCS was created by ``updatewcs`` (``NGOODPIX=0``, ...).

These changes, and additional significant bug fixes, were implemented using
the following github PRs:

- Implemented deblending of segmentation source catalogs ONLY
  for sources larger than the PSF kernel. [#1131]

- Insure SVM processing always generates point-source and
  segmentation (extended) source catalogs, even if empty. [#1129]

- Implemented an efficient single-image identifier of possible
  cosmic-rays/defects, and applied it to help make image
  alignment more reliable.  [#1129]

- Update logic for fitting between source lists to minimize/eliminate
  use of fitting with less than 4 sources. [#1129]

- Implemented model PSF-based point-source identification for SVM
  point-source catalog generation. [#903, #971, #1127]

- Removed dependence on private photutils functions while enabling
  support for all photutils versions >= 1.0.0.
  [#1127, #1117, #1116, #1096]

- Set values for crowding, biggest source, and source
  fraction for use when to use the RickerWavelet kernel and
  when to deblend sources when identifying extended sources
  using segmentation for the segment catalog. [#1115]

- Implemented a more efficient algorithm based on Harris corner
  detection for computing the ``S_REGION`` keyword for pipeline
  and SVM drizzle products. [#1106]

- Fix a memory corruption issue in ``interpolate_bilinear()`` in
  ``cdrizzleblot.c`` which could result in segfault. [#1048]

- Fixed multiprocessing incompatibility with ``Python >= 3.8``. [#1101]

- Add support for environment variable switch, ``PIPELINE_RESET_IDCTAB``,
  to ``runastrodriz`` which will automatically reset ``IDCTAB``
  in FLT/FLC files if different from ``IDCTAB`` in RAW files.  [#1046]

- Update documentation based on revisions to the code.
  [#941, #947, #953]

- Update default astrometry catalogs for alignment to try alignment to
  the ``GAIA eDR3`` catalog first. [#986, #1012]

- Enable user epoch selection when a user requests a GAIA catalog from
  the astrometry catalog web service. [#1006]

- Insure that ``HDRNAME`` is always valid for updated WCS solutions. [#966]

- Revised ``S_REGION`` keyword value to reflect actual outline of chips in
  drizzle products.  [#951]

- Sky Subtraction step will automatically downgrade from ``match`` to ``localmin``,
  and from ``globalmin+match`` to ``globalmin`` when sky matching runs into an
  Exception. [# 1007]

- Changed to insure that ``EXTNAME`` and ``EXTVER`` are always removed from
  simple FITS drizzle product headers. [#954]

- Changed to insure that all the distortion keywords (e.g., ``TDD*``, ``D2IM*``,...)
  are removed from from the output drizzle product headers [#954].

- Set a common active WCS for direct as well as corresponding Grism/Prism images [#929, #946]

- Fix a bug in ``tweakback`` that may cause incorrect "updated" WCS to be
  picked up from the drizzled image. [#913]

- Added ``DRIZPARS`` keyword to final output drizzle product primary header
  to document the name of the associated trailer file. [#934, #1078]

In addition, numerous changes were made to insure this code stayed
compatible with numpy versions > 1.20 and astropy versions > 4.1.

Updates to the ``STWCS`` package version >= 1.6.0 also translated to
the following changes to the Drizzlepac processing:
- Insure HDRNAME keyword is never empty
- Remove duplicate headerlet extensions when running updatewcs
- Compute new a priori WCS solutions for new IDCTAB not already in astrometry database

***API Changes:***

**imageObject.py:**
  - **class imageObject**: Added parameter ``output`` to enable determination
    of rootname for use in processing of each detector.

**adrizzle.py:**
  - **drizSeparate**: Added optional parameter ``logfile`` for specifying
    what file to use for log messages.
  - **drizFinal**: Added optional parameter ``logfile`` for specifying
    what file to use for log messages.

**wcs_functions.py:**
  - Removed ``hdulist`` as parameter from ``get_hstwcs``.

**haputils/analyze.py:**
  - **analyze_data**: Added parameter ``type`` to customize logic for SVM
    processing.

**haputils/astrometric_utils.py:**
  - **retrieve_observation**:  Added parameter ``product_type`` to allow for selection of
    type of products to be returned; pipeline, HAP, or both.

**haputils/make_poller_files.py:**
  - New function ``generate_poller_file`` added to create inputs for SVM processing
    from files on disk.

**haputils/processing_utils.py:**
  - New function ``find_footprint`` added to determine corners of all chips
    in an image for computation of ``S_REGION`` keyword.
  - New function ``interpret_sregion`` added to convert ``S_REGION`` keyword
    value into list of RA/Dec points for visualization.


3.2.1 (16-Feb-2021)
===================

- Fix problems with testing code for this package [#940]


3.2.0 (7-Dec-2020)
==================

This version provides the first operational implementation of the single-visit
mosaic processing used to create the single-visit mosaics products.

- revise naming convention for the StaticMask file so that it has a
  dataset-specific name instead of a generic common name. [#876]

- Update ``runastrodriz`` to work under Windows while adding documentation
  to tell the user to run with ``num_cores`` set to 1.  [#794]

- Fixed a bug in ``TweakReg`` due to which ``TweakReg`` would crash when
  ``updatehdr`` was set to `False`. [#801]


3.1.8 (11-Aug-2020)
===================

A number of changes have been implemented to either correct problems or
improve the processed results.  The most significant of the changes are:

  - rscale only used for alignment.
  - a minimum of 6 sources now gets used for alignment
  - no proper motions used in astrometric (GAIA) catalog when attempting a posteriori fitting
  - chip-to-chip alignment errors were corrected


In addition to a few dozen bug fixes, the following updates to the algorithms
were also implemented.

- Simplified the logic in ``tweakreg`` for deciding how to archive primary WCS
  resulting in a reduction of duplicate WCSes in image headers. [#715]

- Added polynomial look-up table distortion keywords to the list of distortion
  keywords used by ``outputimage.deleteDistortionKeywords`` so that
  distortions can be removed from ACS images that use ``NPOLFILE``.
  This now allows removal of alternate WCS from blotted image headers. [#709]

- Added ``rules_file`` parameter to AstroDrizzle to enable use of custom
  files in pipeline processing. [#674]

- Only apply solutions from the astrometry database which were non-aposteriori
  WCS solutions as the PRIMARY WCS.  This allows the pipeline to compare the
  true apriori WCS solutions (e.g., GSC or HSC WCSs) to aposteriori solutions
  computed using the latest distortion-models and alignment algorithms being
  used at the time of processing. [#669]

- Verification using a similarity index gets reported in the trailer file and
  does not get used as a Pass/Fail criteria for alignment.  [#619]

- If verification fails for either pipeline-default or apriori solution, reset
  cosmic-ray(CR) flag (4096) in DQ arrays.  This will allow subsequent attempt to
  align the images to not be impacted by potentially mis-identified CRs that most
  likely blanked out real sources in the field.  As a result, the image alignment
  process became more robust when computing the aposteriori alignment.  [#614]

- Fix a crash in ``tweakreg`` when finding sources in very large images
  due to a bug in ``scipy.signal.convolve2d``. [#670]

- Fix a bug in ``tweakreg`` due to which the number of matched sources needed to be
  *strictly* greater than ``minobj``. Now the minimum number of matched sources
  maust be *at least* equal or greater than ``minobj``. [#604]

- Fix a crash in ``tweakreg`` when ``2dhist`` is enabled and ``numpy``
  version is ``1.18.1`` and later. [#583, #587]

- Update calibrated (FLC/FLT) files with RMS and NMATCH keywords when it successfully
  aligns the data to GAIA using the a posteriori fit.  Headerlet files for this fit
  which already have these keywords are now retained and provided as the final output
  headerlets as well.  [#555]

- Insure HDRNAME keyword gets added to successfully aligned FLC/FLT files. [#580]

- Fix problem with 'tweakback' task when trying to work with updated WCS names. [#551]

- Fix problems found in processing data with NGOODPIX==0, DRC files not getting
  generated for singletons, alignment trying to use a source too near the chip edge,
  catch the case were all inputs have zero exposure time, lazily remove alignment
  sub-directories, fixed a bug in overlap computation that showed up in oblong mosaics,
  recast an input to histogram2d as int,  defined default values for tables when no
  sources were found. [#593]

- Updated to be compatible with tweakwcs v0.6.0 to correct chip-to-chip alignment issues
  in aposteriori WCS solutions. [#596]

- Correctly define output drizzle product filename during pipeline processing
  for exposures with 'drz' in the rootname. [#523]

- Implement multiple levels of verification for the drizzle products generated
  during pipeline processing (using runastrodriz); including overlapp difference
  computations [#520], and magnitude correlation [#512].

- Replace alignimages module with O-O based align [#512]

- Fix problem with NaNs when looking for sources to use for aligning images [#512]

- Fixed code that selected the brightest sources to use for alignment allowing
  alignment to work (more often) for images with saturated sources. [#512]

- Use logic for defining the PSF extracted from the images to shrink it in each
  axis by one-half for images of crowded fields to allow for more sources to be
  extracted by daofind-like algorithm. This enables source finding and alignment
  to work more reliably on crowded field images. [#512]

- Insure all input files, especially those with zero exposure time or grism
  images, get updated with the latest pipeline calibration for the distortion. [ #495]

This version also relies on updates in the following packages to get correctly
aligned and combined images with correctly specified WCS keywords:

- TWEAKWCS 0.6.4:  This version corrects problems with the chip-to-chip separation
  that arose when applying a single fit solution to the entire observation.

- STWCS 1.5.4:  This version implements a couple of fixes to insure that use of
  headerlets defines the full correct set of keywords from the headerlet for
  the PRIMARY WCS in the science exposure without introducing multiple copies of
  some keywords.

- Numpy 1.18: Changes in numpy data type definitions affected some of the code used
  for computing the offset between images when performing aposteriori alignment
  during pipeline processing and when running the 'tweakreg' task.


3.1.3 (5-Dec-2019)
==================

- Fixed a bug in the ``updatehdr.update_from_shiftfile()`` function that would
  crash while reading shift files. [#448]

- Migration of the HAP portion of the package to an object-oriented
  implemenation. [#427]

- Added support for providing HSTWCS object as input to 'final_refimage'
  or 'single_refimage' parameter. [#426]

- Implementation of grid definition interface to support returning SkyCell
  objects that overlap a mosaic footprint. [#425]

- Complete rewrite of ``runastrodriz`` for pipeline processing to include
  multi-level verification of alignment.  [#440]

3.0.2 (15-Jul-2019)
====================

- Removed deprecated parameter ``coords`` from the parameter list of
  ``pixtopix.tran()`` function. [#406]

- Modified the behavior of the ``verbose`` parameter in ``pixtopix.tran()``
  to not print coordinates when not run as a script and when ``output``
  is `None`. [#406]

- Fixed a compatibility issue in ``tweakutils`` that would result in crash in
  ``skytopix`` when converting coordinates in ``hms`` format. [#385]

- Fixed a bug in the ``astrodrizzle.sky`` module due to which sky matching
  fails with "Keyword 'MDRIZSKY' not found" error when some of the
  input images do not overlap at all with the other images. [#380]

- Fixed a bug in the ``util.WithLogging`` decorator due to which incorrect
  log file was reported when user-supplied log file name does not have ``.log``
  extension. [#365]

- Fixed a bug introduced in #364 returning in ``finally`` block. [#365]

- Improved ``util.WithLogging`` decorator to handle functions that return
  values. [#364]

- Fixed a bug in the automatic computation of the IVM weights when IVM
  was not provided by the user. [#320]

- Fixed a bug in the 2D histogram code used for estimating shifts for
  catalog pre-matching. This may result in better matching. [#286]

- Now ``tolerance`` (in ``tweakreg``) is no longer ignored when ``use2dhist``
  is enabled. [#286]

- Fixed VS compiler errors with pointer artithmetic on void pointers. [#273]

- Fix logic so that code no longer tries to update headers when no valid fit
  could be determined. [#241]

- Fixed a bug in the computation of interpolated large scale flat field
  for STIS data. The bug was inconsequential in practice.
  Removed the dependency on ``stsci.imagemanip`` package. [#227]

- Removed the dependency on ``stsci.ndimage`` (using ``scipy`` routines
  instead). [#225]

- Added ``'Advanced Pipeline Products'`` alignment code to ``drizzlepac``
  package. Enhance ``runastrodriz`` to compute and apply absolute astrometric
  corrections to GAIA (or related) frame to images where possible.
  [#200, #213, #216, #223, #234, #235, #244, #248, #249, #250, #251,
  #259, #260, #268, #271, #283, #294, #302]

- Add computation and reporting of the fit's
  `Root-Mean-Square Error (RMSE) <https://en.wikipedia.org/wiki/Root-mean-square_deviation>`_
  and `Mean Absolute Error (MAE) <https://en.wikipedia.org/wiki/Mean_absolute_error>`_.
  [#210]

- Replaced the use of ``WCS._naxis1`` and ``WCS._naxis2`` with
  ``WCS.pixel_shape`` [#207]

- Removed support for Python 2. Only versions >= 3.5 are supported. [#207]

- Use a more numerically stable ``numpy.linalg.inv`` instead of own matrix
  inversion. [#205]

- The intermediate fit match catalog, with the name ``_catalog_fit.match``
  generated by ``tweakreg`` now has correct RA and DEC values for the sources
  after applying the fit. [#200, #202]

- Simplify logic for determining the chip ID for each source. [#200]


2.2.6 (02-Nov-2018)
===================

- Fix a bug that results in ``tweakreg`` crashing when no sources are found
  with user-specified source-finding parameters and when ``tweakreg`` then
  attempts to find sources using default parameters. [#181]

- Updated unit_tests to use original inputs, rather than updated inputs used by
  nightly regression tests.

- Fix ``numpy`` "floating" deprecation warnings. [#175]

- Fix incorrect units in CR-cleaned images created by ``astrodrizzle``. Now
  CR-cleaned images should have the same units as input images. [#190]


2.2.5 (14-Aug-2018)
===================

- Changed the color scheme of the ``hist2d`` plots to ``viridis``. [#167]

- Refactored test suite

- ``sdist`` now packages C extension source code


2.2.4 (28-June-2018)
====================

- Replace ``pyregion`` with ``stregion``


2.2.3 (13-June-2018)
====================

- Updated links in the documentation to point to latest
  ``drizzlepac`` website and online API documentation.

- Code cleanup.

- Updated C code to be more compatible with latest numpy releases in order
  to reduce numerous compile warnings.

- Updated documentation to eliminate (at this moment) all sphinx documentation
  generation warnings.

- Moved ``'release_notes.rst'`` to ``'CHANGELOG.rst'`` in the top-level
  directory.

- Improved setup to allow documentation build. See
  `drizzlepac PR #142 <https://github.com/spacetelescope/drizzlepac/pull/142>`_
  and `Issue #129 <https://github.com/spacetelescope/drizzlepac/issues/129>`_
  for more details.

- Fixed a bug in a print statement in the create median step due to which
  background values for input images used in this step were not printed.

- Fixed a bug due to which ``TweakReg`` may have effectively ignored
  ``verbose`` setting.

- Fixed a bug in ``drizzlepac.util.WithLogging`` due to which ``astrodrizzle``
  would throw an error trying when to raise another error.
  See `Issue #157 <https://github.com/spacetelescope/drizzlepac/issues/157>`_
  for more details.


2.2.2 (18-April-2018)
=====================

- Fixed a bug in ``TweakReg`` introduced in ``v2.2.0`` due to which, when
  ``TweakReg`` is run from the interpreter, the code may crash when trying to
  interpret input files.


2.2.1 (12-April-2018)
=====================

- Fixed problems with processing WFPC2 data provided by the archive.  User will
  need to make sure they run ``updatewcs`` on all input WFPC2 data before
  combining them with ``astrodrizzle``.


2.2.0 (11-April-2018)
=====================

- Implemented a major refactor of the project directory structure. Building no
  longer requires ``d2to1`` or ``stsci.distutils``. Drizzlepac's release
  information (i.e. version, build date, etc) is now handled by ``relic``.
  See https://github.com/spacetelescope/relic

- Added basic support for compiling Drizzlepac's C extensions under Windows.

- Documentation is now generated during the build process. This ensures the
  end-user always has access to documentation that applies to the version of
  ``drizzlepac`` being used.

- Swapped the effect of setting ``configobj`` to `None` or ``'defaults'`` in
  ``AstroDrizzle`` and ``TweakReg``. When calling one of these tasks with
  ``configobj`` parameter set to `None`, values for the
  not-explicitly-specified parameters should be set to the default values
  for the task. When ``configobj`` is set to ``'defaults'``
  not-explicitly-specified parameters will be loaded from the
  ``~/.teal/astrodrizzle.cfg`` or ``~/.teal/tweakreg.cfg`` files that store
  latest used settings (or from matching configuration files in the current
  directory). See https://github.com/spacetelescope/drizzlepac/pull/115
  for more details.


2.1.22 (15-March-2018)
======================

- Changed the definition of Megabyte used to describe the size of the buffer
  for create median step (``combine_bufsize``). Previously a mixed
  (base-2 and base-10) definition was used with 1MB = 1000x1024B = 1024000B.
  Now 1MB is defined in base-2 (MiB) as 1MB = 1024x1024B = 1048576B.

- Redesigned the logic in ``createMedian`` step used to split large
  ``single_sci`` images into smaller chunks: new logic is more straightforward
  and fixes errors in the old algorithm that resulted in crashes or
  unnecessarily small chunk sizes that slowed down ``createMedian`` step.

- Due to the above mentioned redesign in the logic for splitting large images
  into smaller chunks, now ``overlap`` can be set to 0 if so desired in the
  ``minmed`` combine type. Also, it is automatically ignored (set to 0) for all
  non-``minmed`` combine types. This will result in additional speed-up in the
  Create Median step.

- Both ``AstroDrizzle()`` and ``TweakReg()`` now can be called with
  ``configobj`` parameter set to ``'defaults'`` in order to indicate that
  values for the not-explicitly-specified parameters should be set to
  the default values for the task instead of being loaded from the
  ``~/.teal/astrodrizzle.cfg`` or ``~/.teal/tweakreg.cfg`` files that store
  latest used settings.

- Updated documentation.


2.1.21 (12-January-2018)
========================

- Restore recording of correct ``EXPTIME`` value in the headers of
  single drizzled ("single_sci") images. See
  https://github.com/spacetelescope/drizzlepac/issues/93 for more details.

- Fixed a bug in ``drizzlepac`` due to which user provided ``combine_lthresh`` or
  ``combine_hthresh`` in the ``CREATE MEDIAN IMAGE`` step were not converted
  correctly to electrons (processing unit). This bug affected processing of
  WFPC2, STIS, NICMOS, and WFC3 data. See
  https://github.com/spacetelescope/drizzlepac/issues/94 for more details.

- Modified print format so that scales, skew and rotations are printed with
  10 significant digits while shifts are printed with 4 digits after the
  decimal point.


2.1.20 (07-October-2017)
========================

- Fixed a bug in expanding reference catalog in ``TweakReg`` that would result
  in the code crashing.
  See https://github.com/spacetelescope/drizzlepac/pull/87 for more details.

- Fixed a bug due to which user catalog fluxes would be interpreted as
  magnitudes when ``fluxunits`` was set to ``'cps'``.
  See https://github.com/spacetelescope/drizzlepac/pull/88 for more details.

- Fixed a bug due to which user-supplied flux limits were ignored for
  the reference catalog.
  See https://github.com/spacetelescope/drizzlepac/pull/89 for more details.


2.1.19 (29-September-2017)
==========================

- Fixed a bug in computing optimal order of expanding reference catalog that
  resulted in code crashes.
  See https://github.com/spacetelescope/drizzlepac/pull/86 for more details.


2.1.18 (05-September-2017)
==========================

- Fixed ``astrodrizzle`` lowers the case of the path of output images issue.
  See https://github.com/spacetelescope/drizzlepac/issues/79 for more
  details.

- Fixed ``tweakreg`` ignores user-specified units of image catalogs (provided
  through the ``refcat`` parameter) issue. See https://github.com/spacetelescope/drizzlepac/issues/81 for more details.

- Corrected a message printed by tweakreg about used WCS for alignment. Also
  improved documentation for the ``refimage`` parameter.


2.1.17 (13-June-2017)
=====================

- ``drizzlepac.adrizzle`` updated to work with numpy >=1.12 when they implemented
  more strict array conversion rules for math. Any input which still has INT
  format will be converted to a float before any operations are performed, explicitly
  implementing what was an automatic operation prior to numpy 1.12.


2.1.16 (05-June-2017)
=====================

- Fixed a bug introduced in release v2.1.15 in the logic for merging WCS due to
  which custom WCS scale was being ignored.


2.1.15 (26-May-2017)
====================

- ``fits.io`` operations will no longer use memory mapping in order
  to reduce the number of file handles used when running either
  ``astrodrizzle`` or ``tweakreg``. See
  `issue #39 <https://github.com/spacetelescope/drizzlepac/issues/39>`_
  for more details.

- Fixed bugs and improved the logic for merging WCS that is used to define
  ``astrodrizzle``'s output WCS.

- Added ``crpix1`` and ``crpix2`` parameters to custom WCS.


2.1.14 (28-Apr-2017)
====================

- Supressed info messages related inconsistent WCS - see
  `issue #60 <https://github.com/spacetelescope/drizzlepac/pull/60>`_ and
  `stwcs issue #25 <https://github.com/spacetelescope/stwcs/issues/25>`_
  for more details.


2.1.13 (11-Apr-2017)
====================

- Fixed a bug due to which sky background was subtracted by ``astrodrizzle``
  from the images even though ``skysub`` was set to `False` when
  ``MDRIZSKY`` was already present in input images' headers.


2.1.12 (04-Apr-2017)
====================

- ``astrodrizzle`` now will run ``updatewcs()`` on newly created images
  when necessary, e.g., after converting WAVERED FITS to MEF format
  (``*c0f.fits`` to ``*_c0h.fits``) or after unpacking multi-imset STIS
  ``_flt`` files. See
  `PR #56 <https://github.com/spacetelescope/drizzlepac/pull/56>`_ for
  more details.

- Fixed a bug that was preventing processing STIS image data.

- Fixed a bug in reading user input (see
  `issue #51 <https://github.com/spacetelescope/drizzlepac/issues/51>`_).


2.1.11 (24-Mar-2017)
====================

Bug fix release (a bug was introduced in v2.1.10).


2.1.10 (23-Mar-2017)
====================

Some of the changes introduced in release v2.1.9 were not backward compatible.
This release makes those changes backward compatible.


2.1.9 (22-Mar-2017)
===================

Compatibility improvements with Python 3 and other STScI software packages.


2.1.8 (08-Feb-2017)
===================

- Drizzlepac code will no longer attempt to delete "original" (WCS key 'O')
  resulting in a decreased number of warnings
  (see `issue #35 <https://github.com/spacetelescope/drizzlepac/issues/34>`_ ).

- Negative values are now zeroed in the 'minmed' step before attempting to
  estimate Poisson errors
  (see `issue #22 <https://github.com/spacetelescope/drizzlepac/issues/22>`_).

- Fixed a bug in ``tweakreg`` due to incorrect matrix inversion.

- Improved compatibility with `astropy.io.fits` ('clobber' parameter) and
  `numpy` which has reduced the number of deprecation warnings).

- Existing static masks in the working directory are now overwritten and not
  simply re-used (see
  `issue #23 <https://github.com/spacetelescope/drizzlepac/issues/23>`_).

- Corrected formula for :math:`\sigma` computation in the "create median" step
  to convert background to electrons before computations. This bug was
  producing incorrect :math:`\sigma` for instruments whose gain was different
  from one.

- Improved ``astrodrizzle`` documentation for ``combine_type`` parameter which
  now also documents the formula for :math:`\sigma` computation
  when ``combine_type`` parameter is set to ``'minmed'``.


2.1.6 and 2.1.7rc (15-Aug-2016)
===============================

Package maintenance release.


2.1.5 (09-Aug-2016)
===================

Technical re-release of ``v2.1.4``.


2.1.4 (01-Jul-2016)
===================

The following bug fixes have been implemented:

- ``tweakreg`` crashes when run with a single input image and
  a reference catalog.

- Fixes an issue due to which ``tweakreg``, when updating image headers,
  would not add '-SIP' suffix to CTYPE


2.1.3 (16-Mar-2016)
===================

- Improved ASN input file handling.

- ``astrodrizzle`` does not delete ``d2imfile`` anylonger allowing multiple
  runs of ``updatewcs`` on the same WFPC2 image, see
  `Ticket 1244 <https://trac.stsci.edu/ssb/stsci_python/ticket/1244>`_
  for more details.

- Allow exclusion regions in ``tweakreg`` to be in a different directory and
  allow relative path in exclusion region file name.

- Improved handling of empty input image lists.

- ``tweakreg`` bug fix: use absolute value of polygon area.



2.1.2 (12-Jan-2016)
===================

- ``runastrodriz`` moved to ``drizzlepac`` from ``acstools`` and
  ``wfc3tools`` packages.

- Improved logic for duplicate input detection.

- Improved logic for handling custom WCS parameters in ``astrodrizzle``.

- Compatibility improvements with Python 3.


2.1.1
=====

**Available under SSBX/IRAFX starting:** Nov 17, 2015

This release includes the following bug fixes:

- Resolved order of operation problems when processing WFPC2 data with
  DGEOFILEs.

- The conversion of the WFPC2 ``DGEOFILE`` into ``D2IMFILE`` is now
  incorporated into ``STWCS`` v1.2.3 (r47112, r47113, r47114) rather than a
  part of ``astrodrizzle``. This requires users to run updatewcs first, then
  ``astrodrizzle``/``tweakreg`` will work with that WFPC2 data seamlessly
  (as if they were ACS or WFC3 data).

- Compatibility improvements with Python 3.


2.1.0
=====

**Available under SSBX/IRAFX starting:** Nov 2, 2015

This version builds upon the major set of changes implemented in v2.0.0 by not
only fixing some bugs, but also cleaning up/changing/revising some APIs and
docstrings. The complete list of changes includes:

- [API Change] The 'updatewcs' parameter was removed from both the
  ``astrodrizzle`` and ``tweakreg`` interactive TEAL interfaces.
  The 'updatewcs' parameter can still be used with the Python interface for
  both the ``astrodrizzle``. ``astrodrizzle``() and ``tweakreg``. Call the
  ``stwcs.updatewcs.updatewcs()`` function separately before running
  ``astrodrizzle`` or ``tweakreg``.

- [API Change] The stand-alone interface for the blot routine
  (``ablot.blot()``) has been revised to work seamlessly with
  astrodrizzle-generated products while being more obvious how to call it
  correctly. The help file for this task was also heavily revised to document
  all the input parameters and to provide an example of how to use the task.

- [API Change] Coordinate transformation task
  (``pixtopix``/``pixtosky``/``skytopix``) interfaces changed to be more
  consistent, yet remain backward-compatible for now.

- Both ``astrodrizzle`` and ``tweakreg`` now return an output CD matrix which
  has identical cross-terms indicating the same scale and orientation in each
  axis (an orthogonal CD matrix). This relies on a revision to the
  ``stwcs.distortion.utils.output_wcs()`` function.

- The user interfaces to all 3 coordinate transformation tasks now use
  'coordfile' as the input file of coordinates to transform. The use
  of 'coords' has been deprecated, but still can be used if needed. However,
  use of 'coordfile' will always override any input provided simultaneously
  with 'coords' parameter.  Help files have been updated to document this as
  clearly as possible for users.

- User-provided list of input catalogs no longer needs to be matched exactly
  with input files. As long as all input images are included in input catalog
  list in any order, ``tweakreg`` will apply the correct catalog to the
  correct file.

- ``tweakreg`` has been updated to correctly and fully apply source selection
  criteria for both input source catalogs and reference source catalogs based
  on ``fluxmin``, ``fluxmax`` and ``nbright`` for each.

- All use of keyword deletion has been updated in ``drizzlepac`` (and
  ``fitsblender``) to avoid warnings from astropy.

- All 3 coordinate transformation tasks rely on the input of valid WCS
  information for the calculations. These tasks now warn the user when it
  could not find a valid WCS and instead defaulted to using a unity WCS, so
  that the user can understand what input needs to be checked/revised to get
  the correct results.

- Exclusion/inclusion region files that can be used with ``tweakreg`` can now
  be specified in image coordinates and sky coordinates and will only support
  files written out using DS9-compatible format.

- The filename for 'final_refimage' in ``astrodrizzle`` and 'refimage' in
  ``tweakreg`` can now be specified with OR without an extension, such as
  '[sci,1]' or '[0]'.  If no extension is specified, it will automatically
  look for the first extension with a valid HSTWCS and use that. This makes
  the use of this parameter in both place consistent and more general than
  before.

- The reported fit as written out to a file has been slightly modified to
  report more appropriate numbers of significant digits for the results.

- Use of astrolib.coords was removed from ``drizzlepac`` and replaced by use
  of astropy functions instead. This eliminated one more obsolete dependency
  in our software.

- Code was revised to rely entirely on ``astropy.wcs`` instead of stand-alone
  pywcs.

- Code was revised to rely entirely on ``astropy.io.fits`` instead of
  stand-alone pyfits.

- Added ``photeq`` task to account for inverse sensitivity variations across
  detector chips and/or epochs.

- WFPC2 data from the archive with ``DGEOFILE`` reference files will now need
  to be processed using ``stwcs.updatewcs`` before running them through
  ``astrodrizzle`` or ``tweakreg``.  This update converts the obsolete,
  unsupported ``DGEOFILE`` correction for the WFPC2 data into a ``D2IMFILE``
  specific for each WFPC2 observation, then uses that to convert the WCS based
  on the new conventions used for ACS and WFC3.

This set of changes represents the last major development effort for
``DrizzlePac`` in support of HST.  Support of this code will continue
throughout the lifetime of HST, but will be limited primarily to bug fixes
to keep the code viable as Python libraries used by ``DrizzlePac`` continue
to develop and evolve with the language.


2.0.0
=====

** Available under SSBX/IRAFX starting:** Aug 4, 2014

This version encompasses a large number of updates and revisions to the
``DrizzlePac`` code, including the addition of new tasks and several parameter
name changes. The scope of these changes indicates the level of effort that
went into improving the ``DrizzlePac`` code to make it easier and more
productive for users. The most significant updates to the ``DrizzlePac``
code include:

- The Python code has been updated to work identically (without change) under
  both Python 2.7 and Python 3.x.

- Implementing sky matching, a new algorithm for matching the sky across a set
  of images being combined by ``astrodrizzle``.

- Updating ``tweakreg`` to now align full mosaics where some images may not
  overlap others in the mosaic.

- Added the option to write out single drizzle step images as compressed images
  (to save disk space for large mosaics, and I/O time for single drizzle step).

- Improved ``tweakreg`` residual plots visually while allowing them to be
  written out automatically when ``tweakreg`` gets run in non-interactive mode.

- Renamed parameters in ``tweakreg`` and imagefind to eliminate name clashes.

- Added option to select sources based on sharpness/roundness when ``tweakreg``
  searches for sources.

- Added support for exclusion and inclusion regions arbitrary shape/size when
  ``tweakreg`` searches for sources.

- Added a full set of source detection parameters for reference image to
  support multi-instrument alignment in ``tweakreg``.

- Added support for new (simpler, more robust) ACS calibration of
  time-dependent distortion.

- A full 6-parameter general linear fit can now be performed using
  ``tweakreg``, in addition to shift and rscale.

- Cleaned up logic for sky-subtraction: user can now turn off sky-subtraction
  with skysub=no, and still specify a user-defined sky value as the skyuser
  keyword.  This will reduce(eliminate?) the need to manually set
  ``MDRIZSKY=0``.

In addition to these major updates/changes, numerous smaller bugs were fixed
and other revisions were implemented which affected a small portion of the
use cases, such as:

- headerlet code now accepts lists of files to be updated.

- source sky positions (RA and Dec) now included in match file.

- DQ flags can now be taken into account when performing source finding in
  ``tweakreg``.

- all intermediate files generated by ``astrodrizzle`` will now be removed when
  using 'clean'='yes'.

- a problem was fixed that caused ``createMedian`` to crash where there were no
  good pixels in one of the images (when they did not overlap).

- interpretation of shiftfile now improved to handle arbitrarily-long
  filenames, rather than being limited to 24 character filenames.

- documentation has been updated, sometimes with a lot more extensive
  descriptions.

This version of ``DrizzlePac`` also requires use of the latest release version
of astropy primarily for WCS and FITS I/O support.


1.1.16
======

**Publicly Released through PyPI:** Mar 27, 2014

**Available under SSBX/IRAFX starting:** Mar 13, 2014

- Support for WFPC2 GEIS input images improved to correctly find the associated
  DQ images.

- Static mask files created for all chips in an image now get deleted when
  using the 'group' parameter to only drizzle a single chip or subset of chips.
- Fixed problem caused by changes to ``stsci.tools`` code so that
  ``drizzlepac`` will reference the correct extensions in input images.


1.1.15 (30-Dec-2013)
====================

**Publicly Released through PyPI:** Jan 14, 2014

**Available under SSBX/IRAFX starting:** Jan 6, 2014

Bug fixes
^^^^^^^^^

- Files created or updated by ``drizzlepac``, ``fitsblender``,
  or ``STWCS`` tasks, e.g. ``tweakreg`` or ``apply_headerlet``,
  will now ensure that the ``NEXTEND`` keyword value correctly reflects the
  number of extensions in the FITS file upon completion.


1.1.14dev (21-Oct-2013)
=======================

**Installed in OPUS:** Dec 11, 2013

**Available starting:** Oct 28, 2013

Bug fixes
^^^^^^^^^

- DQ arrays in input images now get updated with cosmic-ray masks
  computed by ``astrodrizzle`` when run with the parameter ``in_memory=True``.
  This restored the cosmic-ray masks detected during pipeline processing.


v1.1.13dev (11-Oct-2013)
========================

**available starting:** Oct 21, 2013

- ``tweakreg`` can now be run in 'batch' mode. This allows the user to generate
  plots and have them saved to disk automatically without stopping processing
  and requiring any user input.


1.1.12dev (05-Sep-2013)
=======================

**available starting:** Sept 9, 2013

This version fixed a couple of bugs in ``astrodrizzle``; namely,

- Logic was updated to support pixfrac = 0.0 without crashing. Ths code will
  now automatically reset the kernel to 'point' in that case.
- ``astrodrizzle`` now forcibly removes all OPUS WCS keywords from drizzle
  product headers.

- Default rules for generating drizzle product headers (as used in the archive)
  were modified to add definitions for 'float_one', 'int_one', 'zero' that
  generate output values of 1.0, 1, and 0 (zero) respectively for use as
  keyword values. This allows the LTM* rules to replace 'first' with
  'float_one' so that the physical and image coordinates for drizzle
  products are consistent.

Additionally, changes were made to ``STWCS`` for reprocessing use:

- Problems with using ``apply_headerlet_as_primary()`` from the ``STWCS``
  package on WFPC2 data have been corrected in this revision.


1.1.11dev (05-Jul-2013)
=======================

**Available starting:** July 15, 2013

- AstroDrizzle now can process all STIS data without crashing.


1.1.10dev (06-Feb-2013)
=======================

**available starting:** May 6, 2013

- The output drizzle image header no longer contains references to D2IM arrays.
  This allows ``tweakreg`` to work with drizzled images as input where 2-D D2IM
  corrections were needed.

- Deprecated references to PyFITS .has_key() methods were also removed from
  the entire package, making it compatible with PyFITS 3.2.x and later.


1.1.8dev (06-Feb-2013)
======================

**available starting:** Feb 11, 2013

- Fixed a bug in ``astrodrizzle`` which caused blot to raise an exception
  when using 'sinc' interpolation.

- Cleaned up the logic for writing out the results from the pixtopix, pixtosky,
  and skytopix tasks to avoid an Exception when a list of inputs are provided
  and no output file is specified.

- A new parameter was added to the tweakback task to allow a user to specify
  the value of ``WCSNAME`` when updating the FLT images with a new solution
  from a DRZ image header.

- Code in tweakback for updating the header with a new WCS will now
  automatically generate a unique ``WCSNAME`` if the there is a WCS solution in
  the FLT headers with the default or user-defined value of ``WCSNAME``.


1.1.7dev (18-Dec-2012)
======================

**available starting:** Feb 4, 2013

- Updated astrodrizzle to work with input images which do not have ``WCSNAME``
  defined. This should make it easier to support non-HST input images in the
  future.

- cleared up confusion between flux parameters in imagefindpars and catalog
  inputs in ``tweakreg``.

- turned of use of fluxes for trimming input source catalogs when no flux
  column can be found in input source catalogs.


1.1.7dev (18-Dec-2012)
======================

**available starting:** Dec 10, 2012

- Update ``tweakreg`` 2d histogram building mode to correctly find the peak
  when all the inputs match with the same offset (no spurious sources in either
  source catalog).

- Fixed a bug so that Ctrl-C does not cause an exception when used while
  ``tweakreg`` is running.

- revised the source finding logic to ignore sources near the image edge,
  a change from how daofind works (daofind expands the image with blanks
  then fits anyway).

- created a new function to apply the nsigma separation criteria to (try to)
  eliminate duplicate entries for the same source from the source list.
  It turns out daofind does have problems with reporting some duplicate sources
  as well. This function does not work perfectly, but works to remove nearly
  all (if not all) duplicates in most cases.


1.1.7dev (8-Jan-2012)
=====================

**available starting:** Jan 14, 2013

- Bug fixed in updatehdr module to allow shiftfiles without RMS columns to work
  as inputs to manually apply shifts to headers of input images.

- Revised ``astrodrizzle`` to update WCS of all input images BEFORE checking
  whether or not they are valid. This ensures that all files provided as input
  to ``astrodrizzle`` in the pipeline have the headers updated with the
  distortion model and new WCS.

- Images with NGOODPIX=0 now identified for WFC3 and WFPC2 inputs, so they
  can be ignored during ``astrodrizzle`` processing.
- Replaced 2d histogram building code originally written in Python with
  a C function that run about 4x faster.


1.1.6dev (5-Dec-2012)
=====================

**available starting:** Dec 10, 2012

- ``tweakreg`` v1.1.0 source finding algorithm now runs many times faster
  (no algorithmic changes). No changes have been made yet to speed
  up the 2d histogram source matching code.

- The 'pixtopix' task was updated to make the 'outimage' parameter optional
  by using the input image as the default. This required no API changes, but
  the help files were updated.

- Very minor update to guard against MDRIZTAB being specified without
  any explicit path.

- Update ``astrodrizzle`` to correctly report the exposure time,
  exposure start, and exposure end for the single drizzle products,
  in addition to insuring the final drizzle values remain correct.

- ``astrodrizzle`` also includes initial changes to safeguard the C code
  from getting improperly cast values from the configObj(TEAL) input.


1.1.5dev (23-Oct-2012)
======================

**available starting:** Oct 29, 2012

- Scaling of sky array for WFC3/IR IVM generation now correct.

- template mask files for WFPC2 no longer generated so that WFPC2 data can now
  be processed using num_cores > 1 (parallel processing).

- interpretation of the 'group' parameter fixed to support a single integer,
  a comma-separated list of integers or a single 'sci,<n>' value. The values
  correspond to the FITS extension number of the extensions that should be
  combined. This fix may also speed up the initialization step as more direct
  use of pyfits was implemented for the interpretation of the 'group'
  parameter.


1.1.1 (31-Aug-2012)
===================

**available starting:** Sept 26, 2012

The HST Archive and operational calibration pipeline started using this
version on Sept 26, 2012.


1.1.4dev (20-Sep-2012)
======================

**available starting:** Sept 24, 2012

- Bug fixed to allow use of final_wht_type=IVM for processing WFPC2 data.

- Revised Initialization processing to speed it up by using more up-to-date,
  direct pyfits calls.


1.1.3 (7-Sep-2012)
==================

**available starting:** Sept 17, 2012

- Fixed the logic so that crclean images always get created regardless of the
  value of the 'clean' parameter.


1.1.2 (5-Sep-2012)
==================

**available starting:** Sept 10, 2012

- Remove the restriction of only being able to process images which have
  ``WCSNAME`` keyword as imposed by r15631. The removal of this restriction
  will now allow for processing of non-updated input files with
  ``updatewcs=False`` for cases where no distortion model exists
  for the data (as required by CADC).

- Added log statements reporting what sky value was actually used in the
  drizzle and blot steps


1.1.1 (30-Aug-2012)
===================

**available starting:** Sept 3, 2012

- Major revision to ``astrodrizzle`` allowing the option to process without
  writing out any intermediate products to disk. The intermediate products
  remain in memory requiring significantly more memory than usual. This
  improves the overall processing time by eliminating as much disk activity
  as possible as long as the OS does not start disk swapping due to lack
  of RAM.

- revised to turn off 'updatewcs' when coeffs=False(no) so that exposures with
  filter combinations not found in the IDCTAB will not cause an error.


1.0.7 (21-Aug-2012)
===================

**available starting:** Aug 27, 2012

- Fixes problems with missing single_sci images.

- Static mask step revised to skip updates to static mask if all pixel data
  falls within a single histogram bin. This avoids problems with masking out
  entire images, which happens if low S/N SBC data is processed with
  ``static_mask=yes``.


1.0.6 (14-Aug-2012)
===================

**available starting:** Aug 20, 2012

Use of IVM for final_wht now correct, as previous code used wrong inputs when
IVM weighting was automatically generated by ``astrodrizzle``.


1.0.5 (8-Aug-2012)
==================

**available starting:** Aug 13, 2012

- Completely removed the use of the TIME arrays for weighting IR drizzle
  products so that the photometry for saturated sources in drizzled products
  now comes out correct.

- Corrected a problem with ``astrodrizzle`` which affected processing of WFPC2
  data where CRPIX2 was not found when creating the output single sci image.


1.0.2 (13-July-2012)
====================

**available starting:** Aug 3, 2012

The complete version of stsci_python can be downloaded from our
`download page <http://www.stsci.edu/institute/software_hardware/pyraf/stsci_python/current/stsci-python-download>`_

- `stsci_python v2.13 Release Notes <http://www.stsci.edu/institute/software_hardware/pyraf/stsci_python/release-notes/releasenotes.2.13>`_

- `Old stsci_python release notes <http://www.stsci.edu/institute/software_hardware/pyraf/stsci_python/release-notes>`_


1.0.1 (20-June-2012)
====================

**Used in archive/pipeline starting:** July 10, 2012

Pipeline and archive started processing ACS data with this version.


1.0.0 (25-May-2012)
===================

**Used in archive/pipeline starting:** June 6, 2012

Pipeline and archive first started using ``astrodrizzle`` by processing WFC3
images.
