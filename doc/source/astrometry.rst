.. _astrometry:

===================================================
Headerlets and You: Astrometry in Drizzle Products
===================================================

The astrometry for any given observation relies upon accurate pointing information from the telescope.   However, HST has evolved over time since it was put into orbit, with the focus changing over time as the telescope desorbs.  The instruments have also changed positions over time relative to the FGS guides, and the coordinates for the guide stars were orginally determined using ground-based information.  All this has limited the calculation of the pointing of any given observation on the sky (absolute astrometry) to no better than 1-2 arc-seconds.

Calibration of the distortion model for each instrument has been done well enough to allow observations to be combined (relative astrometry) with an accuracy of better than 5 milli-arcseconds.  This has made the absolute astrometry inaccuracy stand out even more, making it more difficult to compare observations taken at different times or to compare HST data with observations from other telescopes (like Chandra).

Therefore, multiple efforts have been undertaken to improve the absolute astrometry to match the accuracy of the relative astrometry.  These efforts have resulted in multiple new world coordinate systems (WCS) solutions to be developed for HST observations, starting with ACS, WFC3, and WFPC2.  Complete details of the results of calibrating the distortion models and astrometry for each instrument can be found in each of the instruments web pages; namely,

    * **ACS**: `ACS Distortion <http://www.stsci.edu/hst/acs/analysis/distortion/>`_
    * **WFC3**: WFC3 Data Handbook `Chapter 4: WFC3 Images: Distortion Correction and AstroDrizzle <http://www.stsci.edu/hst/wfc3/documents/handbooks/currentDHB/Chapter4_astrometry1.html#>`_

Each solution has its own advantages and errors, making some good for one use but inadequate for others.  As a result,  **all new WCS solutions which are approved by STScI** are being offered with HST data provided by MAST with headerlets serving as the mechanism for providing and applying all WCS solutions.


Where are all the WCS solutions?
================================

All calibrated HST products delivered by the Barbara A. Mikulski Archive for Space Telescopes (MAST) contain the best available WCS solutions available at the time the data was last processed.  Unfortunately, only 1 WCS can be used at a time to transform the position in the image to an undistorted position on the sky.  The FITS Standard, however, describes how multiple WCS solutions can be defined for an image in `FITS Paper I (Greisen, E. W., and Calabretta, M. R., Astronomy & Astrophysics, 395, 1061-1075, 2002) <http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=2002A%26A...395.1061G&db_key=AST&high=3db47576cf06933>`_.

This standard defines a keyword, WCSNAME, which serves as the label for the WCS specified in the header of the observation.  It also defines additional keywords, WCSNAME[A-Z], that represent additional WCS solutions which are specified as alternate WCS solutions in the same header as the active WCS associated with WCSNAME.  Headerlets rely on the use of the WCSNAME keyword to manage the multiple solutions that may be defined for a given observation, with each WCS having it's own unique WCSNAME value.  Headerlets also build the HDRNAME keyword which specifies the exposure name as well WCS to insure that this WCS will only ever be applied to this exposure and no other.

The use of headerlets extends this standard by allowing each new WCS, specified using the FITS Paper I standards along with the SIP convention, to be stored separately as a new extension at the end of the image's FITS file.  This can be seen in the extensions listed by astropy.fits for a sample ACS/WFC image which includes 6 additional WCS solutions from the astrometry database.::

  No.    Name      Ver    Type      Cards   Dimensions   Format
    0  PRIMARY       1 PrimaryHDU     279   ()
    1  SCI           1 ImageHDU       201   (4096, 2048)   float32
    2  ERR           1 ImageHDU        57   (4096, 2048)   float32
    3  DQ            1 ImageHDU        49   (4096, 2048)   int16
    4  SCI           2 ImageHDU       199   (4096, 2048)   float32
    5  ERR           2 ImageHDU        57   (4096, 2048)   float32
    6  DQ            2 ImageHDU        49   (4096, 2048)   int16
    7  D2IMARR       1 ImageHDU        15   (64, 32)   float32
    8  D2IMARR       2 ImageHDU        15   (64, 32)   float32
    9  D2IMARR       3 ImageHDU        15   (64, 32)   float32
   10  D2IMARR       4 ImageHDU        15   (64, 32)   float32
   11  WCSDVARR      1 ImageHDU        15   (64, 32)   float32
   12  WCSDVARR      2 ImageHDU        15   (64, 32)   float32
   13  WCSDVARR      3 ImageHDU        15   (64, 32)   float32
   14  WCSDVARR      4 ImageHDU        15   (64, 32)   float32
   15  WCSCORR       1 BinTableHDU     59   14R x 24C   [40A, I, A, 24A, 24A, 24A, 24A, D, D, D, D, D, D, D, D, 24A, 24A, D, D, D, D, J, 40A, 128A]
   16  HDRLET        1 HeaderletHDU     22   ()
   17  HDRLET        2 HeaderletHDU     22   ()
   18  HDRLET        3 HeaderletHDU     26   ()
   19  HDRLET        4 HeaderletHDU     26   ()
   20  HDRLET        5 HeaderletHDU     26   ()
   21  HDRLET        6 HeaderletHDU     26   ()

In this example, extensions 16 through 21 contain headerlets each of which specifies different astrometric alignment (each may potentially be based on a different distortion model) and each has been assigned a unique name given by the **WCSNAME** keyword.  A quick listing of the names of all these WCSs can be obtained using the Python library STWCS:

.. code-block:: python

    from stwcs.wcsutil import headerlet
    headerlet.get_headerlet_kw_names("j95y04hq_flc.fits", kw='WCSNAME')
    ['OPUS',
    'IDC_0461802ej',
    'OPUS-GSC240',
    'IDC_0461802ej-GSC240',
    'OPUS-HSC30',
    'IDC_0461802ej-HSC30']


.. _wcsname-conventions:

Interpreting WCS names
-----------------------
The first WCS solution listed, OPUS, corresponds to the default WCS defined by the HST calibration pipeline based on original telemetry and only a first-order distortion model. Solutions which contain **IDC_** are based on the distortion model based on the reference data associated with the **IDCTAB** reference file whose rootname is specified in the WCSNAME.  For example, **IDC_0461802ej** refers to a WCS solution based on the **IDCTAB** reference file 0461802ej_idc.fits along with all the associated secondary distortion reference files given by **DGEOFILE**, **NPOLFILE** and, for WFPC2, **OFFTAB** and **DXYFILE**.

Additional WCS solutions will then be based on either the original **OPUS** WCS or the distortion-corrected **IDC_** WCS solutions.  Two types of solutions can be defined for images; namely, *a priori* and *a posteriori* solutions.

a priori solutions
^^^^^^^^^^^^^^^^^^
The *a priori* solutions have been determined for **ALL HST data** by correcting the coordinates of the guide stars that were used from the originally specified coordinates to the coordinates of those guide stars as determined by GAIA.  The naming convention for these *a priori* solutions are::

  <Starting WCS>-<Astrometric Catalog>

  For example,
  'IDC_0461802ej-GSC240'

where the **Astrometric Catalog** refers the exact astrometric catalog used to correct the guide star positions.  A number of **Astrometric Catalogs** are available through MAST for aligning images.  Solutions generated for the database were initially based on catalogs which were based on the GAIA catalog.  These GAIA-based catalogs include:

  * **GSC240**

    - This catalog contains version 2.4.0 of the *Guide Star Coordinates* (GSC) catalog,
    - All guide stars in the catalog were cross-matched with the GAIA DR1 catalog and corrected to the coordinates reported in GAIA DR1.
    - **APPLIES TO**:  All HST datasets which had a successful guide star acquisition, which is nearly all data in the archive.
    - **NOTE**:  HST Observations taken after September 2017.

  * **HSC30**

    - This catalog contains version 3.0 of the *Hubble Source Catalog* (HSC)
    - Technically, this is an *a posteriori* solution, but it is applied blindly without further verification that the correction fully aligns the image to GAIA;hence, it is included as an *a priori* solution.
    - Sources in the HSC were cross-matched with the GAIA DR1 catalog.
    - Those cross-matched sources were then used to determine a fit to the GAIA catalog.
    - The fit to GAIA was then applied to all remaining sources in the catalog.
    - **APPLIES TO**:  Only datasets which had a sufficient number of sources in the exposure to be aligned to GAIA by the *Hubble Legacy Archive(HLA)* project.

  * **GAIADR1**

    - A MAST-provided version of the first data release (DR1) version of the official GAIA astrometric catalog.
    - This version does not have proper motions for a majority of the sources in the catalog.

  * **GAIADR2**

    - A MAST-provided version of the second data release version of the official GAIA astrometric catalog.
    - This catalog contains initial proper motion measurements (and errors) for most sources in the catalog.

Although all solutions are appended to each FITS file, only 1 WCS (referred to as the **'active' WCS**) can be used at a time to represent the transformation from pixel coordinates to world coordinates. The active WCS is defined by the standard WCS keywords found in the header of the science extension for each chip in the exposure; *e.g.*, CRVAL1, CRVAL2, CRPIX1, CRPIX2, and so on.

The *a priori* solution which gets selected to replace the active WCS solution represents the most accurate solution available in the astrometry database at the time, and will be chosen based on the following hierarchy (as of Summer 2019):

  #. HSC30
  #. GSC240
  #. IDC_<rootname> (distortion-corrected pipeline default WCS)

If the first type of solution is not available, the next solution in the list is selected. As new solutions are added to the astrometry database, these rules will be modified to always try to return the WCS correction which best aligns the data to the GAIA catalog.

a posteriori solutions
^^^^^^^^^^^^^^^^^^^^^^^
The *a posteriori* solutions, on the other hand, get determined from measuring sources in each image, finding overlapping sources from an astrometric catalog, identifying and cross-matching image sources with sources from the astrometric catalog and performing a fit to correct the WCS.  These type of solutions can not be determined for all datasets due to a number of reasons, such as lack of sources in the image and/or lack of overlapping sources from an astrometric catalog.  When these solutions can be determined for an observation, they are given a value for the `WCSNAME` keyword which follows the convention:

  **<Starting WCS>-FIT_<REL|IMG|EVM|SVM>_<Astrometric Catalog>**

For example,

  'IDC_0461802ej-FIT_REL_GAIADR2'

The terms are defined as:

  * **<Starting WCS>**

    - Value of WCSNAME for the exposure prior to applying any astrometric Solutions
    - IDC_<rootname> (like `IDC_041802ej`) refers to a distortion-corrected model based on the IDCTAB reference file `0461802ej_idc.fits`.

  * **`FIT`**

    - This term refers to the fact that sources from the image were identified, cross-matched and fit to sources from an astrometric catalog to create an *a posteriori* WCS solution.  
        
  * **`<REL|IMG|EVM|SVM>`**

    - `REL` : This term denotes the fact that all images were aligned relative (REL) to each other and then aligned to an astrometric catalog.  This attempts to maintain the original relative alignment between the images in a given visit.
    - `IMG` : This term denotes the fact the the images were fit individually to the astrometric catalog.  These solutions are applied only when relative alignment does not yield a viable fit to the astrometric catalog.
    - `EVM` : The cross-match and fit to an astrometric catalog was performed on a single exposure by itself as part of processing the exposures of an entire visit. This will typically only apply to those rare visits which do not have enough valid exposures in the visit for alignment.
    - `SVM` : This term refers to alignment of all the exposures in a single-visit to an astrometric catalog.  The exposures of a visit are aligned to each other (relative alignment), then, as a group, all the exposures are cross-matched and fit to the astrometric catalog specified in the next term in the WCSNAME.


  * **<Astrometric Catalog>**

    - This term describes the astrometric catalog, as listed for use with the *a priori* solutions, which was used for the cross-matching and fitting sources identified in the image(s).  If a value of **NONE** is specified here, it indicates that although the image appears (according to the code) to have been successfully relatively aligned one exposure to another, there were indications that the alignment to an astrometric catalog like **GAIADR2** failed.  The user will need to carefully review the state of alignment of this data when **NONE** is listed in the output WCS.


These separate terms provide as succinct a description of the solution determined for and applied to the exposure as possible. Additional keywords have been written out to the headerlet extension for the *a posteriori* fit which further describe the solution, including:

  * number of sources used in the fit
  * RMS in RA and Dec of the fit
  * parameters determined for the fit
  * and more...

.. note::

    A successfully determined *a posteriori* solution will **always** be used to replace the active WCS (after insuring the previous WCS has been saved as a headerlet extension already) regardless of the original solution.

Pipeline Processing
-------------------
All HST observations get processed in an automated environment using standard
parameters for the calibration code, including the alignment and combination of 
individual exposures into undistorted products.  The standard pipeline processing
to create the undistorted drizzled images (drc.fits or drz.fits) gets performed 
using the 'runastrodriz' task in this package.  This same processing can be 
run at any time using:

.. code-block:: bash

    runastrodriz j8cw03010_asn.fits
    
    runastrodriz j8cw03f6q_raw.fits
    
The files which need to be present are:

    * RAW files (\*raw.fits)
    * FLT files (\*flt.fits)
    * FLC files (\*flc.fits, if any were created by the pipeline)
    * ASN file  (\*asn.fits, if applicable)
    
This processing includes a lot of logic intended to not only apply pre-defined (apriori) 
WCS solutions, but also to try and determine a new aposteriori solution then 
verify which solution (default pipeline, apriori or aposteriori) actually provides
the WCS which comes closest to the GAIA astrometric frame.  
The :ref:`runastrodriz-description` of the runastrodriz task provides 
the full discussion of the logic used to define the
defined 'active' WCS that gets used to create the products which get archived.
 

Choosing a WCS
---------------
The **only** WCS solution that gets used to perform coordinate transformations on the pixel values will be the 'active' or 'primary' WCS associated with the WCSNAME keyword.  The pipeline generated products will include an active WCS which the pipeline specifies as the *best* available WCS given the information used at the time of processing.  However, this default 'active' WCS may not be appropriate for all science, so this WCS may need to be replaced by one of the other WCSs instead to best support the analysis necessary for the research.

Dependent Packages
^^^^^^^^^^^^^^^^^^^^
Working with the WCS solutions and headerlets gets performed using `STWCS package <https://stwcs.readthedocs.io/en/latest/>`_.  Examples of how to work with this package will assume that the user has already installed this package into their working Python environment and has started a python shell.  In addition, the following example relies on the Astropy IO package to work with the FITS headers and extensions.

Finally, the example described here will rely on additional functionality included in the V3.2.0 or later of the Drizzlepac package.  These new functions support the generation of drizzle combined products which have been aligned to an astrometric standard catalog such as GAIA DR2.

Sample Session
^^^^^^^^^^^^^^^
This example will work with 4 exposures taken using ACS/WFC as the association *j95y04010*, with most of the example being performed on the single exposure *j95y04hpq_flc.fits*.

All the necessary libraries for working on this example can be imported using:

.. code-block:: python

    from drizzlepac.haputils import astroquery_utils as aqutils
    from drizzlepac.haputils import astrometric_utils as amutils
    from astropy.io import fits
    from stwcs.wcsutil import headerlet

The data can be obtained from MAST with `astroquery <https://astroquery.readthedocs.io/en/latest/>`_ using a simplified interface developed in **drizzlepac** using the commands:

.. code-block:: python

    filenames = aqutils.retrieve_observation('j95y04010')
    # filenames = ['j95y04hpq_flc.fits', 'j95y04hqq_flc.fits',
    #              'j95y04hsq_flc.fits', 'j95y04huq_flc.fits']

The default 'active' WCS can be determined using:

.. code-block:: python

    default_wcsname = fits.getval(filenames[0], 'wcsname', ext=1)

For this example, we find that **default_wcsname='IDC_0461802ej'**, or as noted earlier, the default distortion correction based WCS provided by the pipeline with no correction for any astrometric catalogs, has been designated by the pipeline as the 'active' WCS.

All available WCSs provided by the astrometry database and attached as headerlet extensions can be queried to find the WCSNAMEs for all the new WCSs using:

.. code-block:: python

    new_wcsnames = headerlet.get_headerlet_kw_names(filenames[0], kw='WCSNAME')

These will be the same ones listed earlier.  For this example, we decide we would like to have this observation aligned using the guide stars corrected to GAIA DR1 through the use of the GSC240-based WCS; specifically, **IDC_0461802ej-GSC240**.  We can replace the 'active' WCS with this new one using:

.. code-block:: python

    new_hdrnames = headerlet.get_headerlet_kw_names(filenames[0])
    # identify hdrname that corresponds with desired WCS with name of IDC_041802ej-GSC240
    new_wcs = new_hdrnames[new_wcsnames.index('IDC_0416802ej-GSC240')]
    headerlet.restore_from_headerlet(filenames[0], hdrname=new_wcs, force=True)
    # confirm new WCS is now 'active'
    fits.getval(filenames[0], 'wcsname', ext=1)

At this point, the exposure has been updated to perform all coordinate transformations with the new GAIA DR1-based WCS as if the guide stars used for taking the observation has GAIA DR1 coordinates in the first place.  This will not mean it will align perfectly with GAIA, but should be within 0.5 arcseconds due to drift in the telescope field-of-view for each of the instruments relative to the FGSs.


Headerlet Primer
=================

The headerlet file itself conforms to FITS standards with the PRIMARY header containing global information about the WCS solution and how it was determined.  Separate extensions in the headerlet then contain the header keywords for specifying the WCS for each chip in the exposure or for the distortion information necessary to correct the pixel positions from the image to the un-distorted position on the sky.  These solutions rely on calibration reference data that describe the distortion observed in each instrument to better than 0.1 pixels in each detector.  Instead of having to retrieve separate files with this distortion information, that distortion information has been folded into the header of each WFC3, ACS and WFPC2 dataset.

Headerlet File Structure
-------------------------
This new object complete with the NPOLFILE and the D2IMFILE extensions
derived from the full FITS file fully describes the WCS of each chip
and serves without further modification as the definition of the
`headerlet`. The listing of the FITS extensions for a `headerlet` for
the sample ACS/WFC exposure after writing it out to a file would then be::

    EXT#  FITSNAME      FILENAME              EXTVE DIMENS       BITPI OBJECT

    0     j8hw27c4q     j8hw27c4q_hdr.fits                       16
    1       IMAGE       D2IMARR               1     4096         -32
    2       IMAGE       WCSDVARR              1     64x32        -32
    3       IMAGE       WCSDVARR              2     64x32        -32
    4       IMAGE       WCSDVARR              3     64x32        -32
    5       IMAGE       WCSDVARR              4     64x32        -32
    6       IMAGE       SIPWCS                1                  8
    7       IMAGE       SIPWCS                2                  8

Detailed Description of headerlet
----------------------------------
The full details on the headerlet, it's required set of keywords, and how the distortion models get described in the headerlet can be found in the `Technical Report on headerlets <https://stwcs.readthedocs.io/en/latest/headerlet_tsr/source/index.html>`_.

Code Interface to headerlets
-----------------------------
The `STWCS package <https://stwcs.readthedocs.io/en/latest/>`_ provides the code used to work with headerlets and WCS solutions.
