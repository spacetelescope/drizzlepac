.. _singlevisit:

==============================
Single-visit Mosaic Processing
==============================
Standard calibration pipeline processing focuses on aligning data taken as 
associations or as single exposures as closely to a standard astrometric coordinate
frame as possible (see :ref:`runastrodriz-description` for full details).  
This involves alignment of very few images all taken under
as identical conditions as possible; namely, detector, filter, guide stars, 
guiding mode, and so on.  These observations were intended by the user to 
represent a single view of the object of interest in a way that allows for 
removal of as many calibration effects as possible. 

Observations taken within a single visit, on the other hand, represent data 
intended to produce a multi-wavelength view of the objects in the desired
field-of-view. Most of the time, those observations are taken using pairs of guide 
stars to provide very stable field-of-view throughout the entire visit resulting 
in images which overlap almost perfectly.  Unfortunately, this is not
always possible due to increasing limitations of the aging telescope systems.
The result is that an increasing number of visits are taken where observations 
drift and/or roll during the course of the visit or re-acquire at slightly 
different pointings from one orbit to the next.  Whatever the reasons, data across
each visit cannot automatically be assumed to align.  Single-visit mosaic (SVM) 
processing attempts to correct these relative
alignment errors and to align the data to an absolute astrometric frame so that
all the data across all the filters used can be drizzled onto the same pixel grid.

Understanding the quality of the SVM products requires knowing what processing
took place to generate that product, and more importantly, what limitations that
processing may have had.  The following sections provide the description of the
single-visit processing.  Definitions of a number of processing-specific terms used 
in this description can be found in the :ref:`hap-glossary`.


Primary User-Interface
=======================
One task has been written to perform the single-visit processing: ``runsinglehap``. 
It gets used by STScI to generate the single-visit products which
can be found in the Mikulsik Archive for Space Telescopes (MAST) archive. This task 
can also be run from the operating system command-line or from within a
Python session to reproduce those results, or with modification of the input 
parameters, perhaps improve on the standard archived results.  Full details on 
how to run this task can be found in the description of the task at :ref:`runsinglehap_api`.


Processing Steps
================
Single-visit processing performed by ``runsinglehap`` 
relies on the results of the standard astrometric 
processing of the individual exposures and associations as the starting point
for alignment. This processing then follows these steps to create the final products:

  * interpret the list of filenames for all exposures taken as part of a single visit
  * copy the pipeline-calibrated (FLT/FLC) files to the current directory for processing
  * rename the input files to conform to the single-visit naming conventions
    
    * This step insures that the original pipeline results remain available in 
      the archive unchanged
  
  * Define what output products can be generated 
  * Align all exposures in a relative sense (all to each other)
  * Create a composite source catalog from all aligned input exposures
  * Cross-match and fit this composite catalog to GAIA to determine new WCS solution
  * Update renamed input exposures with results of alignment to GAIA
  * Create each of the output products using the updated WCS solutions
 

.. _svm_naming_convention:

Single Visit Naming Convention
==============================
All files processed as part of a single visit get renamed from the standard
pipeline filenames into something which describes the data more clearly.  The 
convention used for the names of these input and output files uses these 
components:

  * **<propid>** : the proposal ID for this visit
  * **<obsetid>** : the 2-digit visit ID from this proposal 
  * **<instr>** : 3 letter designation of the instrument used for the observations
  * **<detector>** : name of the detector used for the observations
  * **<filter>** : hyphen-separated list of filter names used for the observations
  * **<ipppssoo>** : standard 8 character rootname for a single exposure defined by the pipeline
  * **<ipppss>** : standard 6 character designation of the **<instr>**/**<propid>**/**<obsetid>** for this visit
  * **dr[cz].fits** : suffix for drizzled products
  * **fl[ct].fits** : suffix for pipeline-calibrated files
  * **asn.fits** : suffix for the association table
  * **hlet.fits** : suffix for headerlet files containing the WCS solution used to create the final drizzle products/mosaics
  
These components get combined to create filenames specific to each type of file being
processed.  The following table provides a complete list of all the products 
created as a result of single-visit processing.

.. list-table:: Single-visit product filenames
  :widths: 8 25 83
  :header-rows: 1
  
  * - Product
    - File Type
    - Filename for Files Produced
  * - Exposure
    - drizzle product
    - hst_<propid>_<obsetid>_<instr>_<detector>_<filter>_<ipppssoo>_dr[cz].fits
  * -
    - flat-field product
    - hst_<propid>_<obsetid>_<instr>_<detector>_<filter>_<ipppssoo>_fl[ct].fits
  * - 
    - headerlet file
    - hst_<propid>_<obsetid>_<instr>_<detector>_<filter>_<ipppssoo>_hlet.fits
  * -
    - trailer file
    - hst_<propid>_<obsetid>_<instr>_<detector>_<filter>_<ipppssoo>_trl.txt
  * -
    - preview (full size)
    - hst_<propid>_<obsetid>_<instr>_<detector>_<filter>_<ipppssoo>_dr[cz].jpg
  * - 
    - preview (thumbnail)
    - hst_<propid>_<obsetid>_<instr>_<detector>_<filter>_<ipppssoo>_dr[cz]_thumb.jpg
  * - Filter
    - drizzle product
    - hst_<propid>_<obsetid>_<instr>_<detector>_<filter>_<ipppss>_dr[cz].fits
  * -
    - point-source catalog
    - hst_<propid>_<obsetid>_<instr>_<detector>_<filter>_<ipppss>_point-cat.ecsv
  * -
    - segment-source catalog
    - hst_<propid>_<obsetid>_<instr>_<detector>_<filter>_<ipppss>_segment-cat.ecsv
  * - 
    - trailer file
    - hst_<propid>_<obsetid>_<instr>_<detector>_<filter>_<ipppss>_trl.txt
  * -
    - preview (full size)
    - hst_<propid>_<obsetid>_<instr>_<detector>_<filter>_<ipppss>_dr[cz].jpg
  * -
    - preview (thumbnail)
    - hst_<propid>_<obsetid>_<instr>_<detector>_<filter>_<ipppss>_dr[cz]_thumb.jpg
  * - Total
    - drizzle product
    - hst_<propid>_<obsetid>_<instr>_<detector>_total_<ipppss>_dr[cz].fits
  * -
    - point-source catalog
    - hst_<propid>_<obsetid>_<instr>_<detector>_total_<ipppss>_point-cat.ecsv
  * -
    - segment-source catalog
    - hst_<propid>_<obsetid>_<instr>_<detector>_total_<ipppss>_segment-cat.ecsv
  * - 
    - trailer file
    - hst_<propid>_<obsetid>_<instr>_<detector>_total_<ipppss>_trl.txt
  * -
    - preview (full size)
    - hst_<propid>_<obsetid>_<instr>_<detector>_total_<ipppss>_dr[cz].jpg
  * - 
    - preview (thumbnail)
    - hst_<propid>_<obsetid>_<instr>_<detector>_total_<ipppss>_dr[cz]_thumb.jpg
  * - 
    - color preview (full size)
    - hst_<propid>_<obsetid>_<instr>_<detector>_total_<ipppss>_<filters>_dr[cz].jpg
  * - 
    - color preview (thumbnail)
    - hst_<propid>_<obsetid>_<instr>_<detector>_total_<ipppss>_<filters>_dr[cz]_thumb.jpg
    

Processing the Input Data
=========================
SVM processing starts with a list of all the single exposures 
which were taken as part of a visit.  Any associations which were defined by the
proposal are ignored, since the visit itself gets treated, in essence, as a new 
association.  The input files can be specified either using the **poller** file format
used by the STScI automated processing or a file with a simple list of filenames
(one filename per line).

Automated poller input file format
----------------------------------
The automated processing performed to populate the MAST archive at 
STScI provides a file with the following format::

    ic0s17h4q_flt.fits,12861,C0S,17,602.937317,F160W,IR,ic0s/ic0s17h4q/ic0s17h4q_flt.fits
    ic0s17h5q_flt.fits,12861,C0S,17,602.937317,F160W,IR,ic0s/ic0s17h5q/ic0s17h5q_flt.fits
    ic0s17h7q_flt.fits,12861,C0S,17,602.937317,F160W,IR,ic0s/ic0s17h7q/ic0s17h7q_flt.fits
    ic0s17hhq_flt.fits,12861,C0S,17,602.937317,F160W,IR,ic0s/ic0s17hhq/ic0s17hhq_flt.fits

This example comes from the 'ic0s1' visit where the columns are:

  #. exposure filename
  #. proposal ID (numeric value)
  #. program ID - ppp value from exposure filename
  #. obset_id - visit number from proposal 
  #. exposure time of the exposure
  #. filters used for the exposure
  #. detector used to take the exposure
  #. location of the exposure in a local cache


Status of Input Data
----------------------
The list of filenames which should be processed as a single-visit provides the
raw science data for creating the new combined output products.  However, these
files need to be properly calibrated prior to SVM processing.  Specifically, the
exposures need to be:

  * fully calibrated using the instruments calibration software, such as 
    ``calacs.e`` for ACS and ``calwf3.e`` for WFC3 data.  This should also
    include CTE-correction for the images whenever possible.
  * processed using ``runastrodriz`` in order to apply the latest distortion
    model calibrations to the astrometry and to align the exposures as closely
    as possible to an external astrometric reference when possible.

These steps insure that the latest calibrations get applied to the data making it
easier for the SVM processing to cross-match the data with minimal interference 
from artifacts in the data.  In addition, the CTE-corrected versions of the data 
get used during pipeline processing in order to allow for better alignment of the 
exposures and to improve the photometry of the data as much as possible.  

These processing steps can be verified in the input data using header keywords from 
the exposures

.. list-table:: Processing keywords 
  :widths: 10 15 40
  :header-rows: 1
  
  * - Header Keyword
    - Valid Values
    - Notes
  * - FLATCORR
    - COMPLETED
    - Completion of basic calibration
  * - DRIZCORR
    - COMPLETED
    - Completion of distortion calibration
  * - WCSNAME
    - \-FIT
    - Successful **a posteriori** alignment
  * -
    - \-HSC30
    - Successful **a priori** alignment
  * -
    - \-GSC240
    - Successful **a priori** alignment

The full set of possibilities for updated WCSs as reported using the **WCSNAME**
keyword can be found in the description of the :ref:`wcsname-conventions`.

As long as the input data meets these requirements, then SVM processing will have
the best chance of success.  Data which has not been able to be aligned successfully
with an **a priori** or **a posteriori** solution can still be processed as part
of a single-visit, however, the alignment may be more difficult to determine due 
to the larger uncertainties for HST pointing prior to October 2017.  


Filtering the input data
--------------------------
Not all HST imaging observations can be aligned using SVM processing.  Observations
taken with the GRISM or in SPATIAL SCAN mode result in sources which can not be 
aligned, for example.  The :ref:`analyze_api` module evaluates all
input exposures using these header keywords for the stated rejection criteria.

.. list-table:: Single-visit product filenames
  :widths: 26 27 60
  :header-rows: 1
  
  * - Header Keyword
    - Values Which Trigger Rejection
    - Explanation
  * - OBSTYPE
    - (not IMAGING)
    - Only Imaging mode data processed
  * - MTFLAG
    - T 
    - No moving targets, WCS and background sources vary
  * - SCAN_TYP
    - C or D (or not N)
    - Can not align streaked sources
  * - FILTER or FILTER1, FILTER2
    - G*, PR*, BLOCK
    - G=Grism and PR=Prism, Can not align streaked sources
  * - EXPTIME
    - 0 
    - no exposure time, no data to align
  * - TARGNAME
    - DARK, TUNGSTEN, BIAS, FLAT, 
    - No alignable external sources in these calibration modes 
  * - 
    - EARTH-CALIB, DEUTERIUM
    - No alignable external sources in these calibration modes 
  * - CHINJECT
    - not NONE
    - No alignable external sources in these calibration modes 


Any observation which meets any of these criteria are flagged to be ignored (not
processed).  In addition, any data taken where the FGSLOCK keyword contains 'COARSE' or 'GY' will be flagged as potentially compromised in the comments generated during
processing.

All observations which are alignable based on these criteria are then
passed along as a table to create the SVM products.  Those inputs which can be
processed are then copied and renamed using the :ref:`svm_naming_convention`.  This 
insures that no SVM processing will affect or otherwise modify the original 
pipeline-processed input files.  Only the SVM named input files will be updated
with new SVM-aligned WCS solutions and then used to produce the drizzle products.  


Defining the Output Products
============================

The table with the set of observations which can be processed now gets interpreted.
The goal is to identify what exposures can be combined to create unique products.  
This grouping will be used to create the **product list**.  
The **product list** is a Python list of 
`drizzlepac/haputils/product/HAPProduct` objects, described in :ref:`product_api` API docs,
which represent each and every output product to be created for the visit.  
Each **Product** instance contains:

  * list of filenames for all input exposures that will contribute to the output drizzlep product
  * WCS for output drizzle product
  * pre-defined names for all output files associated with this **Product** including:

    * drizzle-combined image
    * point-source catalog determined from the drizzle-combined image
    * segmentation-based catalog determined from the drizzle-combined image
    * astrometric catalog used to align the input exposures
      
  * methods for:
    
    * determining average number of images per pixel
    * defining the final WCS
    * aligning the exposures to an astrometric reference (GAIA)
    * applying the selected parameters to ``AstroDrizzle``
    * drizzling the inputs to create the output drizzle product
    * determining the source catalogs fron the drizzle product 

This interpretation of the list of input filenames gets performed using the 
code in :ref:`poller_utils_api` by
grouping similar observations.    The rules used for grouping the inputs into output
products result in outputs which have the same detector and filter.  These output 
products are referred to as **filter products** defined as a ``product/FilterProduct``
instance.  

All exposures for a single detector are also identified and grouped to 
define a **total product** using the ``product/TotalProduct`` class.  
This **total product** drizzle image provides the deepest available 
view of the field-of-view from this visit which will be used to produce the master
catalog of sources for this visit.  The master catalog of source positions will
be used to perform photometry on each exposure, whether the source can be identified
in the exposure at that position or not.  This **forced photometry** results in
limits for the photometry in cases where the sources are not bright enough to be
identified in a given filter.

Two separate source catalogs for each filter are also pre-defined; namely, 

  * a point-source catalog derived using ``photutils`` ``DAOStarFinder`` 
  * a segmentation-based catalog derived using ``photutils`` segmentation code

These two catalogs provide complimentary views of each field-of-view to try to
highlight all types of compact sources found in the exposures. 


Example Visit
--------------
For example, a relatively simple visit of a fairly bright and crowded field with 
6 F555W exposures (two 15-second and four 30-second exposures) and 
6 F814W exposures (two 5-second and four 15-second exposures)
would result in the definition of these output products: 

  * a drizzled image for each separate exposure
  * a single F555W product
  * a single F814W product, and
  * a single **total product**
  * a point-source catalog for the F555W product
  * a segmentation-based source catalog for the F555W product
  * a point-source catalog for the F814W product
  * a segmentation-based source catalog for the F814W product
  * a point-source catalog for the total product
  * a segmentation-based catalog for the total product
  
The function ``haputils.poller_utils.interpret_obset_input`` serves as the sole interface 
for this interpretation. A basic tree gets defined (as a dictionary of dictionaries) 
by this function where the
output exposures are identified along with all the names of the input exposures.
This tree then serves as the basis for organizing the rest of the SVM processing.

In addition to defining what output products need to be generated, all the SVM
products names are defined using the :ref:`svm_naming_convention`.  This insures
that all the output products have filenames which are not only unique but also 
understandable (if a bit long) that are easily grouped on disk.  


Aligning the Input Data
=======================
All input exposures should have already been aligned either individually or by 
association table as close to GAIA as possible during standard pipeline calibration
processing.  However, each exposure or association (of exposures) can be aligned
to slightly different fits or catalogs due to differences in the source objects 
which can be identified in each separate exposure.  The primary goal of SVM 
processing is to refine this alignment so that all exposures in the visit for 
the same detector (those exposures which contribute to each **total product**)
share the same WCS (pixels on the sky).  

Alignment of all the exposures for a **total product** uses the same alignment
code as the standard calibration pipeline.  The basic steps it follows is:

  * generate a source catalog for each exposure (using :ref:`amutils_api`)
  * obtain the WCS from each exposure
  * perform a relative fit between the exposures using ``tweakwcs``
  * obtain an astrometric catalog for the field-of-view 
  * perform a final fit of all the exposures at once to the astrometric catalog
  * update each WCS with the final corrected WCS generated by ``tweakwcs``
  
The limits for performing the relative alignment and absolute fit to the astrometric
catalog (defaults to **GAIADR2**) are lower under the expectation that large 
offsets (> 0.5 arcseconds) have already been removed in the pipeline processing.  
This makes the SVM alignment more robust across a wider range of types of fields-of-view.
The final updated WCS will be provided with a name that reflects this cross-filter
alignment using **-FIT-SVM-<catalog name>** as the final half of the **WCSNAME** 
keyword.  More details on the WCS naming conventions can be found in the
:ref:`wcsname-conventions` section.


Creating the Output Products
============================
Successful alignment of the exposures allows them to be combined into the
pre-defined output products; primarily, the **filter products**  and the **total product**.
These products get created using ``drizzlepac.astrodrizzle.AstroDrizzle``. 

Selecting Drizzle Parameters
-----------------------------
Optimal parameters for creating every possible type of output product or mosaic
would require knowledge of not only the input exposures, but also expert
knowledge of the science.  Parameters optimized for one science goal may not be
optimal for another science goal.  Therefore, automated pipeline processing has
defined a basic set of parameters which will result in a reasonably consistent 
set of products as opposed to trying to optimize for any specific science case.  

The default parameters have been included as part of the ``drizzlepac`` package 
in the ``drizzlepac/pars/hap_pars`` directory.  Index JSON files provide the options
that have been developed for selecting the best available default parameter set
for processing.  The INDEX JSON files point to different parameter files (also in
JSON format) that are also stored in sub-directories the code organized by instrument
and detector.  

Selection criteria are also listed in these Index JSON files for each
step in the SVM processing pipeline; namely, 

  * alignment
  * astrodrizzle
  * catalog generation
  * quality control
  
Initially, only the **astrodrizzle** step defines any selection criteria for use
in processing.  The criteria is based on the number of images being combined for 
the specific instrument and detector of the exposures.  

The SVM processing interprets the input data and verifies what input data can be 
processed.  At that point, the code determines what selection criteria apply to
the data and uses that to obtain the appropriate parameter settings for the processing
steps.  Applying the selection to select the appropriate parameter file simply requires
matching up the key in the JSON file with the selection information. For example,
a **filter product** would end up using the **filter_basic** criteria, while an
8 exposure ACS/WFC association would end up selecting the **acs_wfc_any_n6** entry.


User-customization of Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The parameter configuration file now included in the ``drizzlepac`` package are
designed to be easily customized for manual processing with both ``runastrodriz`` 
(pipeline astrometry processing) and ``runsinglehap`` (SVM processing).  These 
ASCII JSON files can be edited prior to manual reprocessing to include whatever
custom settings would best suit the science needs of the research being performed
with the data.  


Defining the Output WCS
-------------------------
The SVM processing steps through the **product list** to generate each of the 
pre-defined products one at a time after the input exposures have all been 
aligned.  One of the primary goals of SVM processing is to produce combined
images which share the same WCS for all the data from the same detector.  This 
simply requires defining a common WCS which can be used to define the output for
all the **filter products** from the visit.  

The common WCS, or **metawcs**, gets defined by reading in all the WCS definitions
as ``stwcs.wcsutil.HSTWCS`` objects
for all the input exposures taken with the same **instrument** in the visit.  This
list of **HSTWCS** objects then gets fed to ``stwcs.distortion.utils.output_wcs``,
the same function used by ``AstroDrizzle`` to define the default output WCS when
the user does not specify one before-hand.  This results in the definition of a
WCS which spans the entire field-of-view for all the input exposures with the same
plate scale and orientation as the first **HSTWCS** in the input list.  This **metawcs**
then gets used to define the shape, size and WCS pointing for all drizzle products
taken with the same detector in the visit.  


Drizzling
-----------
Each output product gets created using ``AstroDrizzle``.  This step:

  * combines all the input exposures associated with the product 
  * uses the parameters read in from the configuration files 
  * defines the output image using the **metawcs** WCS definition
  * writes out a multi-extension FITS (MEF) file for the drizzled image using
    the pre-defined name 

This drizzled output image has the same structure as the standard pipeline drizzle
products; namely,

  * PRIMARY extension:  all information common to the product such as 
    instrument and detector.
  * SCI extension: the drizzled science image along with header keywords
    describing the combined array such as total exposure time.
  * WHT extension: an array reporting the drizzled weight for each pixel
  * CON extension: an array reporting what input exposures contributed to each output pixel
  
The headers of each extension gets defined as using the ``fitsblender`` software with 
much the same rules used to create the standard pipeline drizzle product headers.  
In short, it uses simple rules files to determine what keywords should be kept in
the output headers from all the input exposures, and how to select or compute the
value from all the input headers for each keyword.  

Unique SVM Keywords
^^^^^^^^^^^^^^^^^^^^^^
A small set of keywords have been added to the standard drizzle headers to reflect
the unique characteristics of the SVM products.  These keywords are:

.. glossary::

  NPIXFRAC
    Fraction of pixels with data
    
  MEANEXPT
    Mean exposure time per pixel with data
  
  MEDEXPT
    Median exposure time per pixel with data
  
  MEANNEXP
    Mean number of exposures per pixel with data
  
  MEDNEXP
    Median number of exposures per pixel with data



Catalog Generation
-------------------
SVM processing does not stop with the creation of the output drizzled images like
the standard calibration pipeline.  Instead, it derives 2 separate source catalogs
from each drizzled **filter product** to provide a standardized measure of each
visit. For more details on how the catalogs are produced, please refer to the :ref:`catalog_generation` documentation
page.



