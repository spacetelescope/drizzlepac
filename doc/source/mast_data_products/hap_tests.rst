Regression Tests
================
.. _regression-tests:


**SVM**

Several HAP specific tests do various checks on different datasets included in the test file name (e.g. test_svm_wfc3ir.py). 
The four groups of tests used are as follows:

A. Simple tests:
   
    test_svm_samewcs
     * Checks that products for both detectors are aligned to the same catalog.

    test_svm_wcs
     * Checks the output primary WCSNAME includes FIT_SVM_GAIA.

    test_svm_manifest_name
     * Ensures the manifest file is created.

   This also includes the general success of SVM alignment and catalog creation code.


B. SVM catalog tests (currently disabled)
    
    test_svm_point_total_cat
     * Tests the number of sources in the point source catalog compared to expected value. 
    
    test_svm_segment_total_cat
     * Tests the number of sources in the segmentation source catalog compared to expected value. 


C. Tests for different detectors
    
    test_svm_wcs_ir, test_svm_wcs_ir_all, test_svm_wcs_uvis, test_svm_wcs_uvis_all
     * SVM alignment and catalog creation for different detectors.


D. Mean magnitude tests

    test_svm_point_cat_meanmag, test_svm_segment_cat_meanmag: 
     * Checks that catalog mean magnitudes are within 0.5% of previous values.


test_svm_je281u.py: A

test_svm_hrcsbc.py: A, B

test_svm_wfc3ir.py: A, B

test_svm_j97e06.py: A, D

test_svm_ibqk07.py: A, B, C, D

test_svm_ibyt50.py: A, B, C

|

**HAP**

test_processing_utils.py

    test_add_skycell_to_header
     * Unit test for function for adding skycell name to SVM headers.

test_pipeline.py

    test_astrometric_singleton
     * A test of runastrodriz.process with varying environment setups. Obtains data using astroquery_utils.retrieve_observations.

test_apriori.py
    
    Tests alignment of all of the available a priori wcs solutions for two datasets for ACS and WFC3. 

test_align.py
    
    Tests alignment of all of the available a posteriori wcs solutions for a variety of datasets and scenarios.

archival_test_run_svmpoller.py (*currently disabled*)

    test_run_svmpoller
     * Tests runsinglehap.perform using poller file as input. 

archival_test_randomlist.py (*currently disabled*)

    test_randomlist
     * Tests SVM alignment (align.perform_align) on a random dataset from "ACSWFC3ListDefault50.csv". Success is marked by a statistical sample (70%) of ACS and WFC3 datasets aligned to within 10mas RMS.

archival_test_alignpipe_randomlist.py (*currently disabled*)

    test_alignpipe_randomlist
     * Similar to test_randomlist but include pipeline processing (runastrodriz.process). 

|

**ACS**

test_acs_narrowband.py
    
    test_acs_narrowband
     * Tests relative fit AstroDrizzle on narrowband association.

test_unit.py
    
    Test do_driz square kernel with point.

test_acs_tweak.py
    
    test_tweak
     * Tests tweakreg and then AstroDrizzle.

    test_pixsky1
     * Tests pixtosky, pixtopix, skytopix on ACS data.

test_acs_kernels.py
    
    test_kernels
     * Tests AstroDrizzle on ACS file over different final combined image kernels.

test_asn_regress.py
    
    test_hrc_asn
     * Relative fit AstroDrizzle of ACS HRC dataset.

|

**WFPC2**

test_wfpc2.py

    test_waiver_single
     * Tests WFPC2 Astrodrizzle association of 1 dataset (c01 and c1f files).

    test_waiver_asn
     * Tests WFPC2 Astrodrizzle association of multiple datasets.

    test_wfpc2_single
     * Tests WFPC2 Astrodrizzle with c01 and c1f with filenames as inputs.

    test_mef_asn
     * Tests WFPC2 Astrodrizzle with data in multi-extension fits file format.

|

**WFC3**

test_wfc3.py

    test_binned_single
     * Tests pipeline processing of WFC3 data with the parameter skysub=False.

    test_uvis_single
     * Tests pipeline processing of a single WFC3/UVIS dataset.

    test_uvis_asn
     * Tests pipeline processing of a WFC3/UVIS association (relative fitting).

    test_wfc3_ir_saturated
     * Tests pipeline processing of a saturated WFC3/IR visit.

test_vary_perf.py

    test_perf
     * Tests Astrodrizzle performence using different numbers of cores.

|

**STIS**

test_stis.py

    test_fuv_mama

    test_nuv_mama
     * Tests for a correctly applied distortion model for STIS NUV MAMA data and the creation of a combined product using AstroDrizzle. 

    test_stis_ccd
     * The same as test_nuv_mama but using CCD data. 

    test_stis_oiii_ccd
     * The same as test_nuv_mama but with STIS F28x50OIII CCD data. 


|

**drizzle algorithm**

test_cdriz.py

    Tests drizzling algorithm for different kernels in small square.

test_kernel.py

    Tests drizzling algorithm for different kernels in larger square.
