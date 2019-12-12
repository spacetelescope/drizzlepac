""" This module runs tests on specific datasets to determine if the datasets can
    be aligned to an astrometric catalog. """
import pytest
from drizzlepac import align as alignimages
from ci_watson.artifactory_helpers import get_bigdata
from .base_test import BaseHLATest

# Nominal acceptable RMS limit for a good solution (IMPROVE THIS)
RMS_LIMIT = 10.0

@pytest.mark.bigdata
class TestAlignMosaic(BaseHLATest):
    """ Tests which validate whether mosaics can be aligned to an astrometric standard.

        Characteristics of these tests:
          * A reference WCS is generated based upon all the input images for
            the field.
          * A source astrometric catalog is created using the Photutils
            package to detect explicitly sources in the images.
          * An astrometric catalog is created to extract astrometric positions
            for the found sources in the input images' field-of-view using
            GAIADR2 (preferred) or GAIADR1.
          * Cross matching/fitting is done between found sources and catalog
            coordinates with the Tweakwcs package.
          * The quality of the fit is evaluated against a minimum threshold and
            potentially another fit algorithm is invoked or an alternative
            catalog is used in an effort to obtain a better quality fit.
          * If the option is set, the WCS information is updated for the
            input exposures. The default is False.
          * No mosaic is generated.
          * An output table containing characterizations of the process and
            associated fit is generated.

        Success Criteria:
          * Success criteria hard-coded for this test represents 10mas RMS for the ACS and
            WFC3 images based on the fit of source positions to the astrometric catalog
            source positions.
              * RMS values are extracted from the table output from `perform_align`

        The environment variable needs to be set in the following manner:
            export TEST_BIGDATA=https://bytesalad.stsci.edu/artifactory/
            OR
            export TEST_BIGDATA=/Users/YourNameHere/TestDataDirectory/

        For this test, the TEST_BIGDATA defines the root of the location where the data
        is stored.  The full path is TEST_BIGDATA plus the path components provided
        in the get_bigdata() invocation.  This test file can be executed in the following
        manner:
            $ pytest -s --bigdata test_align.py >& test_align_output.txt &
            $ tail -f test_align_output.txt

    """

    @pytest.mark.slow
    def test_align_ngc188(self):
        """ Verify whether NGC188 exposures can be aligned to an astrometric standard.

        Characteristics of this test:
            of NGC188 suitable for creating a combined mosaic using both instruments.
        """
        total_rms = 0.0
        input_filenames = ['iaal01hxq_flc.fits', 'iaala3btq_flc.fits',
                           'iaal01hyq_flc.fits', 'iaala3bsq_flc.fits',
                           'j8boa1m8q_flc.fits', 'j8boa1m4q_flc.fits',
                           'j8boa1maq_flc.fits', 'j8boa1m6q_flc.fits']

        # Since these are full file names (*_flc.fits) which cannot be obtained via astroquery from
        # MAST, get the data now using ci_watson.
        for input_file in input_filenames:
            get_bigdata('hst-hla-pipeline', 'dev', 'mosaic_ngc188', input_file)

        dataset_table = alignimages.perform_align(input_filenames, archive=False, clobber=False,
                                                  debug=False, update_hdr_wcs=False,
                                                  print_fit_parameters=True, print_git_info=False,
                                                  output=False)

        # Examine the output table to extract the RMS for the entire fit and the compromised
        # information
        if dataset_table:
            total_rms = dataset_table['total_rms'][0]

        assert 0.0 < total_rms <= RMS_LIMIT

    @pytest.mark.slow
    def test_align_47tuc(self):
        """ Verify whether 47Tuc exposures can be aligned to an astrometric standard.

        Characteristics of this test:
          * Input exposures include both ACS and WFC3 images of the same general field-of-view
            of 47Tuc suitable for creating a combined mosaic using both instruments.
        """
        total_rms = 0.0
        input_filenames = ['ib6v06c4q_flc.fits', 'ib6v06c7q_flc.fits',
                           'ib6v25aqq_flc.fits', 'ib6v25atq_flc.fits',
                           'jddh02gjq_flc.fits', 'jddh02glq_flc.fits',
                           'jddh02goq_flc.fits']

        # Since these are full file names (*_flc.fits) which cannot be obtained via astroquery from
        # MAST, get the data now using ci_watson.
        for input_file in input_filenames:
            get_bigdata('hst-hla-pipeline', 'dev', 'mosaic_47tuc', input_file)

        dataset_table = alignimages.perform_align(input_filenames, archive=False, clobber=False,
                                                  debug=False, update_hdr_wcs=False,
                                                  print_fit_parameters=True, print_git_info=False,
                                                  output=False)

        # Examine the output table to extract the RMS for the entire fit and the compromised
        # information
        if dataset_table:
            total_rms = dataset_table['total_rms'][0]

        assert 0.0 < total_rms <= RMS_LIMIT


    @pytest.mark.parametrize("input_filenames", [['j8ura1j1q_flt.fits', 'j8ura1j2q_flt.fits',
                                                  'j8ura1j4q_flt.fits', 'j8ura1j6q_flt.fits',
                                                  'j8ura1j7q_flt.fits', 'j8ura1j8q_flt.fits',
                                                  'j8ura1j9q_flt.fits', 'j8ura1jaq_flt.fits',
                                                  'j8ura1jbq_flt.fits', 'j8ura1jcq_flt.fits',
                                                  'j8ura1jdq_flt.fits', 'j8ura1jeq_flt.fits',
                                                  'j8ura1jfq_flt.fits', 'j8ura1jgq_flt.fits',
                                                  'j8ura1jhq_flt.fits', 'j8ura1jiq_flt.fits',
                                                  'j8ura1jjq_flt.fits', 'j8ura1jkq_flt.fits'],
                                                 ['j92c01b4q_flc.fits', 'j92c01b5q_flc.fits',
                                                  'j92c01b7q_flc.fits', 'j92c01b9q_flc.fits'],
                                                 ['jbqf02gzq_flc.fits', 'jbqf02h5q_flc.fits',
                                                  'jbqf02h7q_flc.fits', 'jbqf02hdq_flc.fits',
                                                  'jbqf02hjq_flc.fits', 'jbqf02hoq_flc.fits',
                                                  'jbqf02hqq_flc.fits', 'jbqf02hxq_flc.fits',
                                                  'jbqf02i3q_flc.fits', 'jbqf02i8q_flc.fits',
                                                  'jbqf02iaq_flc.fits'],
                                                 ['ib2u12kaq_flt.fits', 'ib2u12keq_flt.fits',
                                                  'ib2u12kiq_flt.fits', 'ib2u12klq_flt.fits'],
                                                 ['ibnh02coq_flc.fits', 'ibnh02cmq_flc.fits',
                                                  'ibnh02c7q_flc.fits', 'ibnh02c5q_flc.fits',
                                                  'ibnh02cpq_flc.fits', 'ibnh02c9q_flc.fits',
                                                  'ibnh02bfq_flc.fits', 'ibnh02beq_flc.fits']])
    @pytest.mark.slow
    def test_align_single_visits(self, input_filenames):
        """ Verify whether single-visit exposures can be aligned to an astrometric standard.

        Characteristics of these tests:
          * Input exposures include exposures from a number of single visit datasets to explore what
            impact differing observing modes (differing instruments, detectors, filters, subarray
            size, etc.) have on astrometry.

        The following datasets are used in these tests:

            * ACS dataset 10048_a1: 2x F344N, 1x F435W, 1x F475W, 2x F502N, 2x F550M, 1x F555W,
              1x F606W, 1x F625W, 2x F658N, 1x F775W, 1x F814W, 1x F850LP, and
              2x F892N ACS/HRC images
            * ACS dataset 10265_01: 4x F606W full-frame ACS/WFC images
            * ACS dataset 12580_02: 5x F475W & 6x F814W ACS/WFC images
            * WFC3 dataset 11663_12: 4x F160W full-frame WFC3/IR images
            * WFC3 dataset 12379_02: 4X F606W, 4x F502N full-frame WFC3/UVIS images

        """
        total_rms = 0.0

        # Since these are full file names (*_flc.fits) which cannot be obtained via astroquery from
        # MAST, get the data now using ci_watson.
        for input_file in input_filenames:
            get_bigdata('hst-hla-pipeline', 'dev', 'base_tests', input_file)

        dataset_table = alignimages.perform_align(input_filenames, archive=False, clobber=False,
                                                  debug=False, update_hdr_wcs=False,
                                                  print_fit_parameters=True, print_git_info=False,
                                                  output=False)

        # Examine the output table to extract the RMS for the entire fit and the compromised
        # information
        if dataset_table:
            total_rms = dataset_table['total_rms'][0]

        assert 0.0 < total_rms <= RMS_LIMIT

    @pytest.mark.xfail
    @pytest.mark.slow
    def test_align_fail_single_visit(self):
        """ Verify whether single-visit exposures can be aligned to an astrometric standard.

        Characteristics of this test:
          * Input exposures include exposures from a number of single visit datasets to explore what
            impact differing observing modes (differing instruments, detectors, filters, subarray
            size, etc.) have on astrometry. This test is known to fail due to "RuntimeError: Number
            of output coordinates exceeded allocation (475)a". It will exercise the code using both
            catalogs for each of the three fitting algorithms at this time. Nans will be present
            in the output table.


        The following datasets are used in these tests:
            * WFC3 dataset 12219_01: 8x F160W full-frame WFC3/IR images, 9x F336W full-frame
              WFC3/UVIS images

        """
        total_rms = 0.0
        input_filenames = ['ibjt01a1q_flc.fits', 'ibjt01a8q_flc.fits',
                           'ibjt01aiq_flt.fits', 'ibjt01amq_flt.fits',
                           'ibjt01aqq_flt.fits', 'ibjt01auq_flt.fits',
                           'ibjt01yqq_flc.fits', 'ibjt01z0q_flc.fits',
                           'ibjt01zwq_flc.fits', 'ibjt01a4q_flc.fits',
                           'ibjt01acq_flc.fits', 'ibjt01akq_flt.fits',
                           'ibjt01aoq_flt.fits', 'ibjt01asq_flt.fits',
                           'ibjt01avq_flt.fits', 'ibjt01yuq_flc.fits',
                           'ibjt01ztq_flc.fits'],

        # Since these are full file names (*_flc.fits) which cannot be obtained via astroquery from
        # MAST, get the data now using ci_watson.
        for input_file in input_filenames:
            get_bigdata('hst-hla-pipeline', 'dev', 'base_tests', input_file)

        dataset_table = alignimages.perform_align(input_filenames, archive=False, clobber=False,
                                                  debug=False, update_hdr_wcs=False,
                                                  print_fit_parameters=True, print_git_info=False,
                                                  output=False)

        # Examine the output table to extract the RMS for the entire fit and the compromised
        # information
        if dataset_table:
            total_rms = dataset_table['total_rms'][0]

        assert 0.0 < total_rms <= RMS_LIMIT

    def test_astroquery(self):
        """Verify that new astroquery interface will work"""

        total_rms = 0.0

        dataset_table = alignimages.perform_align(['IB6V06060'], archive=False, clobber=True,
                                                  debug=False, update_hdr_wcs=False,
                                                  print_fit_parameters=True, print_git_info=False,
                                                  output=False)

        # Examine the output table to extract the RMS for the entire fit and the compromised
        # information
        if dataset_table:
            total_rms = dataset_table['total_rms'][0]

        assert 0.0 < total_rms <= RMS_LIMIT
