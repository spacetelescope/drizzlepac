""" This module is the high-level wrapper for running a test on a list of input
    datasets in order to collect statistics on alignment of each dataset to an
    astrometric catalog."""
import pytest
from .process_randomlist import TestAlignMosaic

@pytest.mark.bigdata
@pytest.mark.xfail
@pytest.mark.slow
@pytest.mark.unit
def test_randomlist(start_row, num_rows):
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
          * Success criterion hard-coded for this test represents whether a
            statistical sample (70%) of ACS and WFC3 datasets were able to be
            aligned to within 10mas RMS.
              * RMS values are extracted from the table output from `perform_align`

        The environment variable needs to be set in the following manner:
            export TEST_BIGDATA=https://bytesalad.stsci.edu/artifactory/
            OR
            export TEST_BIGDATA=/Users/YourNameHere/TestDataDirectory/

        For this test, the TEST_BIGDATA defines the root of the location where
        the CSV file is stored.  The full path is
        TEST_BIGDATA/self.input_repo/self.tree/self.input_loc.  The CSV is
        output from a database and lists associations and singletons for ACS
        and WFC3 instruments randomly sorted.  The actual data files are
        downloaded from MAST via astroquery.

        This test file can be executed in the following manner:
            $ pytest -s --bigdata test_randomlist.py >& test_random_output.txt &
            $ tail -f test_random_output.txt
    """
    run_random_test = TestAlignMosaic()

    assert run_random_test.test_align_randomfields(start_row, num_rows) > 0.7
