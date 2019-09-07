""" This module tests full pipeline use of the drizzlepac package.

"""
import os
import pytest

from astropy.utils.data import conf

from ci_watson.artifactory_helpers import get_bigdata
from ci_watson.hst_helpers import download_crds, ref_from_image

from drizzlepac.hlautils import astroquery_utils as aqutils
from drizzlepac import runastrodriz
from astropy.io import fits


class BasePipeline:
    prevdir = os.getcwd()
    use_ftp_crds = True
    timeout = 30  # seconds
    tree = ''

    # Numpy default for allclose comparison
    rtol = 1e-6
    atol = 1e-5

    # To be defined by instrument
    refstr = ''
    prevref = ''
    input_loc = ''
    ref_loc = ''
    ignore_keywords = []

    # To be defined by individual test
    subdir = ''

    @pytest.fixture(autouse=True)
    def setup_class(self, tmpdir, envopt, pytestconfig):
        """
        Run test in own dir so we can keep results separate from
        other tests.
        """
        if not tmpdir.ensure(self.subdir, dir=True):
            p = tmpdir.mkdir(self.subdir).strpath
        else:
            p = tmpdir.join(self.subdir).strpath
        os.chdir(p)

        # NOTE: This could be explicitly controlled using pytest fixture
        #       but too many ways to do the same thing would be confusing.
        #       Refine this logic if using pytest fixture.
        # HSTCAL cannot open remote CRDS on FTP but central storage is okay.
        # So use central storage if available to avoid FTP.
        if self.prevref is None or self.prevref.startswith(('ftp', 'http')):
            os.environ[self.refstr] = p + os.sep
            self.use_ftp_crds = True
        os.environ['TEST_BIGDATA'] = p

        # This controls astropy.io.fits timeout
        conf.remote_timeout = self.timeout

        # Update tree to point to correct environment
        self.tree = envopt

        # Collect pytest configuration values specified in setup.cfg or pytest.ini
        self.inputs_root = pytestconfig.getini('inputs_root')[0]
        self.results_root = pytestconfig.getini('results_root')[0]

    def teardown_class(self):
        """Reset path and variables."""
        conf.reset('remote_timeout')
        os.chdir(self.prevdir)
        if self.use_ftp_crds and self.prevref is not None:
            os.environ[self.refstr] = self.prevref

    def get_data(self, *args, docopy=True):
        """
        Download `filename` into working directory using
        `get_bigdata`.  This will then return the full path to
        the local copy of the file.
        """
        local_file = get_bigdata(*args, docopy=docopy)

        return local_file

    def get_input_file(self, *args, refsep='$', docopy=True):
        """
        Download or copy input file (e.g., RAW) into the working directory.
        The associated CRDS reference files in ``refstr`` are also
        downloaded, if necessary.
        """
        filename = self.get_data(*args, docopy=docopy)
        ref_files = ref_from_image(filename, ['IDCTAB', 'OFFTAB', 'NPOLFILE', 'D2IMFILE',
                                              'DGEOFILE', 'MDRIZTAB'])
        print("Looking for REF_FILES: {}".format(ref_files))

        for ref_file in ref_files:
            if ref_file.strip() == '':
                continue
            if refsep not in ref_file:  # Local file
                refname = self.get_data('customRef', ref_file)
            else:  # Download from FTP, if applicable
                refname = os.path.join(ref_file)
                if self.use_ftp_crds:
                    download_crds(refname, self.timeout)
        return filename


class BaseWFC3Pipeline(BasePipeline):
    refstr = 'iref'
    input_loc = ''
    ref_loc = 'wfc3/ref'
    prevref = os.environ.get(refstr)
    ignore_keywords = ['origin', 'filename', 'date', 'iraf-tlm', 'fitsdate',
                       'upwtim', 'wcscdate', 'upwcsver', 'pywcsver',
                       'history', 'prod_ver', 'rulefile']



class TestSingleton(BaseWFC3Pipeline):

    @pytest.mark.parametrize(
        'dataset_names', ['iaaua1n4q', 'iacs01t4q']
    )

    def test_astrometric_singleton(self, dataset_names):
        """ Tests pipeline-style processing of a singleton exposure using runastrodriz.
        """
        # Get sample data through astroquery
        flcfile = aqutils.retrieve_observation(dataset_names, suffix=['FLC'])[0]
        fltfile = aqutils.retrieve_observation(dataset_names, suffix=['FLT'])[0]
        rawfile = aqutils.retrieve_observation(dataset_names, suffix=['RAW'])[0]

        # Retrieve reference files for these as well
        self.get_input_file('', fltfile, docopy=False)

        # Insure environment variables are set for full processing
        os.environ['ASTROMETRY_STEP_CONTROL'] = 'on'
        os.environ['ASTROMETRY_COMPUTE_APOSTERIORI'] = 'on'
        os.environ['ASTROMETRY_APPLY_APRIORI'] = 'on'

        # Run pipeline processing using
        runastrodriz.process(rawfile, force=True, inmemory=True)

        # compare WCSNAMEs from flt and flc files
        flc_wcsname = fits.getval(flcfile, 'wcsname', ext=1)
        flt_wcsname = fits.getval(fltfile, 'wcsname', ext=1)

        # Perform comparisons:
        #   - WCSNAME values should contain '-' from either a priori or a posteriori solution
        #   - WCSNAME value should be the same for FLT and FLC images
        assert('-' in flc_wcsname)
        assert('-' in flt_wcsname)
        assert(flc_wcsname == flt_wcsname)
