""" This module evaluates the Concentration Index in a Point source
    catalog to ensure the value stays close to a nominal, CI_REF,
    value. """
import glob
import numpy as np
import os
import pytest

from astropy.io import ascii
from astropy.utils.data import conf

from ci_watson.artifactory_helpers import get_bigdata
from ci_watson.hst_helpers import download_crds, ref_from_image

from drizzlepac.haputils import astroquery_utils as aqutils
from drizzlepac import hapsequencer, runastrodriz

# Nominal acceptable CI value for a good solution
CI_REF = 1.209
CI_LIMIT = 0.05

"""
    This module processes a dataset with both runastrodriz and runsinglehap
    to generate a Point source catalog.  The Point catalog is checked,
    in particular, as it is a diagnostic for degradations in the
    alignment (Photutils deprecations associated with IterativePSFPhotometry)
    that may cause the Concentration Index (CI) to grow.  As the
    IterativePSFPhotometry class is used for both runastrodriz and
    runsinglehap, it is cleanest to use freshly made FLTs for this
    test which is why calwf3 is invoked.

    This test file can be executed in the following manner:
        $ pytest -s test_svm_ci.py --basetemp=/Users/yourNameHere/PYTEST >& test_svm_ci.log &
        $ tail -f test_svm_ci.log

        The actual testing and files can be found in /Users/yourNameHere/PYTEST.
"""

class BaseWFC3Pipeline:
    prevdir = os.getcwd()
    use_ftp_crds = True
    timeout = 30  # seconds
    tree = 'dev'

    # To be defined by instrument
    refstr = 'iref'
    input_loc = ''
    ref_loc = 'wfc3/ref'
    prevref = os.environ.get(refstr)

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

    def get_input_file(self, *args, refsep='$', docopy=True):
        """
        Download or copy input file (e.g., RAW) into the working directory.
        The associated CRDS reference files in ``refstr`` are also
        downloaded, if necessary.
        """
        filename = args[1]
        ref_files = ref_from_image(filename, ['IDCTAB', 'OFFTAB', 'NPOLFILE', 'D2IMFILE',
                                              'DGEOFILE', 'MDRIZTAB'])
        print("Looking for REF_FILES: {}".format(ref_files))

        for ref_file in ref_files:
            if ref_file.strip() == '':
                continue
            refname = os.path.join(ref_file)
            if self.use_ftp_crds:
                if isinstance(refname, str):
                    download_crds(refname)
                else:
                    raise TypeError("Expected string, got {}".format(type(refname)))

        return filename


class TestSVMCI(BaseWFC3Pipeline):

    def test_svm_ci(self):
        """ Tests pipeline-style processing of a ASN using runastrodriz.
        """
        # Get test data through astroquery - retrieve the ASN and FLT files
        _ = aqutils.retrieve_observation("idxo15030", suffix=["ASN", "FLT"], product_type="pipeline")

        # Generate the list of FLT files
        inputFiles = glob.glob("idxo*_flt.fits") 

        # Retrieve reference files for these as well
        self.get_input_file('', inputFiles[0], docopy=False)

        # Insure environment variables are set for full processing
        os.environ['ASTROMETRY_STEP_CONTROL'] = 'on'
        os.environ['ASTROMETRY_COMPUTE_APOSTERIORI'] = 'on'
        os.environ['ASTROMETRY_APPLY_APRIORI'] = 'on'

        # Perform the standard pipeline drizzle step (runastrodriz)
        try:
            runastrodriz.process("idxo15030_asn.fits")
        except Exception as x_cept:
            print("")
            print("Exception encountered executing runastrodriz on file idxo15030_asn.fits: {}".format(x_cept))
            assert False

        # Write out the files as a list of names
        outputFilename = "idxo_flt.lst"
        with open(outputFilename, "w") as file:
            for name in inputFiles:
                file.write(name + '\n')

        # Run the SVM processing to produce the source catalogs (runsinglehap/hapsequencer) 
        path = os.path.join(os.path.dirname(__file__), outputFilename)

        # Perform the HAP SVM processing
        try:
            returnValue = hapsequencer.run_hap_processing(outputFilename)
        except Exception as x_cept:
            print("")
            print("Exception encountered executing runsinglehap: {}".format(x_cept))
            assert False

        # Examine the output Point ECSV file and compute the median of all the CI values
        if not returnValue: 
            srcTab = ascii.read("hst_15652_15_wfc3_ir_total_idxo15_point-cat.ecsv", format="ecsv")
            median = np.median(srcTab["CI_f160w"].data)
            ci_diff = abs(median - CI_REF)
        else:
            assert False

        assert ci_diff <= CI_LIMIT

