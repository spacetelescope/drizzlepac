""" This module evaluates the Concentration Index in a Point source
    catalog to ensure the value stays close to a nominal, CI_REF,
    value. """
from astropy.io import ascii
from drizzlepac import hapsequencer, runastrodriz
from drizzlepac.haputils import astroquery_utils as aqutils
import glob
import numpy as np
import os
import subprocess

# Nominal acceptable CI value for a good solution
CI_REF = 1.209
CI_LIMIT = 0.035

"""
    This module processes a dataset through the full standard pipeline
    (calwf3 + runastrodriz) plus HAP SVM processing (runsinglehap) to
    generate a Point source catalog.  The Point catalog is checked,
    in particular, as it is a diagnostic for degradations in the
    alignment (Photutils deprecations associated with IterativePSFPhotometry)
    that may cause the Concentration Index (CI) to grow.  As the
    IterativePSFPhotometry class is used for both runastrodriz and
    runsinglehap, it is cleanest to use freshly made FLTs for this
    test which is why calwf3 is invoked.

    This test file can be executed in the following manner:
        $ pytest -s test_svm_ci.py --basetemp-/Users/yourNameHere/PYTEST >& test_svm_ci.log &
        $ tail -f test_svm_ci.log

        The actual testing and files can be found in /Users/yourNameHere/PYTEST.
        
        Note: The "iref" environment variable needs to be set.
        $ export iref="/grp/crds/cache/references/hst/"
"""

def test_svm_ci(tmp_path_factory):
    """Verify the CI is within an acceptable limit.  """

    # Create working directory specified for the test
    curdir = tmp_path_factory.mktemp(os.path.basename(__file__)) 
    os.chdir(curdir)

    # Get test data through astroquery - retrieve the ASN and RAW files
    _ = aqutils.retrieve_observation("idxo15030", suffix=["ASN", "RAW"], product_type="pipeline")

    # Process the data through CALWF3 to get fresh FLT files
    try:
        subprocess.call(["calwf3.e", "idxo15030_asn.fits", "-vt"])
    except Exception as x_cept:
        print("")
        print("Exception encountered executing calwf3 on file idxo15030_asn.fits: {}.".format(x_cept))
        assert False

    # Perform the standard pipeline drizzle step (runastrodriz)
    try:
        runastrodriz.process("idxo15030_asn.fits")
    except Exception as x_cept:
        print("")
        print("Exception encountered executing runastrodriz on file idxo15030_asn.fits: {}".format(x_cept))
        assert False
    
    # Generate the list of FLT files
    inputFiles = glob.glob("idxo*_flt.fits") 
    
    # Write out the files as a list of names
    outputFilename = "idxo_flt.lst"
    with open(outputFilename, "w") as file:
        for name in inputFiles:
            file.write(name + '\n')

    # Run the SVM processing to produce the source catalogs (runsinglehap/hapsequencer) 
    path = os.path.join(os.path.dirname(__file__), outputFilename)
    returnValue = hapsequencer.run_hap_processing(outputFilename)

    # Examine the output Point ECSV file and compute the median of all the CI values
    if not returnValue: 
        srcTab = ascii.read("hst_15652_15_wfc3_ir_total_idxo15_point-cat.ecsv", format="ecsv")
        median = np.median(srcTab["CI_f160w"].data)
        ci_diff = abs(median - CI_REF)
    else:
        assert False

    assert ci_diff <= CI_LIMIT
