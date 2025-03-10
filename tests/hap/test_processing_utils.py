import os
import pytest
import shutil
from astropy.io import fits
from drizzlepac.haputils import processing_utils

test_path = "tests/hap/"

@pytest.mark.parametrize(
    "filename",
    [
        "sample_ipppssoot_flc.fits", # WFC3/UVIS
        "sample_ipppssoot_flt.fits", # WFC3/IR
        "sample_svm_flc.fits", # WFC3/UVIS SVM file
        "sample_svm_flt.fits", # WFC3/IR SVM file
        "sample_ipppssoot_flt_w_skycell_header_keyword.fits"
    ],
)
def test_add_skycell_to_header(filename, tmpdir):
    temp_path = os.path.join(tmpdir, f"temp_{filename}")
    # creating copy as to not alter the sample file. 
    shutil.copy2(test_path+filename, temp_path)
    processing_utils.add_skycell_to_header(temp_path)
    hdu = fits.open(temp_path)
    assert hdu[1].header["SKYCELL"] == "p0121x12y16"

def test_add_skycell_to_header_invalid_filename():
    with pytest.raises(Exception):
        processing_utils.add_skycell_to_header('invalid_filename.fits')

