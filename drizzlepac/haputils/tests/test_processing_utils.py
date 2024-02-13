import os
import pytest
import shutil
import tempfile
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
from drizzlepac.haputils import processing_utils

abs_path = "drizzlepac/haputils/tests/"


@pytest.mark.parametrize(
    "filename",
    [
        "sample_ipsoot_flc.fits",
        "sample_ipsoot_flt.fits",
        "sample_svm_flc.fits",
        "sample_svm_flt.fits",
    ],
)
def test_add_skycell_to_header(filename, tmpdir):
    temp_path = os.path.join(tmpdir, f"temp_{filename}")
    shutil.copy2(abs_path+filename, temp_path)
    processing_utils.add_skycell_to_header(temp_path)
    hdu = fits.open(temp_path)
    assert hdu[1].header["SKYCELL"] == "p0121x12y16"
