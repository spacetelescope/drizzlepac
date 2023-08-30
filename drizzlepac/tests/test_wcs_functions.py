import os
import pytest
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
from drizzlepac.wcs_functions import make_mosaic_wcs

abs_path = 'drizzlepac/tests/'

# checks if function works with single filenames and lists of filenames
# does not check if wcs is different, but this was confirmed independently. 
@pytest.mark.parametrize(
    "filename_array",
    [
        # (abs_path+"sample_drc_header.fits"), # doesn't work for drc files yet
        (abs_path+"sample_flc_header.fits"),
        ([abs_path+"sample_flc_header.fits", abs_path+"sample_flc_header.fits"]),
        pytest.param([abs_path+"sample_flc_header.fits", abs_path+"sample_drc_header.fits"], marks=pytest.mark.xfail),
        pytest.param(["bad_filename.fits"], marks=pytest.mark.xfail),
    ],
)
def test_make_mosaic_wcs(filename_array):
    new_wcs = make_mosaic_wcs(filename_array)
