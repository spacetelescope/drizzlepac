import os
import pytest
import shutil
from astropy.io import fits
from drizzlepac.haputils import processing_utils
from unittest import mock

test_path = "tests/hap/"


@pytest.mark.parametrize(
    "filename",
    [
        "sample_ipppssoot_flc.fits",  # WFC3/UVIS
        "sample_ipppssoot_flt.fits",  # WFC3/IR
        "sample_svm_flc.fits",  # WFC3/UVIS SVM file
        "sample_svm_flt.fits",  # WFC3/IR SVM file
        "sample_ipppssoot_flt_w_skycell_header_keyword.fits",
    ],
)
def test_add_skycell_to_header(filename, tmpdir):
    temp_path = os.path.join(tmpdir, f"temp_{filename}")
    # creating copy as to not alter the sample file.
    shutil.copy2(test_path + filename, temp_path)
    processing_utils.add_skycell_to_header(temp_path)
    hdu = fits.open(temp_path)
    assert hdu[1].header["SKYCELL"] == "p0121x12y16"


def test_add_skycell_to_header_invalid_filename():
    """Similar to test_add_skycell_to_header
    except that the filename is invalid.
    """
    with pytest.raises(Exception):
        processing_utils.add_skycell_to_header("invalid_filename.fits")


@mock.patch("drizzlepac.hapmultisequencer.run_mvm_processing")
def test_add_svm_inputs_to_mvm_header(mock_filter_product, tmpdir):
    """This test checks if the function add_svm_inputs_to_mvm_header.
    It checks if the values are being updated in the HDRTABLE of the
    MVM file correctly

    Parameters
    ----------
    mock_filter_product : mock.Mock
        Mock of a single filter_product for the sample_flc.fits file.
    tmpdir : str
        temporary directory for storing files.
    """

    sample_flcs = [
        "sample_svm_flc.fits",
        "sample_mvm_drc.fits",
    ]

    path_to_temp_flcs = []

    # saves the local example files to temprary location as to not alter them
    for filename in sample_flcs:
        temp_path = os.path.join(tmpdir, f"{filename}")
        shutil.copy2(test_path + filename, temp_path)
        path_to_temp_flcs.append(temp_path)

    # path to mvm file that we are updating
    file_to_update = path_to_temp_flcs[-1]

    # mock the filter product
    mock_filter_product = mock.Mock(
        name="add_svm_inputs_to_mvm_header",
        **{
            "drizzle_filename": file_to_update,
            "svm_input_filename": path_to_temp_flcs[0],
        },
    )

    hdu = processing_utils.add_svm_inputs_to_mvm_header(
        mock_filter_product, return_hdu=True
    )

    assert isinstance(hdu[4].data["SVMROOTNAME"][0], str)
    assert isinstance(hdu[4].data["GENDATE"][0], str)
    assert len(hdu[4].data["GENDATE"]) == 12
    assert len(hdu[4].data["GENDATE"][0]) == 10
    assert len(hdu[4].data["SVMROOTNAME"]) == 12
    assert len(hdu[4].data["SVMROOTNAME"][0]) == 128


@mock.patch("drizzlepac.hapmultisequencer.run_mvm_processing")
def test_add_svm_inputs_to_mvm_header_invalid_svm_filename(mock_filter_product, tmpdir):
    """Similar to test_add_svm_inputs_to_mvm_header but with invalid svm_filename.

    Parameters
    ----------
    mock_filter_product : mock.Mock
        Mock of a single filter_product for the sample_flc.fits file.
    tmpdir : str
        temporary directory for storing files.
    """

    sample_flcs = [
        "sample_svm_flc.fits",
        "sample_mvm_drc.fits",
    ]

    path_to_temp_flcs = []

    # saves the local example files to temprary location as to not alter them
    for filename in sample_flcs:
        temp_path = os.path.join(tmpdir, f"{filename}")
        shutil.copy2(test_path + filename, temp_path)
        path_to_temp_flcs.append(temp_path)

    # path to mvm file that we are updating
    file_to_update = path_to_temp_flcs[-1]

    # mock the filter product
    mock_filter_product = mock.Mock(
        name="add_svm_inputs_to_mvm_header",
        **{
            "drizzle_filename": file_to_update,
            "svm_input_filename": 'invalid_svm_filename.fits',
        },
    )

    hdu = processing_utils.add_svm_inputs_to_mvm_header(
        mock_filter_product, return_hdu=True
    )
