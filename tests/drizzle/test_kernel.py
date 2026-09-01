import pytest
import numpy as np
import cdriz_setup


@pytest.fixture
def kernel_pars(request):
    kernel = request.getfixturevalue("kernel")

    # The point (nearest-neighbour) kernel places each input pixel's flux in
    # a single output pixel. With the default 50x60 input the input and
    # output grids differ in parity along x (50 vs 51), so a half-pixel crpix
    # offset maps the source exactly onto an output pixel boundary; which
    # side it lands on is then decided by sub-ULP round-off in the WCS
    # transform, which varies between WCSLIB versions.

    # Apply a small irrational offset to the output WCS for the "point" kernel
    # to avoid input pixel falling on an output pixel boundary.

    offset = -1.0e-5 * np.pi * np.ones(2, dtype=float) if kernel == "point" else None
    _params = cdriz_setup.Get_Grid(inx=50, iny=60, outx=51, outy=66, offset=offset, background=0.0)
    return _params


@pytest.mark.parametrize("kernel", ["square", "point", "turbo", "gaussian", "lanczos3"])
def test_point_data_kernel(kernel_pars, kernel):
    """Function tests different c code kernels (inputs already created on instantiation).

    Parameters
    ----------
    kernel : str
        String associated with one of the c code point kernel options.
    kernel_pars : Class
        The Class initialized in Get_Grid which includes all of the inputs needed to run cdriz.tdriz.
    """

    # truth filename
    output_name = f"{kernel}_truth"
    relative_path = "truth_files"
    output_fullpath = cdriz_setup.get_output_fullpath(relative_path, output_name)

    # add missing/flagged pixels in inwht
    kernel_pars.insci[20:22, 21:22] = 100

    # resample:
    cdriz_setup.cdriz_call(kernel_pars, kernel)

    # save truth file as figure
    cdriz_setup.generate_png(kernel_pars, f"{output_fullpath}.png")

    try:
        truth_array = np.genfromtxt(f"{output_fullpath}.csv", delimiter=",")
    except:
        cdriz_setup.save_array(kernel_pars.outsci, f"{output_fullpath}.csv")
        truth_array = kernel_pars.outsci.copy()

    assert np.allclose(
        kernel_pars.outsci,
        truth_array,
        atol=1e-4,
        rtol=1e-5,
    ), cdriz_setup.error_message(kernel_pars.outsci, f"{output_fullpath}_new.csv")

@pytest.mark.parametrize("kernel", ["gaussian"])
def test_cdriz_edge(kernel_pars, kernel):
    """Similar to test_point_kernel but looking at bright pixels at edge of field."""

    output_name = f"edge_{kernel}_truth"
    relative_path = "truth_files"
    output_fullpath = cdriz_setup.get_output_fullpath(relative_path, output_name)

    kernel_pars.insci[0, 21] = 100
    cdriz_setup.cdriz_call(kernel_pars, kernel)

    cdriz_setup.generate_png(kernel_pars, f"{output_fullpath}.png")

    try:
        truth_array = np.genfromtxt(f"{output_fullpath}.csv", delimiter=",")
    except:
        cdriz_setup.save_array(kernel_pars.outsci, f"{output_fullpath}.csv")
    assert np.allclose(
        kernel_pars.outsci, truth_array, atol=1e-4
    ), cdriz_setup.error_message(kernel_pars.outsci, f"{output_fullpath}_new.csv")


@pytest.mark.parametrize("kernel", ["gaussian"])
def test_cdriz_large(kernel_pars, kernel):
    """Similar to test_point_kernel but looking at large pixel."""

    output_name = f"large_square_{kernel}_truth"
    relative_path = "truth_files"
    output_fullpath = cdriz_setup.get_output_fullpath(relative_path, output_name)

    kernel_pars.insci[21:25, 22:26] = 100
    cdriz_setup.cdriz_call(kernel_pars, kernel)

    cdriz_setup.generate_png(kernel_pars, f"{output_fullpath}.png")

    try:
        truth_array = np.genfromtxt(f"{output_fullpath}.csv", delimiter=",")
    except:
        cdriz_setup.save_array(kernel_pars.outsci, f"{output_fullpath}.csv")
    assert np.allclose(
        kernel_pars.outsci, truth_array, atol=1e-4
    ), cdriz_setup.error_message(kernel_pars.outsci, f"{output_fullpath}_new.csv")


@pytest.mark.parametrize("kernel", ["gaussian"])
def test_cdriz_non_symmetrical(kernel_pars, kernel):
    """Similar to test_point_kernel but looking at non-symmetrical pixel."""

    output_name = f"nonsymmetrical_{kernel}_truth"
    relative_path = "truth_files"
    output_fullpath = cdriz_setup.get_output_fullpath(relative_path, output_name)

    kernel_pars.insci[21:25, 22:23] = 100
    cdriz_setup.cdriz_call(kernel_pars, kernel)

    cdriz_setup.generate_png(kernel_pars, f"{output_fullpath}.png")

    try:
        truth_array = np.genfromtxt(f"{output_fullpath}.csv", delimiter=",")
    except:
        cdriz_setup.save_array(kernel_pars.outsci, f"{output_fullpath}.csv")

    assert np.allclose(
        kernel_pars.outsci, truth_array, atol=1e-4
    ), cdriz_setup.error_message(kernel_pars.outsci, f"{output_fullpath}_new.csv")


@pytest.mark.parametrize("kernel", ["square", "point", "turbo", "gaussian", "lanczos3"])
def test_zero_input_weight(kernel_pars, kernel):
    """Tests that do_driz ignores bad pixels, or those that have an input weight (inwht) of 0.

    Parameters
    ----------
    kernel : str
        String associated with one of the c code point kernel options.
    kernel_pars : Class
        The Class initialized in Get_Class which includes all of the inputs need to run cdriz.tdriz.
    """
    # add bad bright pixels in insci
    kernel_pars.insci[0:4, 0:4] = 1e8
    kernel_pars.inwht[0:4, 0:4] = 0

    # adding two additional bright "sources"
    kernel_pars.insci[6, 7] = 1000
    kernel_pars.insci[9, 6] = 1000

    # resample
    cdriz_setup.cdriz_call(kernel_pars, kernel)

    # check that original bad pixel flux still coming through:
    assert np.sum(kernel_pars.outsci) > 1e8
