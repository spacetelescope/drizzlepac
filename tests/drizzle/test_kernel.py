import pytest
import numpy as np
from astropy import wcs
from drizzlepac import cdriz
import test_setup


@pytest.fixture
def kernel_pars():
    return test_setup.Get_Grid(inx=50, iny=60, outx=51, outy=66)


@pytest.mark.parametrize("kernel", ["square", "point", "turbo", "gaussian", "lanczos3"])
def test_point_kernel(kernel, kernel_pars, new_truth=False, return_png=False):
    """Function tests different c code point kernels (inputs already created on instantiation).

    Parameters
    ----------
    kernel : str
        String associated with one of the c code point kernel options.
    kernel_pars : Class
        The Class inintialized in Get_Grid which includes all of the inputs need to run cdriz.tdriz.
    new_truth : bool, optional
        Flag, whether to save a new csv truth file from the output data, by default False
    return_png: bool, optional
        Flag, whether to create an output image map (png).
    """
    # truth filename
    truth_filename = f"./tests/drizzle/truth_files/{kernel}_truth"

    # add missing/flagged pixels in inwht
    kernel_pars.insci[20:22, 21:22] = 100

    # resample:
    test_setup.cdriz_call(kernel_pars, kernel)

    if new_truth:
        # save new truth file
        test_setup.save_array(kernel_pars.outsci, f"{truth_filename}.csv")

    if return_png:
        # save truth file as figure
        test_setup.generate_png(kernel_pars, f"{truth_filename}.png")

    truth_array = np.genfromtxt(f"{truth_filename}.csv", delimiter=",")
    assert np.allclose(kernel_pars.outsci, truth_array, atol=1e-4)


def test_cdriz_edge(kernel_pars, kernel="gaussian", new_truth=False, return_png=False):
    """Similar to test_point_kernel but looking at bright pixels at edge of field."""

    truth_filename = f"./tests/drizzle/truth_files/edge_{kernel}_truth"
    kernel_pars.insci[0, 21] = 100
    test_setup.cdriz_call(kernel_pars, kernel)
    if new_truth:
        test_setup.save_array(kernel_pars.outsci, f"{truth_filename}.csv")
    if return_png:
        test_setup.generate_png(kernel_pars, f"{truth_filename}.png")

    truth_array = np.genfromtxt(f"{truth_filename}.csv", delimiter=",")
    assert np.allclose(kernel_pars.outsci, truth_array, atol=1e-4)


def test_cdriz_large(kernel_pars, kernel="gaussian", new_truth=False, return_png=False):
    """Similar to test_point_kernel but looking at large pixel."""

    truth_filename = f"./tests/drizzle/truth_files/large_sqaure_{kernel}_truth"
    kernel_pars.insci[21:25, 22:26] = 100
    test_setup.cdriz_call(kernel_pars, kernel)
    if new_truth:
        test_setup.save_array(kernel_pars.outsci, f"{truth_filename}.csv")
    if return_png:
        test_setup.generate_png(kernel_pars, f"{truth_filename}.png")

    truth_array = np.genfromtxt(f"{truth_filename}.csv", delimiter=",")
    assert np.allclose(kernel_pars.outsci, truth_array, atol=1e-4)


def test_cdriz_non_symmetrical(
    kernel_pars, kernel="gaussian", new_truth=False, return_png=False
):
    """Similar to test_point_kernel but looking at non-symmetrical pixel."""

    truth_filename = f"./tests/drizzle/truth_files/nonsymmetrical_{kernel}_truth"
    kernel_pars.insci[21:25, 22:23] = 100
    test_setup.cdriz_call(kernel_pars, kernel)
    if new_truth:
        test_setup.save_array(kernel_pars.outsci, f"{truth_filename}.csv")
    if return_png:
        test_setup.generate_png(kernel_pars, f"{truth_filename}.png")

    truth_array = np.genfromtxt(f"{truth_filename}.csv", delimiter=",")
    assert np.allclose(
        kernel_pars.outsci, truth_array, atol=1e-4
    ), test_setup.error_message(kernel_pars.outsci, f"{truth_filename}_new.csv")


@pytest.mark.parametrize("kernel", ["square", "point", "turbo", "gaussian", "lanczos3"])
def test_zero_input_weight(kernel, kernel_pars):
    """Tests that do_driz ignores bad pixels, or those that have an input weight (inwht) of 0.

    Parameters
    ----------
    kernel : str
        String associated with one of the c code point kernel options.
    kernel_pars : Class
        The Class inintialized in Get_Class which includes all of the inputs need to run cdriz.tdriz.
    """

    # zero for all insci
    kernel_pars.zero_background()

    # add bad bright pixels in insci
    kernel_pars.insci[0:4, 0:4] = 1e8
    kernel_pars.inwht[0:4, 0:4] = 0

    # adding two additinoal bright "sources"
    kernel_pars.insci[6, 7] = 1000
    kernel_pars.insci[9, 6] = 1000

    # resample
    test_setup.cdriz_call(kernel_pars, kernel)

    # check that any pixel with 0 weight has any counts:
    assert np.sum(kernel_pars.outsci) < 1e5
