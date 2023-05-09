import pytest
import os
import numpy as np
import cdriz_setup


@pytest.fixture
def kernel_pars():
    # get cdriz.triz inputs from cdriz_setup.py
    _params = cdriz_setup.Get_Grid(inx=10, iny=10, outx=13, outy=13)
    _params.zero_background()
    return _params


# "square", "point", "turbo", "gaussian", "lanczos3"
@pytest.mark.parametrize("kernel", ["square", "point", "turbo"])
def test_square_kernel(kernel_pars, kernel, return_png=True):
    """Function tests different c code point kernels (inputs already created on instantiation).

    Parameters
    ----------
    kernel_pars : Class
        The Class inintialized in Get_Grid which includes all of the inputs need to run cdriz.tdriz.
    kernel: str [argument passed with parameterize]
        The name of the kernel being used.
    return_png:
        Whether to return a png of the outputs.
    """

    output_name = f"cdriz_{kernel}.png"
    relative_path = "truth_files"
    output_fullpath = cdriz_setup.get_output_fullpath(relative_path, output_name)

    # add bright pixel at center
    kernel_pars.insci[3, 3] = 1e4

    # resample:
    cdriz_setup.cdriz_call(kernel_pars, kernel)

    if return_png:
        # save truth file as figure
        cdriz_setup.generate_png(kernel_pars, output_fullpath)

    assert np.allclose(np.sum(kernel_pars.outsci), 10000, 1e-7)


def test_gaussian_kernel(kernel_pars, return_png=True):
    """Same as above test but with gaussian kernel."""

    output_name = "cdriz_gaussian.png"
    relative_path = "truth_files"
    output_fullpath = cdriz_setup.get_output_fullpath(relative_path, output_name)

    kernel_pars.insci[3, 3] = 1e4
    cdriz_setup.cdriz_call(kernel_pars, "gaussian")
    if return_png:
        cdriz_setup.generate_png(kernel_pars, output_fullpath)
    assert np.allclose(np.sum(kernel_pars.outsci), 10000, 1e-3)


def test_lanczos3_kernel(kernel_pars, return_png=True):
    """Same as above test but with lanczos3 kernel."""

    output_name = "cdriz_lanczos3.png"
    relative_path = "truth_files"
    output_fullpath = cdriz_setup.get_output_fullpath(relative_path, output_name)

    kernel_pars.insci[3, 3] = 1e4
    cdriz_setup.cdriz_call(kernel_pars, "lanczos3")
    if return_png:
        cdriz_setup.generate_png(kernel_pars, output_fullpath)
    assert np.allclose(np.sum(kernel_pars.outsci), 9882.103, 1e-3)
