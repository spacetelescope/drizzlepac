import pytest
import numpy as np
import test_setup


@pytest.fixture
def kernel_pars():
    # get cdriz.triz inputs from test_setup.py
    _params = test_setup.Get_Grid(inx=3, iny=3, outx=4, outy=4)
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
    # add bright pixel at center
    kernel_pars.insci[1:2, 1:2] = 1e4

    # resample:
    test_setup.cdriz_call(kernel_pars, kernel)

    if return_png:
        # save truth file as figure
        test_setup.generate_png(kernel_pars, "cdriz_square.png")

    assert np.allclose(np.sum(kernel_pars.outsci), 5000, 1e-7)


def test_gaussian_kernel(kernel_pars, return_png=True):
    """Same as above test but with gaussian kernel."""
    kernel_pars.insci[1:2, 1:2] = 1e4
    test_setup.cdriz_call(kernel_pars, "gaussian")
    if return_png:
        test_setup.generate_png(kernel_pars, "cdriz_gaussian.png")
    assert np.allclose(np.sum(kernel_pars.outsci), 5314.6846, 1e-3)


def test_lanczos3_kernel(kernel_pars, return_png=True):
    """Same as above test but with lanczos3 kernel."""
    kernel_pars.insci[1:2, 1:2] = 1e4
    test_setup.cdriz_call(kernel_pars, "lanczos3")
    if return_png:
        test_setup.generate_png(kernel_pars, "cdriz_lanczos3.png")
    assert np.allclose(np.sum(kernel_pars.outsci), 4960.454, 1e-3)
