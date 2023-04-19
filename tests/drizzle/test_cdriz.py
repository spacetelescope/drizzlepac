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
@pytest.mark.parametrize("kernel", ["square", "point", 'turbo'])
def test_square_kernel(kernel_pars, kernel):
    """Function tests different c code point kernels (inputs already created on instantiation).

    Parameters
    ----------
    kernel_pars : Class
        The Class inintialized in Get_Grid which includes all of the inputs need to run cdriz.tdriz.
    """
    # add bright pixel at center
    kernel_pars.insci[1:2, 1:2] = 1E4

    # resample:
    test_setup.cdriz_call(kernel_pars, kernel)
    assert np.allclose(np.sum(kernel_pars.outsci), 5000, 1E-7)

def test_gaussian_kernel(kernel_pars):
    """Same as above test but with gaussian kernel."""
    kernel_pars.insci[1:2, 1:2] = 1E4
    test_setup.cdriz_call(kernel_pars, "gaussian")
    assert np.allclose(np.sum(kernel_pars.outsci), 5314.6846, 1E-3)

def test_lanczos3_kernel(kernel_pars):
    """Same as above test but with lanczos3 kernel."""
    kernel_pars.insci[1:2, 1:2] = 1E4
    test_setup.cdriz_call(kernel_pars, "lanczos3")
    assert np.allclose(np.sum(kernel_pars.outsci), 4960.454, 1E-3)
