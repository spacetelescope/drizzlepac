import pytest
import numpy as np
from astropy import wcs
from drizzlepac import cdriz
import test_setup

@pytest.fixture
def kernel_pars():
    return test_setup.Get_Grid(inx=4,iny=4, outx=5, outy=5)


# "square", "point", "turbo", "gaussian", "lanczos3"

def test_square_kernel(kernel_pars):
    """Function tests different c code point kernels (inputs already created on instantiation).

    Parameters
    ----------
    kernel : str
        String associated with one of the c code point kernel options.
    kernel_pars : Class
        The Class inintialized in Get_Grid which includes all of the inputs need to run cdriz.tdriz.
    """

    kernel='square'

    # add missing/flagged pixels in inwht
    kernel_pars.insci[2:3, 2:3] = 10

    # resample:
    test_setup.cdriz_call(kernel_pars, kernel)

    assert np.allclose(kernel_pars.outsci, truth_array, atol=1e-4)


