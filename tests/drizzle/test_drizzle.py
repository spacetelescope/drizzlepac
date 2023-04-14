import pytest
import numpy as np
from astropy import wcs
from drizzlepac import cdriz


# testing set up
class Get_Grid:
    def __init__(self):
        def get_wcs(_grid):
            w = wcs.WCS()
            w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
            w.wcs.crpix = [_grid[0], _grid[1]]
            w.wcs.crval = [10, 10]
            w.wcs.cdelt = [1e-3, 1e-3]
            w.wcs.set()
            return w

        np.random.seed(0)  # keep same random across each instance
        in_grid = (50, 60)
        out_grid = (51, 66)
        self.insci = np.random.randn(in_grid[0], in_grid[1]).astype("float32")
        self.inwht = np.ones(in_grid, dtype=np.float32)
        self.dny = out_grid[1]  # output y_grid
        self.outsci = np.zeros(out_grid, dtype=np.float32)
        self.outwht = np.ones(out_grid, dtype=np.float32)
        self.outctx = np.ones(out_grid, dtype=np.int32)
        self.w1 = get_wcs(in_grid)
        self.w2 = get_wcs(out_grid)
        self.mapping = cdriz.DefaultWCSMapping(
            self.w1, self.w2, out_grid[0], out_grid[1], 1
        )


def cdriz_call(_set_point_function, kernel):
    cdriz.tdriz(
        _set_point_function.insci,
        _set_point_function.inwht,
        _set_point_function.outsci,
        _set_point_function.outwht,
        _set_point_function.outctx,
        1,  # uniqid
        0,  # ystart
        1,
        1,
        _set_point_function.dny,  # _dny
        1.0,  # pix_ratio
        1.0,
        1.0,
        "center",
        1.0,  # pixfrac
        kernel,
        "cps",  # units
        1.0,  # expscale
        1.0,  # wt_scale
        "INDEF",  # fillval
        0,  # nmiss
        0,  # nskip
        1,
        _set_point_function.mapping,
    )


def save_array(_data, _name):
    np.savetxt(
        _name,
        X=_data,
        fmt="%1.8f",
        delimiter=",",
    )


# ########### TESTS
@pytest.fixture
def point_function():
    return Get_Grid()


@pytest.mark.parametrize("kernel", ["square", "point", "turbo", "gaussian", "lanczos3"])
def test_point_kernel(kernel, point_function, new_truth=False):
    """Function tests different c code point kernels (inputs already created on instantiation).

    Parameters
    ----------
    kernel : str
        String associated with one of the c code point kernel options.
    point_function : Class
        The Class inintialized in Get_Class which includes all of the inputs need to run cdriz.tdriz.
    new_truth : bool, optional
        Flag, whether to save a new csv truth file from the output data, by default False
    """

    # add missing/flagged pixels in inwht
    point_function.insci[20:22, 21:22] = 100

    # resample:
    cdriz_call(point_function, kernel)

    if new_truth:
        # save new truth file
        save_array(point_function.outsci, f"./tests/drizzle/{kernel}_truth.csv")

    truth_array = np.genfromtxt(f"./tests/drizzle/{kernel}_truth.csv", delimiter=",")
    assert np.allclose(point_function.outsci, truth_array, atol=1e-4)


@pytest.mark.parametrize("kernel", ["square", "point", "turbo", "gaussian", "lanczos3"])
def test_zero_input_weight(kernel, point_function):
    """Tests that do_driz ignores bad pixels, or those that have an input weight (inwht) of 0.

    Parameters
    ----------
    kernel : str
        String associated with one of the c code point kernel options.
    point_function : Class
        The Class inintialized in Get_Class which includes all of the inputs need to run cdriz.tdriz.
    """
    # add missing/flagged pixels in inwht
    point_function.inwht[20:21, 22:23] = 0

    # resample
    cdriz_call(point_function, kernel)

    # check that no pixel with 0 weight has any counts:
    assert np.allclose(
        np.sum(np.abs(point_function.outsci[(point_function.outwht == 0)])), 0
    )
