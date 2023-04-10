import pytest
import numpy as np
from astropy import wcs
from drizzlepac import cdriz


class get_grid:
    def __init__(self):
        def get_wcs(crpix1, crpix2):
            w = wcs.WCS()
            w.wcs.ctype = ["RA---CAR", "DEC--CAR"]
            w.wcs.crpix = [crpix1, crpix2]
            w.wcs.crval = [10, 10]
            w.wcs.cdelt = [1e-3, 1e-3]
            w.wcs.set()
            return w

        self.insci = np.ones((50, 60), dtype=np.float32)
        self.inwht = np.ones((50, 60), dtype=np.float32)
        self.outsci = np.zeros((51, 66), dtype=np.float32)
        self.outwht = np.ones((51, 66), dtype=np.float32)
        self.outctx = np.ones((51, 66), dtype=np.int32)
        self.w1 = get_wcs(50, 60)
        self.w2 = get_wcs(51, 66)
        self.mapping = cdriz.DefaultWCSMapping(self.w1, self.w2, 51, 66, 1)


@pytest.fixture
def point_function():
    return get_grid()


@pytest.mark.parametrize("kernel", ["square", "point", "turbo", "gaussian", "lanczos3"])
def test_zero_input_weight(kernel, point_function):
    """
    Test do_driz with missing pixels
    """
    # add missing/flagged pixels in inwht
    point_function.inwht[20:21, 22:23] = 0

    # resample:
    cdriz.tdriz(
        point_function.insci,
        point_function.inwht,
        point_function.outsci,
        point_function.outwht,
        point_function.outctx,
        1, # uniqid
        0, # ystart
        1,
        1,
        200, #_dny
        1.0, # pix_ratio
        1.0,
        1.0,
        "center",
        1.0, # pixfrac
        kernel,
        "cps", # units
        1.0, # expscale
        1.0, # wt_scale
        "INDEF", # fillval
        0, # nmiss
        0, # nskip
        1,
        point_function.mapping,
    )             

    # check that no pixel with 0 weight has any counts:
    assert np.allclose(
        np.sum(np.abs(point_function.outsci[(point_function.outwht == 0)])), 0
    )

#test kernel response function
@pytest.mark.parametrize("kernel", ["square", "point", "turbo", "gaussian", "lanczos3"])
def test_point_kernel(kernel, point_function):
    """
    Test do_driz kernel with point
    """
    # add missing/flagged pixels in inwht
    point_function.insci[20:21, 21:22] = 50
    # resample:
    cdriz.tdriz(
        point_function.insci,
        point_function.inwht,
        point_function.outsci,
        point_function.outwht,
        point_function.outctx,
        1,
        0,
        1,
        1,
        200,
        1.0,
        1.0,
        1.0,
        "center",
        1.0,
        kernel,
        "cps",
        1.0,
        1.0,
        "INDEF",
        0,
        0,
        1,
        point_function.mapping,
    )
    # for generating truth files
    np.savetxt('./tests/drizzle/'+kernel+'_truth.csv', X=point_function.outsci, delimiter=',')
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(4,2))
    ax1 = fig.add_subplot(111, projection=point_function.w1)
    ax1.imshow(point_function.outsci, origin='lower', cmap='Greys')
    ax1.set_ylabel(' ')
    ax1.set_xlabel(' ')
    fig.savefig('./tests/drizzle/'+kernel+'_truth.png') 

    #for generating input file
    np.savetxt('./tests/drizzle/input_truth.csv', X=point_function.insci, fmt='%.32e',delimiter=',')
    fig = plt.figure(figsize=(4,2))
    ax1 = fig.add_subplot(111, projection=point_function.w1)
    ax1.imshow(point_function.insci, origin='lower', cmap='Greys')
    ax1.set_ylabel(' ')
    ax1.set_xlabel(' ')
    fig.savefig('./tests/drizzle/input_truth.png')    

    truth_array = np.genfromtxt('./tests/drizzle/'+kernel+'_truth.csv', delimiter=',')

    assert np.allclose(point_function.outsci, truth_array, atol=1e-4)
