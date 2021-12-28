import pytest

import numpy as np
from astropy import wcs

from drizzlepac import cdriz

@pytest.mark.parametrize(
    'kernel', ['square', 'point', 'turbo', 'gaussian', 'lanczos3']
)
def test_zero_input_weight(kernel):
    """
    Test do_driz square kernel with grid
    """
    # initialize input:
    insci = np.ones((200, 400), dtype=np.float32)
    inwht = np.ones((200, 400), dtype=np.float32)
    inwht[:, 150:155] = 0

    # initialize output:
    outsci = np.zeros((210, 410), dtype=np.float32)
    outwht = np.zeros((210, 410), dtype=np.float32)
    outctx = np.zeros((210, 410), dtype=np.int32)

    # define coordinate mapping:
    w1 = wcs.WCS()
    w1.wcs.ctype = ['RA---CAR', 'DEC--CAR']
    w1.wcs.crpix = [201, 101]
    w1.wcs.crval = [10, 10]
    w1.wcs.cdelt = [1e-3, 1e-3]
    w1.wcs.set()

    w2 = wcs.WCS()
    w2.wcs.ctype = ['RA---CAR', 'DEC--CAR']
    w2.wcs.crpix = [206, 106]
    w2.wcs.crval = [10, 10]
    w2.wcs.cdelt = [1e-3, 1e-3]
    w2.wcs.set()

    mapping = cdriz.DefaultWCSMapping(w1, w2, 400, 200, 1)

    # resample:
    cdriz.tdriz(
        insci, inwht, outsci, outwht,
        outctx, 1, 0, 1, 1, 200,
        1.0, 1.0, 1.0, 'center', 1.0,
        kernel, 'cps', 1.0, 1.0,
        'INDEF', 0, 0, 1, mapping
    )

    # check that no pixel with 0 weight has any counts:
    assert np.allclose(np.sum(np.abs(outsci[(outwht == 0)])), 0)
