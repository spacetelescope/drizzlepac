import numpy as np
from astropy import wcs

from drizzlepac import cdriz

def test_zero_input_weight():
    """
    Test do_driz square kernel with grid
    """
    insci = np.ones((200, 400), dtype=np.float32)
    inwht = np.ones((200, 400), dtype=np.float32)
    inwht[:, 150:155] = 0
    tflux = np.sum(inwht)

    for kernel in ['square', 'point', 'turbo']:
        # initialize output:
        outsci = np.zeros((210, 410), dtype=np.float32)
        outwht = np.zeros((210, 410), dtype=np.float32)
        outctx = np.zeros((210, 410), dtype=np.int32)

        w = wcs.WCS()
        w.wcs.ctype = ['RA---CAR', 'DEC--CAR']
        w.wcs.crpix = [201, 101]
        w.wcs.crval = [10, 10]
        w.wcs.cdelt = [1e-3, 1e-3]
        w.wcs.set()

        mapping = cdriz.DefaultWCSMapping(w, w, 400, 200, 1)

        cdriz.tdriz(
            insci, inwht, outsci, outwht,
            outctx, 1, 0, 1, 1, insci.shape[0],
            1.0, 1.0, 1.0, 'center', 1.0,
            kernel, 'cps', 1.0, 1.0,
            'INDEF', 0, 0, 1, mapping)

        assert np.allclose(np.sum(outwht), tflux)
