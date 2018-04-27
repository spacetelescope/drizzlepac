import os
import shutil

from .utils import check_url, download

__all__ = ['BigdataError', 'get_bigdata']


BIGDATA_PATHS = [
    os.environ.get('TEST_BIGDATA', '/data4/jwst_test_data'),
    'https://bytesalad.stsci.edu/artifactory/rt-drizzlepac-dev'
]


class BigdataError(Exception):
    pass


def _select_bigdata():
    """ Find and returns the path to the nearest big datasets
    """
    for path in BIGDATA_PATHS:
        if os.path.exists(path) or check_url(path):
            return path

    return None


def get_bigdata(*args):
    """ Acquire requested data from a managed resource

    Usage:
        filename = get_bigdata('abc', '123', 'sample.fits')
        with open(filename, 'rb') as data:
            example = data.read()

    Returns:
        Absolute path to local copy of data (i.e. /path/to/example.fits)
    """

    src = os.path.join(_select_bigdata(), *args)
    filename = os.path.basename(src)
    dest = os.path.abspath(os.path.join(os.curdir, filename))

    if os.path.exists(src):
        if src == dest:
            raise BigdataError('Source and destination paths are identical: '
                               '{}'.format(src))
        shutil.copy2(src, dest)

    elif check_url(src):
        download(src, dest)

    else:
        raise BigdataError('Failed to retrieve data: {}'.format(src))

    return dest
