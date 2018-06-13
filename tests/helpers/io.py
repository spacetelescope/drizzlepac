import copy
import json
import os
import shutil

from .utils import check_url, download

UPLOAD_SCHEMA = {"files": [
                    {"pattern": "",
                     "target": "",
                     "props": None,
                     "recursive": "false",
                     "flat": "true",
                     "regexp": "false",
                     "explode": "false",
                     "excludePatterns": []
                    }
                  ]
                }

__all__ = ['BigdataError', 'get_bigdata', 'upload_results']


BIGDATA_PATHS = [
    os.environ.get('TEST_BIGDATA', '/srv/rt/betadrizzle'),
    'https://bytesalad.stsci.edu/artifactory/drizzlepac'
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


def upload_results(**kwargs):
    """Write out JSON file to upload results from test to storage area.

    This function relies on the JFROG JSON schema for uploading data into
    artifactory using the Jenkins plugin.  Docs can be found at::

        https://www.jfrog.com/confluence/display/RTF/Using+File+Specs

    Parameters
    ----------
    pattern : str or list of strings
        Specifies the local file system path to test results which should be
        uploaded to Artifactory. You can specify multiple artifacts by using
        wildcards or a regular expression as designated by the regexp property.

    target : str
        Specifies the target path in Artifactory in the following format:
            [repository_name]/[repository_path]

    testname : str
        Name of test that generate the results.  This will be used to create the
        name of the JSON file to enable these results to be uploaded to Artifactory.

    recursive : bool, optional
        Specify whether or not to identify files listed in sub-directories
        for uploading.  Default: False

    """
    # Interpret mandatory inputs
    pattern = kwargs.get("pattern")
    target = kwargs.get("target")
    testname = kwargs.get("testname")

    # Finish interpreting inputs
    jsonfile = "{}_results.json".format(testname)
    recursive = repr(kwargs.get("recursive", False)).lower()

    if isinstance(pattern, list):
        # Populate schema for this test's data
        upload_schema = {"files": []}

        for p in pattern:
            temp_schema = copy.deepcopy(UPLOAD_SCHEMA["files"][0])
            temp_schema.update({"pattern": p, "target": target, "recursive": recursive})
            upload_schema["files"].append(temp_schema)

    else:
        # Populate schema for this test's data
        upload_schema = copy.deepcopy(UPLOAD_SCHEMA)
        upload_schema["files"][0].update({"pattern": pattern, "target": target, "recursive": recursive})

    # Write out JSON file with description of test results
    with open(jsonfile, 'w') as outfile:
        json.dump(upload_schema, outfile)
