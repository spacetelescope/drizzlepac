"""HSTCAL regression test helpers."""
from __future__ import absolute_import, division, print_function
from astropy.extern.six.moves import urllib

import os
import shutil
from distutils.spawn import find_executable

import pytest
import requests
from astropy.io import fits
from astropy.io.fits import FITSDiff
from astropy.table import Table
from astropy.utils.data import conf

__all__ = ['slow', 'download_crds',
           'download_file_cgi', 'ref_from_image', 'raw_from_asn', 'BaseACS',
           'BaseSTIS', 'BaseWFC3IR', 'BaseWFC3UVIS', 'BaseWFPC2']

# pytest marker to mark resource-intensive tests that should not be
# executed with every commit.
slow = pytest.mark.slow


def _download_file(url, filename, filemode='wb', timeout=None):
    """Generic remote data download."""
    if url.startswith('http'):
        r = requests.get(url, timeout=timeout)
        with open(filename, filemode) as fout:
            fout.write(r.content)
    elif url.startswith('ftp'):  # TODO: Support filemode and timeout.
        urllib.request.urlretrieve(url, filename=filename)
    else:  # pragma: no cover
        raise ValueError('Unsupported protocol for {}'.format(url))


def download_crds(refdir, refname, timeout=None):
    """Download a CRDS file from HTTP or FTP to current directory."""
    # CRDS file for given name never changes, so no need to re-download.
    if os.path.exists(refname):
        return

    try:
        url = 'http://ssb.stsci.edu/cdbs/{}/{}'.format(refdir, refname)
        _download_file(url, refname, timeout=timeout)
    except Exception:  # Fall back to FTP
        url = 'ftp://ftp.stsci.edu/cdbs/{}/{}'.format(refdir, refname)
        _download_file(url, refname, timeout=timeout)


def download_file_cgi(tree, project, filename, filemode='wb', timeout=30,
                      allow_remote_ref=False):
    """
    Download remote data to current directory from Central Storage
    or using CGI script interface.

    Parameters
    ----------
    tree : {'rt', 'rtx', 'null'}
        Test tree:

        * rt = dev
        * rtx = public
        * null = neither, grab first match (e.g., for STAK)

    project : str
        Path containing data in the test tree.
        For example, ``hstcal/acs/calacs_e``.

    filename : str
        Filename of the data.

    filemode : str
        Mode of the output file writer. By default, it saves as binary
        and overwrites existing file.

    timeout : number or `None`
        This is not a time limit on the entire response download;
        rather, an exception is raised if the server has not issued
        a response for timeout seconds (more precisely, if no bytes
        have been received on the underlying socket for timeout seconds).
        If no timeout is specified explicitly, requests do not time out.

    allow_remote_ref : bool
        If set to `True`, instead of copying or downloading to current
        directory, return the full path. This should only be used for
        "truth" images.

    """
    # NOTE: This could be explicitly controlled using pytest fixture
    #       but too many ways to do the same thing would be confusing.
    #       Refine this logic if using pytest fixture.
    # cs_file = os.path.join('/eng/ssb2/tests', tree, project, filename)
    cs_file = os.path.join('/user/hack/data/mdtng/RT', project, filename)

    if os.path.isfile(cs_file):
        if allow_remote_ref:
            return cs_file
        shutil.copyfile(cs_file, filename)
    else:
        url = ('http://ssb.stsci.edu/cgi-bin/remote_testing.cgi?'
               'tree={}&project={}&name={}'.format(tree, project, filename))
        if allow_remote_ref:
            return url
        _download_file(url, filename, filemode=filemode, timeout=timeout)


def _get_reffile(hdr, key):
    """Get ref file from given key in given FITS header."""
    ref_file = None
    if key in hdr:  # Keyword might not exist
        ref_file = hdr[key].strip()
        if ref_file.upper() == 'N/A':  # Not all ref file is defined
            ref_file = None
    return ref_file


def ref_from_image(input_image):
    """
    Return a list of reference filenames, as defined in the primary
    header of the given input image, necessary for calibration; i.e.,
    only those associated with ``*CORR`` set to ``PERFORM`` will be
    considered.
    """
    # NOTE: Add additional mapping as needed.
    # Map mandatory CRDS reference file for instrument/detector combo.

    reffile_lookup = ['IDCTAB', 'OFFTAB', 'NPOLFILE', 'D2IMFILE']
    
    ref_files = []
    hdr = fits.getheader(input_image, ext=0)

    for reffile in reffile_lookup:
        s = _get_reffile(hdr, reffile)
        if s is not None:
            ref_files.append(s)

    return ref_files


def raw_from_asn(asn_file, suffix='_raw.fits'):
    """Return a list of RAW input files in a given ASN."""
    raw_files = []
    tab = Table.read(asn_file, format='fits')

    for row in tab:
        if row['MEMTYPE'].startswith('PROD'):
            continue
        pfx = row['MEMNAME'].lower().strip().replace('\x00', '')
        raw_files.append(pfx + suffix)

    return raw_files


# Base classes for actual tests.
# NOTE: Named in a way so pytest will not pick them up here.

@pytest.mark.remote_data
class BaseCal(object):
    prevdir = os.getcwd()
    use_ftp_crds = False
    timeout = 30  # seconds
    tree = 'rt'  # Use dev for now

    # Numpy default for allclose comparison
    rtol = 1e-7
    atol = 0

    # To be defined by instrument
    refstr = ''
    prevref = ''
    input_loc = ''
    ref_loc = ''
    ignore_keywords = []

    # To be defined by individual test
    subdir = ''

    @pytest.fixture(autouse=True)
    def setup_class(self, tmpdir):
        """
        Run test in own dir so we can keep results separate from
        other tests.
        """        
        p = tmpdir.mkdir(self.subdir).strpath
        os.chdir(p)

        # NOTE: This could be explicitly controlled using pytest fixture
        #       but too many ways to do the same thing would be confusing.
        #       Refine this logic if using pytest fixture.
        # HSTCAL cannot open remote CRDS on FTP but central storage is okay.
        # So use central storage if available to avoid FTP.
        if self.prevref is None or self.prevref.startswith(('ftp', 'http')):
            os.environ[self.refstr] = p + os.sep
            self.use_ftp_crds = True

        # This controls astropy.io.fits timeout
        conf.remote_timeout = self.timeout

    def teardown_class(self):
        """Reset path and variables."""
        conf.reset('remote_timeout')
        os.chdir(self.prevdir)
        if self.use_ftp_crds and self.prevref is not None:
            os.environ[self.refstr] = self.prevref

    def get_input_file(self, filename, refsep='$'):
        """
        Download or copy input file (e.g., RAW) into the working directory.
        The associated CRDS reference files in ``refstr`` are also
        downloaded, if necessary.
        """
        download_file_cgi(self.tree, self.input_loc, filename,
                          timeout=self.timeout)
        ref_files = ref_from_image(filename)

        for ref_file in ref_files:
            if refsep not in ref_file:  # Local file
                download_file_cgi(self.tree, self.input_loc, ref_file,
                                  timeout=self.timeout)
            else:  # Download from FTP, if applicable
                s = ref_file.split(refsep)
                refdir = s[0]
                refname = s[1]
                if self.use_ftp_crds:
                    download_crds(refdir, refname, timeout=self.timeout)

    def compare_outputs(self, outputs, raise_error=True):
        """
        Compare CALXXX output with "truth" using ``fitsdiff``.

        Parameters
        ----------
        outputs : list of tuple
            A list of tuples, each containing filename (without path)
            of CALXXX output and truth, in that order.

        raise_error : bool
            Raise ``AssertionError`` if difference is found.

        Returns
        -------
        report : str
            Report from ``fitsdiff``.
            This is part of error message if ``raise_error=True``.

        """
        all_okay = True
        creature_report = ''

        for actual, desired in outputs:
            # Get "truth" image
            s = download_file_cgi(self.tree, self.ref_loc, desired,
                                  allow_remote_ref=True, timeout=self.timeout)
            if s is not None:
                desired = s

            fdiff = FITSDiff(actual, desired, rtol=self.rtol, atol=self.atol,
                             ignore_keywords=self.ignore_keywords)
            creature_report += fdiff.report()

            if not fdiff.identical and all_okay:
                all_okay = False

        if not all_okay and raise_error:
            raise AssertionError(os.linesep + creature_report)

        return creature_report


class BaseACS(BaseCal):
    refstr = 'jref'
    prevref = os.environ.get(refstr)
    input_loc = 'drizzlepac/acs'
    ref_loc = 'drizzlepac/acs'
    ignore_keywords = ['origin', 'filename', 'date', 'iraf-tlm', 'fitsdate',
                       'upwtim', 'wcscdate', 'upwcsver', 'pywcsver',
                       'history']


class BaseACSHRC(BaseACS):
    input_loc = 'drizzlepac/acs/hrc'
    ref_loc = 'drizzlepac/acs/hrc/ref'


class BaseACSWFC(BaseACS):
    input_loc = 'drizzlepac/acs/wfc'
    ref_loc = 'drizzlepac/acs/wfc/ref'


class BaseWFC3(BaseCal):
    refstr = 'iref'
    input_loc = 'drizzlepac/wf3'
    ref_loc = 'drizzlepac/wf3/ref'
    prevref = os.environ.get(refstr)
    ignore_keywords = ['origin', 'filename', 'date', 'iraf-tlm', 'fitsdate',
                       'upwtim', 'wcscdate', 'upwcsver', 'pywcsver',
                       'history']


class BaseSTIS(BaseCal):
    refstr = 'oref'
    prevref = os.environ.get(refstr)
    input_loc = 'drizzlepac/stis'
    ref_loc = 'drizzlepac/stis/ref'
    ignore_keywords = ['origin', 'filename', 'date', 'iraf-tlm', 'fitsdate',
                       'upwtim', 'wcscdate', 'upwcsver', 'pywcsver',
                       'history']

class BaseWFPC2(BaseCal):
    refstr = 'uref'
    prevref = os.environ.get(refstr)
    input_loc = 'drizzlepac/wfpc2'
    ref_loc = 'drizzlepac/wfpc2/ref'
    ignore_keywords = ['origin', 'filename', 'date', 'iraf-tlm', 'fitsdate',
                       'upwtim', 'wcscdate', 'upwcsver', 'pywcsver',
                       'history']

