import getpass
import os
import sys
import math
from io import StringIO
import shutil
import datetime
from os.path import splitext
from difflib import unified_diff

import pytest
from astropy.io import fits
from astropy.io.fits import FITSDiff
from astropy.utils.data import conf

import numpy as np
import stwcs
from stsci.tools import fileutil

from ci_watson.artifactory_helpers import get_bigdata, generate_upload_schema
from ci_watson.hst_helpers import download_crds, ref_from_image


# Base classes for actual tests.
# NOTE: Named in a way so pytest will not pick them up here.
@pytest.mark.bigdata
class BaseCal:
    prevdir = os.getcwd()
    use_ftp_crds = True
    timeout = 30  # seconds
    tree = 'dev'

    # Numpy default for allclose comparison
    rtol = 1e-6
    atol = 1e-5

    # To be defined by instrument
    refstr = ''
    prevref = ''
    input_loc = ''
    ref_loc = ''
    ignore_keywords = []

    # To be defined by individual test
    subdir = ''

    @pytest.fixture(autouse=True)
    def setup_class(self, tmpdir, envopt, pytestconfig):
        """
        Run test in own dir so we can keep results separate from
        other tests.
        """
        if not tmpdir.ensure(self.subdir, dir=True):
            p = tmpdir.mkdir(self.subdir).strpath
        else:
            p = tmpdir.join(self.subdir).strpath
        os.chdir(p)

        # NOTE: This could be explicitly controlled using pytest fixture
        #       but too many ways to do the same thing would be confusing.
        #       Refine this logic if using pytest fixture.
        # HSTCAL cannot open remote CRDS on FTP but central storage is okay.
        # So use central storage if available to avoid FTP.
        if self.prevref is None or self.prevref.startswith(('ftp', 'http')):
            os.environ[self.refstr] = p + os.sep
            self.use_ftp_crds = True

        # Turn off Astrometry updates
        os.environ['ASTROMETRY_STEP_CONTROL'] = 'OFF'

        # This controls astropy.io.fits timeout
        conf.remote_timeout = self.timeout

        # Update tree to point to correct environment
        self.tree = envopt

        # Collect pytest configuration values specified in setup.cfg or pytest.ini
        self.inputs_root = pytestconfig.getini('inputs_root')[0]
        self.results_root = pytestconfig.getini('results_root')[0]

    def teardown_class(self):
        """Reset path and variables."""
        conf.reset('remote_timeout')
        os.chdir(self.prevdir)
        if self.use_ftp_crds and self.prevref is not None:
            os.environ[self.refstr] = self.prevref

    def get_data(self, *args):
        """
        Download `filename` into working directory using
        `get_bigdata`.  This will then return the full path to
        the local copy of the file.
        """
        local_file = get_bigdata(self.inputs_root, self.tree, self.input_loc, *args)

        return local_file

    def get_input_file(self, *args, refsep='$'):
        """
        Download or copy input file (e.g., RAW) into the working directory.
        The associated CRDS reference files in ``refstr`` are also
        downloaded, if necessary.
        """
        filename = self.get_data(*args)
        ref_files = ref_from_image(filename, ['IDCTAB', 'OFFTAB', 'NPOLFILE', 'D2IMFILE', 'DGEOFILE'])
        print("Looking for REF_FILES: {}".format(ref_files))

        for ref_file in ref_files:
            if ref_file.strip() == '':
                continue
            if refsep not in ref_file:  # Local file
                refname = self.get_data('customRef', ref_file)
            else:  # Download from FTP, if applicable
                refname = os.path.join(ref_file)
                if self.use_ftp_crds:
                    download_crds(refname, self.timeout)
        return filename

    def compare_outputs(self, outputs, raise_error=True):
        """
        Compare output with "truth" using appropriate
        diff routine; namely,
            ``fitsdiff`` for FITS file comparisons
            ``unified_diff`` for ASCII products.

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
        # Create instructions for uploading results to artifactory for use
        # as new comparison/truth files
        testpath, testname = os.path.split(os.path.abspath(os.curdir))
        # organize results by day test was run...could replace with git-hash
        whoami = getpass.getuser() or 'nobody'
        dt = datetime.datetime.now().strftime("%d%b%YT")
        ttime = datetime.datetime.now().strftime("%H_%M_%S")
        user_tag = 'NOT_CI_{}_{}'.format(whoami, ttime)
        build_tag = os.environ.get('BUILD_TAG',  user_tag)
        build_suffix = os.environ.get('BUILD_MATRIX_SUFFIX', 'standalone')
        testdir = "{}_{}_{}".format(testname, build_tag, build_suffix)
        tree = os.path.join(self.results_root, self.input_loc,
                            dt, testdir) + os.sep

        updated_outputs = []
        for actual, desired in outputs:
            # Get "truth" image
            s = self.get_data('truth', desired)
            if s is not None:
                desired = s

            if actual.endswith('fits'):
                # Working with FITS files...
                fdiff = FITSDiff(actual, desired, rtol=self.rtol, atol=self.atol,
                                 ignore_keywords=self.ignore_keywords)
                creature_report += fdiff.report()
                if not fdiff.identical:
                    # Only keep track of failed results which need to
                    # be used to replace the truth files (if OK).
                    updated_outputs.append((actual, desired))
                if not fdiff.identical and all_okay:
                    all_okay = False
            else:
                # Try ASCII-based diff
                with open(actual) as afile:
                    actual_lines = afile.readlines()
                with open(desired) as dfile:
                    desired_lines = dfile.readlines()
                udiff = unified_diff(actual_lines, desired_lines,
                                     fromfile=actual, tofile=desired)

                old_stdout = sys.stdout
                udiffIO = StringIO()
                sys.stdout = udiffIO
                sys.stdout.writelines(udiff)
                sys.stdout = old_stdout
                udiff_report = udiffIO.getvalue()
                creature_report += udiff_report
                if len(udiff_report) > 2 and all_okay:
                    all_okay = False
                if len(udiff_report) > 2:
                    # Only keep track of failed results which need to
                    # be used to replace the truth files (if OK).
                    updated_outputs.append((actual, desired))

        if not all_okay:
            # Write out JSON file to enable retention of different results
            new_truths = [os.path.abspath(i[1]) for i in updated_outputs]
            for files in updated_outputs:
                print("Renaming {} as new 'truth' file: {}".format(
                      files[0], files[1]))
                shutil.move(files[0], files[1])
            log_pattern = [os.path.join(os.path.dirname(x), '*.log') for x in new_truths]
            generate_upload_schema(pattern=new_truths + log_pattern,
                           testname=testname,
                           target= tree)

        if not all_okay and raise_error:
            raise AssertionError(os.linesep + creature_report)


        return creature_report


class BaseACS(BaseCal):
    refstr = 'jref'
    prevref = os.environ.get(refstr)
    input_loc = 'acs'
    ref_loc = 'acs'
    ignore_keywords = ['origin', 'filename', 'date', 'iraf-tlm', 'fitsdate',
                       'upwtim', 'wcscdate', 'upwcsver', 'pywcsver',
                       'history', 'prod_ver', 'rulefile']


class BaseACSHRC(BaseACS):
    input_loc = 'acs/hrc'
    ref_loc = 'acs/hrc/ref'


class BaseACSWFC(BaseACS):
    input_loc = 'acs/wfc'
    ref_loc = 'acs/wfc/ref'


class BaseWFC3(BaseCal):
    refstr = 'iref'
    input_loc = 'wfc3'
    ref_loc = 'wfc3/ref'
    prevref = os.environ.get(refstr)
    ignore_keywords = ['origin', 'filename', 'date', 'iraf-tlm', 'fitsdate',
                       'upwtim', 'wcscdate', 'upwcsver', 'pywcsver',
                       'history', 'prod_ver', 'rulefile']


class BaseSTIS(BaseCal):
    refstr = 'oref'
    prevref = os.environ.get(refstr)
    input_loc = 'stis'
    ref_loc = 'stis/ref'
    ignore_keywords = ['origin', 'filename', 'date', 'iraf-tlm', 'fitsdate',
                       'upwtim', 'wcscdate', 'upwcsver', 'pywcsver',
                       'history', 'prod_ver', 'rulefile']


class BaseWFPC2(BaseCal):
    refstr = 'uref'
    prevref = os.environ.get(refstr)
    input_loc = 'wfpc2'
    ref_loc = 'wfpc2/ref'
    ignore_keywords = ['origin', 'filename', 'date', 'iraf-tlm', 'fitsdate',
                       'upwtim', 'wcscdate', 'upwcsver', 'pywcsver',
                       'history', 'prod_ver', 'rulefile']

def centroid_compare(centroid):
    return centroid[1]

class BaseUnit(BaseCal):
    buff = 0
    refstr = 'jref'
    prevref = os.environ.get(refstr)
    input_loc = 'acs'
    ref_loc = 'acs'
    ignore_keywords = ['origin', 'filename', 'date', 'iraf-tlm', 'fitsdate',
                       'upwtim', 'wcscdate', 'upwcsver', 'pywcsver',
                       'history', 'prod_ver', 'rulefile']
    atol = 1.0e-5

    def bound_image(self, image):
        """
        Compute region where image is non-zero
        """
        coords = np.nonzero(image)
        ymin = coords[0].min()
        ymax = coords[0].max()
        xmin = coords[1].min()
        xmax = coords[1].max()
        return (ymin, ymax, xmin, xmax)

    def centroid(self, image, size, center):
        """
        Compute the centroid of a rectangular area
        """
        ylo = int(center[0]) - size // 2
        yhi = min(ylo + size, image.shape[0])
        xlo = int(center[1]) - size // 2
        xhi = min(xlo + size, image.shape[1])

        center = [0.0, 0.0, 0.0]
        for y in range(ylo, yhi):
            for x in range(xlo, xhi):
                center[0] += y * image[y,x]
                center[1] += x * image[y,x]
                center[2] += image[y,x]

        if center[2] == 0.0: return None

        center[0] /= center[2]
        center[1] /= center[2]
        return center

    def centroid_close(self, list_of_centroids, size, point):
        """
        Find if any centroid is close to a point
        """
        for i in range(len(list_of_centroids)-1, -1, -1):
            if (abs(list_of_centroids[i][0] - point[0]) < size / 2 and
                abs(list_of_centroids[i][1] - point[1]) < size / 2):
                return 1

        return 0

    def centroid_distances(self, image1, image2, amp, size):
        """
        Compute a list of centroids and the distances between them in two images
        """
        distances = []
        list_of_centroids, lst_pts = self.centroid_list(image2, amp, size)

        for center2, pt in zip(list_of_centroids, lst_pts):
            center1 = self.centroid(image1, size, pt)
            if center1 is None: continue

            disty = center2[0] - center1[0]
            distx = center2[1] - center1[1]
            dist = math.sqrt(disty * disty + distx * distx)
            dflux = abs(center2[2] - center1[2])
            distances.append([dist, dflux, center1, center2])

        distances.sort(key=centroid_compare)
        return distances

    def centroid_list(self, image, amp, size):
        """
        Find the next centroid
        """
        list_of_centroids = []
        list_of_points = []
        points = np.transpose(np.nonzero(image > amp))

        for point in points:
            if not self.centroid_close(list_of_centroids, size, point):
                center = self.centroid(image, size, point)
                list_of_centroids.append(center)
                list_of_points.append(point)

        return list_of_centroids, list_of_points

    def centroid_statistics(self, title, fname, image1, image2, amp, size):
        """
        write centroid statistics to compare differences btw two images
        """
        stats = ("minimum", "median", "maximum")
        images = (None, None, image1, image2)
        im_type = ("", "", "test", "reference")

        diff = []
        distances = self.centroid_distances(image1, image2, amp, size)
        indexes = (0, len(distances)//2, len(distances)-1)
        fd = open(fname, 'w')
        fd.write("*** %s ***\n" % title)

        if len(distances) == 0:
            diff = [0.0, 0.0, 0.0]
            fd.write("No matches!!\n")

        elif len(distances) == 1:
            diff = [distances[0][0], distances[0][0], distances[0][0]]

            fd.write("1 match\n")
            fd.write("distance = %f flux difference = %f\n" % (distances[0][0], distances[0][1]))

            for j in range(2, 4):
                ylo = int(distances[0][j][0]) - (1+self.buff)
                yhi = int(distances[0][j][0]) + (2+self.buff)
                xlo = int(distances[0][j][1]) - (1+self.buff)
                xhi = int(distances[0][j][1]) + (2+self.buff)
                subimage = images[j][ylo:yhi,xlo:xhi]
                fd.write("\n%s image centroid = (%f,%f) image flux = %f\n" %
                         (im_type[j], distances[0][j][0], distances[0][j][1], distances[0][j][2]))
                fd.write(str(subimage) + "\n")

        else:
            fd.write("%d matches\n" % len(distances))

            for k in range(0,3):
                i = indexes[k]
                diff.append(distances[i][0])
                fd.write("\n%s distance = %f flux difference = %f\n" % (stats[k], distances[i][0], distances[i][1]))

                for j in range(2, 4):
                    ylo = int(distances[i][j][0]) - (1+self.buff)
                    yhi = int(distances[i][j][0]) + (2+self.buff)
                    xlo = int(distances[i][j][1]) - (1+self.buff)
                    xhi = int(distances[i][j][1]) + (2+self.buff)
                    subimage = images[j][ylo:yhi,xlo:xhi]
                    fd.write("\n%s %s image centroid = (%f,%f) image flux = %f\n" %
                             (stats[k], im_type[j], distances[i][j][0], distances[i][j][1], distances[i][j][2]))
                    fd.write(str(subimage) + "\n")

        fd.close()
        return tuple(diff)

    def make_point_image(self, input_image, point, value):
        """
        Create an image with a single point set
        """
        output_image = np.zeros(input_image.shape, dtype=input_image.dtype)
        output_image[point] = value
        return output_image

    def make_grid_image(self, input_image, spacing, value):
        """
        Create an image with points on a grid set
        """
        output_image = np.zeros(input_image.shape, dtype=input_image.dtype)

        shape = output_image.shape
        for y in range(spacing//2, shape[0], spacing):
            for x in range(spacing//2, shape[1], spacing):
                output_image[y,x] = value

        return output_image

    def print_wcs(self, title, wcs):
        """
        Print the wcs header cards
        """
        print("=== %s ===" % title)
        print(wcs.to_header_string())


    def read_image(self, filename):
        """
        Read the image from a fits file
        """
        hdu = fits.open(filename)

        image = hdu[1].data
        hdu.close()
        return image

    def read_wcs(self, filename):
        """
        Read the wcs of a fits file
        """
        hdu = fits.open(filename)

        wcs = stwcs.wcsutil.HSTWCS(hdu, 1)
        hdu.close()
        return wcs

    def write_wcs(self, hdu, image_wcs):
        """
        Update header with WCS keywords
        """
        hdu.header['ORIENTAT'] = image_wcs.orientat
        hdu.header['CD1_1'] = image_wcs.wcs.cd[0][0]
        hdu.header['CD1_2'] = image_wcs.wcs.cd[0][1]
        hdu.header['CD2_1'] = image_wcs.wcs.cd[1][0]
        hdu.header['CD2_2'] = image_wcs.wcs.cd[1][1]
        hdu.header['CRVAL1'] = image_wcs.wcs.crval[0]
        hdu.header['CRVAL2'] = image_wcs.wcs.crval[1]
        hdu.header['CRPIX1'] = image_wcs.wcs.crpix[0]
        hdu.header['CRPIX2'] = image_wcs.wcs.crpix[1]
        hdu.header['CTYPE1'] = image_wcs.wcs.ctype[0]
        hdu.header['CTYPE2'] = image_wcs.wcs.ctype[1]
        hdu.header['VAFACTOR'] = 1.0

    def write_image(self, filename, wcs, *args):
        """
        Read the image from a fits file
        """
        extarray = ['SCI', 'WHT', 'CTX']

        pimg = fits.HDUList()
        phdu = fits.PrimaryHDU()
        phdu.header['NDRIZIM'] = 1
        phdu.header['ROOTNAME'] = filename
        pimg.append(phdu)

        for img in args:
            # Create a MEF file with the specified extname
            extn = extarray.pop(0)
            extname = fileutil.parseExtn(extn)

            ehdu = fits.ImageHDU(data=img)
            ehdu.header['EXTNAME'] = extname[0]
            ehdu.header['EXTVER'] = extname[1]
            self.write_wcs(ehdu, wcs)
            pimg.append(ehdu)

        pimg.writeto(filename)
        del pimg


def add_suffix(fname, suffix, range=None):
    """Add suffix to file name

    Parameters
    ----------
    fname: str
        The file name to add the suffix to

    suffix: str
        The suffix to add_suffix

    range: range
        If specified, the set of indexes will be added to the
        outputs.

    Returns
    -------
    fname, fname_with_suffix
        2-tuple of the original file name and name with suffix.
        If `range` is defined, `fname_with_suffix` will be a list.

    """
    fname_root, fname_ext = splitext(fname)
    if range is None:
        with_suffix = ''.join([
            fname_root,
            '_',
            suffix,
            fname_ext
        ])
    else:
        with_suffix = []
        for idx in range:
            with_suffix.append(''.join([
                fname_root,
                '_',
                str(idx),
                '_',
                suffix,
                fname_ext
            ]))

    return fname, with_suffix
