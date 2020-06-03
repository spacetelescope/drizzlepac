import os
import pytest
import math

from astropy.io import fits
from astropy.table import Table

import numpy as np
import stwcs
from stwcs import updatewcs
from stsci.tools import fileutil

from ci_watson.artifactory_helpers import get_bigdata_root
from ci_watson.hst_helpers import raw_from_asn, ref_from_image, download_crds
try:
    from ci_watson.artifactory_helpers import check_url
except ImportError:
    from ci_watson.artifactory_helpers import _is_url as check_url

from .base_classes import BaseTest

__all__ = ['BaseHLATest', 'BaseHLAParTest', 'centroid_compare', 'BaseUnit']

@pytest.mark.usefixtures('_jail')
class BaseHLATest(BaseTest):
    ignore_hdus = []
    input_repo = 'hst-hla-pipeline'
    results_root = 'hst-hla-pipeline-results'
    output_shift_file = None
    fit_limit = 0.010 # 10 milli-arcseconds
    
    docopy = False  # Do not make additional copy by default
    rtol = 1e-6

    refstr = 'jref'
    prevref = os.environ.get(refstr)

    ignore_keywords = ['origin', 'filename', 'date', 'iraf-tlm', 'fitsdate',
                       'upwtim', 'wcscdate', 'upwcsver', 'pywcsver',
                       'history', 'prod_ver', 'rulefile']

    reffile_lookup = ['IDCTAB', 'OFFTAB', 'NPOLFILE', 'D2IMFILE', 'DGEOFILE']

    def set_environ(self):
        # Enforce copies of data when TEST_BIGDATA is URL
        input_dir = get_bigdata_root()

        if input_dir and check_url(input_dir):
            self.docopy = True

        # NOTE: This could be explicitly controlled using pytest fixture
        #       but too many ways to do the same thing would be confusing.
        #       Refine this logic if using pytest fixture.
        # HSTCAL cannot open remote CRDS on FTP but central storage is okay.
        # So use central storage if available to avoid FTP.
        if self.prevref is None or self.prevref.startswith(('ftp', 'http')):
            os.environ[self.refstr] = self.curdir + os.sep
            self.use_ftp_crds = True

        # Turn off Astrometry updates
        os.environ['ASTROMETRY_STEP_CONTROL'] = 'OFF'

    def raw_from_asn(self, asn_file, suffix='_flt.fits'):
        return raw_from_asn(asn_file, suffix='_flt.fits')

    def get_input_file(self, *args, refsep='$', **kwargs):

        # If user has specified action for docopy, apply it with
        # default behavior being whatever was defined in the base class.
        docopy = kwargs.get('docopy', self.docopy)

#        Download or copy input file (e.g., RAW) into the working directory.
#        The associated CRDS reference files in ``refstr`` are also
#        downloaded, if necessary.
        curdir = os.getcwd()
        filenames = self.get_data(*args, docopy=docopy)
        for filename in filenames:
            ref_files = ref_from_image(filename, reffile_lookup=self.reffile_lookup)
            print("Looking for {} REF_FILES: {}".format(filename, ref_files))

            for ref_file in ref_files:
                if ref_file.strip() == '':
                    continue
                if refsep not in ref_file:  # Local file
                    self.get_data('customRef', ref_file, docopy=docopy)
                else:
                    # Start by checking to see whether IRAF variable *ref/*tab
                    # has been added to os.environ
                    refdir, refname = ref_file.split(refsep)
                    refdir_parent = os.path.split(refdir)[0]
                    # Define refdir to point to current directory if:
                    #   i. refdir is not defined in environment already
                    #  ii. refdir in os.environ points to another test directory
                    # This logic should leave refdir unchanged if it already
                    # points to a globally defined directory.
                    if refdir not in os.environ or refdir_parent in curdir:
                        os.environ[refdir] = curdir + os.sep

                    # Download from FTP, if applicable
                    if self.use_ftp_crds:
                        download_crds(ref_file, timeout=self.timeout)
        return filenames

# Pytest function to support the parameterization of these classes
def pytest_generate_tests(metafunc):
    # called once per each test function
    funcarglist = metafunc.cls.params[metafunc.function.__name__]
    argnames = sorted(funcarglist[0])
    idlist = [funcargs['id'] for funcargs in funcarglist]
    del argnames[argnames.index('id')]
    metafunc.parametrize(argnames, [[funcargs[name] for name in argnames]
            for funcargs in funcarglist], ids=idlist)


@pytest.mark.usefixtures('_jail')
class BaseHLAParTest(BaseHLATest):

    params = {'test_modes':[dict(input="",
                                 test_dir=None,
                                 step_class=None,
                                 step_pars=dict(),
                                 output_truth="",
                                 output_hdus=[])
                            ]
             }

    def test_modes(self, input, test_dir, step_class, step_pars,
                   output_truth, output_hdus):
        """
        Template method for parameterizing some tests based on JWST code.
        """
        if test_dir is None:
            return

        self.test_dir = test_dir
        self.ref_loc = [self.test_dir, 'truth']

        # can be removed once all truth files have been updated
        self.ignore_keywords += ['FILENAME']

        input_file = self.get_data(self.test_dir, input)

        result = step_class.call(input_file, **step_pars)

        output_file = result.meta.filename
        result.save(output_file)
        result.close()

        output_pars = None
        if isinstance(output_truth, tuple):
            output_pars = output_truth[1]
            output_truth = output_truth[0]

        if not output_pars:
            if output_hdus:
                output_spec = (output_file, output_truth, output_hdus)
            else:
                output_spec = (output_file, output_truth)
        else:
            output_spec = {'files':(output_file, output_truth),
                           'pars':output_pars}
        outputs = [output_spec]
        self.compare_outputs(outputs)


def centroid_compare(centroid):
    return centroid[1]

class BaseUnit(BaseHLATest):
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
