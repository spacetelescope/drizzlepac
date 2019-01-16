import os
import pytest

import numpy as np
from stwcs import distortion
from ..resources import BaseUnit

import drizzlepac.adrizzle as adrizzle
import drizzlepac.ablot as ablot


class TestDriz(BaseUnit):

    def test_square_with_point(self):
        """
        Test do_driz square kernel with point
        """
        input = os.path.basename(self.get_input_file('input','j8bt06nyq_unit.fits'))
        output = 'output_square_point.fits'
        output_difference = 'difference_square_point.txt'
        output_template = os.path.basename(self.get_data('truth',
                                           'reference_square_point.fits'))

        insci = self.read_image(input)
        input_wcs = self.read_wcs(input)
        insci = self.make_point_image(insci, (500, 200), 100.0)
        inwht = np.ones(insci.shape,dtype=insci.dtype)
        output_wcs = self.read_wcs(output_template)

        naxis1, naxis2 = output_wcs.pixel_shape
        outsci = np.zeros((naxis2, naxis1), dtype='float32')
        outwht = np.zeros((naxis2, naxis1), dtype='float32')
        outcon = np.zeros((1, naxis2, naxis1), dtype='i4')

        expin = 1.0
        wt_scl = expin
        in_units = 'cps'
        wcslin = distortion.utils.output_wcs([input_wcs],undistort=False)

        adrizzle.do_driz(insci, input_wcs, inwht,
                         output_wcs, outsci, outwht, outcon,
                         expin, in_units, wt_scl, wcslin_pscale=wcslin.pscale)

        output_bounds = self.bound_image(outsci)
        self.write_image(output, output_wcs, outsci, outwht, outcon[0])

        template_data = self.read_image(output_template)
        template_bounds = self.bound_image(template_data)

        (min_diff, med_diff, max_diff) = self.centroid_statistics("square with point", output_difference,
                                                                  outsci, template_data, 20.0, 8)

        assert(med_diff < 1.0e-6)
        assert(max_diff < 1.0e-5)

    def test_square_with_grid(self):
        """
        Test do_driz square kernel with grid
        """
        input = os.path.basename(self.get_input_file('input','j8bt06nyq_unit.fits'))
        output = 'output_square_grid.fits'
        output_difference = 'difference_square_grid.txt'
        output_template = os.path.basename(self.get_data('truth',
                                           'reference_square_grid.fits'))

        insci = self.read_image(input)
        input_wcs = self.read_wcs(input)
        insci = self.make_grid_image(insci, 64, 100.0)
        inwht = np.ones(insci.shape,dtype=insci.dtype)
        output_wcs = self.read_wcs(output_template)

        naxis1, naxis2 = output_wcs.pixel_shape
        outsci = np.zeros((naxis2, naxis1), dtype='float32')
        outwht = np.zeros((naxis2, naxis1), dtype='float32')
        outcon = np.zeros((1, naxis2, naxis1), dtype='i4')

        expin = 1.0
        wt_scl = expin
        in_units = 'cps'
        wcslin = distortion.utils.output_wcs([input_wcs],undistort=False)

        adrizzle.do_driz(insci, input_wcs, inwht,
                         output_wcs, outsci, outwht, outcon,
                         expin, in_units, wt_scl, wcslin_pscale=wcslin.pscale)

        self.write_image(output, output_wcs, outsci, outwht, outcon[0])

        template_data = self.read_image(output_template)

        (min_diff, med_diff, max_diff) = self.centroid_statistics("square with grid", output_difference,
                                                                  outsci, template_data, 20.0, 8)

        assert(med_diff < 1.0e-6)
        assert(max_diff < 1.0e-5)

    def test_turbo_with_grid(self):
        """
        Test do_driz turbo kernel with grid
        """
        input = os.path.basename(self.get_input_file('input','j8bt06nyq_unit.fits'))
        output = 'output_turbo_grid.fits'
        output_difference = os.path.basename(self.get_data('truth',
                                             'difference_turbo_grid.txt'))
        output_template = os.path.basename(self.get_data('truth',
                                           'reference_turbo_grid.fits'))

        insci = self.read_image(input)
        input_wcs = self.read_wcs(input)
        insci = self.make_grid_image(insci, 64, 100.0)
        inwht = np.ones(insci.shape,dtype=insci.dtype)
        output_wcs = self.read_wcs(output_template)

        naxis1, naxis2 = output_wcs.pixel_shape
        outsci = np.zeros((naxis2, naxis1), dtype='float32')
        outwht = np.zeros((naxis2, naxis1), dtype='float32')
        outcon = np.zeros((1, naxis2, naxis1), dtype='i4')

        expin = 1.0
        wt_scl = expin
        in_units = 'cps'
        wcslin = distortion.utils.output_wcs([input_wcs],undistort=False)

        adrizzle.do_driz(insci, input_wcs, inwht,
                         output_wcs, outsci, outwht, outcon,
                         expin, in_units, wt_scl,
                         kernel='turbo', wcslin_pscale=wcslin.pscale)

        self.write_image(output, output_wcs, outsci, outwht, outcon[0])

        template_data = self.read_image(output_template)

        (min_diff, med_diff, max_diff) = self.centroid_statistics("turbo with grid", output_difference,
                                                                  outsci, template_data, 20.0, 8)

        assert(med_diff < 1.0e-6)
        assert(max_diff < 1.0e-5)

    def test_gaussian_with_grid(self):
        """
        Test do_driz gaussian kernel with grid
        """
        input = os.path.basename(self.get_input_file('input','j8bt06nyq_unit.fits'))
        output = 'output_gaussian_grid.fits'
        output_difference = os.path.basename(self.get_data('truth',
                                             'difference_gaussian_grid.txt'))
        output_template = os.path.basename(self.get_data('truth',
                                           'reference_gaussian_grid.fits'))

        insci = self.read_image(input)
        input_wcs = self.read_wcs(input)
        insci = self.make_grid_image(insci, 64, 100.0)
        inwht = np.ones(insci.shape,dtype=insci.dtype)
        output_wcs = self.read_wcs(output_template)

        naxis1, naxis2 = output_wcs.pixel_shape
        outsci = np.zeros((naxis2, naxis1), dtype='float32')
        outwht = np.zeros((naxis2, naxis1), dtype='float32')
        outcon = np.zeros((1, naxis2, naxis1), dtype='i4')

        expin = 1.0
        wt_scl = expin
        in_units = 'cps'
        wcslin = distortion.utils.output_wcs([input_wcs],undistort=False)

        adrizzle.do_driz(insci, input_wcs, inwht,
                         output_wcs, outsci, outwht, outcon,
                         expin, in_units, wt_scl,
                         kernel='gaussian', wcslin_pscale=wcslin.pscale)

        self.write_image(output, output_wcs, outsci, outwht, outcon[0])

        template_data = self.read_image(output_template)

        (min_diff, med_diff, max_diff) = self.centroid_statistics("gaussian with grid", output_difference,
                                                                  outsci, template_data, 20.0, 8)

        assert(med_diff < 1.0e-6)
        assert(max_diff < 2.0e-5)

    def test_lanczos_with_grid(self):
        """
        Test do_driz lanczos kernel with grid
        """
        input = os.path.basename(self.get_input_file('input', 'j8bt06nyq_unit.fits'))
        output = 'output_lanczos_grid.fits'
        output_difference = os.path.basename(self.get_data('truth',
                                             'difference_lanczos_grid.txt'))
        output_template = os.path.basename(self.get_data('truth',
                                           'reference_lanczos_grid.fits'))

        insci = self.read_image(input)
        input_wcs = self.read_wcs(input)
        insci = self.make_grid_image(insci, 64, 100.0)
        inwht = np.ones(insci.shape,dtype=insci.dtype)
        output_wcs = self.read_wcs(output_template)

        naxis1, naxis2 = output_wcs.pixel_shape
        outsci = np.zeros((naxis2, naxis1), dtype='float32')
        outwht = np.zeros((naxis2, naxis1), dtype='float32')
        outcon = np.zeros((1, naxis2, naxis1), dtype='i4')

        expin = 1.0
        wt_scl = expin
        in_units = 'cps'
        wcslin = distortion.utils.output_wcs([input_wcs],undistort=False)

        adrizzle.do_driz(insci, input_wcs, inwht,
                         output_wcs, outsci, outwht, outcon,
                         expin, in_units, wt_scl,
                         kernel='lanczos3', wcslin_pscale=wcslin.pscale)

        self.write_image(output, output_wcs, outsci, outwht, outcon[0])

        template_data = self.read_image(output_template)

        (min_diff, med_diff, max_diff) = self.centroid_statistics("lanczos with grid", output_difference,
                                                                  outsci, template_data, 20.0, 8)

        assert(med_diff < 1.0e-6)
        assert(max_diff < 1.0e-5)

    def test_tophat_with_grid(self):
        """
        Test do_driz tophat kernel with grid
        """
        input = os.path.basename(self.get_input_file('input', 'j8bt06nyq_unit.fits'))
        output = 'output_tophat_grid.fits'
        output_difference = os.path.basename(self.get_data('truth',
                                             'difference_tophat_grid.txt'))
        output_template = os.path.basename(self.get_data('truth',
                                           'reference_tophat_grid.fits'))

        insci = self.read_image(input)
        input_wcs = self.read_wcs(input)
        insci = self.make_grid_image(insci, 64, 100.0)
        inwht = np.ones(insci.shape,dtype=insci.dtype)
        output_wcs = self.read_wcs(output_template)

        naxis1, naxis2 = output_wcs.pixel_shape
        outsci = np.zeros((naxis2, naxis1), dtype='float32')
        outwht = np.zeros((naxis2, naxis1), dtype='float32')
        outcon = np.zeros((1, naxis2, naxis1), dtype='i4')

        expin = 1.0
        wt_scl = expin
        in_units = 'cps'
        wcslin = distortion.utils.output_wcs([input_wcs],undistort=False)

        adrizzle.do_driz(insci, input_wcs, inwht,
                         output_wcs, outsci, outwht, outcon,
                         expin, in_units, wt_scl,
                         kernel='tophat', wcslin_pscale=wcslin.pscale)

        self.write_image(output, output_wcs, outsci, outwht, outcon[0])

        template_data = self.read_image(output_template)

        (min_diff, med_diff, max_diff) = self.centroid_statistics("tophat with grid", output_difference,
                                                                  outsci, template_data, 20.0, 8)

        assert(med_diff < 1.0e-6)
        assert(max_diff < 1.0e-5)

    def test_point_with_grid(self):
        """
        Test do_driz point kernel with grid
        """
        input = os.path.basename(self.get_input_file('input', 'j8bt06nyq_unit.fits'))
        output = 'output_point_grid.fits'
        output_difference = os.path.basename(self.get_data('truth',
                                         'difference_point_grid.txt'))
        output_template = os.path.basename(self.get_data('truth',
                                           'reference_point_grid.fits'))

        insci = self.read_image(input)
        input_wcs = self.read_wcs(input)
        insci = self.make_grid_image(insci, 64, 100.0)
        inwht = np.ones(insci.shape,dtype=insci.dtype)
        output_wcs = self.read_wcs(output_template)

        naxis1, naxis2 = output_wcs.pixel_shape
        outsci = np.zeros((naxis2, naxis1), dtype='float32')
        outwht = np.zeros((naxis2, naxis1), dtype='float32')
        outcon = np.zeros((1, naxis2, naxis1), dtype='i4')

        expin = 1.0
        wt_scl = expin
        in_units = 'cps'
        wcslin = distortion.utils.output_wcs([input_wcs],undistort=False)

        adrizzle.do_driz(insci, input_wcs, inwht,
                         output_wcs, outsci, outwht, outcon,
                         expin, in_units, wt_scl,
                         kernel='point', wcslin_pscale=wcslin.pscale)

        self.write_image(output, output_wcs, outsci, outwht, outcon[0])

        template_data = self.read_image(output_template)

        (min_diff, med_diff, max_diff) = self.centroid_statistics("point with grid", output_difference,
                                                                  outsci, template_data, 20.0, 8)

        assert(med_diff < 1.0e-6)
        assert(max_diff < 1.0e-5)

    def test_square_with_image(self):
        """
        Test do_driz square kernel
        """
        input = os.path.basename(self.get_input_file('input', 'j8bt06nyq_unit.fits'))
        output = 'output_square_image.fits'
        output_template = os.path.basename(self.get_data('truth',
                                           'reference_square_image.fits'))

        insci = self.read_image(input)
        input_wcs = self.read_wcs(input)
        inwht = np.ones(insci.shape,dtype=insci.dtype)

        output_wcs = self.read_wcs(output_template)
        naxis1, naxis2 = output_wcs.pixel_shape
        outsci = np.zeros((naxis2, naxis1), dtype='float32')
        outwht = np.zeros((naxis2, naxis1), dtype='float32')
        outcon = np.zeros((1, naxis2, naxis1), dtype='i4')

        expin = 1.0
        wt_scl = expin
        in_units = 'cps'
        wcslin = distortion.utils.output_wcs([input_wcs],undistort=False)

        adrizzle.do_driz(insci, input_wcs, inwht,
                         output_wcs, outsci, outwht, outcon,
                         expin, in_units, wt_scl, wcslin_pscale=wcslin.pscale)

        self.write_image(output, output_wcs, outsci, outwht, outcon[0])

        template_data = self.read_image(output_template)

        self.ignore_keywords += ['rootname']
        self.compare_outputs([(output, output_template)])

    def test_turbo_with_image(self):
        """
        Test do_driz turbo kernel
        """
        input = os.path.basename(self.get_input_file('input', 'j8bt06nyq_unit.fits'))
        output = 'output_turbo_image.fits'
        output_template = os.path.basename(self.get_data('truth',
                                           'reference_turbo_image.fits'))

        insci = self.read_image(input)
        input_wcs = self.read_wcs(input)
        inwht = np.ones(insci.shape,dtype=insci.dtype)

        output_wcs = self.read_wcs(output_template)
        naxis1, naxis2 = output_wcs.pixel_shape
        outsci = np.zeros((naxis2, naxis1), dtype='float32')
        outwht = np.zeros((naxis2, naxis1), dtype='float32')
        outcon = np.zeros((1, naxis2, naxis1), dtype='i4')

        expin = 1.0
        wt_scl = expin
        in_units = 'cps'
        wcslin = distortion.utils.output_wcs([input_wcs],undistort=False)

        adrizzle.do_driz(insci, input_wcs, inwht, output_wcs, outsci,
                         outwht, outcon, expin, in_units, wt_scl,
                         kernel='turbo', wcslin_pscale=wcslin.pscale)

        self.write_image(output, output_wcs, outsci, outwht, outcon[0])

        template_data = self.read_image(output_template)

        self.ignore_keywords += ['rootname']
        self.compare_outputs([(output, output_template)])

class TestBlot(BaseUnit):
    buff = 1

    def test_blot_with_point(self):
        """
        Test do_blot with point image
        """
        input = os.path.basename(self.get_input_file('input', 'j8bt06nyq_unit.fits'))
        output = 'output_blot_point.fits'
        output_difference = os.path.basename(self.get_data('truth',
                                         'difference_blot_point.txt'))
        output_template = os.path.basename(self.get_data('truth',
                                           'reference_blot_point.fits'))

        insci = self.read_image(input)
        input_wcs = self.read_wcs(input)
        insci = self.make_point_image(insci, (500, 200), 40.0)
        output_wcs = self.read_wcs(output_template)

        expin = 1.0
        outsci = ablot.do_blot(insci, input_wcs, output_wcs, expin, coeffs = False)
        output_bounds = self.bound_image(outsci)

        self.write_image(output, output_wcs, outsci)

        template_data = self.read_image(output_template)
        template_bounds = self.bound_image(template_data)

        (min_diff, med_diff, max_diff) = self.centroid_statistics("blot with point", output_difference,
                                                                  outsci, template_data, 20.0, 16)

        assert(med_diff < 1.0e-6)
        assert(max_diff < 1.0e-5)

    def test_blot_with_default(self):
        """
        Test do_blot with default grid image
        """
        input = os.path.basename(self.get_input_file('input', 'j8bt06nyq_unit.fits'))
        output = 'output_blot_default.fits'
        output_difference = os.path.basename(self.get_data('truth',
                                         'difference_blot_default.txt'))
        output_template = os.path.basename(self.get_data('truth',
                                           'reference_blot_default.fits'))

        insci = self.read_image(input)
        insci = self.make_grid_image(insci, 64, 100.0)
        input_wcs = self.read_wcs(input)
        output_wcs = self.read_wcs(output_template)

        expin = 1.0
        outsci = ablot.do_blot(insci, input_wcs, output_wcs, expin, coeffs = False)

        self.write_image(output, output_wcs, outsci)

        template_data = self.read_image(output_template)

        (min_diff, med_diff, max_diff) = self.centroid_statistics("blot with defaults", output_difference,
                                                                  outsci, template_data, 20.0, 16)

        assert(med_diff < 1.0e-6)
        assert(max_diff < 1.0e-5)

    def test_blot_with_lan3(self):
        """
        Test do_blot with lan3 grid image
        """
        input = os.path.basename(self.get_input_file('input', 'j8bt06nyq_unit.fits'))
        output = 'output_blot_lan3.fits'
        output_difference = os.path.basename(self.get_data('truth',
                                             'difference_blot_lan3.txt'))
        output_template = os.path.basename(self.get_data('truth',
                                           'reference_blot_lan3.fits'))

        insci = self.read_image(input)
        insci = self.make_grid_image(insci, 64, 100.0)
        input_wcs = self.read_wcs(input)
        output_wcs = self.read_wcs(output_template)

        expin = 1.0
        outsci = ablot.do_blot(insci, input_wcs, output_wcs, expin, interp="lan3", coeffs = False)

        self.write_image(output, output_wcs, outsci)

        template_data = self.read_image(output_template)

        (min_diff, med_diff, max_diff) = self.centroid_statistics("blot with lan3", output_difference,
                                                                  outsci, template_data, 20.0, 16)

        assert(med_diff < 1.0e-6)
        assert(max_diff < 1.0e-5)

    def test_blot_with_lan5(self):
        """
        Test do_blot with lan5 grid image
        """
        input = os.path.basename(self.get_input_file('input', 'j8bt06nyq_unit.fits'))
        output = 'output_blot_lan5.fits'
        output_difference = os.path.basename(self.get_data('truth',
                                         'difference_blot_lan5.txt'))
        output_template = os.path.basename(self.get_data('truth',
                                           'reference_blot_lan5.fits'))

        insci = self.read_image(input)
        insci = self.make_grid_image(insci, 64, 100.0)
        input_wcs = self.read_wcs(input)
        output_wcs = self.read_wcs(output_template)

        expin = 1.0
        outsci = ablot.do_blot(insci, input_wcs, output_wcs, expin, interp="lan5", coeffs = False)

        self.write_image(output, output_wcs, outsci)

        template_data = self.read_image(output_template)

        (min_diff, med_diff, max_diff) = self.centroid_statistics("blot with lan5", output_difference,
                                                                  outsci, template_data, 20.0, 16)

        assert(med_diff < 1.0e-6)
        assert(max_diff < 1.0e-5)

    @pytest.mark.xfail(reason='Input with different distortion')
    def test_blot_with_image(self):
        """
        Test do_blot with full image
        """
        input = os.path.basename(self.get_input_file('input', 'j8bt06nyq_unit.fits'))
        output = 'output_blot_image.fits'
        output_template = os.path.basename(self.get_data('truth',
                                           'reference_blot_image.fits'))

        insci = self.read_image(input)
        input_wcs = self.read_wcs(input)
        output_wcs = self.read_wcs(output_template)

        expin = 1.0
        outsci = ablot.do_blot(insci, input_wcs, output_wcs, expin, coeffs = False)
        self.write_image(output, output_wcs, outsci)

        self.compare_outputs([(output, output_template)])
