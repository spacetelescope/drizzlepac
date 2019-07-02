#!/usr/bin/env python

"""This script contains code to support creation of photometric sourcelists using two techniques: aperture photometry
segmentation-map based photometry.
"""
import pdb
import sys


from astropy.io import fits
from astropy.stats import mad_std,gaussian_fwhm_to_sigma, gaussian_sigma_to_fwhm
import numpy as np
from photutils import aperture_photometry, CircularAperture, DAOStarFinder
from photutils import Background2D, MedianBackground
from photutils import detect_sources, source_properties
from stsci.tools import logutil

from drizzlepac import util
from drizzlepac.hlautils import astrometric_utils
from drizzlepac.hlautils import se_source_generation

__taskname__ = 'sourcelist_generation_oo'

log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout)


# ----------------------------------------------------------------------------------------------------------------------

class build_catalogs(object):
    """Using aperture photometry, generate photometric sourcelist for specified image(s).
    """
    def __init__(self,fitsfile):
        self.label="build_catalogs"
        self.description="A set of routines to generate photometric sourcelists using aperture photometry"
        self.imgname = fitsfile
        self.point_sourcelist_filename = self.imgname.replace(self.imgname[-9:],"_point-cat.ecsv")
        self.seg_sourcelist_filename = self.imgname.replace(self.imgname[-9:], "_segment-cat.ecsv")
        self.param_dict = {
        "ACS HRC": {
            "astrodrizzle": {
                "SCALE": 0.025,
                "PIXFRAC": 1.0,
                "KERNEL": "square",
                "OUTNX": None,
                "OUTNY": None,
                "ROT": 0.0,
                "BITS": 256},
            "ci filter": {
                "ci_daolower_limit": 0.9,
                "ci_daoupper_limit": 1.6,
                "ci_selower_limit": 0.9,
                "ci_seupper_limit": 1.6},
            "dao": {
                "TWEAK_FWHMPSF": 0.073,
                "TWEAK_THRESHOLD": 3.0,
                "aperture_1": 0.03,
                "aperture_2": 0.125,
                "bthresh": 5.0},
            "sourcex": {
                "fwhm": 0.073,
                "thresh": 1.4,
                "bthresh": 5.0,
                "source_box": 7},
            "swarm filter": {
                "upper_epp_limit": 70000.,
                "lower_epp_limit": 2000.,
                "eppsky_limit": 1000.,
                "swarm_thresh": 1.,
                "clip_radius_list": [120.0, 100.0, 80.0, 60.0, 40.0, 20.0, 10.0, 5.0, 2.0, 0.0],
                "scale_factor_list": [0.0, 1.778106e-05, 3.821292e-05, 9.017166e-05, 2.725184e-04, 1.269197e-03, 7.007126e-03, 3.839166e-02, 2.553349e-01, 1.000000e+00],
                "proximity_binary": "no"}},
        "ACS SBC": {
            "astrodrizzle": {
                "SCALE": 0.03,
                "PIXFRAC": 1.0,
                "KERNEL": "square",
                "OUTNX": None,
                "OUTNY": None,
                "ROT": 0.0,
                "BITS": 256},
            "ci filter": {
                "ci_daolower_limit": 0.15,
                "ci_daoupper_limit": 0.45,
                "ci_selower_limit": 0.15,
                "ci_seupper_limit": 0.45},
            "dao": {
                "TWEAK_FWHMPSF": 0.065,
                "TWEAK_THRESHOLD": 3.0,
                "aperture_1": 0.07,
                "aperture_2": 0.125,
                "bthresh": 5.0},
            "sourcex": {
                "fwhm": 0.065,
                "thresh": 1.4,
                "bthresh": 5.0,
                "source_box": 7},
            "swarm filter": {
                "upper_epp_limit": 70000.,
                "lower_epp_limit": 2000.,
                "eppsky_limit": 1000.,
                "swarm_thresh": 1.,
                "clip_radius_list": [120.0, 100.0, 80.0, 60.0, 40.0, 20.0, 10.0, 5.0, 2.0, 0.0],
                "scale_factor_list": [0.0, 1.778106e-05, 3.821292e-05, 9.017166e-05, 2.725184e-04, 1.269197e-03, 7.007126e-03, 3.839166e-02, 2.553349e-01, 1.000000e+00],
                "proximity_binary": "no"}},
        "ACS WFC": {
            "astrodrizzle": {
                "SCALE": 0.05,
                "PIXFRAC": 1.0,
                "KERNEL": "square",
                "OUTNX": None,
                "OUTNY": None,
                "ROT": 0.0,
                "BITS": 256},
            "ci filter": {
                "ci_daolower_limit": 0.9,
                "ci_daoupper_limit": 1.23,
                "ci_selower_limit": 0.9,
                "ci_seupper_limit": 1.23},
            "dao": {
                "TWEAK_FWHMPSF": 0.076,
                "TWEAK_THRESHOLD": None,
                "aperture_1": 0.05,  # update from 0.15
                "aperture_2": 0.15,  # update from 0.25
                "bthresh": 5.0},
            "sourcex": {
                "fwhm": 0.13,
                "thresh": None,
                "bthresh": 5.0,
                "source_box": 5},
            "swarm filter": {
                "upper_epp_limit": 70000.,
                "lower_epp_limit": 2000.,
                "eppsky_limit": 1000.,
                "swarm_thresh": 1.,
                "clip_radius_list": [120., 100., 80., 60., 40., 30., 20., 10., 5., 2., 0.],
                "scale_factor_list": [0.0, 0.000000e+00, 6.498530e-06, 3.687270e-05, 1.412972e-04, 3.151877e-04, 1.023391e-03, 3.134859e-03, 2.602436e-02, 1.820539e-01, 1.000000e+00],
                "proximity_binary": "no"}},
        "WFC3 IR": {
            "astrodrizzle": {
                "SCALE": 0.09,
                "PIXFRAC": 1.0,
                "KERNEL": "square",
                "OUTNX": None,
                "OUTNY": None,
                "ROT": 0.0,
                "BITS": 768},
            "ci filter": {
                "ci_daolower_limit": 0.25,
                "ci_daoupper_limit": 0.55,
                "ci_selower_limit": 0.25,
                "ci_seupper_limit": 0.55},
            "dao": {
                "TWEAK_FWHMPSF": 0.14,
                "TWEAK_THRESHOLD": 3.0,
                "aperture_1": 0.15,
                "aperture_2": 0.45,
                "bthresh": 5.0},
            "sourcex": {
                "fwhm": 0.14,
                "thresh": 1.4,
                "bthresh": 5.0,
                "source_box": 7},
            "swarm filter": {
                "upper_epp_limit": 70000.,
                "lower_epp_limit": 2000.,
                "eppsky_limit": 100.,
                "swarm_thresh": 1.,
                "clip_radius_list": [140., 120., 100., 80., 60., 40., 20., 10., 5., 2., 0.],
                #                   x10    x10    x10   x10   x10   x10    x10   x10  x10  x2,
                "scale_factor_list": [1.5e-5, 2.3e-5, 4.e-5, 8.e-5, 2.e-4, 0.0006, 0.015, 0.05, 0.15, 0.9, 1.],
                # "scale_factor_list_orig": [1.5e-5, 2.3e-5, 4.e-5, 8.e-5, 2.e-4, 0.0006, 0.005, 0.05, 0.15, 0.9, 1.],
                "proximity_binary": "yes"}},
        "WFC3 UVIS": {
            "astrodrizzle": {
                "SCALE": 0.04,
                "PIXFRAC": 1.0,
                "KERNEL": "square",
                "OUTNX": None,
                "OUTNY": None,
                "ROT": 0.0,
                "BITS": 256},
            "ci filter": {
                "ci_daolower_limit": 0.75,
                "ci_daoupper_limit": 1.0,
                "ci_selower_limit": 0.75,
                "ci_seupper_limit": 1.0},
            "dao": {
                "TWEAK_FWHMPSF": 0.076,
                "TWEAK_THRESHOLD": 3.0,
                "aperture_1": 0.05,
                "aperture_2": 0.15,
                "bthresh": 5.0},
            "sourcex": {
                "fwhm": 0.076,
                "thresh": 1.4,
                "bthresh": 5.0,
                "source_box": 7},
            "swarm filter": {
                "upper_epp_limit": 70000.,
                "lower_epp_limit": 2000.,
                "eppsky_limit": 1000.,
                "swarm_thresh": 1.,
                "clip_radius_list": [120., 100., 80., 60., 40., 20., 10., 5., 2., 0.],
                "scale_factor_list": [2.3e-6, 4.e-6, 8.e-6, 2.e-5, 0.0005, 0.005, 0.005, 0.015, 0.45, 1.],
                # "scale_factor_list_orig": [2.3e-6, 4.e-6, 8.e-6, 2.e-5, 6.e-5, 0.0005, 0.005, 0.015, 0.45, 1.],
                "proximity_binary": "yes"}}} # TODO: remove para_dict definition once we have fleshed out the config object
        self.inst_det = "{} {}".format(self.imgname.split("_")[3].upper(),self.imgname.split("_")[4].upper())
        self.param_dict = self.param_dict[self.inst_det] # TODO: remove para_dict redefinition once we have fleshed out the config object


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def identify_point_sources(self,fitsfile,param_dict,bkgsig_sf=4.,dao_ratio=0.8):
        """Create a master coordinate list of sources identified in the specified total detection product image

        Parameters
        ----------
        fitsfile : string
            Name of the drizzle-combined filter product to used to generate photometric sourcelists.

        param_dict : dictionary
            Dictionary of instrument/detector - specific drizzle, source finding and photometric parameters

        dao_fwhm : float
            (photutils.DAOstarfinder param 'fwhm') The full-width half-maximum (FWHM) of the major axis of the
            Gaussian kernel in units of pixels. Default value = 3.5.

        bkgsig_sf : float
            multiplictive scale factor applied to background sigma value to compute DAOfind input parameter
            'threshold'. Default value = 2.

        dao_ratio : float
            The ratio of the minor to major axis standard deviations of the Gaussian kernel.

        Returns
        -------
        sources : astropy table
            Table containing x, y coordinates of identified sources
        """
        # read in sci, wht extensions of drizzled product

        hdulist = fits.open(fitsfile)
        image = hdulist['SCI'].data
        image -= np.nanmedian(image)
        wht_image = hdulist['WHT'].data

        bkg_sigma = mad_std(image, ignore_nan=True)

        detect_sources_thresh = bkgsig_sf * bkg_sigma
        default_fwhm = param_dict['dao']['TWEAK_FWHMPSF'] / param_dict['astrodrizzle']['SCALE']

        # Estimate background for DaoStarfinder 'threshold' input.
        bkg_estimator = MedianBackground()
        bkg = None
        threshold = param_dict['dao']['TWEAK_THRESHOLD']
        exclude_percentiles = [10, 25, 50, 75]
        for percentile in exclude_percentiles:
            try:
                bkg = Background2D(image, (50, 50), filter_size=(3, 3),
                                   bkg_estimator=bkg_estimator,
                                   exclude_percentile=percentile)
            except Exception:
                bkg = None
                continue
            if bkg is not None:
                # If it succeeds, stop and use that value
                bkg_rms = (5. * bkg.background_rms)
                bkg_rms_mean = bkg.background.mean() + 5. * bkg_rms.std()
                default_threshold = bkg.background + bkg_rms
                if threshold is None:
                    threshold = default_threshold
                elif threshold < 0:
                    threshold = -1 * threshold * default_threshold
                    log.info("{} based on {}".format(threshold.max(), default_threshold.max()))
                    bkg_rms_mean = threshold.max()
                else:
                    bkg_rms_mean = 3. * threshold

                if bkg_rms_mean < 0:
                    bkg_rms_mean = 0.
                break

        # If Background2D does not work at all, define default scalar values for
        # the background to be used in source identification
        if bkg is None:
            bkg_rms_mean = max(0.01, imgarr.min())
            bkg_rms = bkg_rms_mean * 5



        # Estimate FWHM from image sources
        kernel = astrometric_utils.build_auto_kernel(image, wht_image, threshold=bkg_rms, fwhm=default_fwhm)
        segm = detect_sources(image, detect_sources_thresh, npixels=param_dict["sourcex"]["source_box"], filter_kernel=kernel)
        cat = source_properties(image, segm)
        source_table = cat.to_table()
        smajor_sigma = source_table['semimajor_axis_sigma'].mean().value
        source_fwhm = smajor_sigma * gaussian_sigma_to_fwhm

        log.info("DAOStarFinder(fwhm={}, threshold={}, ratio={})".format(source_fwhm,bkg_rms_mean,bkg_rms_mean))
        daofind = DAOStarFinder(fwhm=source_fwhm, threshold=bkg_rms_mean, ratio=dao_ratio)
        sources = daofind(image)
        hdulist.close()

        for col in sources.colnames:
            sources[col].info.format = '%.8g'  # for consistent table output

        return(sources)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def perform_point_photometry(self,fitsfile,sources,aper_radius=4.):
        """Perform aperture photometry on identified sources

        Parameters
        ----------
        fitsfile : string
            Name of the drizzle-combined filter product to used to generate photometric sourcelists.

        sources : astropy table
            Table containing x, y coordinates of identified sources

        aper_radius : float
            Aperture radius (in pixels) used for photometry. Default value = 4.

        Returns
        -------
        phot_table : astropy table
            Table containing photometric information for specified sources based on image data in the specified image.
        """
        # Open and background subtract image
        hdulist = fits.open(fitsfile)
        image = hdulist['SCI'].data
        image -= np.nanmedian(image)


        # Aperture Photometry
        positions = (sources['xcentroid'], sources['ycentroid'])
        apertures = CircularAperture(positions, r=aper_radius)
        phot_table = aperture_photometry(image, apertures)

        for col in phot_table.colnames: phot_table[col].info.format = '%.8g'  # for consistent table output
        hdulist.close()
        return(phot_table)
        # # Write out sourcelist
        # tbl_length = len(phot_table)
        # phot_table.write(sl_filename, format="ascii.ecsv")
        # log.info("Created sourcelist file '{}' with {} sources".format(sl_filename, tbl_length))
        #
        # # Write out ds9-compatable .reg file
        # if make_region_file:
        #     reg_filename = sl_filename.replace(".ecsv",".reg")
        #     out_table = phot_table.copy()
        #     out_table['xcenter'].data = out_table['xcenter'].data + np.float64(1.0)
        #     out_table['ycenter'].data = out_table['ycenter'].data + np.float64(1.0)
        #     out_table.remove_column('id')
        #     out_table.write(reg_filename, format="ascii")
        #     log.info("Created region file '{}' with {} sources".format(reg_filename, tbl_length))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def write_catalog_to_file(self,catalog,catalog_filename,write_region_file=False):
        """Write specified catalog to file on disk

        Parameters
        ----------
        catalog : astropy table
            table data to write to disk

        catalog_filename : string
            output filename that catalog data will be written to.

        write_region_file : Boolean
           Write ds9-compatible region file along with the catalog file? Default value = False

        Returns
        -------
        Nothing!

        """
        # Write out catalog to ecsv file
        catalog.write(catalog_filename, format="ascii.ecsv")
        log.info("Wrote catalog file '{}' containing {} sources".format(catalog_filename, len(catalog)))

        # Write out region file if input 'write_region_file' is turned on.
        if write_region_file:
            out_table = catalog.copy()
            if 'xcentroid' in out_table.keys(): # for point-source source catalogs
                # Remove all other columns besides xcentroid and ycentroid
                out_table.keep_columns(['xcentroid','ycentroid'])
                # Add offset of 1.0 in X and Y to line up sources in region file with image displayed in ds9.
                out_table['xcentroid'].data[:] += np.float64(1.0)
                out_table['ycentroid'].data[:] += np.float64(1.0)
            elif 'xcenter' in out_table.keys(): # for point-source photometric catalogs
                # Remove all other columns besides xcenter and ycenter
                out_table.keep_columns(['xcenter', 'ycenter'])
                # Add offset of 1.0 in X and Y to line up sources in region file with image displayed in ds9.
                out_table['xcenter'].data = out_table['xcenter'].data + np.float64(1.0)
                out_table['ycenter'].data = out_table['ycenter'].data + np.float64(1.0)
            else: # Bail out if anything else is encountered.
                log.info("Error: unrecognized catalog format. Skipping region file generation.")
                return()
            reg_filename = catalog_filename.replace(".ecsv",".reg")
            out_table.write(reg_filename, format="ascii")
            log.info("Wrote region file '{}' containing {} sources".format(reg_filename, len(out_table)))


# ----------------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    """Super simple testing interface for the above code."""
    import argparse
    parser = argparse.ArgumentParser(description='test interface for sourcelist_generation')
    parser.add_argument('total_product_name',help="total product filename")
    parser.add_argument('-f', '--filter_product_list',nargs='+',
                        help="Space-seperated list of one or more total filter products")
    args = parser.parse_args()

    total_product = build_catalogs(args.total_product_name)
    total_product.ps_source_cat = total_product.identify_point_sources(total_product.imgname,total_product.param_dict)
    total_product.write_catalog_to_file(total_product.ps_source_cat,
                                        total_product.point_sourcelist_filename,
                                        write_region_file=True)

    for filter_img_name in args.filter_product_list:
        filter_product = build_catalogs(filter_img_name)
        filter_product.ps_phot_cat = filter_product.perform_point_photometry(filter_product.imgname,
                                                                             total_product.ps_source_cat)
        filter_product.write_catalog_to_file(filter_product.ps_phot_cat,
                                             filter_product.point_sourcelist_filename,
                                             write_region_file=True)
