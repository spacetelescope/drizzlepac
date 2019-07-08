#!/usr/bin/env python

"""This script contains code to support creation of photometric sourcelists using two techniques: aperture photometry
segmentation-map based photometry.
"""
import argparse
import datetime
import os
import pdb
import sys

import astropy.units as u
from astropy.io import ascii
from astropy.io import fits as fits
from astropy.convolution import Gaussian2DKernel, MexicanHat2DKernel
from astropy.stats import mad_std, gaussian_fwhm_to_sigma, gaussian_sigma_to_fwhm
from astropy.table import Column, Table
import numpy as np
import photutils
from photutils import aperture_photometry, CircularAperture, DAOStarFinder
from photutils import Background2D, MedianBackground, SExtractorBackground, StdBackgroundRMS
from photutils import detect_sources, source_properties, deblend_sources
from stsci.tools import logutil
from stwcs.wcsutil import HSTWCS

from drizzlepac import util
from drizzlepac.hlautils import astrometric_utils


try:
    from matplotlib import pyplot as plt
except Exception:
    plt = None

__taskname__ = 'catalog_utils'

log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout)


# ----------------------------------------------------------------------------------------------------------------------

class build_catalogs(object):
    """Using aperture photometry, generate photometric sourcelist for specified image(s).
    """
    def __init__(self,fitsfile):
        self.label="build_catalogs"
        self.description="A set of routines to generate photometric sourcelists using aperture photometry"
        
        # Filename stuff
        self.imgname = fitsfile
        self.point_sourcelist_filename = self.imgname.replace(self.imgname[-9:],"_point-cat.ecsv")
        self.seg_sourcelist_filename = self.imgname.replace(self.imgname[-9:], "_segment-cat.ecsv")
        
        # Fits file read
        self.imghdu = fits.open(self.imgname)
        
        # Parameter dictionary definition
        self.inst_det = "{} {}".format(self.imgname.split("_")[3].upper(), self.imgname.split("_")[4].upper())
        self.full_param_dict = {
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
        self.param_dict=self.full_param_dict[self.inst_det].copy() # TODO: remove para_dict redefinition once we have fleshed out the config object


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def identify_point_sources(self,bkgsig_sf=4.,dao_ratio=0.8):
        """Create a master coordinate list of sources identified in the specified total detection product image

        Parameters
        ----------
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
        image = self.imghdu['SCI'].data
        image -= np.nanmedian(image)
        wht_image = self.imghdu['WHT'].data

        bkg_sigma = mad_std(image, ignore_nan=True)

        detect_sources_thresh = bkgsig_sf * bkg_sigma
        default_fwhm = self.param_dict['dao']['TWEAK_FWHMPSF'] / self.param_dict['astrodrizzle']['SCALE']

        # Estimate background for DaoStarfinder 'threshold' input.
        bkg_estimator = MedianBackground()
        bkg = None
        threshold = self.param_dict['dao']['TWEAK_THRESHOLD']
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
        segm = detect_sources(image, detect_sources_thresh, npixels=self.param_dict["sourcex"]["source_box"],
                              filter_kernel=kernel)
        cat = source_properties(image, segm)
        source_table = cat.to_table()
        smajor_sigma = source_table['semimajor_axis_sigma'].mean().value
        source_fwhm = smajor_sigma * gaussian_sigma_to_fwhm

        log.info("DAOStarFinder(fwhm={}, threshold={}, ratio={})".format(source_fwhm,bkg_rms_mean,bkg_rms_mean))
        daofind = DAOStarFinder(fwhm=source_fwhm, threshold=bkg_rms_mean, ratio=dao_ratio)
        sources = daofind(image)

        for col in sources.colnames:
            sources[col].info.format = '%.8g'  # for consistent table output

        return(sources)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def perform_point_photometry(self,sources,aper_radius=4.):
        """Perform aperture photometry on identified sources

        Parameters
        ----------
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
        image = self.imghdu['SCI'].data
        image -= np.nanmedian(image)


        # Aperture Photometry
        positions = (sources['xcentroid'], sources['ycentroid'])
        apertures = CircularAperture(positions, r=aper_radius)
        phot_table = aperture_photometry(image, apertures)

        for col in phot_table.colnames: phot_table[col].info.format = '%.8g'  # for consistent table output
        return(phot_table)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def write_catalog_to_file(self,catalog,write_region_file=False):
        """Write specified catalog to file on disk

        Parameters
        ----------
        catalog : astropy table
            table data to write to disk

        write_region_file : Boolean
           Write ds9-compatible region file along with the catalog file? Default value = False

        Returns
        -------
        Nothing!

        """
        # Write out catalog to ecsv file
        catalog.write(self.point_sourcelist_filename, format="ascii.ecsv")
        log.info("Wrote catalog file '{}' containing {} sources".format(self.point_sourcelist_filename, len(catalog)))

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
            reg_filename = self.point_sourcelist_filename.replace(".ecsv",".reg")
            out_table.write(reg_filename, format="ascii")
            log.info("Wrote region file '{}' containing {} sources".format(reg_filename, len(out_table)))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def check_param_dict(self):
        print("="*100)
        print("=" * 100)
        print("                                 param_dict check!")
        overall_status = "ALL PARAMS OK"
        for key1 in self.param_dict.keys():
            for key2 in self.param_dict[key1].keys():
                if self.param_dict[key1][key2] == self.full_param_dict["ACS WFC"][key1][key2]:
                    status = "OK  "
                else:
                    status = "BAD!"
                    overall_status = "PROBLEMS FOUND\a"
                print(status,key1,key2,self.param_dict[key1][key2],self.full_param_dict["ACS WFC"][key1][key2])

            print("\n")
        print(overall_status)
        if overall_status.startswith("PROBLEMS"):
            pdb.set_trace()
        print("=" * 100)
        print("=" * 100)

# ----------------------------------------------------------------------------------------------------------------------
#       Contents of Michele's se_source_generation.py, as of commit b2db3ec9c918188cea2d3b0e4b64e39cc79c4146
# ----------------------------------------------------------------------------------------------------------------------
    def create_sextractor_like_sourcelists(self,catalog_filename, se_debug=False):
        """Use photutils to find sources in image based on segmentation.

        Parameters
        ----------
        catalog_filename : string
            Name of the output source catalog for the total detection product

        self.param_dict : dictionary
            dictionary of drizzle, source finding, and photometric parameters

        se_debug : bool, optional
            Specify whether or not to plot the image and segmentation image for
            visualization and debugging purposes

        Returns
        -------
        segm : `photutils.segmentation.SegmentationImage`
            Two-dimensional segmentation image where found source regions are labeled with
            unique, non-zero positive integers.

        kernel :

        bkg : `~photutils.background.Background2D` or None
            A background map based upon the `~photutils.background.SExtractorBackground`
            estimator

        bkg_rms_mean : float
            Mean bkg.background FIX

        """
        # get the TDP SCI image data

        imgarr = self.imghdu['sci', 1].data

        # Get the HSTWCS object from the first extension
        imgwcs = HSTWCS(self.imghdu, 1)

        # Get header information to annotate the output catalogs
        keyword_dict = self._get_header_data()

        # Get the instrument/detector-specific values from the self.param_dict
        fwhm = self.param_dict["sourcex"]["fwhm"]
        size_source_box = self.param_dict["sourcex"]["source_box"]
        threshold_flag = self.param_dict["sourcex"]["thresh"]

        # Report configuration values to log
        log.info("{}".format("=" * 80))
        log.info("")
        log.info("SExtractor-like source finding settings for Photutils segmentation")
        log.info("Total Detection Product - Input Parameters")
        log.info("FWHM: {}".format(fwhm))
        log.info("size_source_box: {}".format(size_source_box))
        log.info("threshold_flag: {}".format(threshold_flag))
        log.info("")
        log.info("{}".format("=" * 80))

        # Only use a single kernel for now
        kernel_list = [Gaussian2DKernel, MexicanHat2DKernel]
        kernel_in_use = kernel_list[0]

        bkg, bkg_dao_rms, threshold = self._compute_background(imgarr, nsigma=5., threshold_flag=threshold_flag)

        # FIX imgarr should be background subtracted, sextractor uses the filtered_data image
        imgarr_bkgsub = imgarr - bkg.background

        # *** FIX: should size_source_box size be used in all these places? ***
        # Create a 2D filter kernel - this will be used to smooth the input
        # image prior to thresholding in detect_sources().
        sigma = fwhm * gaussian_fwhm_to_sigma
        kernel = kernel_in_use(sigma, x_size=size_source_box, y_size=size_source_box)
        kernel.normalize()

        # Source segmentation/extraction
        # If the threshold includes the background level, then the input image
        # should NOT be background subtracted.
        # Note: SExtractor has "connectivity=8" which is the default for this function
        segm = detect_sources(imgarr, threshold, npixels=size_source_box, filter_kernel=kernel)

        # For debugging purposes...
        if se_debug:
            # Write out a catalog which can be used as an overlay for image in ds9
            cat = source_properties(imgarr_bkgsub, segm, background=bkg.background, filter_kernel=kernel, wcs=imgwcs)
            table = cat.to_table()

            # Copy out only the X and Y coordinates to a "debug table" and
            # cast as an Astropy Table
            tbl = Table(table["xcentroid", "ycentroid"])

            # Construct the debug output filename and write the catalog
            indx = catalog_filename.find("ecsv")
            outname = catalog_filename[0:indx] + "reg"

            tbl["xcentroid"].info.format = ".10f"  # optional format
            tbl["ycentroid"].info.format = ".10f"

            # Add one to the X and Y table values to put the data onto a one-based system,
            # particularly for display with DS9
            tbl["xcentroid"] = tbl["xcentroid"] + 1
            tbl["ycentroid"] = tbl["ycentroid"] + 1
            tbl.write(outname, format="ascii.commented_header")
            log.info("Wrote debug source catalog: {}".format(outname))

            """
            # Generate a graphic of the image and the segmented image
            norm = ImageNormalize(stretch=SqrtStretch())
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))
            ax1.imshow(imgarr, origin="lower", cmap="Greys_r", norm=norm)
            ax1.set_title("Data")
            ax2.imshow(segm, origin="lower", cmap=segm.cmap(random_state=12345))
            ax2.set_title("Segmentation Image")
            plt.show()
            """

        # TROUBLESOME at this time
        # Deblending is a combination of multi-thresholding and watershed
        # segmentation. Sextractor uses a multi-thresholding technique.
        # npixels = number of connected pixels in source
        # npixels and filter_kernel should match those used by detect_sources()
        # Note: SExtractor has "connectivity=8" which is the default for this function
        """
        segm = deblend_sources(imgarr, segm, npixels=size_source_box,
                               filter_kernel=kernel, nlevels=32,
                               contrast=0.005)
        print("after deblend. ", segm)
        """

        """
        if se_debug:
            # Generate a graphic of the image and the segmented image
            norm = ImageNormalize(stretch=SqrtStretch())
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))
            ax1.imshow(imgarr, origin="lower", cmap="Greys_r", norm=norm)
            ax1.set_title("Data")
            ax2.imshow(segm, origin="lower", cmap=segm.cmap(random_state=12345))
            ax2.set_title("Segmentation Image")
            plt.show()
        """

        # Regenerate the source catalog with presumably now only good sources
        seg_cat = source_properties(imgarr_bkgsub, segm, background=bkg.background, filter_kernel=kernel, wcs=imgwcs)

        self._write_catalog(seg_cat, keyword_dict, catalog_filename)

        return segm, kernel, bkg_dao_rms

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def measure_source_properties(self,segm, kernel, catalog_filename, param_dict):
        """Use the positions of the sources identified in the white light image to
        measure properties of these sources in the filter images

        An instrument/detector combination may have multiple filter-level products.
        This routine is called for each filter image which is then measured to generate
        a filter-level source catalog based on object positions measured in the total
        detection product image.

        Parameters
        ----------
        segm : `~astropy.photutils.segmentation` Segmentation image
            Two-dimensional image of labeled source regions based on the "white light" drizzed product

        kernel : `~astropy.convolution`
            Two dimensional function of a specified FWHM used to smooth the image and
            used in the detection of sources as well as for the determination of the
            source properties (this routine)

        catalog_filename : string
            Name of the output source catalog for the filter detection product

        Returns
        -------

        """

        # get filter-level science data

        imgarr = self.imghdu['sci', 1].data

        # Get the HSTWCS object from the first extension
        imgwcs = HSTWCS(self.imghdu, 1)

        # Get header information to annotate the output catalogs
        keyword_dict = self._get_header_data(product="fdp")

        # Get the instrument/detector-specific values from the param_dict
        fwhm = param_dict["sourcex"]["fwhm"]
        size_source_box = param_dict["sourcex"]["source_box"]
        threshold_flag = param_dict["sourcex"]["thresh"]

        # Report configuration values to log
        log.info("{}".format("=" * 80))
        log.info("")
        log.info("SExtractor-like source property measurements based on Photutils segmentation")
        log.info("Filter Level Product - Input Parameters")
        log.info("FWHM: {}".format(fwhm))
        log.info("size_source_box: {}".format(size_source_box))
        log.info("threshold_flag: {}".format(threshold_flag))
        log.info("")
        log.info("{}".format("=" * 80))

        # The data needs to be background subtracted when computing the source properties
        bkg, _, _ = self._compute_background(imgarr, nsigma=5., threshold_flag=threshold_flag)

        imgarr_bkgsub = imgarr - bkg.background

        # Compute source properties...
        seg_cat = source_properties(imgarr_bkgsub, segm, background=bkg.background, filter_kernel=kernel, wcs=imgwcs)

        # Write the source catalog
        self._write_catalog(seg_cat, keyword_dict, catalog_filename, product="fdp")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def _compute_background(self,image, box_size=50, win_size=3, nsigma=5., threshold_flag=None):
        """Use Background2D to determine the background of the input image.

        Parameters
        ----------
        image : ndarray
            Numpy array of the science extension from the observations FITS file.

        box_size : int
            Size of box along each axis

        win_size : int
            Size of 2D filter to apply to the background image

        nsigma : float
            Number of sigma above background

        threshold_flag : float or None
            Value from the image which serves as the limit for determining sources.
            If None, compute a default value of (background+5*rms(background)).
            If threshold < 0.0, use absolute value as scaling factor for default value.


        Returns
        -------
        bkg : 2D ndarray
            Background image

        bkg_dao_rms : ndarry
            Background RMS image

        threshold : ndarray
            Numpy array representing the background plus RMS

        """
        # Report configuration values to log
        log.info("")
        log.info("Computation of white light image background - Input Parameters")
        log.info("Box size: {}".format(box_size))
        log.info("Window size: {}".format(win_size))
        log.info("NSigma: {}".format(nsigma))

        # SExtractorBackground ans StdBackgroundRMS are the defaults
        bkg_estimator = SExtractorBackground()
        bkgrms_estimator = StdBackgroundRMS()
        bkg = None
        bkg_dao_rms = None

        exclude_percentiles = [10, 25, 50, 75]
        for percentile in exclude_percentiles:
            log.info("")
            log.info("Percentile in use: {}".format(percentile))
            try:
                bkg = Background2D(image, box_size, filter_size=win_size, bkg_estimator=bkg_estimator,
                                   bkgrms_estimator=bkgrms_estimator, exclude_percentile=percentile, edge_method="pad")
            except Exception:
                bkg = None
                continue

            if bkg is not None:
                # Set the bkg_rms at "nsigma" sigma above background
                bkg_rms = nsigma * bkg.background_rms
                default_threshold = bkg.background + bkg_rms
                bkg_rms_mean = bkg.background.mean() + nsigma * bkg_rms.std()
                bkg_dao_rms = bkg.background_rms
                if threshold_flag is None:
                    threshold = default_threshold
                elif threshold_flag < 0:
                    threshold = -1 * threshold_flag * default_threshold
                    log.info("Background threshold set to {} based on {}".format(threshold.max(), default_threshold.max()))
                    bkg_rms_mean = threshold.max()
                else:
                    bkg_rms_mean = 3. * threshold_flag
                    threshold = bkg_rms_mean

                if bkg_rms_mean < 0:
                    bkg_rms_mean = 0.
                break

        # If Background2D does not work at all, define default scalar values for
        # the background to be used in source identification
        if bkg is None:
            bkg_rms_mean = max(0.01, image.min())
            bkg_rms = nsigma * bkg_rms_mean
            bkg_dao_rms = bkg_rms_mean
            threshold = bkg_rms_mean + bkg_rms

        # *** FIX: Need to do something for bkg if bkg is None ***

        # Report other useful quantities
        log.info("")
        log.info("Mean background: {}".format(bkg.background.mean()))
        log.info("Mean threshold: {}".format(bkg_rms_mean))
        log.info("")
        log.info("{}".format("=" * 80))

        return bkg, bkg_dao_rms, threshold


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def _write_catalog(self,seg_cat, keyword_dict, catalog_filename, product="tdp"):
        """Actually write the specified source catalog out to disk

        Parameters
        ----------
        seg_cat : list of `~photutils.SourceProperties` objects
            List of SourceProperties objects, one for each source found in the
            specified detection product

        keyword_dict : dict
            Dictionary containing FITS keyword values pertaining to the data product

        catalog_filename : str
            Official generated name for the output catalog

        product : str, optional
            Identification string for the catalog product being written.  This
            controls the data being put into the catalog product
        """

        # Convert the list of SourceProperties objects to a QTable and
        # document in column metadata Photutils columns which map to SExtractor columns
        seg_table = Table(seg_cat.to_table())
        radec_data = seg_table["sky_centroid_icrs"]
        ra_icrs = radec_data.ra.degree
        dec_icrs = radec_data.dec.degree

        num_sources = len(seg_table)

        # If the output is for the total detection product, then only
        # a subset of the full catalog is needed.
        if product.lower() == "tdp":

            # [x|y]centroid are in pixels, physical data coordinates
            seg_subset_table = seg_table["xcentroid", "ycentroid"]

            # Add metadata to the output subset table
            seg_subset_table = self._annotate_table(seg_subset_table, keyword_dict, num_sources, product=product)

            seg_subset_table["xcentroid"].description = "SExtractor Column x_image"
            seg_subset_table["ycentroid"].description = "SExtractor Column y_image"
            seg_subset_table["RA_icrs"] = ra_icrs
            seg_subset_table["Dec_icrs"] = dec_icrs
            seg_subset_table["RA_icrs"].description = "SExtractor Column RA"
            seg_subset_table["Dec_icrs"].description = "SExtractor Column Dec"
            seg_subset_table["RA_icrs"].unit = u.deg
            seg_subset_table["Dec_icrs"].unit = u.deg

            # Write out the official total detection product source catalog
            seg_subset_table["xcentroid"].info.format = ".10f"
            seg_subset_table["ycentroid"].info.format = ".10f"
            seg_subset_table["RA_icrs"].info.format = ".10f"
            seg_subset_table["Dec_icrs"].info.format = ".10f"
            print("seg_subset_table (white light image): ", seg_subset_table)

            seg_subset_table.write(catalog_filename, format="ascii.ecsv")
            log.info("Wrote source catalog: {}".format(catalog_filename))

        # else the product is the "filter detection product"
        else:

            seg_table = self._annotate_table(seg_table, keyword_dict, num_sources, product=product)

            # Rework the current table for output
            del seg_table["id"]
            del seg_table["sky_centroid"]
            del seg_table["sky_centroid_icrs"]
            rr = Column(ra_icrs, name="RA_icrs", description="SExtractor Column RA", unit=u.deg)
            dd = Column(dec_icrs, name="Dec_icrs", description="SExtractor Column Dec", unit=u.deg)
            seg_table.add_column(dd, index=2)
            seg_table.add_column(rr, index=2)

            # Add a description for columns which map to SExtractor catalog columns
            seg_table["xcentroid"].description = "SExtractor Column x_image"
            seg_table["ycentroid"].description = "SExtractor Column y_image"
            seg_table["background_at_centroid"].description = "SExtractor Column background"
            seg_table["source_sum"].description = "SExtractor Column flux_iso"
            seg_table["source_sum_err"].description = "SExtractor Column fluxerr_iso"
            # FIX: is mapping to _image or _world?  _image
            seg_table["cxx"].description = "SExtractor Column cxx_image, ellipse parameter"
            seg_table["cyy"].description = "SExtractor Column cyy_image, ellipse parameter"
            seg_table["cxy"].description = "SExtractor Column cxy_image, ellipse parameter"
            # FIX: is the mapping to _image or _world?
            seg_table["covar_sigx2"].description = "SExtractor Column x2_image, (0,0) element of covariance matrix"
            seg_table["covar_sigy2"].description = "SExtractor Column y2_image, (1,1) element of covariance matrix"
            seg_table[
                "covar_sigxy"].description = "SExtractor Column xy_image, (0,1) and (1,0) elements of covariance matrix"

            seg_table["xmin"].description = "SExtractor Column xmin_image"
            seg_table["xmax"].description = "SExtractor Column xmax_image"
            seg_table["ymin"].description = "SExtractor Column ymin_image"
            seg_table["ymin"].description = "SExtractor Column ymax_image"

            # Write out the official filter detection product source catalog
            seg_table["xcentroid"].info.format = ".10f"
            seg_table["ycentroid"].info.format = ".10f"
            seg_table["RA_icrs"].info.format = ".10f"
            seg_table["Dec_icrs"].info.format = ".10f"
            print("seg_table (filter): {}".format(seg_table))

            seg_table.write(catalog_filename, format="ascii.ecsv")
            log.info("Wrote filter source catalog: {}".format(catalog_filename))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def _get_header_data(self,product="tdp"):
        """Read FITS keywords from the primary or extension header and store the
        information in a dictionary

        Parameters
        ----------
        None.
        """

        keyword_dict = {}

        keyword_dict["proposal_id"] = self.imghdu[0].header["PROPOSID"]
        keyword_dict["image_file_name"] = self.imghdu[0].header['FILENAME'].upper()
        keyword_dict["target_name"] = self.imghdu[0].header["TARGNAME"].upper()
        keyword_dict["date_obs"] = self.imghdu[0].header["DATE-OBS"]
        keyword_dict["instrument"] = self.imghdu[0].header["INSTRUME"].upper()
        keyword_dict["detector"] = self.imghdu[0].header["DETECTOR"].upper()
        keyword_dict["target_ra"] = self.imghdu[0].header["RA_TARG"]
        keyword_dict["target_dec"] = self.imghdu[0].header["DEC_TARG"]
        keyword_dict["expo_start"] = self.imghdu[0].header["EXPSTART"]
        keyword_dict["texpo_time"] = self.imghdu[0].header["TEXPTIME"]
        keyword_dict["ccd_gain"] = self.imghdu[0].header["CCDGAIN"]
        keyword_dict["aperture_pa"] = self.imghdu[0].header["PA_V3"]

        # The total detection product has the FILTER keyword in
        # the primary header - read it for any instrument.
        #
        # For the filter detection product:
        # WFC3 only has FILTER, but ACS has FILTER1 and FILTER2
        # in the primary header.
        if product.lower() == "tdp":
            keyword_dict["filter"] = self.imghdu[0].header["FILTER"]
        # The filter detection product...
        else:
            if keyword_dict["instrument"] == "ACS":
                keyword_dict["filter1"] = self.imghdu[0].header["FILTER1"]
                keyword_dict["filter2"] = self.imghdu[0].header["FILTER2"]
            else:
                keyword_dict["filter1"] = self.imghdu[0].header["FILTER"]
                keyword_dict["filter2"] = ""

        # Get the HSTWCS object from the first extension
        keyword_dict["wcs_name"] = self.imghdu[1].header["WCSNAME"]
        keyword_dict["wcs_type"] = self.imghdu[1].header["WCSTYPE"]
        print('WCSTYPE: {}'.format(keyword_dict["wcs_type"]))
        keyword_dict["orientation"] = self.imghdu[1].header["ORIENTAT"]
        keyword_dict["aperture_ra"] = self.imghdu[1].header["RA_APER"]
        keyword_dict["aperture_dec"] = self.imghdu[1].header["DEC_APER"]

        return (keyword_dict)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def _annotate_table(self,data_table, keyword_dict, num_sources, product="tdp"):
        """Add state metadata to the output source catalog

        Parameters
        ----------
        data_table : QTable
            Table of source properties

        keyword_dict : dict
            Dictionary containing FITS keyword values pertaining to the data product

        num_sources : int
            Number of sources (items) in table

        product : str, optional
            Identification string for the catalog product being written.  This
            controls the data being put into the catalog product

        Returns
        -------
        data_table : QTable
            Table of source properties updatd to contain state metadata

        """

        data_table.meta["WCSNAME"] = keyword_dict["wcs_name"]
        data_table.meta["WCSTYPE"] = keyword_dict["wcs_type"]
        data_table.meta["Proposal ID"] = keyword_dict["proposal_id"]
        data_table.meta["Image File Name"] = keyword_dict['image_file_name']
        data_table.meta["Target Name"] = keyword_dict["target_name"]
        data_table.meta["Date Observed"] = keyword_dict["date_obs"]
        # FIX
        if product.lower() == "tdp":
            data_table.meta["Time Observed"] = " "
            data_table.meta["Filter"] = keyword_dict["filter"]
        else:
            data_table.meta["Time Observed"] = "FIX ME"
            data_table.meta["Filter 1"] = keyword_dict["filter1"]
            data_table.meta["Filter 2"] = keyword_dict["filter2"]
        data_table.meta["Instrument"] = keyword_dict["instrument"]
        data_table.meta["Detector"] = keyword_dict["detector"]
        data_table.meta["Target RA"] = keyword_dict["target_ra"]
        data_table.meta["Target DEC"] = keyword_dict["target_dec"]
        data_table.meta["Orientation"] = keyword_dict["orientation"]
        data_table.meta["Aperture RA"] = keyword_dict["aperture_ra"]
        data_table.meta["Aperture DEC"] = keyword_dict["aperture_dec"]
        data_table.meta["Aperture PA"] = keyword_dict["aperture_pa"]
        data_table.meta["Exposure Start"] = keyword_dict["expo_start"]
        data_table.meta["Total Exposure Time"] = keyword_dict["texpo_time"]
        data_table.meta["CCD Gain"] = keyword_dict["ccd_gain"]
        data_table.meta["Number of sources"] = num_sources
        data_table.meta[""] = " "
        data_table.meta[""] = "Absolute coordinates are in a zero-based coordinate system."

        return (data_table)


# ======================================================================================================================


if __name__ == '__main__':
    """Super simple testing interface for the above code."""

    starting_dt = datetime.datetime.now()


    parser = argparse.ArgumentParser(description='test interface for sourcelist_generation')
    parser.add_argument('total_product_name',help="total product filename")
    parser.add_argument('-f', '--filter_product_list',nargs='+',required=True,
                        help="Space-separated list of one or more total filter products")
    parser.add_argument('-d', '--debug',required=False,choices=['True','False'],default='False',help='debug mode on? (generate region files?)')
    parser.add_argument('-m', '--phot_mode',required=False,choices=['point','seg','both'],default='both',help="which photometry mode should be run? 'point' for point-soruce only; 'seg' for segment only, and 'both' for both point-source and segment photometry. ")
    args = parser.parse_args()
    if args.debug == "True":
        args.debug = True
    else:
        args.debug = False

    run_catalog_utils(args, starting_dt)

@util.with_logging
def run_catalog_utils(args,starting_dt):
    log.info("Run start time: {}".format(str(starting_dt)))
    log.info("python {} {} -f {} -d {} -m {}".format(os.path.realpath(__file__),
                                               args.total_product_name,
                                               " ".join(args.filter_product_list),
                                               args.debug,args.phot_mode))

    total_product = build_catalogs(args.total_product_name)

    total_product.check_param_dict()
    if args.phot_mode in ['seg', 'both']:
        total_product.segmap, \
        total_product.kernel, \
        total_product.bkg_dao_rms = \
            total_product.create_sextractor_like_sourcelists(total_product.seg_sourcelist_filename,se_debug=args.debug)


    if args.phot_mode in ['point', 'both']:
        total_product.ps_source_cat = total_product.identify_point_sources()
        total_product.write_catalog_to_file(total_product.ps_source_cat, write_region_file=args.debug)

        print("\a\a")
        sys.exit()

    for filter_img_name in args.filter_product_list:
        filter_product = build_catalogs(filter_img_name)
        if args.phot_mode in ['point', 'both']:
            filter_product.ps_phot_cat = filter_product.perform_point_photometry(total_product.ps_source_cat)
            filter_product.write_catalog_to_file(filter_product.ps_phot_cat,write_region_file=args.debug)

        if args.phot_mode in ['seg', 'both']:
            filter_product.measure_source_properties(total_product.segmap,
                                                     total_product.kernel,
                                                     filter_product.seg_sourcelist_filename,
                                                     filter_product.param_dict)

    log.info('Total processing time: {} sec'.format((datetime.datetime.now() - starting_dt).total_seconds()))