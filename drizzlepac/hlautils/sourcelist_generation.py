#!/usr/bin/env python

"""This script contains code to support creation of source extractor-like and daophot-like sourcelists.

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

__taskname__ = 'sourcelist_generation'

log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout)

# ----------------------------------------------------------------------------------------------------------------------


def create_dao_like_coordlists(fitsfile,sourcelist_filename,param_dict,make_region_file=False,bkgsig_sf=4.,
                               dao_ratio=0.8):
    """Make daofind-like coordinate lists

    Parameters
    ----------
    fitsfile : string
        Name of the drizzle-combined filter product to used to generate photometric sourcelists.

    sourcelist_filename : string
        Name of optionally generated ds9-compatible region file

    param_dict : dictionary
        Dictionary of instrument/detector - specific drizzle, source finding and photometric parameters

    dao_fwhm : float
        (photutils.DAOstarfinder param 'fwhm') The full-width half-maximum (FWHM) of the major axis of the
        Gaussian kernel in units of pixels. Default value = 3.5.

    make_region_file : Boolean
        Generate ds9-compatible region file? Default value = True

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
    bad_srcs = np.where(astrometric_utils.classify_sources(cat) == 0)[0] + 1
    segm.remove_labels(bad_srcs)

    source_table = cat.to_table()
    smajor_sigma = source_table['semimajor_axis_sigma'].mean().value
    source_fwhm = smajor_sigma * gaussian_sigma_to_fwhm

    log.info("DAOStarFinder(fwhm={}, threshold={}, ratio={})".format(source_fwhm,bkg_rms_mean,bkg_rms_mean))
    daofind = DAOStarFinder(fwhm=source_fwhm, threshold=bkg_rms_mean, ratio=dao_ratio)
    sources = daofind(image)
    hdulist.close()

    for col in sources.colnames:
        sources[col].info.format = '%.8g'  # for consistent table output

    # Write out ecsv file
    tbl_length = len(sources)
    sources.write(sourcelist_filename, format="ascii.ecsv")
    log.info("Created coord list  file '{}' with {} sources".format(sourcelist_filename, tbl_length))

    if make_region_file:
        out_table = sources.copy()
        # Remove all other columns besides xcentroid and ycentroid
        out_table.keep_columns(['xcentroid','ycentroid'])

        # Add offset of 1.0 in X and Y to line up sources in region file with image displayed in ds9.
        out_table['xcentroid'].data[:] += np.float64(1.0)
        out_table['ycentroid'].data[:] += np.float64(1.0)

        reg_filename = sourcelist_filename.replace(".ecsv",".reg")
        out_table.write(reg_filename, format="ascii")
        log.info("Created region file '{}' with {} sources".format(reg_filename, len(out_table)))

    return(sources)


# ----------------------------------------------------------------------------------------------------------------------


def create_dao_like_sourcelists(fitsfile,sl_filename,sources,aper_radius=4.,make_region_file=False):
    """Make DAOphot-like sourcelists

    Parameters
    ----------
    fitsfile : string
        Name of the drizzle-combined filter product to used to generate photometric sourcelists.


    sl_filename : string
        Name of the sourcelist file that will be generated by this subroutine

    sources : astropy table
        Table containing x, y coordinates of identified sources

    aper_radius : float
        Aperture radius (in pixels) used for photometry. Default value = 4.

    make_region_file : Boolean
        Generate ds9-compatible region file(s) along with the sourcelist? Default value = False

    Returns
    -------
    Nothing.
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

    # Write out sourcelist
    tbl_length = len(phot_table)
    phot_table.write(sl_filename, format="ascii.ecsv")
    log.info("Created sourcelist file '{}' with {} sources".format(sl_filename, tbl_length))

    # Write out ds9-compatable .reg file
    if make_region_file:
        reg_filename = sl_filename.replace(".ecsv",".reg")
        out_table = phot_table.copy()
        out_table['xcenter'].data = out_table['xcenter'].data + np.float64(1.0)
        out_table['ycenter'].data = out_table['ycenter'].data + np.float64(1.0)
        out_table.remove_column('id')
        out_table.write(reg_filename, format="ascii")
        log.info("Created region file '{}' with {} sources".format(reg_filename, tbl_length))


# ----------------------------------------------------------------------------------------------------------------------


def create_sourcelists(obs_info_dict, param_dict):
    """Main calling code. Make sourcelists

    Parameters
    ----------
    obs_info_dict : dictionary
        Dictionary containing all information about the images being processed

    param_dict : dictionary
        Dictionary of instrument/detector - specific drizzle, source finding and photometric parameters

    Returns
    -------
    """

    log.info("-" * 118)
    for product_type in obs_info_dict:
        for item_type in obs_info_dict[product_type]:
            log.info("obs_info_dict[{}][{}]: {}".format(product_type,item_type,obs_info_dict[product_type][item_type]))  # TODO: REMOVE THIS SECTION BEFORE ACTUAL USE
    log.info("-"*118)

    for tdp_keyname in [oid_key for oid_key in list(obs_info_dict.keys()) if
                        oid_key.startswith('total detection product')]:  # loop over total filtered products
        log.info("=====> {} <======".format(tdp_keyname))
        parse_tdp_info = obs_info_dict[tdp_keyname]['info'].split()
        inst_det = "{} {}".format(parse_tdp_info[2].upper(),parse_tdp_info[3].upper())

        detection_image = obs_info_dict[tdp_keyname]['product filenames']['image']
        tdp_seg_catalog_filename = obs_info_dict[tdp_keyname]['product filenames']['segment source catalog']
        tdp_ps_catalog_filename = obs_info_dict[tdp_keyname]['product filenames']['point source catalog']

        # segmap, kernel, bkg_dao_rms = se_source_generation.create_sextractor_like_sourcelists(
        #     detection_image, tdp_seg_catalog_filename, param_dict[inst_det], se_debug=False) # TODO: UNCOMMENT PRIOR TO DEPLOYMENT

        dao_coord_list = create_dao_like_coordlists(detection_image,tdp_ps_catalog_filename,param_dict[inst_det])


        for fp_keyname in obs_info_dict[tdp_keyname]['associated filter products']:
            filter_combined_imagename = obs_info_dict[fp_keyname]['product filenames']['image']

            point_source_catalog_name = obs_info_dict[fp_keyname]['product filenames']['point source catalog']
            seg_source_catalog_name = obs_info_dict[fp_keyname]['product filenames']['segment source catalog']

            log.info("Filter combined image... {}".format(filter_combined_imagename))
            log.info("Point source catalog.... {}".format(point_source_catalog_name))
            log.info("Segment source catalog.. {}".format(seg_source_catalog_name))

            # se_source_generation.measure_source_properties(segmap, kernel, filter_combined_imagename,
            #                                                seg_source_catalog_name, param_dict[inst_det]) # TODO: UNCOMMENT PRIOR TO DEPLOYMENT

            create_dao_like_sourcelists(filter_combined_imagename, point_source_catalog_name, dao_coord_list)
            

# ----------------------------------------------------------------------------------------------------------------------


@util.with_logging
def run_create_sourcelists(obs_info_dict, param_dict):
    """ subroutine to run create_sourcelists and produce log file when not run sourcelist_generation is not run from
    hlaprocessing.py.

    Parameters
    ----------
    obs_info_dict : dictionary
        Dictionary containing all information about the images being processed

    param_dict : dictionary
        Dictionary of instrument/detector - specific drizzle, source finding and photometric parameters

    Returns
    -------
    """

    create_sourcelists(obs_info_dict, param_dict)


