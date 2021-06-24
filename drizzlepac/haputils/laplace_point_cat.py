"""Module to test using enhanced edge detection for point-source identification

    Steps:
        1. Drizzle all inputs to final output frame as 'drz_ref'
        2. Drizzle all inputs to final output frame using 'kernel=turbo'
        3. Add offset (0.25 or 0.5 or user-specified value) to CRPIX1 and/or CRPIX2 of all input exposures as 'turbo_drz'
        4. Drizzle all updated inputs to final output frame (using "reference_wcs=drz_ref[1]") using 'kernel=turbo' as 'offset_drz'
        5. Create 'delta_drz' as 'turbo_drz - offset_drz'
        6. Apply 'np.nan_to_num' to 'delta_drz'  [not necessary if 'final_fillval=0.0']
        7. Smooth 'delta_drz' with gaussian_filter(sigma=1.5)
        8. Compute background stats 'bkg' using 'astropy.stats.sigma_clipped_stats(niters=1)'
        9. Perform segmentation using 'photutils.segmentation.detect_sources(threshold=bkg[2]*nsigma, npixels=npix), where nsigma=3 or user-supplied and npix is user-supplied(default=8? 24?)
       10. Create source mask 'rawmask' using 'np.clip(0,1)'
       11. Perform binary dilation (niterations=1) to close holes/edges created by source detection of pos and neg regions for each source
       12. Perform segmentation on mask
       13. For each slice from mask segmentation image, find pixel with peak value in 'drz_ref'
       14. Feed list of initial source positions to 'UserStarFinder' for final catalog measurements.

"""
import os
import sys
from distutils.version import LooseVersion
from datetime import datetime

import numpy as np  # nan_to_num, clip
from scipy import ndimage  # gaussian_filter, binary_dilation
from scipy.spatial import distance

from astropy.stats import sigma_clipped_stats
from astropy.io import fits
from astropy.table import Table, Column, vstack

import photutils
from photutils.detection import DAOStarFinder, IRAFStarFinder
from photutils.segmentation import (detect_sources, source_properties,
                                    deblend_sources)

from stwcs.wcsutil import HSTWCS

from stsci.tools import logutil

from .deconvolve_utils import UserStarFinder
from .cell_utils import SkyFootprint
from . import astrometric_utils as amutils

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)

def find_peaks(segmap, imgarr):
    """Return list of (X,Y) values for all peaks of segmap"""
    peaks = np.zeros((len(segmap.slices), 2))
    saturated = np.zeros(len(segmap.slices), dtype=np.bool)

    for i,segslice in enumerate(segmap):
        s = segslice.slices
        maxv = np.where(imgarr[s] == imgarr[s].max())
        zp = [s[0].start, s[1].start]
        cen = [s[0].start + (s[0].stop - s[0].start) // 2,
               s[1].start + (s[1].stop - s[1].start) // 2]
        # save peak in (X,Y) order
        peaks[i] = [maxv[1][0] + zp[1], maxv[0][0] + zp[0]]
        peak_dist = np.sqrt((peaks[i][0] - cen[1]) ** 2 + (peaks[i][1] - cen[0]) ** 2)
        saturated[i] = True if peak_dist > 3 else False

    peak_tab = Table(data=[Column(data=peaks[:,0].astype(np.int32), name='x_peak'),
                           Column(data=peaks[:,1].astype(np.int32), name='y_peak'),
                           Column(data=saturated, name='sat_flag'),
                           segmap.labels
                           ])
    return peak_tab


def build_seg_mask(refimg, ref_arr, fwhm, footprint_mask,
                   nsigma=3, npixels=8,
                   invert=False,
                   classify=True,
                   diagnostic_mode=False):

    sf = int(fwhm * 100)

    inv_factor = -1 if invert else 1

    files = []
    # Create delta image from base drizzle images
    dt = datetime.ctime(datetime.now())
    log.info('Computing Laplacian source edges now')
    delta_drz = inv_factor * ndimage.gaussian_laplace(ref_arr, fwhm)

    if diagnostic_mode:
        invstr = '_inv' if invert else ''
        offset_hdu = fits.PrimaryHDU(data=delta_drz, header=fits.getheader(refimg, ext=1))
        tmpname = refimg.replace('_dr', '{}_laplace{:03d}_dr'.format(invstr, sf))
        offset_hdu.writeto(tmpname, overwrite=True)
        files.append(tmpname)
        del offset_hdu

    # Compute background stats
    log.info("Computing background stats")
    bkg_stats = sigma_clipped_stats(delta_drz, maxiters=5)

    if diagnostic_mode:
       log.debug('Background computed as: {}'.format(bkg_stats))

    # Perform segmentation on smoothed delta image
    threshold = max(0., bkg_stats[0]) + bkg_stats[2] * nsigma

    dt = datetime.ctime(datetime.now())
    if diagnostic_mode:
        log.debug("Using detection threshold of: {}, {}".format(threshold, nsigma))
    log.info("Identifying sources from Laplacian edge map.")
    offset_segm = detect_sources(delta_drz,
                                 threshold=threshold,
                                 npixels=npixels,
                                 mask=footprint_mask)

    if classify:
        dt = datetime.ctime(datetime.now())
        log.info("Identifying potential cosmic-ray sources.")
        offset_segm = amutils.find_crs(ref_arr, offset_segm, fwhm)

    if diagnostic_mode:
        offset_hdu = fits.PrimaryHDU(data=offset_segm, header=fits.getheader(refimg, ext=1))
        tmpname = refimg.replace('_dr', '{}_laplace{:03d}_dr'.format(invstr, sf))
        offset_hdu.writeto(tmpname, overwrite=True)
        files.append(tmpname)
        del offset_hdu

    # create raw mask
    offset_mask = np.clip(offset_segm.data, 0, 1).astype(np.bool)

    if not invert:
        # dilate mask to close up plus/neg delta regions for same source
        log.debug("Filling holes in laplacian mask...")
        offset_mask = ndimage.binary_erosion(ndimage.binary_dilation(offset_mask, iterations=2), iterations=2)
        offset_mask = ndimage.binary_fill_holes(offset_mask)

    return offset_mask, files



def find_saturated_sources(peaks, point_cat, tolerance=1.0):

    # extract arrays of all detected peaks and of measured positions
    peak_arr = np.array([peaks['x_peak'], peaks['y_peak']]).T
    cat_arr = np.array([point_cat['xcentroid'], point_cat['ycentroid']]).T
    # look for any peak that does not correspond (within a pixel or two?) to
    # a position in the measured catalog
    peak_seps = distance.cdist(peak_arr, cat_arr)
    sat_indx = np.where(np.min(peak_seps, axis=1) > tolerance)[0]

    return sat_indx


def measure_saturated_source(img, segment, seg_table, diagnostic_mode=False):

    seg_slice = segment.slices

    # Define raw data from this slice
    detection_img = img[seg_slice]
    # zero out any pixels which do not have this segments label
    detection_img[segment.data == 0] = 0

    # This logic will flag/measure sources which have more than 3 pixels
    # within 10% of the max value in the source segment, a situation
    # which would indicate the presence of a saturated source
    if (detection_img > detection_img.max() * 0.9).sum() > 3:
        # Revert to segmentation photometry for sat. source posns
        segment_properties = source_properties(detection_img, segment.data)
        sat_table = segment_properties.to_table()
        # Add new row to table, then populate it with the newly measured values
        seg_table.add_row([v for v in seg_table[-1].values()])
        xpeak = sat_table['xcentroid'][0].value + seg_slice[1].start
        ypeak = sat_table['ycentroid'][0].value + seg_slice[0].start
        seg_table['flux'][-1] = sat_table['source_sum'][0]
        seg_table['peak'][-1] = sat_table['max_value'][0]
        seg_table['xcentroid'][-1] = xpeak
        seg_table['ycentroid'][-1] = ypeak
        seg_table['npix'][-1] = sat_table['area'][0].value
        sky = sat_table['background_mean'][0]
        seg_table['sky'][-1] = sky.value if sky is not None and not np.isnan(sky) else 0.0
        seg_table['mag'][-1] = -2.5 * np.log10(sat_table['source_sum'][0])
        seg_table['sat_flag'] = True

        log.debug("Use segmentation to measure a saturated source at {},{}".format(xpeak, ypeak))

    return seg_table


#
#
#  Primary User Interface for this module
#
#
def generate_catalog(refimg, fwhm,
                     threshold=None,
                     mask=None, imglist=None,
                     log_level=logutil.logging.NOTSET,
                     diagnostic_mode=False,
                     **pars):
    """

    Parameters
    ----------
    refimg : str
        Filename of final drizzle product to generate catalog for
    fwhm : float
        FWHM of PSF for the image
    mask : ndarray or None, optional
        Array of footprint where 1 corresponds to areas which were exposed.
        If `None`, either use images from `imglist` to
        create one **if `imglist` is specified**, or simply use the entire image.
    imglist : list, optional
        If `mask` is None, use this list of input images to be drizzled
        together to create the mask.
    log_level : int, optional
        Specify the level of log messages to generate during processing
    diagnostic_mode : bool, optional
        Specify whether or not to generate additional intermediate outputs
        to try to understand the computations better.
    pars : dict
        Additional parameters for processing.  Optional parameters that can
        specified include:
        * edge_distance : int, optional [Default: 10]
        * classify : bool, optional [Default: True]
        * nsigma : float, optional [Default : 3.0]
        * npixels : int, optional [Default: 8]
        * nlevels : int, optional [Default: 32]
        * contrast : float, optional  [Default: 0.005]
        * deblend_mode : str, optional [Default: 'exponential']
        * sharplo : float, optional [Default: None]
        * sharphi : float, optional [Default: None]

    Returns
    -------
    point_cat : `~astropy.table.Table`
        Table of identified and measured point sources derived from `refimg`
    files : list
        List of filenames of products written to disk during processing

    """
    log.setLevel(log_level)

    # detection parameters
    nsigma = pars.pop('nsigma') if 'nsigma' in pars else 3.0
    npixels = pars.pop('npixels') if 'npixels' in pars else 8
    edge_distance = pars.pop('edge_distance') if 'edge_distance' in pars else 10

    classify = pars.pop('classify') if 'classify' in pars else True

    # deblending parameters
    nlevels = pars.pop('nlevels') if 'nlevels' in pars else 32
    contrast = pars.pop('contrast') if 'contrast' in pars else 0.001
    deblend_mode = pars.pop('deblend_mode') if 'deblend_mode' in pars else 'exponential'

    mask_segm, mfiles = detect_LoG_segments(refimg, sciext=1, imglist=imglist,
                                            edge_distance=edge_distance,
                                            nsigma=nsigma,
                                            npixels=npixels,
                                            classify=classify,
                                            diagnostic_mode=diagnostic_mode)
    files += mfiles
    offset_mask = ~np.clip(mask_segm.data, 0, 1).astype(np.bool)

    if threshold is None:
        masked_stats = sigma_clipped_stats(ref_arr, mask=offset_mask, maxiters=3)
        low_cut = masked_stats[0] + nsigma * masked_stats[2]
    else:
        low_cut = threshold
    daofind = DAOStarFinder(fwhm=fwhm, threshold=low_cut)
    dao_sources = daofind(ref_arr, mask=offset_mask)

    if diagnostic_mode:
        offset_hdu = fits.PrimaryHDU(data=offset_mask.astype(np.int16), header=fits.getheader(refimg, ext=1))
        tmpname = refimg.replace('_dr', '_mask_filled_dr')
        offset_hdu.writeto(tmpname, overwrite=True)
        files.append(tmpname)
        del offset_hdu


    # deblend segmentation map
    log.info('Deblending segments from laplacian mask')
    mask_deblend = deblend_sources(ref_arr, mask_segm,
                                   npixels=npixels,
                                   nlevels=nlevels,
                                   mode=deblend_mode,
                                   contrast=contrast)

    # find peaks from original DRZ image 'refimg' based on mask segmentation map
    dt = datetime.ctime(datetime.now())
    log.info('Looking for the peak value in each detected segment')
    drz_peaks = find_peaks(mask_deblend, ref_arr)

    if diagnostic_mode:
        offset_hdu = fits.PrimaryHDU(data=mask_deblend, header=fits.getheader(refimg, ext=1))
        tmpname = refimg.replace('_dr', '_mask_segmap_dr')
        offset_hdu.writeto(tmpname, overwrite=True)
        files.append(tmpname)
        del offset_hdu
        indx = refimg.find('_dr')
        drz_peaks.write('{}_raw_peaks.ecsv'.format(refimg[:indx]), format='ascii.ecsv')

    # Generate final measured point source catalog
    log.info("Generating catalog based on initial peak positions")
    usf = UserStarFinder(threshold=0, fwhm=fwhm, coords=drz_peaks, **pars)
    point_cat = usf.find_stars(ref_arr)

    # Add measurements for saturated sources
    # sat_indices = find_saturated_sources(drz_peaks, point_cat)
    sat_indices = np.where(drz_peaks['sat_flag'])[0]
    if diagnostic_mode:
        log.debug('Processing {} rejected stars'.format(len(sat_indices)))

    # add new column to record whether any sources were
    # flagged and measured as saturated sources
    sat_col = Column(data=np.zeros(len(point_cat), dtype=np.bool), name='sat_flag')
    point_cat.add_column(sat_col)

    for segment_indx in sat_indices:
        segment = mask_deblend.segments[segment_indx]
        log.debug('Measuring sat source with label={}'.format(segment.label))
        point_cat = measure_saturated_source(ref_arr, segment, point_cat)

    # concatentate dao catalog with laplace catalog
    point_cat = vstack([point_cat, dao_sources])

    if 'sharplo' in pars:
        indx = np.where(point_cat['sharpness'] <= pars['sharplo'])[0][::-1]
        for i in indx: del point_cat[i]
    if 'sharphi' in pars:
        indx = np.where(point_cat['sharpness'] >= pars['sharphi'])[0][::-1]
        for i in indx: del point_cat[i]
    if 'roundlo' in pars:
        indx = np.where(point_cat['roundness1'] <= pars['roundlo'])[0][::-1]
        for i in indx: del point_cat[i]
        indx = np.where(point_cat['roundness2'] <= pars['roundlo'])[0][::-1]
        for i in indx: del point_cat[i]
    if 'roundhi' in pars:
        indx = np.where(point_cat['roundness1'] >= pars['roundhi'])[0][::-1]
        for i in indx: del point_cat[i]
        indx = np.where(point_cat['roundness2'] >= pars['roundhi'])[0][::-1]
        for i in indx: del point_cat[i]

    return point_cat, files


def detect_LoG_segments(img, sciext=1, imglist=None, edge_distance=10,
                        nsigma=3.0, npixels=8, classify=False,
                        diagnostic_mode=False):
    """Core function to create segmentation map of LoG sources from img

    Returns
    -------
    mask_segm : `photutils.segmentation.SegmentationImage`
        Segmentation image for all identified sources based on Laplace-of-Gaussian
        filtering of the input image.
    files : list
        List of any intermediate files that were generated during processing.
    """
    if not isinstance(img, np.ndarray):
        img_arr = np.nan_to_num(fits.getdata(img, ext=('sci', sciext)), 0)
    else:
        img_arr = img.copy()

    # find peaks from input image 'img' based on mask segmentation map
    img_arr = np.clip(img_arr, 0, img_arr.max())

    files = []
    # determine footprint mask, if not provided by user
    if mask is None:
        if imglist:
            imgwcs = HSTWCS(img, ext=1)
            footprint = SkyFootprint(imgwcs)
            footprint.build(imglist)
            # shrink footprint by 'edge_distance' number of pixels
            footprint_mask = ~ndimage.binary_erosion(footprint.total_mask, iterations=edge_distance)
        else:
            footprint_mask = ndimage.binary_erosion((img_arr != 0), iterations=edge_distance)
            # footprint_mask = np.zeros_like(img_arr, dtype=np.bool)
    else:
        footprint_mask = mask

    if diagnostic_mode:
        offset_hdu = fits.PrimaryHDU(data=footprint_mask.astype(np.int16), header=fits.getheader(refimg, ext=1))
        tmpname = refimg.replace('_dr', '_footprint_mask_dr')
        offset_hdu.writeto(tmpname, overwrite=True)
        files.append(tmpname)
        del offset_hdu

    # create detection masks using gaussian_laplace filter
    # Identify smaller, fainter sources by applying smoothing with slightly smaller sigma
    # pull out brighter sources more readily by smoothing with larger sigma
    offset_mask, dfiles = build_seg_mask(img, img_arr, 1.0, footprint_mask,
                                          nsigma=nsigma*2, npixels=npixels,
                                          classify=classify,
                                          diagnostic_mode=diagnostic_mode)
    files += dfiles

    # get final segmentation map
    log.info('Detecting segments in laplacian mask')
    mask_segm = detect_sources(img_arr,
                               threshold=0,
                               npixels=npixels,
                               mask=~offset_mask)
    #                          mask=footprint_mask)

    return mask_segm, files

def clean_files(files):
    """Remove intermediate files created by generate_catalog"""

    for f in files:
        try:
            if f.strip() not in ['', None]:
                os.remove(f)
        except Exception:
            print("WARNING: Could not remove:\n    {}\n".format(f))
