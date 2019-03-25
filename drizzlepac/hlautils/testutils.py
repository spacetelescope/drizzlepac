import sys
import os

from astropy.io import fits

from stsci.tools import logutil
from stwcs import updatewcs
from stwcs.wcsutil import headerlet

from .. import alignimages

__taskname__ = 'testutils'

log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout)

def compare_wcs_alignment(dataset):
    """Return results from aligning dataset using all available WCS solutions.

        ASSUMPTIONS:
            - All images in dataset have the same set of a priori solutions
    """
    # Setup
    #   Insure that database will be queried for new WCS solutions
    control = os.environ['ASTROMETRY_STEP_CONTROL']
    os.environ['ASTROMETRY_STEP_CONTROL'] = 'ON'

    # Step 1:
    #   Determine alignment for pipeline-defined WCS
    results = alignimages.perform_align([dataset], clobber=False, debug=True)
    imglist = results['imageName'].astype(str).tolist()

    # Step 2:
    #   Create results output organized by WCSNAME
    img0 = imglist[0]
    default_wcsname = fits.getval(img0, 'wcsname', ext=1)
    log.info("Default WCSNAME: {}".format(default_wcsname))
    alignment = {default_wcsname:extract_results(results)}

    # Step 3:
    #   Update inputs with latest distortion model and pull in solutions from dB
    imglist = updatewcs.updatewcs(imglist)
    img0 = imglist[0]
    # Step 4:
    #   Loop over each WCS solution and perform alignment to GAIA
    wcsnames = headerlet.get_headerlet_kw_names(img0, kw='WCSNAME')
    if len(wcsnames) == 0:
        msg = "No a priori solutions found for {}".format(img0)
        log.error(msg)
        raise ValueError(msg)

    for wcs in wcsnames:
        log.info("Starting with {}".format(wcs))
        if wcs in [default_wcsname, 'OPUS']:
            continue # skip default pipeline solutions, since we have already aligned it
        # apply WCS from headerlet
        for img in imglist:
            wnames = headerlet.get_headerlet_kw_names(img, kw='WCSNAME')
            hnames = headerlet.get_headerlet_kw_names(img)
            for w,h in zip(wnames, hnames):
                if w == wcs:
                    hdrlet = h
                    break
            log.info("Applying WCS {} to {}".format(hdrlet, img))
            headerlet.restore_from_headerlet(img, hdrname=hdrlet,
                                             archive=False, force=True)

        results = alignimages.perform_align([dataset], clobber=False, debug=True)
        alignment[wcs] = extract_results(results)

    # Restore user environment to original state
    os.environ['ASTROMETRY_STEP_CONTROL'] = control

    return alignment

def extract_results(results):
    """Return dict with select columns from alignment results Table."""
    results_dict = {'images':results['imageName'].astype(str).tolist(),
                                  'offset_x':results['offset_x'],
                                  'offset_y':results['offset_y'],
                                  'rotation':results['rotation'],
                                  'rms_x': results['rms_x'], # RMS in pixels
                                  'rms_y':results['rms_y'],
                                  'fit_rms':results['fit_rms'], # RMS in arcsec
                                  'total_rms':results['total_rms'],
                                  'status': results['status'],
                                  'fit_qual': results['fit_qual'],
                                  'matched_sources': results['matchSources']}
    return results_dict
