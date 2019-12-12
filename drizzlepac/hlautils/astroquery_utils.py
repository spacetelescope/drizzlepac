"""Wrappers for astroquery-related functionality"""
import shutil
import os

try:
    from astroquery.mast import Observations
except FileExistsError:
    Observations = None

import sys
from stsci.tools import logutil

__taskname__ = 'astroquery_utils'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout, 
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)

def retrieve_observation(obsid, suffix=['FLC'], archive=False, clobber=False):
    """Simple interface for retrieving an observation from the MAST archive

    If the input obsid is for an association, it will request all members with
    the specified suffixes.

    Parameters
    -----------
    obsid : string
        ID for observation to be retrieved from the MAST archive.  Only the
        IPPSSOOT (rootname) of exposure or ASN needs to be provided; eg.,
        ib6v06060.

    suffix : list, optional
        List containing suffixes of files which should be requested from MAST.
        Default value  "['FLC']".

    archive : Boolean, optional
        Retain copies of the downloaded files in the astroquery created
        sub-directories? Default is "False".

    clobber : Boolean, optional
        Download and Overwrite existing files? Default is "False".

    Returns
    -------
    local_files : list
        List of filenames
    """
    local_files = []

    if Observations is None:
        log.warning("The astroquery package was not found.  No files retrieved!")
        return local_files

    # Query MAST for the data with an observation type of either "science" or
    # "calibration"
    obs_table = Observations.query_criteria(obs_id=obsid, obstype='all')
    # Catch the case where no files are found for download
    if not obs_table:
        log.info("WARNING: Query for {} returned NO RESULTS!".format(obsid))
        return local_files

    dpobs = Observations.get_product_list(obs_table)
    data_products_by_id = Observations.filter_products(dpobs,
                                                       productSubGroupDescription=suffix,
                                                       extension='fits',
                                                       mrp_only=False)

    # After the filtering has been done, ensure there is still data in the
    # table for download. If the table is empty, look for FLT images in lieu
    # of FLC images. Only want one or the other (not both!), so just do the
    # filtering again.
    if not data_products_by_id:
        log.info("WARNING: No FLC files found for {} - will look for FLT "
                 "files instead.".format(obsid))
        suffix = ['FLT']
        data_products_by_id = Observations.filter_products(dpobs,
                                                           productSubGroupDescription=suffix,
                                                           extension='fits',
                                                           mrp_only=False)

        # If still no data, then return.  An exception will eventually be
        # thrown in the higher level code.
        if not data_products_by_id:
            log.info(
                "WARNING: No FLC or FLT files found for {}.".format(obsid))
            return local_files
    all_images = data_products_by_id['productFilename'].tolist()
    log.info(all_images)
    if not clobber:
        rows_to_remove = []
        for row_idx, row in enumerate(data_products_by_id):
            fname = row['productFilename']
            if os.path.isfile(fname):
                log.info(fname + " already exists. File download skipped.")
                rows_to_remove.append(row_idx)
        data_products_by_id.remove_rows(rows_to_remove)

    manifest = Observations.download_products(data_products_by_id,
                                              mrp_only=False)

    if not clobber:
        for rownum in rows_to_remove[::-1]:
            if manifest:
                manifest.insert_row(rownum,
                                    vals=[all_images[rownum], "LOCAL", "None", "None"])
            else:
                return all_images

    download_dir = None
    for file, file_status in zip(manifest['Local Path'], manifest['Status']):
        if file_status != "LOCAL":
            # Identify what sub-directory was created by astroquery for the
            # download
            if download_dir is None:
                download_dir = os.path.dirname(os.path.abspath(file))
            # Move or copy downloaded file to current directory
            local_file = os.path.abspath(os.path.basename(file))
            if archive:
                shutil.copy(file, local_file)
            else:
                shutil.move(file, local_file)
            # Record what files were downloaded and their current location
            local_files.append(os.path.basename(local_file))
        else:
            local_files.append(file)
    if not archive:
        # Remove astroquery created sub-directories
        shutil.rmtree('mastDownload')
    return local_files
