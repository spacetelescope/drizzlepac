import os
import sys
import shutil
import datetime
import traceback
import glob

import numpy as np
import astropy
from astropy.io import fits
from astropy.table import Table

from stwcs.wcsutil import HSTWCS
from stsci.tools import logutil

from . import astrometric_utils as amutils
from . import astroquery_utils as aqutils

from .. import wcs_functions

__taskname__ = 'align_utils'

log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout)


def check_and_get_data(input_list, **pars):
    """Verify that all specified files are present. If not, retrieve them from MAST.

    Parameters
    ----------
    input_list : list
        List of one or more calibrated fits images that will be used for catalog generation.

    Returns
    =======
    total_input_list: list
        list of full filenames

    """
    empty_list = []
    retrieve_list = []    # Actual files retrieved via astroquery and resident on disk
    candidate_list = []   # File names gathered from *_asn.fits file
    ipppssoot_list = []   # ipppssoot names used to avoid duplicate downloads
    total_input_list = []  # Output full filename list of data on disk
    member_suffix = '_flc.fits'

    # Loop over the input_list to determine if the item in the input_list is a full association file
    # (*_asn.fits), a full individual image file (aka singleton, *_flt.fits), or a root name specification
    # (association or singleton, ipppssoot).
    for input_item in input_list:
        log.info('Input item: {}'.format(input_item))
        indx = input_item.find('_')

        # Input with a suffix (_xxx.fits)
        if indx != -1:
            lc_input_item = input_item.lower()
            suffix = lc_input_item[indx + 1:indx + 4]
            log.info('file: {}'.format(lc_input_item))
            # For an association, need to open the table and read the image names as this could
            # be a custom association.  The assumption is this file is on local disk when specified
            # in this manner (vs just the ipppssoot of the association).
            # This "if" block just collects the wanted full file names.
            if suffix == 'asn':
                try:
                    asntab = Table.read(input_item, format='fits')
                except FileNotFoundError:
                    log.error('File {} not found.'.format(input_item))
                    return(empty_list)
                for row in asntab:
                    if row['MEMTYPE'].startswith('PROD'):
                        continue
                    memname = row['MEMNAME'].lower().strip()
                    # Need to check if the MEMNAME is a full filename or an ipppssoot
                    if memname.find('_') != -1:
                        candidate_list.append(memname)
                    else:
                        # Define suffix for all members based on what files are present
                        if not os.path.exists(memname + member_suffix):
                            member_suffix = '_flt.fits'

                        candidate_list.append(memname + member_suffix)
            elif suffix in ['flc', 'flt']:
                if lc_input_item not in candidate_list:
                    candidate_list.append(lc_input_item)
            else:
                log.error(
                    'Inappropriate file suffix: {}.  Looking for "asn.fits", '
                    '"flc.fits", or "flt.fits".'.format(
                        suffix))
                return (empty_list)

        # Input is an ipppssoot (association or singleton), nine characters by definition.
        # This "else" block actually downloads the data specified as ipppssoot.
        elif len(input_item) == 9:
            try:
                if input_item not in ipppssoot_list:
                    # An ipppssoot of an individual file which is part of an association cannot be
                    # retrieved from MAST
                    retrieve_list = aqutils.retrieve_observation(input_item, **pars)

                    # If the retrieved list is not empty, add filename(s) to the total_input_list.
                    # Also, update the ipppssoot_list so we do not try to download the data again.  Need
                    # to do this since retrieve_list can be empty because (1) data cannot be acquired (error)
                    # or (2) data is already on disk (ok).
                    if retrieve_list:
                        total_input_list += retrieve_list
                        ipppssoot_list.append(input_item)
                    else:
                        log.error('File {} cannot be retrieved from MAST.'.format(input_item))
                        return(empty_list)
            except Exception:
                exc_type, exc_value, exc_tb = sys.exc_info()
                traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)

    # Only the retrieve_list files via astroquery have been put into the total_input_list thus far.
    # Now check candidate_list to detect or acquire the requested files from MAST via astroquery.
    for file in candidate_list:
        # If the file is found on disk, add it to the total_input_list and continue
        if glob.glob(file):
            total_input_list.append(file)
            continue
        else:
            log.error('File {} cannot be found on the local disk.'.format(file))
            return(empty_list)

    log.info("TOTAL INPUT LIST: {}".format(total_input_list))
    return(total_input_list)



class AlignmentTable:
    def __init__(self, input_list, clobber=False):
        log.info("{} STEP 1: Get data {}".format("-" * 20, "-" * 66))
        self.zero_dt = starting_dt = datetime.datetime.now()
        log.info(str(starting_dt))
        imglist = check_and_get_data(input_list, archive=False, clobber=clobber)
        log.info("SUCCESS")

        current_dt = datetime.datetime.now()
        delta_dt = (current_dt - starting_dt).total_seconds()
        log.info('Processing time of [STEP 1]: {} sec'.format(delta_dt))
        starting_dt = current_dt
        # 2: Apply filter to input observations to insure that they meet minimum criteria for being able to be aligned
        log.info(
            "{} STEP 2: Filter data {}".format("-" * 20, "-" * 63))
        self.filtered_table = filter.analyze_data(imglist)

        if self.filtered_table['doProcess'].sum() == 0:
            log.warning("No viable images in filtered table - no processing done.\n")
            current_dt = datetime.datetime.now()
            delta_dt = (current_dt - starting_dt).total_seconds()
            log.info('Processing time of [STEP 2]: {} sec'.format(delta_dt))
            return

        # Get the list of all "good" files to use for the alignment
        process_list = self.filtered_table['imageName'][np.where(self.filtered_table['doProcess'])]
        self.process_list = list(process_list)  # Convert process_list from numpy list to regular python list
        log.info("SUCCESS")

        # Convert input images to tweakwcs-compatible FITSWCS objects and
        # attach source catalogs to them.
        self.imglist = []
        for group_id, image in enumerate(self.process_list):
            img = amutils.build_wcscat(image, group_id,
                                       self.extracted_sources[image]['catalog_table'])
            # add the name of the image to the imglist object
            for im in img:
            #    im.meta['name'] = image
                log.info('im.meta[name] = {}'.format(im.meta['name']))
            self.imglist.extend(img)

        self.group_id_dict = {}
        for image in self.imglist:
            self.group_id_dict["{}_{}".format(image.meta["filename"], image.meta["chip"])] = image.meta["group_id"]

    def set_catalog(self, catalog_names):
        """Define the reference catalog(s) to be used for alignment

        Parameters
        ----------
        catalog_names : list
            List of astrometric catalog names to use for alignment.
            Options would include (but not be limited to):
            "GAIADR1" and "GAIADR2".
            If a filename is provided, file would be used instead
            of deriving catalogs from standard astrometric catalogs
            like GAIADR2.
        """
        pass

    def set_method(self, method):
        """Define what alignment method(s) are to be used

        Parameters
        ----------
        method : list
            List of alignment method names to be used.
            Supported options: relative, 2dhist, threshold.
        """
        pass

    def find_sources(self):
        """Find observable sources in each input exposure."""
        self.extracted_sources = None

    def get_reference_catalog(self, catalog):
        """Return the desired reference catalog to be used for alignment"""
        pass

    def reset_group_id(self):
        for image in self.imglist:
            image.meta["group_id"] = self.group_id_dict["{}_{}".format(image.meta["filename"], image.meta["chip"])]

    def perform_fit(self, method):
        """Perform fit using specified method, then determine fit quality"""
        pass
