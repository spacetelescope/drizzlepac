import os
import shutil
import datetime

import numpy as np
import astropy
from astropy.io import fits
from astropy.table import Table

from stwcs.wcsutil import HSTWCS
from stsci.tools import logutil

from . import astrometric_utils as amutils

from .. import wcs_functions

__taskname__ = 'align_utils'

log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout)


class AlignmentTable:
    def __init__(self, input_list):
        log.info("{} STEP 1: Get data {}".format("-" * 20, "-" * 66))
        zero_dt = starting_dt = datetime.datetime.now()
        log.info(str(starting_dt))
        self.imglist = check_and_get_data(input_list, archive=archive, clobber=clobber)
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
        pass
    
