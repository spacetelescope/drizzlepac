#!/usr/bin/env python

"""Code related to the creation and modifaction of diagnostic .json files"""

import collections
from datetime import datetime
import json
import os
import pdb
import random
import string
import sys

import numpy as np

from stsci.tools import logutil

__taskname__ = 'diagnostic_utils'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)
# ======================================================================================================================

class HapDiagnosticObj(object):
    def __init__(self,prop_id,obset_id,telescope,instrument,detector,filter,data_source,description,log_level=logutil.logging.NOTSET):
        """base class used to set up a HapDiagnostic object.

        Parameters
        ----------
        prop_id : string
            proposal ID

        obset_id : string
            obset ID

        telescope: string
            telescope name (e.g. hst, jwst, etc)

        instrument : string
            instrument name

        detector : string
            detector name

        filter : string
            filter name

        data_source : string
            name of the script that generated the data that will be stored in the "data" section

        description : string
            brief description of what the data is, and how it should be used.

        log_level : int, optional
            The desired level of verboseness in the log statements displayed on the screen and written to the .log file.
            Default value is 'NOTSET'.

        Returns
        -------
        Nothing.
        """
        self.prop_id = prop_id
        self.obset_id = obset_id
        self.telescope = telescope
        self.instrument = instrument
        self.detector =  detector
        self.filter = filter
        self.data_source = data_source
        self.description = description
        log.setLevel(log_level)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def instantiate(self):
        """Creates a new diagnostic dictionary using the standard format. The standard header values are as follows:

        - Generation date
        - Generation time (local 24-hour format)
        - Proposal ID
        - Obset ID
        - Telescope name
        - Instrument name
        - Detector name
        - Filter name
        - Data source (name of the piece of code that produced the data)
        - Description (brief description of what the data is, and how it should be used)

        Parameters
        ----------
        Nothing.

        Updates
        -------
        self.out_dict : Ordered dictionary
            dictionary that will ultimately be written to a json file
        """
        # summon nested orderedDict into existence
        self.out_dict = collections.OrderedDict()
        self.out_dict['header'] = collections.OrderedDict()
        self.out_dict['data'] = collections.OrderedDict()

        # Populate standard header fields
        timestamp = datetime.now().strftime("%m/%d/%YT%H:%M:%S")
        self.out_dict['header']['generation date'] = timestamp.split("T")[0]
        self.out_dict['header']['generation time'] = timestamp.split("T")[1]
        header_item_list = ["prop_id", "obset_id", "telescope", "instrument", "detector", "filter", "data_source", "description"]
        for header_item in header_item_list:
            self.out_dict['header'][header_item] = self.__dict__[header_item]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def addDataItem(self,dataset,title):
        """main subroutine for adding data to self.out_table.

        Parameters
        ----------
        dataset : varies
            data to add to self.out_dict

        title : str
            Name of the dictionary key that will be used to store dataset in self.out_dict

        Updates
        -------
        self.out_dict : Ordered dictionary
            dictionary that will ultimately be written to a json file
        """
        dataset_type = str(type(dataset))
        if dataset_type == "<class 'numpy.ndarray'>": #For numpy arrays
            self.out_dict['data'][title]=collections.OrderedDict()
            self.out_dict['data'][title]["original format"] = dataset_type
            self.out_dict['data'][title]["data"] = dataset.tolist()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def addheaderItem(self):
        pass
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def readJsonFile(self):
        pass

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def writeJsonFile(self,json_filename,clobber=False):
        """Writes self.out_dict to user-specified filename.

        Parameters
        ----------
        json_filename : string
            name of the json file to write

        clobber : bool, optional
            Overwrite file with same name? Default value = False

        Returns
        -------
        Nothing.
        """
        file_exists = os.path.exists(json_filename)
        print(file_exists,clobber)
        if clobber:
            if file_exists:
                os.remove(json_filename)
        else:
            if file_exists:
                random_string = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(4))
                json_filename = json_filename.replace(".json","_{}.json".format(random_string))
        with open(json_filename,"w") as json_file:
            json.dump(self.out_dict, json_file, indent=4)
        log.info("Wrote json file {}".format(json_filename))


# ======================================================================================================================
if __name__ == "__main__":
    """
    Preliminary data section item format:
    Each numpy table column would be converted into a nested dictionary with the fillowing items
    - title (same as the key to this dictionary(?)
    - units
    - dtype
    - data
    - masking information
    """

    blarg = HapDiagnosticObj(telescope="hst",
                             instrument = "wfc3",
                             detector = "ir",
                             filter = "f160w",
                             prop_id = "11979",
                             obset_id = "01",
                             data_source = "hla_flag_filter",
                             description = "test item please ignore",
                             log_level=10)
    blarg.instantiate()

    blarg.writeJsonFile("diag_test.json", clobber=True)

