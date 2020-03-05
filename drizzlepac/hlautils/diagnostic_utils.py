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

    def addUpdateHeaderItem(self,element_name, new_element_value, clobber=True, addnew=True):
        """add or update a single user-specified header item

        Parameters
        ----------
        element_name : str
            Name of the header element to add or update

        new_element_value : varies
            desired new value for the header element

        clobber : bool, optional
            overwrite existing value if specified header element exists? Default = True

        addnew : bool, optional
            if specified header element does not already exist, add it as a new header element? Default = True

        Updates
        -------
        self.out_dict : Ordered dictionary
            dictionary that will ultimately be written to a json file
        """

        if element_name in self.out_dict['header'].keys():
            if clobber:
                log.info("{}: {} -> {} value update successful.".format(element_name,
                                                                        self.out_dict['header'][element_name],
                                                                        new_element_value))
                self.out_dict['header'][element_name] = new_element_value

            else:
                log.warning("Header element '{}' already exists. Update NOT performed.".format(element_name))
        else:
            if addnew:
                self.out_dict['header'][element_name] = new_element_value
                log.info("New element {} = {} successfully added to header".format(element_name,new_element_value))
            else:
                log.warning("Unable to add new element {} = {} to header".format(element_name,new_element_value))

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
        header_item_list = ["prop_id", "obset_id", "telescope", "instrument", "detector", "filter", "data_source",
                            "description"]
        for header_item in header_item_list:
            self.out_dict['header'][header_item] = self.__dict__[header_item]


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



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def readJsonFile(json_filename):
    """extracts header and data sections from specified json file.

    Parameters
    ----------
    json_filename : str
        Name of the json file to extract data from

    Returns
    -------
    json_data : dictionary
        dictionary structured similarly to self.out_dict with seperate 'header' and 'data' sections. The
        information stored in the 'data' section will be in the same format that it was in before it was serialized
        and stored as a json file.
    """
    if os.path.exists(json_filename):
        with open(json_filename) as f:
            json_data = json.load(f)

    else:
        errmsg = "json file {} not found!".format(json_filename)
        log.error(errmsg)
        raise Exception(errmsg)
    return(json_data)
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
    
    blarg.addUpdateHeaderItem("filter3",None,clobber=False,addnew=True)
    blarg.writeJsonFile("diag_test.json", clobber=True)

