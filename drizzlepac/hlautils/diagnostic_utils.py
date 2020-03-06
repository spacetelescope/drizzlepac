#!/usr/bin/env python

"""Code related to the creation and modification of diagnostic .json files"""

import collections
from datetime import datetime
import json
import os
import pdb
import random
import string
import sys

from astropy.table import Table
import numpy as np

from drizzlepac.hlautils import get_git_rev_info
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
        # gobble up inputs
        self.prop_id = prop_id
        self.obset_id = obset_id
        self.telescope = telescope
        self.instrument = instrument
        self.detector =  detector
        self.filter = filter
        self.data_source = data_source
        self.description = description

        #set logging level
        log.setLevel(log_level)

        #instantiate data storage dictionary
        self._instantiate()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def _astropy_table_to_dict(self, table):
        """Convert Astropy Table to Python dict.

        Numpy arrays are converted to lists, so that
        the output is JSON serialisable.

        Can work with multi-dimensional array columns,
        by representing them as list of list.

        Method based on code written by Christoph Deil.
        URL = https://github.com/astropy/astropy/issues/4604#issuecomment-184551578

        Parameters
        ----------
        table : astropy.table.table.Table
            Astropy table to be converted to python dictionary format.
        """
        total_data = collections.OrderedDict()
        for colname in table.colnames:
            total_data[colname] = collections.OrderedDict()
            total_data[colname]['dtype'] = str(table[colname].dtype)
            total_data[colname]['unit'] = str(table[colname].unit)
            total_data[colname]['format'] = table[colname].format
            total_data[colname]['data'] = table[colname].tolist()
        return total_data

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def _instantiate(self):
        """Creates a new diagnostic dictionary using the standard format. The standard header values are as follows:

        - Generation date
        - Generation time (local 24-hour format)
        - git commit ID
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
        # Add generation date/time
        timestamp = datetime.now().strftime("%m/%d/%YT%H:%M:%S")
        self.out_dict['header']['generation date'] = timestamp.split("T")[0]
        self.out_dict['header']['generation time'] = timestamp.split("T")[1]
        # add git commit id
        reporootpath = "/"
        for item in __file__.split("/")[0:-3]:
            reporootpath = os.path.join(reporootpath, item)
        self.out_dict['header']['commit id'] = get_git_rev_info.get_rev_id(reporootpath)
        del reporootpath
        header_item_list = ["prop_id", "obset_id", "telescope", "instrument", "detector", "filter", "data_source",
                            "description"]
        for header_item in header_item_list:
            self.out_dict['header'][header_item] = self.__dict__[header_item]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def addDataItem(self,dataset,title):
        """main subroutine for adding data to self.out_table.

        Supported data types:

        - all basic single-value python datatypes (float, int, string, Boolean, etc.)
        - lists
        - simple key-value dictionaries/ordered dictionaries
        - multi-layer nested dictionaries and ordered dictionaries
        - tuples
        - numpy arrays
        - astropy tables

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
        self.out_dict['data'][title] = collections.OrderedDict()
        self.out_dict['data'][title]["original format"] = dataset_type
        if dataset_type == "<class 'numpy.ndarray'>": #For numpy arrays
            self.out_dict['data'][title]["dtype"] = str(dataset.dtype)
            self.out_dict['data'][title]["data"] = dataset.tolist()
        elif dataset_type =="<class 'astropy.table.table.Table'>": # for astropy tables
            self.out_dict['data'][title]["data"] = self._astropy_table_to_dict(dataset)
        else: # For everything else. Add more types!
            self.out_dict['data'][title]["original format"] = dataset_type
            self.out_dict['data'][title]["data"] = dataset

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


#----------------------------------------------------------------------------------------------------------------------
def dict_to_astropy_table(in_dict):
    """Converts a dictionary back to astropy table format.

    Parameters
    ----------
    in_dict : dictionary
        dataset to convert back to astropy table format

    Returns
    -------
    out_table : astropy.table.table.Table
        astropy table generated in_data
    """
    colname_list = []
    dtype_list = []
    data_list = []
    for colname in in_dict.keys(): # load up lists by dictionary item type in preparation for table generation
        colname_list.append(colname)
        dtype_list.append(in_dict[colname]['dtype'])
        data_list.append(in_dict[colname]['data'])
    out_table = Table(data_list, names=colname_list, dtype=dtype_list) # generate table, but without units or format details
    for colname in colname_list: # add units and format details
        out_table[colname].unit = in_dict[colname]['unit']
        out_table[colname].format = in_dict[colname]['format']
    return out_table


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def readJsonFile(json_filename):
    """extracts header and data sections from specified json file and returns the header and data (in it's original
    pre-json format) as a nested ordered dictionary

    Supported output data types:

    - all basic single-value python datatypes (float, int, string, Boolean, etc.)
    - lists
    - simple key-value dictionaries and ordered dictionaries
    - multi-layer nested dictionaries and ordered dictionaries
    - tuples
    - numpy arrays
    - astropy tables

    Parameters
    ----------
    json_filename : str
        Name of the json file to extract data from

    Returns
    -------
    out_dict : dictionary
        dictionary structured similarly to self.out_dict with separate 'header' and 'data' keys. The
        information stored in the 'data' section will be in the same format that it was in before it was serialized
        and stored as a json file.
    """
    if os.path.exists(json_filename):
        out_dict = collections.OrderedDict()
        with open(json_filename) as f:
            json_data = json.load(f)
        out_dict['header'] = json_data['header'] # copy over the 'header' section directly.
        out_dict['data'] = collections.OrderedDict() # set up blank data section
        for datakey in json_data['data'].keys():
            if json_data['data'][datakey]['original format'] == "<class 'numpy.ndarray'>": # Extract numpy array
                log.info("Converting dataset '{}' back to format '{}', dtype = {}".format(datakey,
                                                                                          json_data['data'][datakey]['original format'],
                                                                                          json_data['data'][datakey]['dtype']))
                out_dict['data'][datakey] = np.asarray(json_data['data'][datakey]['data'],
                                                       dtype=json_data['data'][datakey]['dtype'])
            elif json_data['data'][datakey]['original format'] == "<class 'astropy.table.table.Table'>": #Extract astropy tables
                log.info("Converting dataset '{}' back to format '{}'".format(datakey,
                                                                              json_data['data'][datakey]['original format']))
                out_dict['data'][datakey] = dict_to_astropy_table(json_data['data'][datakey]['data'])
            elif json_data['data'][datakey]['original format'] == "<class 'tuple'>": # Extract tuples
                out_dict['data'][datakey] = tuple(json_data['data'][datakey]['data'])
            else: # Catchall for everything else
                out_dict['data'][datakey] = json_data['data'][datakey]['data']

    else:
        errmsg = "json file {} not found!".format(json_filename)
        log.error(errmsg)
        raise Exception(errmsg)
    return(out_dict)
# ======================================================================================================================
if __name__ == "__main__":
    blarg = HapDiagnosticObj(telescope = "hst",
                             instrument = "wfc3",
                             detector = "ir",
                             filter = "f160w",
                             prop_id = "11979",
                             obset_id = "01",
                             data_source = "hla_flag_filter",
                             description = "test item please ignore",
                             log_level=10)
    catfile = "hst_11665_06_wfc3_ir_f160w_ib4606_point-cat.ecsv"
    catdata = Table.read(catfile, format='ascii.ecsv')
    blarg.addDataItem(catdata,"CATALOG")
    test_tuple = (True,None,"A",4,5,6,7,8,9,10)
    blarg.addDataItem(test_tuple, "test_tuple")
    test_nested_dict = {}
    test_nested_dict["a"] = "AA"
    test_nested_dict["b"] = {}
    test_nested_dict["b"]["b0"] = "BA"
    test_nested_dict["b"]["b1"] = "BB"
    test_nested_dict["b"]["b2"] = {}
    test_nested_dict["b"]["b2"]["BB0"] = "BBB"
    blarg.addDataItem(test_nested_dict, "test_nested_dict")
    blarg.writeJsonFile("diag_test.json", clobber=True)

    foo = readJsonFile("diag_test.json")
    pdb.set_trace()