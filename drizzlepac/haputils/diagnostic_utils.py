#!/usr/bin/env python

"""Code related to the creation and modification of diagnostic .json files"""

import collections
from datetime import datetime
import glob
import json
import os
import pdb
import random
import re
import string
import sys
import time

from astropy.io.fits import getheader
from astropy.table import Table
import numpy as np

from drizzlepac.haputils import get_git_rev_info
from drizzlepac.haputils import poller_utils
from stsci.tools import logutil


__taskname__ = 'diagnostic_utils'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)
# ======================================================================================================================


class HapDiagnostic(object):
    def __init__(self, log_level=logutil.logging.NOTSET):
        """base class used to set up a HapDiagnostic object.

        Parameters
        ----------
        log_level : int, optional
            The desired level of verboseness in the log statements displayed on the screen and written to the
            .log file. Default value is 'NOTSET'.

        Returns
        -------
        Nothing.
        """
        log.setLevel(log_level)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def _astropy_table_to_dict(self, table):
        """Convert Astropy Table to Python dict.

        Numpy arrays are converted to lists, so that
        the output is JSON serializable.

        Can work with multi-dimensional array columns,
        by representing them as list of list.

        Method based on code written by Christoph Deil.
        URL = https://github.com/astropy/astropy/issues/4604#issuecomment-184551578

        Parameters
        ----------
        table : astropy.table.table.Table
            Astropy table to be converted to python dictionary format.

        Returns
        -------
        total_data : dict
            input astropy table 'table' converted to dictionary form.
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
        """Creates a new diagnostic dictionary using the following standardized format:

        - 'header' section: Contains the primary fits header of the relevant image
        - 'general information' section: Kind of like the 'header' section of distilled down to just the most
        important pieces of information, some other additional information. Included fields are as follows:
            - telescope name
            - proposal ID
            - visit number
            - instrument name
            - detector name
            - filter name
            - dataset name
            - pandas DataFrame index title
            - fits image name
            - generation date (local)
            - generation time (local)
            - seconds since epoch (UTC)
            - git commit ID
            - data source (name of the piece of code that produced the data)
            - description (brief description of what the data is, and how it should be used)

        - 'data' section: this section will contain the test results and depending on the test, may also
        contain relevant test data as well. It should be noted that depending on the test being run,
        additional 'data' sections may be appended to the diagnostic dictionary after instantiation.

        Parameters
        ----------
        Nothing.

        Updates
        -------
        self.out_dict : Ordered dictionary
            dictionary that will ultimately be written to a json file
        """
        # Trim off select items from self.header
        header_items_to_remove = ['',
                                  'HISTORY',
                                  'COMMENT']
        if self.header['INSTRUME'] == 'ACS':
            header_items_to_remove.append('FILTER1')
            header_items_to_remove.append('FILTER2')
        if self.header['INSTRUME'] == 'WFC3':
            header_items_to_remove.append('FILTER')

        header_patterns_to_remove = ["D\d+VER",
                                     "D\d+GEOM",
                                     "D\d+DATA",
                                     "D\d+DEXP",
                                     "D\d+OUDA",
                                     "D\d+OUWE",
                                     "D\d+OUCO",
                                     "D\d+MASK",
                                     "D\d+WTSC",
                                     "D\d+KERN",
                                     "D\d+PIXF",
                                     "D\d+COEF",
                                     "D\d+OUUN",
                                     "D\d+FVAL",
                                     "D\d+WKEY",
                                     "D\d+SCAL",
                                     "D\d+ISCL"]

        # summon nested orderedDict into existence
        self.out_dict = collections.OrderedDict()
        self.out_dict['header'] = collections.OrderedDict()
        self.out_dict['general information'] = collections.OrderedDict()
        self.out_dict['data'] = collections.OrderedDict()

        # add trimmed fits header from self.header
        # and also generate the 'general information' section.
        for header_item in self.header.keys():
            self.out_dict['header'][header_item] = self.header[header_item]

        # use regex to find additional drizzle header keywords to remove
        for pattern in header_patterns_to_remove:
            r = re.compile(pattern)
            more_header_items_to_remove = list(filter(r.match, self.out_dict['header'].keys()))
            header_items_to_remove += more_header_items_to_remove

        for header_item_to_remove in header_items_to_remove:
            if header_item_to_remove in self.out_dict['header'].keys():
                del(self.out_dict['header'][header_item_to_remove])
        
        # Now populate the general information section
        dict_keys = {"TELESCOP": "telescope", 
                     "PROPOSID": "proposal_id", 
                     "INSTRUME": "instrument", 
                     "DETECTOR": "detector"}
        for key in dict_keys:
            self.out_dict['general information'][dict_keys[key]] = self.header[key]
        # Now, add items which require more interpretation
        try:
            self.out_dict['general information']['visit'] = self.header['linenum'].split(".")[0]
        except:
            self.out_dict['general information']['visit'] = self.header['filename'].split("_")[2]
        # determine filter...
        filter_names = ';'.join([self.header[f] for f in self.header['filter*']])
        self.out_dict['general information']['filter'] = poller_utils.determine_filter_name(filter_names)

        rootname = self.header['rootname'].split('_')
        if len(rootname) > 2:
            # This case is the SVM-compatible filename format
            dataset = rootname[-1]
        else:
            # Pipeline default filename format 
            dataset = rootname[0]
        self.out_dict['general information']['dataset'] = dataset
        
        self.out_dict['general information']["dataframe_index"] = self.out_dict['header']['FILENAME'][:-9]
        self.out_dict['general information']["imgname"] = self.out_dict['header']['FILENAME']
        # Add generation date/time
        if self.timestamp:
            timestamp = self.timestamp
        else:
            timestamp = datetime.now().strftime("%m/%d/%YT%H:%M:%S")
        self.out_dict['general information']['generation date'] = timestamp.split("T")[0]  # TODO: is 'generation date' too generic? should this be renamed something more descriptive?
        self.out_dict['general information']['generation time'] = timestamp.split("T")[1]  # TODO: is 'generation date' too generic? should this be renamed something more descriptive?
        # Add time since epoch (January 1, 1970, 00:00:00 UTC)
        if self.time_since_epoch:
            time_since_epoch = self.time_since_epoch
        else:
            time_since_epoch = time.time()
        self.out_dict['general information']['seconds since epoch'] = time_since_epoch
        # add git commit id
        reporootpath = "/"
        for item in __file__.split("/")[0:-3]:
            reporootpath = os.path.join(reporootpath, item)
        self.out_dict['general information']['commit id'] = get_git_rev_info.get_rev_id(reporootpath)
        del reporootpath
        # add data_source and description # TODO: THESE MAY BE REMOVED LATER ON
        header_item_list = ["data_source", "description"]
        for header_item in header_item_list:
            self.out_dict['general information'][header_item] = self.__dict__[header_item]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def instantiate_from_fitsfile(self, filename, data_source=None, description=None, timestamp=None, time_since_epoch=None):
        """Get necessary information for execution of _instantiate() from user-specified hap product, and
        execute _instantiate()

        Parameters
        ----------
        filename : str
            fits file to pull header information and filter information from to populate the json
            "header" section.

        data_source : str, optional
            name of the script that generated the data that will be stored in the "data" section.  If not
            specified, default value is logical 'None'

        description : str, optional
            brief description of what the data is, and how it should be used.  If not specified, default
            value is logical 'None'

        timestamp: str, optional
            .json file generation date and time (local timezone). Format: MM/DD/YYYYTHH:MM:SS
            (Example: 05/04/2020T13:46:35). If not specified, default value is logical 'None'

        time_since_epoch : float
            .json file generation time. Format: Time (in seconds) elapsed since
            January 1, 1970, 00:00:00 (UTC). If not specified, default value is logical 'None'

        Returns
        -------
        Nothing.
        """
        if os.path.exists(filename):
            self.header = getheader(filename)
            if self.header['INSTRUME'].lower() == "acs":
                self.filter = poller_utils.determine_filter_name("{};{}".format(self.header['FILTER1'], self.header['FILTER2']))
            elif self.header['INSTRUME'].lower() == "wfc3":
                self.filter = poller_utils.determine_filter_name(self.header['FILTER'])
            else:
                errmsg = "Invalid instrument."
                log.error(errmsg)
                raise Exception(errmsg)
        else:
            errmsg = "Invalid input. File {} does not exist.".format(filename)
            log.error(errmsg)
            raise Exception(errmsg)

        # gobble up other inputs
        self.data_source = data_source
        self.description = description
        self.timestamp = timestamp
        self.time_since_epoch = time_since_epoch

        # instantiate data storage dictionary
        self._instantiate()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def instantiate_from_hap_obj(self, hap_obj, data_source=None, description=None, timestamp=None, time_since_epoch=None):
        """Get necessary information for execution of _instantiate() from user-specified hap product, and
        execute _instantiate()

        Parameters
        ----------
        header_source : drizzlepac.haputils.Product.TotalProduct, drizzlepac.haputils.Product.FilterProduct, or
        drizzlepac.haputils.Product.ExposureProduct, depending on input.
            hap product object to pull header information and filter information from to populate the json
        "header" section.

        data_source : str, optional
            name of the script that generated the data that will be stored in the "data" section. If not
            specified, default value is logical 'None'

        description : str, optional
            brief description of what the data is, and how it should be used. If not specified,
            default value is logical 'None'

        timestamp: str, optional
            .json file generation date and time (local timezone). Format: MM/DD/YYYYTHH:MM:SS
            (Example: 05/04/2020T13:46:35). If not specified, default value is logical 'None'

        time_since_epoch : float
            .json file generation time. Format: Time (in seconds) elapsed since
            January 1, 1970, 00:00:00 (UTC). If not specified, default value is logical 'None'

        Returns
        -------
        Nothing.
        """
        self.header = hap_obj.primary_header.copy()
        if hasattr(hap_obj, "filters"):
            self.filter = hap_obj.filters
        else:
            if self.header['INSTRUME'].lower() == "acs":
                self.filter = poller_utils.determine_filter_name(
                    "{};{}".format(self.header['FILTER1'], self.header['FILTER2']))
            elif self.header['INSTRUME'].lower() == "wfc3":
                self.filter = poller_utils.determine_filter_name(self.header['FILTER'])
            else:
                errmsg = "Invalid instrument."
                log.error(errmsg)
                raise Exception(errmsg)

        # gobble up other inputs
        self.data_source = data_source
        self.description = description
        self.timestamp = timestamp
        self.time_since_epoch = time_since_epoch

        # instantiate data storage dictionary
        self._instantiate()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def add_data_item(self, dataset, title, item_description="", descriptions=None, units=None):
        """main subroutine for adding data to self.out_table.

        Supported data types:

        - all basic single-value python data types (float, int, string, Boolean, etc.)
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
            
        item_description : str
            Single string description for this item as a whole

        descriptions : dict, optional
            dictionary containing description strings for each element of the dataset stored in the 'data'
            section

        units : dict, optional
            dictionary containing units for each element of the dataset stored in the 'data' section

        Updates
        -------
        self.out_dict : Ordered dictionary
            dictionary that will ultimately be written to a json file
        """
        dataset_type = str(type(dataset))
        self.out_dict['data'][title] = collections.OrderedDict()
        self.out_dict['data'][title]["original format"] = dataset_type
        self.out_dict['data'][title]["description"] = item_description
        if dataset_type == "<class 'numpy.ndarray'>":  # For numpy arrays
            self.out_dict['data'][title]["dtype"] = str(dataset.dtype)
            self.out_dict['data'][title]["data"] = dataset.tolist()
            self.out_dict['data'][title]["descriptions"] = descriptions
            self.out_dict['data'][title]["units"] = units
        elif dataset_type == "<class 'astropy.table.table.Table'>":  # for astropy tables
            self.out_dict['data'][title]["data"] = self._astropy_table_to_dict(dataset)
            self.out_dict['data'][title]["descriptions"] = descriptions
            self.out_dict['data'][title]["units"] = units
        else:  # For everything else. Add more types!
            self.out_dict['data'][title]["original format"] = dataset_type
            self.out_dict['data'][title]["data"] = dataset
            self.out_dict['data'][title]["descriptions"] = descriptions
            self.out_dict['data'][title]["units"] = units

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def add_update_info_item(self, section_name, element_name, new_element_value, clobber=True, addnew=True):
        """add or update a single user-specified 'header' or 'general information' item

        Parameters
        ----------
        section_name : str
            dictionary section to update. Choices are either 'header' or 'general information'.

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
        if section_name in ['header', 'general information']:
            if element_name in self.out_dict[section_name].keys():
                if clobber:
                    log.info("{}: {} -> {} value update successful.".format(element_name,
                                                                            self.out_dict[section_name][element_name],
                                                                            new_element_value))
                    self.out_dict['header'][element_name] = new_element_value
                else:
                    log.warning("{} element '{}' already exists. Update NOT performed.".format(section_name, element_name))
            else:
                if addnew:
                    self.out_dict[section_name][element_name] = new_element_value
                    log.info("New element {} = {} successfully added to {}".format(element_name, new_element_value, section_name))
                else:
                    log.warning("Unable to add new element {} = {} to {}".format(element_name, new_element_value, section_name))
        else:
            log.warning("**NO UPDATES PERFORMED** Only the 'header' or 'general information' sections that can be updated.")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def write_json_file(self, json_filename, clobber=False):
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
        # TODO: JSON filename should have a common suffix; perhaps, '*diag.json'. (much like cat.ecsv)
        file_exists = os.path.exists(json_filename)
        if clobber:
            if file_exists:
                os.remove(json_filename)
        else:
            if file_exists:
                random_string = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(4))
                json_filename = json_filename.replace(".json", "_{}.json".format(random_string))
        with open(json_filename, "w") as json_file:
            json.dump(self.out_dict, json_file, indent=4)
        log.info("Wrote json file {}".format(json_filename))


def dict_to_astropy_table(in_dict):
    """Converts an astropy table stored as a dictionary back to astropy table format.

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
    for colname in in_dict.keys():  # load up lists by dictionary item type in preparation for table generation
        colname_list.append(colname)
        dtype_list.append(in_dict[colname]['dtype'])
        data_list.append(in_dict[colname]['data'])
    out_table = Table(data_list, names=colname_list, dtype=dtype_list)  # generate table, but without units or format details
    for colname in colname_list:  # add units and format details
        out_table[colname].unit = in_dict[colname]['unit']
        out_table[colname].format = in_dict[colname]['format']
    return out_table


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def read_json_file(json_filename):
    """extracts header and data sections from specified json file and returns the header and data (in it's original
    pre-json format) as a nested ordered dictionary

    Supported output data types:

    - all basic single-value python data types (float, int, string, Boolean, etc.)
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

        out_dict['header'] = json_data['header']  # copy over the 'header' section directly.
        out_dict['general information'] = json_data['general information']
        out_dict['data'] = collections.OrderedDict()  # set up blank data section
        out_dict['descriptions'] = collections.OrderedDict()
        out_dict['units'] = collections.OrderedDict()
        for datakey in json_data['data'].keys():
            out_dict['descriptions'][datakey] = json_data['data'][datakey]['descriptions']
            out_dict['units'][datakey] = json_data['data'][datakey]['units']
            if json_data['data'][datakey]['original format'] == "<class 'numpy.ndarray'>":  # Extract numpy array
                log.info("Converting dataset '{}' back to format '{}', dtype = {}".format(datakey,
                                                                                          json_data['data'][datakey]['original format'],
                                                                                          json_data['data'][datakey]['dtype']))
                out_dict['data'][datakey] = np.asarray(json_data['data'][datakey]['data'],
                                                       dtype=json_data['data'][datakey]['dtype'])
            elif json_data['data'][datakey]['original format'] == "<class 'astropy.table.table.Table'>":  # Extract astropy tables
                log.info("Converting dataset '{}' back to format '{}'".format(datakey,
                                                                              json_data['data'][datakey]['original format']))
                out_dict['data'][datakey] = dict_to_astropy_table(json_data['data'][datakey]['data'])
            elif json_data['data'][datakey]['original format'] == "<class 'tuple'>":  # Extract tuples
                out_dict['data'][datakey] = tuple(json_data['data'][datakey]['data'])
            else:  # Catchall for everything else
                out_dict['data'][datakey] = json_data['data'][datakey]['data']

    else:
        errmsg = "json file {} not found!".format(json_filename)
        log.error(errmsg)
        raise Exception(errmsg)
    return(out_dict)
# ======================================================================================================================


if __name__ == "__main__":
    # Testing
    header_fits_filename = "hst_10265_01_acs_wfc_f606w_j92c01_drc.fits"
    print(header_fits_filename)
    blarg = HapDiagnostic(log_level=10)
    blarg.instantiate_from_fitsfile(header_fits_filename, data_source="hla_flag_filter", description="test item please ignore",)
    catfile = glob.glob("*_point-cat.ecsv")[0]
    catdata = Table.read(catfile, format='ascii.ecsv')
    blarg.add_data_item(catdata, "CATALOG")
    test_tuple = (True, None, "A", 4, 5, 6, 7, 8, 9, 10)
    blarg.add_data_item(test_tuple, "test_tuple")
    test_nested_dict = {}
    test_nested_dict["a"] = "AA"
    test_nested_dict["b"] = {}
    test_nested_dict["b"]["b0"] = "BA"
    test_nested_dict["b"]["b1"] = "BB"
    test_nested_dict["b"]["b2"] = {}
    test_nested_dict["b"]["b2"]["BB0"] = "BBB"
    blarg.add_data_item(test_nested_dict, "test_nested_dict")
    blarg.write_json_file("diag_test.json", clobber=True)

    foo = read_json_file("diag_test.json")
    pdb.set_trace()


