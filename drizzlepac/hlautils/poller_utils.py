"""Utilities to interpret the pipeline poller obset information and generate product filenames

The function, interpret_obset_input, parses the file generated by the pipeline
poller, and produces a tree listing of the output products.  The function,
parse_obset_tree, converts the tree into product catagories.

"""
import os
import sys
from collections import OrderedDict
import numpy as np

from sklearn.cluster import KMeans

from stsci.tools import logutil

from astropy.io import fits
from astropy.io import ascii
from astropy.io.fits import getheader
from astropy.table import Table, Column
from drizzlepac.hlautils.product import ExposureProduct, FilterProduct, TotalProduct
from . import analyze
from . import astroquery_utils as aqutils
from . import processing_utils

# Define information/formatted strings to be included in output dict
SEP_STR = 'single exposure product {:02d}'
FP_STR = 'filter product {:02d}'
TDP_STR = 'total detection product {:02d}'

# Define the mapping between the first character of the filename and the associated instrument
INSTRUMENT_DICT = {'i': 'WFC3', 'j': 'ACS', 'o': 'STIS', 'u': 'WFPC2', 'x': 'FOC', 'w': 'WFPC'}
POLLER_COLNAMES = ['filename', 'proposal_id', 'program_id', 'obset_id',
                   'exptime', 'filters', 'detector', 'pathname']

EXP_LABELS = {2: 'short', 1: 'med', 0: 'long', None: 'all'}

__taskname__ = 'poller_utils'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)

# -----------------------------------------------------------------------------
# Single Visit Processing Functions
#
def interpret_obset_input(results, log_level):
    """

    Parameters
    -----------
    results : str or list
        Input poller file name or Python list of dataset names to be processed as a single visit.
        Dataset names have to be either the filename of a singleton (non-associated exposure) or the
        name of an ASN file (e.g., jabc01010_asn.fits).

    log_level : int
        The desired level of verboseness in the log statements displayed on the screen and written to the .log file.


    Notes
    -------
    Interpret the database query for a given obset to prepare the returned
    values for use in generating the names of all the expected output products.

    Input will have format of:
        ib4606c5q_flc.fits,11665,B46,06,1.0,F555W,UVIS,/ifs/archive/ops/hst/public/ib46/ib4606c5q/ib4606c5q_flc.fits
        which are
        filename, proposal_id, program_id, obset_id, exptime, filters, detector, pathname

    Output dict will have format (as needed by further code for creating the
        product filenames) of:

        obs_info_dict["single exposure product 00": {"info": '11665 06 wfc3 uvis ib4606c5q f555w drc',
                                                     "files": ['ib4606c5q_flc.fits']}
        .
        .
        .
        obs_info_dict["single exposure product 08": {"info": '11665 06 wfc3 ir ib4606clq f110w drz',
                                                     "files": ['ib4606clq_flt.fits']}

        obs_info_dict["filter product 00": {"info": '11665 06 wfc3 uvis ib4606c5q f555w drc',
                                            "files": ['ib4606c5q_flc.fits', 'ib4606c6q_flc.fits']},
        .
        .
        .
        obs_info_dict["filter product 01": {"info": '11665 06 wfc3 ir ib4606cmq f160w drz',
                                            "files": ['ib4606cmq_flt.fits', 'ib4606crq_flt.fits']},


        obs_info_dict["total detection product 00": {"info": '11665 06 wfc3 uvis ib4606c5q f555w drc',
                                                     "files": ['ib4606c5q_flc.fits', 'ib4606c6q_flc.fits']}
        .
        .
        .
        obs_info_dict["total detection product 01": {"info": '11665 06 wfc3 ir ib4606cmq f160w drz',
                                                     "files": ['ib4606cmq_flt.fits', 'ib4606crq_flt.fits']}

    """
    # set logging level to user-specified level
    log.setLevel(log_level)

    log.debug("Interpret the poller file for the observation set.")
    obset_table = build_poller_table(results, log_level)
    # Add INSTRUMENT column
    instr = INSTRUMENT_DICT[obset_table['filename'][0][0]]
    # convert input to an Astropy Table for parsing
    obset_table.add_column(Column([instr] * len(obset_table)), name='instrument')

    # Sort the rows of the table in an effort to optimize the number of quality sources found in the initial images
    obset_table = sort_poller_table(obset_table)
    log.debug("Sorted input:")
    log.debug(obset_table)
    log.debug("\n\n")
    # parse Table into a tree-like dict
    log.debug("Build the observation set tree.")
    obset_tree = build_obset_tree(obset_table)
    # Now create the output product objects
    log.debug("Parse the observation set tree and create the exposure, filter, and total detection objects.")
    obset_dict, tdp_list = parse_obset_tree(obset_tree, log_level)

    # This little bit of code adds an attribute to single exposure objects that is True 
    # if a given filter only contains one input (e.g. n_exp = 1)
    for tot_obj in tdp_list:
        for filt_obj in tot_obj.fdp_list:
            if len(filt_obj.edp_list) == 1:
                is_singleton = True
            else:
                is_singleton = False
            for edp_obj in filt_obj.edp_list:
                edp_obj.is_singleton = is_singleton

    return obset_dict, tdp_list

# Translate the database query on an obset into actionable lists of filenames
def build_obset_tree(obset_table):
    """Convert obset table into a tree listing all products to be created."""

    # Each product will consist of the appropriate string as the key and
    # a dict of 'info' and 'files' information

    # Start interpreting the obset table
    obset_tree = {}
    for row in obset_table:
        # Get some basic information from the first row - no need to check
        # for multiple instruments as data from different instruments will
        # not be combined.
        det = row['detector']
        orig_filt = row['filters']
        # Potentially need to manipulate the 'filters' string for instruments
        # with two filter wheels
        filt = determine_filter_name(orig_filt)
        row['filters'] = filt
        row_info, filename = create_row_info(row)
        # Initial population of the obset tree for this detector
        if det not in obset_tree:
            obset_tree[det] = {}
            obset_tree[det][filt] = [(row_info, filename)]
        else:
            det_node = obset_tree[det]
            if filt not in det_node:
                det_node[filt] = [(row_info, filename)]
            else:
                det_node[filt].append((row_info, filename))

    return obset_tree


def create_row_info(row):
    """Build info string for a row from the obset table"""
    info_list = [str(row['proposal_id']), "{}".format(row['obset_id']), row['instrument'],
                 row['detector'], row['filename'][:row['filename'].find('_')], row['filters']]
    return ' '.join(map(str.upper, info_list)), row['filename']


# -----------------------------------------------------------------------------
# Multi-Visit Processing Functions
#
def interpret_mvm_input(results, skycells, log_level, exp_limit=2.0):
    """

    Parameters
    -----------
    results : str or list
        Input poller file name or Python list of dataset names to be processed as a single visit.
        Dataset names have to be either the filename of a singleton (non-associated exposure) or the
        name of an ASN file (e.g., jabc01010_asn.fits).

    skycells : str
        Input file name containing the list of existing sky cell layers and the
        exposures names that are members of each layer. This will be used to
        reprocess any layer which has new observations from the `results` input
        that contribute to that layer. If there is no existing layer for that
        skycell that matches the new input exposures, a new layer will be defined.

    log_level : int
        The desired level of verboseness in the log statements displayed on the screen
        and written to the .log file.

    exp_limit : float, optional
        This specifies the ratio between long and short exposures that should
        be used to define each exposure time bin (short, med, long).  If None,
        all exposure times will be combined into a single bin (all). Splitting
        of the bins gets computed using Kmeans clustering for 3 output bins.


    Notes
    -------
    Interpret the database query for a given obset to prepare the returned
    values for use in generating the names of all the expected output products.

    Input will have format of:

        ib4606c5q_flc.fits,11665,B46,06,1.0,F555W,UVIS,/ifs/archive/ops/hst/public/ib46/ib4606c5q/ib4606c5q_flc.fits

        which are
        filename, proposal_id, program_id, obset_id, exptime, filters, detector, pathname
    
    The SkyCell ID will be included in this input information to allow for 
    grouping of exposures into the same SkyCell layer based on filter, exptime and year.

    Output dict will have only have definitions for each defined sky cell layer to
    be either created or updated based on the input exposures being processed.
    There will be a layer defined for each unique combination of filter/exp class/year.
    An example would be:

        obs_info_dict["layer 00": {"info": '11665 06 wfc3 uvis f555w short 2011 drc',
                                   "files": ['ib4606c5q_flc.fits', 'ib4606cgq_flc.fits']}
        .
        .
        .

        obs_info_dict["layer 01": {"info": '11665 06 wfc3 uvis f555w med 2011 drc',
                                   "files": ['ib4606c6q_flc.fits', 'ib4606cdq_flc.fits']}
        .
        .
        .
        obs_info_dict["layer 01": {"info": '11665 06 wfc3 uvis f555w long 2011 drc',
                                   "files": ['ib4606c9q_flc.fits', 'ib4606ceq_flc.fits']}
        .
        .
        .

    """
    # set logging level to user-specified level
    log.setLevel(log_level)

    log.debug("Interpret the poller file for the observation set.")
    obset_table = build_poller_table(results, log_level)
    # Add INSTRUMENT column
    instr = INSTRUMENT_DICT[obset_table['filename'][0][0]]
    # convert input to an Astropy Table for parsing
    obset_table.add_column(Column([instr] * len(obset_table)), name='instrument')

    # Add Date column    
    years = [int(its.getval(f, 'date-obs').split('-')[0]) for f in obset_table['filename']]
    obset_table.add_column(Column(years), name='year')

    # Sort the rows of the table in an effort to optimize the number of quality 
    # sources found in the initial images
    obset_table = sort_poller_table(obset_table)
    exp_obset_table = define_exp_layers(year_obset_table, exp_limit=exp_limit)
        
    # parse Table into a tree-like dict
    log.debug("Build the multi-visit layers tree.")
    obset_tree = build_mvm_tree(obset_table)
    # Now create the output product objects
    log.debug("Parse the observation set tree and create the exposure, filter, and total detection objects.")
    obset_dict, tdp_list = parse_mvm_tree(obset_tree, log_level)
    
    # Now we need to merge any pre-existing layer products with the new products
    # This will add the exposure products from the pre-existing definitions of
    # the overlapping sky cell layers with those defined for the new exposures 
    # 
    # obset_dict, tdp_list = merge_mvm_trees(obset_tree, layer_tree, log_level)
    

    # This little bit of code adds an attribute to single exposure objects that is True 
    # if a given filter only contains one input (e.g. n_exp = 1)
    for tot_obj in tdp_list:
        for filt_obj in tot_obj.fdp_list:
            if len(filt_obj.edp_list) == 1:
                is_singleton = True
            else:
                is_singleton = False
            for edp_obj in filt_obj.edp_list:
                edp_obj.is_singleton = is_singleton

    return obset_dict, tdp_list

# Translate the database query on an obset into actionable lists of filenames
def build_mvm_tree(obset_table):
    """Convert obset table into a tree listing all products to be created."""

    # Each product will consist of the appropriate string as the key and
    # a dict of 'info' and 'files' information

    # Start interpreting the obset table
    obset_tree = {}
    for row in obset_table:
        # Get some basic information from the first row - no need to check
        # for multiple instruments as data from different instruments will
        # not be combined.
        det = row['detector']
        orig_filt = row['filters']
        exp_layer = row['exp_layer']
        year_layer = row['year_layer']

        # Potentially need to manipulate the 'filters' string for instruments
        # with two filter wheels
        filt = determine_filter_name(orig_filt)
        row['filters'] = filt
        layer = (filt, exp_layer, year_layer)
        row_info, filename = create_mvm_info(row)
        # Initial population of the obset tree for this detector
        if det not in obset_tree:
            obset_tree[det] = {}
            obset_tree[det][layer] = [(row_info, filename)]
        else:
            det_node = obset_tree[det]
            if layer not in det_node:
                det_node[layer] = [(row_info, filename)]
            else:
                det_node[layer].append((row_info, filename))

    return obset_tree


def create_mvm_info(row):
    """Build info string for a row from the obset table"""
    info_list = [str(row['proposal_id']), "{}".format(row['obset_id']), row['instrument'],
                 row['detector'], row['filters'], row['exp_layer'], row['year_layer']]
    return ' '.join(map(str.upper, info_list)), row['filename']


def parse_mvm_tree(det_tree, log_level):
    """Convert tree into products for sky cell layers

    Tree generated by `create_mvm_info()` will always have the following
    levels:
          * detector
              * filters
                  * layer
    Each exposure will have an entry dict with keys 'info' and 'filename'.

    Products created will be:
      * a single sky cell layer product per 'detector/filter/layer' combination

    """
    log.setLevel(log_level)

    # Initialize products dict
    obset_products = {}

    # For each level, define a product, starting with the detector used...
    prev_det_indx = 0
    det_indx = 0
    filt_indx = 0
    sep_indx = 0

    tdp_list = []

    # Determine if the individual files being processed are flt or flc and
    # set the filetype accordingly (flt->drz or flc->drc).
    filetype = ''
    # Setup products for each detector used
    #
    # det_tree = {'UVIS': {'f200lp': [('11515 01 WFC3 UVIS IACS01T9Q F200LP', 'iacs01t9q_flt.fits'), ('11515 01 WFC3 UVIS IACS01TBQ F200LP', 'iacs01tbq_flt.fits')]}, 'IR': {'f160w': [('11515 01 WFC3 IR IACS01T4Q F160W', 'iacs01t4q_flt.fits')]}}
    
    for filt_tree in det_tree.values():
        totprod = TDP_STR.format(det_indx)
        obset_products[totprod] = {'info': "", 'files': []}
        det_indx += 1
        # Find all filters used...
        # filt_tree = {'f200lp': [('11515 01 WFC3 UVIS IACS01T9Q F200LP', 'iacs01t9q_flt.fits'), 
        #            ('11515 01 WFC3 UVIS IACS01TBQ F200LP', 'iacs01tbq_flt.fits')]}

        for filter_files in filt_tree.values():
            # Use this to create and populate filter product dictionary entry
            fprod = FP_STR.format(filt_indx)
            obset_products[fprod] = {'info': "", 'files': []}
            filt_indx += 1
            # Populate single exposure dictionary entry now as well
            # filename = ('11515 01 WFC3 UVIS IACS01T9Q F200LP', 'iacs01t9q_flt.fits')
            #
            for filename in filter_files:
                # Parse the first filename[1] to determine if the products are flt or flc
                if det_indx != prev_det_indx:
                    filetype = "drc"
                    if filename[1][10:13].lower().endswith("flt"):
                        filetype = "drz"
                    prev_det_indx = det_indx

                # Generate the full product dictionary information string:
                #
                # [row['proposal_id'], row['obset_id'], row['instrument'],
                # row['detector'], row['filters'], row['exp_layer'], row['year_layer']]
                #
                prod_info = (filename[0] + " " + filetype).lower()
                #
                # svm prod_info = '11515 01 wfc3 uvis f200lp all 2009 drz'
                #
                prod_list = prod_info.split(" ")
                layer = (prod_list[4], prod_list[5], prod_list[6])
                ftype = prod_list[-1]
                
                # Set up the single exposure product dictionary
                sep = SEP_STR.format(sep_indx)
                obset_products[sep] = {'info': prod_info,
                                       'files': [filename[1]]}

                # Create a single exposure product object
                sep_obj = ExposureProduct(prod_list[0], prod_list[1], prod_list[2], prod_list[3],
                                          filename[1], prod_list[4], ftype, log_level)

                # Set up the filter product dictionary and create a filter product object
                # Initialize `info` key for this filter product dictionary
                if not obset_products[fprod]['info']:
                    obset_products[fprod]['info'] = prod_info

                    # Create a filter product object for this instrument/detector
                    # FilterProduct(prop_id, obset_id, instrument, detector, 
                    #               filename, filters, filetype, log_level)
                    filt_obj = SkyCellProduct(prod_list[0], prod_list[1], prod_list[2], prod_list[3],
                                             'layerwcs', layer, ftype, log_level)
                # Append exposure object to the list of exposure objects for this specific filter product object
                filt_obj.add_member(sep_obj)
                # Populate filter product dictionary with input filename
                obset_products[fprod]['files'].append(filename[1])

                # Increment single exposure index
                sep_indx += 1

            # Append filter object to the list of filter objects for this specific total product object
            log.debug("Attach the sky cell layer object {} to its associated total product object {}/{}.".format(filt_obj.filters,
                                                                                                                   tdp_obj.instrument,
                                                                                                                   tdp_obj.detector))
            tdp_obj.add_product(filt_obj)

        # Add the total product object to the list of TotalProducts
        tdp_list.append(tdp_obj)

    # Done... return dict and object product list
    return obset_products, tdp_list


# -----------------------------------------------------------------------------
# Utility Functions
#
def parse_obset_tree(det_tree, log_level):
    """Convert tree into products

    Tree generated by `create_row_info()` will always have the following
    levels:
          * detector
              * filters
                  * exposure
    Each exposure will have an entry dict with keys 'info' and 'filename'.

    Products created will be:
      * total detection product per detector
      * filter products per detector
      * single exposure product
    """
    log.setLevel(log_level)

    # Initialize products dict
    obset_products = {}

    # For each level, define a product, starting with the detector used...
    prev_det_indx = 0
    det_indx = 0
    filt_indx = 0
    sep_indx = 0

    tdp_list = []

    # Determine if the individual files being processed are flt or flc and
    # set the filetype accordingly (flt->drz or flc->drc).
    filetype = ''
    # Setup products for each detector used
    for filt_tree in det_tree.values():
        totprod = TDP_STR.format(det_indx)
        obset_products[totprod] = {'info': "", 'files': []}
        det_indx += 1
        # Find all filters used...
        for filter_files in filt_tree.values():
            # Use this to create and populate filter product dictionary entry
            fprod = FP_STR.format(filt_indx)
            obset_products[fprod] = {'info': "", 'files': []}
            filt_indx += 1
            # Populate single exposure dictionary entry now as well
            for filename in filter_files:
                # Parse the first filename[1] to determine if the products are flt or flc
                if det_indx != prev_det_indx:
                    filetype = "drc"
                    if filename[1][10:13].lower().endswith("flt"):
                        filetype = "drz"
                    prev_det_indx = det_indx

                # Generate the full product dictionary information string:
                # proposal_id, obset_id, instrument, detector, ipppssoot, filter, and filetype
                # filename = ('11515 01 WFC3 UVIS IACS01T9Q F200LP', 'iacs01t9q_flt.fits')
                prod_info = (filename[0] + " " + filetype).lower()

                # Set up the single exposure product dictionary
                sep = SEP_STR.format(sep_indx)
                obset_products[sep] = {'info': prod_info,
                                       'files': [filename[1]]}

                # Create a single exposure product object
                prod_list = prod_info.split(" ")
                sep_obj = ExposureProduct(prod_list[0], prod_list[1], prod_list[2], prod_list[3],
                                          filename[1], prod_list[5], prod_list[6], log_level)
                # Set up the filter product dictionary and create a filter product object
                # Initialize `info` key for this filter product dictionary
                if not obset_products[fprod]['info']:
                    obset_products[fprod]['info'] = prod_info

                    # Create a filter product object for this instrument/detector
                    filt_obj = FilterProduct(prod_list[0], prod_list[1], prod_list[2], prod_list[3],
                                             prod_list[4], prod_list[5], prod_list[6], log_level)
                # Append exposure object to the list of exposure objects for this specific filter product object
                filt_obj.add_member(sep_obj)
                # Populate filter product dictionary with input filename
                obset_products[fprod]['files'].append(filename[1])

                # Set up the total detection product dictionary and create a total detection product object
                # Initialize `info` key for total detection product
                if not obset_products[totprod]['info']:
                    obset_products[totprod]['info'] = prod_info

                    # Create a total detection product object for this instrument/detector
                    tdp_obj = TotalProduct(prod_list[0], prod_list[1], prod_list[2], prod_list[3],
                                           prod_list[4], prod_list[6], log_level)

                # Append exposure object to the list of exposure objects for this specific total detection product
                tdp_obj.add_member(sep_obj)
                # Populate total detection product dictionary with input filename
                obset_products[totprod]['files'].append(filename[1])

                # Increment single exposure index
                sep_indx += 1

            # Append filter object to the list of filter objects for this specific total product object
            log.debug("Attach the filter object {} to its associated total detection product object {}/{}.".format(filt_obj.filters,
                                                                                                                   tdp_obj.instrument,
                                                                                                                   tdp_obj.detector))
            tdp_obj.add_product(filt_obj)

        # Add the total product object to the list of TotalProducts
        tdp_list.append(tdp_obj)

    # Done... return dict and object product list
    return obset_products, tdp_list

# ------------------------------------------------------------------------------


def determine_filter_name(raw_filter):
    """
    Generate the final filter name to be used for an observation.

    Parameters
    ----------
    raw_filter : string
        filters component one exposure from an input observation visit

    Returns
    -------
    filter_name : string
        final filter name

    If the raw_filter is actually a combination of two filter names, as
    can be true for instruments with two filter wheels, then generate
    a new filter string according the following rules:
    - If one filter name is 'clear*', then use the other filter name.
    - If both filter names are 'clear*', then use 'clear'.
    - If there are two filters in use, then use 'filter1-filter2'.
    - If one filter is a polarizer ('pol*'), then always put the polarizer
      name second (e.g., 'f606w-pol60').
    - If one filter is 'n/a' (as is the case for SBC), then use the other filter
      name. NOTE: At this time (February 2020) the filter names for the SBC
      detector of ACS have malformed names (e.g., F140LP;N/A;F140LP) where the
      final token delimited by the ";" should not be present. Remove the final
      delimiter and entry.  Unlike the other ACS detectors, SBC only uses a
      single filter wheel.

    - NOTE: There should always be at least one filter name provided to
      this routine or this input is invalid.
    """

    raw_filter = raw_filter.lower()

    # There might be multiple filters, so split the filter names into a list
    # and only retain the first two entries.  SBC has a bogus third entry.
    filter_list = raw_filter.split(';')[0:2]
    output_filter_list = []

    for filt in filter_list:
        if not any(x in filt for x in ['clear', 'n/a']):
            output_filter_list.append(filt)

    if not output_filter_list:
        output_filter_list = ['clear']
    else:
        if output_filter_list[0].startswith('pol'):
            output_filter_list.reverse()

    delimiter = '-'
    filter_name = delimiter.join(output_filter_list)

    return filter_name

# ------------------------------------------------------------------------------


def build_poller_table(input, log_level):
    """Create a poller file from dataset names.

    Parameters
    -----------
    input : str, list
        Filename with list of dataset names, or just a Python list of dataset names, provided by the user.

    Returns
    --------
    poller_table : Table
        Astropy table object with the same columns as a poller file.

    """
    log.setLevel(log_level)

    # Check the input file is not empty
    if not os.path.getsize(input):
        log.error('Input poller manifest file, {}, is empty - processing is exiting.'.format(input))
        sys.exit(0)

    datasets = []
    is_poller_file = False
    obs_converters = {'col4': [ascii.convert_numpy(np.str)]}
    if isinstance(input, str):
        input_table = ascii.read(input, format='no_header', converters=obs_converters)
        if len(input_table.columns) == len(POLLER_COLNAMES):
            # We were provided a poller file
            # Now assign column names to table
            for i, colname in enumerate(POLLER_COLNAMES):
                input_table.columns[i].name = colname

            # Convert to a string column, instead of int64
            input_table['obset_id'] = input_table['obset_id'].astype(np.str)
            is_poller_file = True

        elif len(input_table.columns) == 1:
            input_table.columns[0].name = 'filename'
            is_poller_file = False

        # Since a poller file was the input, it is assumed all the input
        # data is in the locale directory so just collect the filenames.
        filenames = list(input_table.columns[0])

    elif isinstance(input, list):
        filenames = input

    else:
        id = '[poller_utils.build_poller_table] '
        log.error("{}: Input {} not supported as input for processing.".format(id, input))
        raise ValueError

    # At this point, we have a poller file or a list of filenames.  If the latter, then any individual
    # filename can be a singleton or an association name.  We need to get the full list of actual
    # filenames from the association name.
    if not is_poller_file:
        for filename in filenames:
            # Look for dataset in local directory.
            if "asn" in filename or not os.path.exists(filename):
                # This retrieval will NOT overwrite any ASN members already on local disk
                # Return value will still be list of all members
                files = aqutils.retrieve_observation([filename[:9]], suffix=['FLC'], clobber=False)
                if len(files) == 0:
                    log.error("Filename {} not found in archive!!".format(filename))
                    log.error("Please provide ASN filename instead!")
                    raise ValueError
            else:
                files = [filename]
            datasets += files
    else:
        datasets = filenames
        
    # Each image, whether from a poller file or from an input list needs to be
    # analyzed to ensure it is viable for drizzle processing.  If the image is not
    # viable, it should not be included in the output "poller" table.
    usable_datasets = analyze.analyze_wrapper(datasets)
    if not usable_datasets:
        log.warning("No usable images in poller file or input list for drizzling. The processing of this data is ending.")
        sys.exit(0)

    cols = OrderedDict()
    for cname in POLLER_COLNAMES:
        cols[cname] = []
    cols['filename'] = usable_datasets

    # If processing a list of files, evaluate each input dataset for the information needed
    # for the poller file
    if not is_poller_file:
        for d in usable_datasets:
            with fits.open(d) as dhdu:
                hdr = dhdu[0].header
                cols['program_id'].append(d[1:4].upper())
                cols['obset_id'].append(str(d[4:6]))
                cols['proposal_id'].append(hdr['proposid'])
                cols['exptime'].append(hdr['exptime'])
                cols['detector'].append(hdr['detector'])
                cols['pathname'].append(os.path.abspath(d))
                # process filter names
                if d[0] == 'j':  # ACS data
                    filters = processing_utils.get_acs_filters(dhdu, all=True)
                elif d[0] == 'i':
                    filters = hdr['filter']
                cols['filters'].append(filters)

        # Build output table
        poller_data = [col for col in cols.values()]
        poller_table = Table(data=poller_data,
                             dtype=('str', 'int', 'str', 'str', 'float', 'object', 'str', 'str'))

        # Now assign column names to obset_table
        for i, colname in enumerate(POLLER_COLNAMES):
            poller_table.columns[i].name = colname
    # The input was a poller file, so just keep the viable data rows for output
    else:
        good_rows = []
        for d in usable_datasets:
            for i, old_row in enumerate(input_table):
                if d == input_table['filename'][i]:
                    good_rows.append(old_row)
            poller_table = Table(rows=good_rows, names=input_table.colnames,
                                 dtype=('str', 'int', 'str', 'str', 'float', 'object', 'str', 'str'))

    return poller_table


# ------------------------------------------------------------------------------

def sort_poller_table(obset_table):
    """Sort the input table by photflam and exposure time.

    The alignment of data within a single visit is first done in a relative way,
    successively adding images to the stack within the visit.  The next step is
    to achieve absolute alignment of the stack of images to the GAIA catalog.
    One scheme to obtain the largest number of quality (unsaturated) fiducials/sources,
    is to sort the images first according to their inverse sensitivity (photflam)
    and then according to their exposure time.

    The obset_table has the following columns:
    filename, proposal_id, program_id, obset_id, exptime, filters, detector, and pathname

    Parameters
    ----------
    obset_table: Astropy table
        Astropy table object with columns which map to the original poller file

    Returns
    -------
    updated_obset_table: Astropy table
        The sorted version of the input Astropy table
    """

    # Create a copy of the input table and add the photflam column with a filler value
    expanded_obset_table = Table(obset_table)
    expanded_obset_table['flam'] = -999999.0

    for row in expanded_obset_table:
        input_file = row[expanded_obset_table.colnames[0]]

        # Open the specified FITS file
        h0 = getheader(input_file, 0)

        # Need to get the instrument and detector keywords in order to determine
        # where to look for the various necessary keywords (i.e., primary or
        # extension)
        instrument = h0['instrume'].upper()
        detector = h0['detector'].upper()

        # HST IMAGE
        # photflam: inverse sensitivity, ergs/s-/cm2-/A-1 for 1 electron/s
        # PHOTFLAM keyword is found in each science extension for the currently
        # supported instruments (as of 01 Jan 2020), except for WFC3/IR where
        # PHOTFLAM is in the Primary.
        #
        # Although the PHOTFLAM keyword is science extension-dependent,
        # the differences in values is so small as to not be relevant in
        # this particular context.
        if instrument == 'WFC3' and detector == 'IR':
            photflam = h0['photflam']
        else:
            h1 = getheader(input_file, 'sci', 1)
            photflam = h1['photflam']

        row['flam'] = photflam

    # Determine the rank order the data with a primary key of photflam and a secondary key
    # of exposure time (in seconds).  The primary and secondary keys both need 
    # to be sorted in descending order.  Use the rank to sort
    # the original input table for output.
    # Unfortunately, there are cases where the default works better; namely, 
    # fields with few point sources which are bright that are saturated in the
    # longer exposure time images (like j8cw03 and j8cw53).  
    # rank = np.lexsort((-expanded_obset_table['flam'], -expanded_obset_table['exptime']))
    # Original implementation:
    # rank = np.lexsort((expanded_obset_table['exptime'], -expanded_obset_table['flam']))
    rank = np.lexsort((-expanded_obset_table['flam'], expanded_obset_table['exptime']))
    updated_obset_table = obset_table[rank]

    return updated_obset_table

# ------------------------------------------------------------------------------


def add_primary_fits_header_as_attr(hap_obj, log_level=logutil.logging.NOTSET):
    """create object attribute containing image primary header info

    Parameters
    ----------
    hap_obj : drizzlepac.hlautils.Product.TotalProduct, drizzlepac.hlautils.Product.FilterProduct, or
    drizzlepac.hlautils.Product.ExposureProduct, depending on input
        object to update

    log_level : int, optional.
        The desired level of verboseness in the log statements displayed on the screen and written to the .log file.

    Returns
    -------
    hap_obj : drizzlepac.hlautils.Product.TotalProduct, drizzlepac.hlautils.Product.FilterProduct, or
    drizzlepac.hlautils.Product.ExposureProduct, depending on input
        updated version of input object
    """
    if os.path.exists(hap_obj.drizzle_filename):
        file_name = hap_obj.drizzle_filename
    else:
        file_name = hap_obj.full_filename

    hap_obj.primary_header = fits.getheader(file_name)
    log.info("added primary header info from file {} to {}".format(file_name,hap_obj))

    return hap_obj
