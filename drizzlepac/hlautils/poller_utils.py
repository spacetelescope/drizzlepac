"""Utilities to interpret the pipeline poller obset information and generate product filenames

The function, interpret_obset_input, parses the file generated by the pipeline
poller, and produces a tree listing of the output products.  The function,
parse_obset_tree, converts the tree into product catagories.

"""
import os
import sys
from collections import OrderedDict
import numpy as np

from stsci.tools import logutil


from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table, Column
from drizzlepac.hlautils.product import ExposureProduct, FilterProduct, TotalProduct
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

__taskname__ = 'poller_utils'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout, 
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)

def interpret_obset_input(results):
    """

    Parameters
    -----------
    results : str or list
        Input poller file name or Python list of dataset names to be processed as a single visit.
        Dataset names have to be either the filename of a singleton (non-associated exposure) or the
        name of an ASN file (e.g., jabc01010_asn.fits).

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
    log.info("Interpret the poller file for the observation set.")
    obset_table = build_poller_table(results)
    print(obset_table['obset_id'])
    # Add INSTRUMENT column
    instr = INSTRUMENT_DICT[obset_table['filename'][0][0]]
    # convert input to an Astropy Table for parsing
    obset_table.add_column(Column([instr] * len(obset_table)), name='instrument')
    # parse Table into a tree-like dict
    log.info("Build the observation set tree.")
    obset_tree = build_obset_tree(obset_table)
    # Now create the output product objects
    log.info("Parse the observation set tree and create the exposure, filter, and total detection objects.")
    obset_dict, tdp_list = parse_obset_tree(obset_tree)

    # This little bit of code adds an attribute to single exposure objects that is True if a given filter only contains
    # one input (e.g. n_exp = 1)
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

def parse_obset_tree(det_tree):
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
                prod_info = (filename[0] + " " + filetype).lower()

                # Set up the single exposure product dictionary
                sep = SEP_STR.format(sep_indx)
                obset_products[sep] = {'info': prod_info,
                                       'files': [filename[1]]}

                # Create a single exposure product object
                prod_list = prod_info.split(" ")
                sep_obj = ExposureProduct(prod_list[0], prod_list[1], prod_list[2], prod_list[3], filename[1], prod_list[5], prod_list[6])
                # Set up the filter product dictionary and create a filter product object
                # Initialize `info` key for this filter product dictionary
                if not obset_products[fprod]['info']:
                    obset_products[fprod]['info'] = prod_info

                    # Create a filter product object for this instrument/detector
                    filt_obj = FilterProduct(prod_list[0], prod_list[1], prod_list[2], prod_list[3], prod_list[4], prod_list[5], prod_list[6])
                # Append exposure object to the list of exposure objects for this specific filter product object
                filt_obj.add_member(sep_obj)
                # Populate filter product dictionary with input filename
                obset_products[fprod]['files'].append(filename[1])

                # Set up the total detection product dictionary and create a total detection product object
                # Initialize `info` key for total detection product
                if not obset_products[totprod]['info']:
                    obset_products[totprod]['info'] = prod_info

                    # Create a total detection product object for this instrument/detector
                    tdp_obj = TotalProduct(prod_list[0], prod_list[1], prod_list[2], prod_list[3], prod_list[4], prod_list[6])

                # Append exposure object to the list of exposure objects for this specific total detection product
                tdp_obj.add_member(sep_obj)
                # Populate total detection product dictionary with input filename
                obset_products[totprod]['files'].append(filename[1])

                # Increment single exposure index
                sep_indx += 1

            # Append filter object to the list of filter objects for this specific total product object
            log.info("Attach the filter object {} to its associated total detection product object {}/{}.".format(filt_obj.filters, tdp_obj.instrument, tdp_obj.detector))
            tdp_obj.add_product(filt_obj)

        # Add the total product object to the list of TotalProducts
        tdp_list.append(tdp_obj)

    # Done... return dict and object product list
    return obset_products, tdp_list

# ----------------------------------------------------------------------------------------------------------

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
    - NOTE: There should always be at least one filter name provided to
      this routine or this input is invalid.
    """

    raw_filter = raw_filter.lower()

    # There might be two filters, so split the filter names into a list
    filter_list = raw_filter.split(';')
    output_filter_list = []

    for filt in filter_list:
        # Get the names of the non-clear filters
        if 'clear' not in filt:
            output_filter_list.append(filt)

    if not output_filter_list:
        filter_name = 'clear'
    else:
        if output_filter_list[0].startswith('pol'):
            output_filter_list.reverse()

        delimiter = '-'
        filter_name = delimiter.join(output_filter_list)

    return filter_name

# ----------------------------------------------------------------------------------------------------------

def build_poller_table(input):
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
    obs_converters = {'col4': [ascii.convert_numpy(np.str)]}
    if isinstance(input, str):
        input = ascii.read(input, format='no_header', converters=obs_converters)
        if len(input.columns) == len(POLLER_COLNAMES):
            # We were provided a poller file, so use as-is
            # Now assign column names to obset_table
            for i, colname in enumerate(POLLER_COLNAMES):
                input.columns[i].name = colname

            # Convert to a string column, instead of int64
            input['obset_id'] = input['obset_id'].astype(np.str)
            return input

        # Return first column
        filenames = input[input.colnames[0]].tolist()
    elif isinstance(input, list):
        filenames = input

    else:
        id = '[poller_utils.build_poller_table] '
        log.error("{}: Input {} not supported as input for processing.".format(id, input))
        raise ValueError

    # At this point, we should have a list of dataset names
    datasets = []  # final set of filenames locally on disk for processing
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

    cols = OrderedDict()
    for cname in POLLER_COLNAMES:
        cols[cname] = []
    cols['filename'] = datasets

    # Now, evaluate each input dataset for the information needed for the poller file
    for d in datasets:
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
                filters = dhdu['filter']
            cols['filters'].append(filters)
    # Build output table
    poller_data = [col for col in cols.values()]
    poller_table = Table(data=poller_data)

    # Now assign column names to obset_table
    for i, colname in enumerate(POLLER_COLNAMES):
        poller_table.columns[i].name = colname

    return poller_table
