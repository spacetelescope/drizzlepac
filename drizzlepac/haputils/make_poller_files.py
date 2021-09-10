#!/usr/bin/env python

"""Makes .out files used as input to runsinglehap.py, runmultihap.py based on the files or rootnames listed
user-specified list file."""

import argparse
import os
import sys

from astropy.io import fits

from drizzlepac.haputils import poller_utils

__taskname__ = 'make_poller_files'


def generate_poller_file(input_list, poller_file_type='svm', output_poller_filename="poller_file.out",
                         skycell_name=None):
    """Creates a properly formatted SVM or MVM poller file.

    Parameters
    ----------
    input_list : str
        Name of the file containing the list of rootnames to process

    poller_file_type : str, optional
        Type of poller file to create. 'svm' for single visit mosaic, 'mvm' for multi-visit mosaic. Default
        value is 'smv'.

    output_poller_filename : str, optional
        Name of the output poller file that will be created. Default value is 'poller_file.out'.

    skycell_name : str, optional
        Name of the skycell to use when creating a MVM poller file. skycell_name is REQUIRED for the creation
        of a MVM poller file, but completely unnecessary for the creation of a SVM poller file. Default value
        is logical 'None'.

    Returns
    -------
    Nothing.
    """
    # Open rootname list file
    f = open(input_list, 'r')
    rootname_list = f.readlines()
    f.close()
    output_list = []
    for rootname in rootname_list:
        rootname = rootname.strip()
        fullfilepath = locate_fitsfile(rootname)
        if len(fullfilepath) > 0:
            print("Rootname {}: Found fits file {}".format(rootname, fullfilepath))
            imgname = fullfilepath.split("/")[-1]
        else:
            # Warn user if no fits file can be located for a given rootname, and skip processing of the file.
            print("WARNING: No fits file found for rootname '{}'. This rootname will be omitted from poller "
                  "file generation.".format(rootname))
            continue
        # Build each individual poller file line
        linelist = []
        linelist.append(imgname)
        imghdu = fits.open(fullfilepath)
        imghdr = imghdu[0].header
        linelist.append("{}".format(imghdr['proposid']))
        linelist.append(imgname.split("_")[-2][1:4].upper())
        linelist.append(imghdr['linenum'].split(".")[0])
        linelist.append("{}".format(imghdr['exptime']))
        if imghdr['INSTRUME'].lower() == "acs":
            filter = poller_utils.determine_filter_name("{};{}".format(imghdr['FILTER1'], imghdr['FILTER2']))
        elif imghdr['INSTRUME'].lower() == "wfc3":
            filter = poller_utils.determine_filter_name(imghdr['FILTER'])
        linelist.append(filter.upper())
        linelist.append(imghdr['detector'].upper())
        if poller_file_type == 'mvm':  # Additional stuff to add to MVM poller files
            if skycell_name.startswith("skycell-"):
                linelist.append("{}".format(skycell_name))
            if skycell_name.startswith("p"):
                linelist.append("skycell-{}".format(skycell_name))
            linelist.append("NEW")
        linelist.append(fullfilepath)
        imghdu.close()
        # Append newly created poller file line to the list of lines to be written to the output file.
        output_list.append(",".join(linelist))
    # adding carriage returns to all but the very last line in the output file.
    list_size = len(output_list)
    for ctr in range(0, list_size):
        if ctr != list_size-1:
            trailing_char = "\n"
        else:
            trailing_char = ""
        output_list[ctr] = output_list[ctr]+trailing_char

    # write output poller file
    with open(output_poller_filename, 'w') as f:
        f.writelines(output_list)
    print("wrote {} poller file '{}'.".format(poller_file_type.upper(), output_poller_filename))


# ============================================================================================================

def locate_fitsfile(search_string):
    """returns full file name (fullpath + filename) for a specified rootname or filename. The search
    algorithm looks for the file in the following order:

    - Search for a _flc.fits file in the current working directory
    - Search for a _flt.fits file in the current working directory
    - Search for a _flc.fits file in subdirectory in the path specified in $DATA_PATH
    - Search for a _flt.fits file in subdirectory in the path specified in $DATA_PATH

    Parameters
    ----------
    search_string : str
        rootname or filename to locate

    Returns
    -------
    fullfilepath : str
         full file path + image name of specified search_string.
    """
    if search_string.endswith("_flt.fits") or search_string.endswith("_flc.fits"):  # Process search_string as a full filename
        # Look for files in CWD first
        if os.path.exists(search_string):
            return os.getcwd()+"/"+search_string
        # If not found in CWD, look elsewhere...
        if not os.getenv("DATA_PATH"):
            sys.exit("ERROR: Undefined online cache data root path. Please set environment variable 'DATA_PATH'")
        fullfilepath = "{}/{}/{}/{}".format(os.getenv("DATA_PATH"), search_string[:4], search_string, search_string)
        if os.path.exists(search_string):
            return fullfilepath
            
    else:  # Process search_string as a rootname
        # Look for files in CWD first
        if os.path.exists("{}_flc.fits".format(search_string)):
            return "{}/{}_flc.fits".format(os.getcwd(), search_string)
        if os.path.exists("{}_flt.fits".format(search_string)):
            return "{}/{}_flt.fits".format(os.getcwd(), search_string)
        # If not found in CWD, look elsewhere...
        if not os.getenv("DATA_PATH"):
            sys.exit("ERROR: Undefined online cache data root path. Please set environment variable 'DATA_PATH'")
        filenamestub = "{}/{}/{}/{}".format(os.getenv("DATA_PATH"), search_string[:4], search_string, search_string)
        if os.path.exists("{}_flc.fits".format(filenamestub)):
            fullfilepath = "{}_flc.fits".format(filenamestub)
        else:
            fullfilepath = "{}_flt.fits".format(filenamestub)
        return fullfilepath


# ============================================================================================================

if __name__ == '__main__':
    # Parse input arguments
    parser = argparse.ArgumentParser(description='Create a HAP SVM or MVM poller file')

    parser.add_argument('input_list',
                        help='Name of a file containing a list of calibrated fits files (ending with '
                             '"_flt.fits" or "_flc.fits") or rootnames (9 characters, usually ending '
                             'with a "q" to process. The corresponding flc.fits or flt.fits files must '
                             'exist in the online cache')
    parser.add_argument('-o', '--output_poller_filename', required=False, default="poller_file.out",
                        help='Name of an output poller file that will be created. If not explicitly '
                             'specified, the poller file will be named "poller_file.out".')
    parser.add_argument('-s', '--skycell_name', required=False, default="None",
                        help='Name of the skycell. NOTE: this input argument is *REQUIRED* only for MVM '
                             'poller file creation. ')
    parser.add_argument('-t', '--poller_file_type', required=False, choices=['svm', 'mvm'], default='svm',
                        help='Type of poller file to be created. "smv" to create a poller file for use with '
                             'the single-visit mosaics pipeline and "mvm" to create a poller file for use '
                             'with the multiple-visit mosaics pipeline. If not explicitly '
                             'specified, the default value is "svm". NOTE: if creating a MVM poller file, '
                             'one must specify the skycell name using the "-s" input argument.')
    in_args = parser.parse_args()

    # reformat input args
    if in_args.skycell_name == 'None':
        in_args.skycell_name = None

    # logic to make sure user has specified the skycell name if a MVM poller file is to be created.
    if in_args.poller_file_type == "mvm" and in_args.skycell_name is None:
        parser.error("ERROR: To create a MVM poller file, a skycell name must be specified with the '-s' argument.")

    generate_poller_file(in_args.input_list,
                         poller_file_type=in_args.poller_file_type,
                         output_poller_filename=in_args.output_poller_filename,
                         skycell_name=in_args.skycell_name)
