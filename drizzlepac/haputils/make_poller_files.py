#!/usr/bin/env python

"""Generates a poller file that will be used as input to runsinglehap.py, hapsequencer.py, runmultihap.py or
hapmultisequencer.py based on the files or rootnames listed user-specified list file.

USAGE
    >>> python drizzlepac/haputils/make_poller_files.py <input filename> -[ost]
    - input filename: Name of a file containing a list of calibrated fits files (ending with "_flt.fits" or
      "_flc.fits") or rootnames (9 characters, usually ending with a "q" to process. The corresponding
      flc.fits or flt.fits files must exist in the user-specified path, the current working directory or the
      online cache

    - The '-o' optional input allows users to input the name of an output poller file that will be created.
      If not explicitly specified, the poller file will be named "poller_file.out".

    - The '-s' optional input allows users to input the Name of the skycell. The correct syntax for skycell
      names is "skycell-pNNNNxXXyXX", where NNNN is the 4-digit projection cell number, and XX and YY are the
      two-digit X and Y skycell indices, respectively. NOTE: this input argument is not needed for SVM poller
      file creation, but *REQUIRED* for MVM poller file creation. Users can determine the skycell(s) that
      their observations occupy using the ``haputils.which_skycell`` script.

    - The '-t' optional input allows users to specify the type of poller file that will be created. The
      valid input options are "svm" to create a poller file for use with the single-visit mosaics pipeline
      or "mvm" to create a poller file for use with the multiple-visit mosaics pipeline. If not explicitly
      specified, the default value is "svm". NOTE: if creating a MVM poller file, one must specify the
      skycell name using the "-s" input argument.

Python USAGE:
    >>> python
    >>> from drizzlepac.haputils import make_poller_files
    >>> make_poller_files.generate_poller_file(input_list, poller_file_type='svm', output_poller_filename="poller_file.out", skycell_name=None):
"""

import argparse
import os
import re
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
        Name of the text file containing the list of filenames or rootnames to process

    poller_file_type : str, optional
        Type of poller file to create. 'svm' for single visit mosaic, 'mvm' for multi-visit mosaic. Default
        value is 'svm'.

    output_poller_filename : str, optional
        Name of the output poller file that will be created. Default value is 'poller_file.out'.

    skycell_name : str, optional
        Name of the skycell to use when creating a MVM poller file. skycell_name is REQUIRED for the creation
        of a MVM poller file, but completely unnecessary for the creation of a SVM poller file. The correct
        syntax for skycell names is 'skycell-pNNNNxXXyXX', where NNNN is the 4-digit projection cell number,
        and XX and YY are the two-digit X and Y skycell indices, respectively. Default value is logical
        'None'. NOTE: this input argument is not needed for SVM poller file creation, but *REQUIRED* for MVM
        poller file creation. Users can determine the skycell(s) that their observations occupy using the
        ``haputils.which_skycell`` script.

    Returns
    -------
    Nothing.
    """
    if poller_file_type == 'svm' and skycell_name:
        print("PROTIP: Users only need to provide a skycell name for the creation of MVM poller files, not SVM poller files.")
    # Open rootname list file
    f = open(input_list, 'r')
    rootname_list = f.readlines()
    f.close()
    output_list = []
    for rootname in rootname_list:
        rootname = rootname.strip()
        fullfilepath = locate_fitsfile(rootname)
        if len(fullfilepath) > 0:
            if rootname.endswith(".fits"):
                print("Found fits file {}".format(fullfilepath))
            else:
                print("Rootname {}: Found fits file {}".format(rootname, fullfilepath))
            imgname = fullfilepath.split(os.sep)[-1]
        else:
            # Warn user if no fits file can be located for a given rootname, and skip processing of the file.
            if rootname.endswith(".fits"):
                item_type = "filename"
            else:
                item_type = "rootname"
            print("WARNING: No fits file found for {} '{}'. This {} will be omitted from the poller file.".format(item_type, rootname, item_type))
            continue
        # Build each individual poller file line
        linelist = []
        linelist.append(imgname)
        imghdu = fits.open(fullfilepath)
        imghdr = imghdu[0].header
        linelist.append("{}".format(imghdr['proposid']))
        linelist.append(imgname.split("_")[-2][1:4].upper())
        if imghdr['primesi'].lower() == imghdr['instrume']:
            linelist.append(imghdr['linenum'].split(".")[0])
        else:
            linelist.append(imghdr['rootname'][-5:-3].upper())
        linelist.append("{}".format(imghdr['exptime']))
        if imghdr['INSTRUME'].lower() == "acs":
            filter = poller_utils.determine_filter_name("{};{}".format(imghdr['FILTER1'], imghdr['FILTER2']))
        elif imghdr['INSTRUME'].lower() == "wfc3":
            filter = poller_utils.determine_filter_name(imghdr['FILTER'])
        linelist.append(filter.upper())
        linelist.append(imghdr['detector'].upper())
        if poller_file_type == 'mvm':  # Additional stuff to add to MVM poller files
            if skycell_name:
                pattern = re.compile("(skycell-p\d{4}x\d{2}y\d{2})")
                skycell_name_format_check = pattern.match(skycell_name)
                if skycell_name_format_check:
                    linelist.append("{}".format(skycell_name))
                else:
                    raise ValueError("'{}' is an improperly formatted skycell name. Please refer to documentation for information regarding correct skycell name syntax.".format(skycell_name))
            else:
                raise Exception("No skycell name was provided. The name of the skycell that the observations occupy is required for MVM poller file creation.")
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
        # Look in user-provided path (assuming they provided one)
        if os.path.exists(search_string) and os.sep in search_string:
            return search_string
        # Look for files in CWD
        if os.path.exists(search_string) and os.sep not in search_string:
            return os.getcwd() + os.sep + search_string
        # If not found in CWD, look elsewhere...
        if not os.getenv("DATA_PATH"):
            sys.exit("ERROR: Undefined online cache data root path. Please set environment variable 'DATA_PATH'")
        fullfilepath = "{}{}{}{}{}{}{}".format(os.getenv("DATA_PATH"), os.sep, search_string[:4],
                                               os.sep, search_string[:-9], os.sep, search_string)
        if os.path.exists(search_string):
            return fullfilepath
        else:
            return ""  # Return a null string if no file is found
            
    else:  # Process search_string as a rootname
        # Look for files in CWD first
        for fits_ext in ["flc", "flt"]:
            if os.path.exists("{}_{}.fits".format(search_string, fits_ext)):
                return "{}{}{}_{}.fits".format(os.getcwd(), os.sep, search_string, fits_ext)
        # If not found in CWD, look elsewhere...
        if not os.getenv("DATA_PATH"):
            sys.exit("ERROR: Undefined online cache data root path. Please set environment variable 'DATA_PATH'")
        filenamestub = "{}{}{}{}{}{}{}".format(os.getenv("DATA_PATH"), os.sep, search_string[:4],
                                               os.sep, search_string, os.sep, search_string)
        for fits_ext in ["flc", "flt"]:
            if os.path.exists("{}_{}.fits".format(filenamestub, fits_ext)):
                return "{}_{}.fits".format(filenamestub, fits_ext)
        # it should never get here unless no file was found either locally or elsewhere in $DATA_PATH.
        return ""  # Return a null string if no file is found

# ============================================================================================================


if __name__ == '__main__':
    # Parse input arguments
    parser = argparse.ArgumentParser(description='Create a HAP SVM or MVM poller file')

    parser.add_argument('input_list',
                        help='Name of a file containing a list of calibrated fits files (ending with '
                             '"_flt.fits" or "_flc.fits") or rootnames (9 characters, usually ending '
                             'with a "q" to process. The corresponding flc.fits or flt.fits files must '
                             'exist in the user-specified path, the current working directory or the online '
                             'cache')
    parser.add_argument('-o', '--output_poller_filename', required=False, default="poller_file.out",
                        help='Name of an output poller file that will be created. If not explicitly '
                             'specified, the poller file will be named "poller_file.out".')
    parser.add_argument('-s', '--skycell_name', required=False, default="None",
                        help='Name of the skycell. The correct syntax for skycell names is '
                             '"skycell-pNNNNxXXyXX", where NNNN is the 4-digit projection cell number, and '
                             'XX and YY are the two-digit X and Y skycell indices, respectively. NOTE: this '
                             'input argument is not needed for SVM poller file creation, but *REQUIRED* for '
                             'MVM poller file creation. Users can determine the skycell(s) that their '
                             'observations occupy using the haputils.which_skycell.py script.')
    parser.add_argument('-t', '--poller_file_type', required=False, choices=['svm', 'mvm'], default='svm',
                        help='Type of poller file to be created. "svm" to create a poller file for use with '
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
