#!/usr/bin/env python

"""Makes .out files used as input to runsinglehap, runmultihap.py based on the files found in the current
working dir"""

import argparse
import glob
import pdb
import os

from astropy.io import fits

from drizzlepac.haputils import poller_utils


def generate_poller_file(poller_file_type='svm', input_list=None, output_poller_filename="poller_file.out", skycell=None):
    # Open rootname list file
    f = open(input_list, 'r')
    rootname_list = f.readlines()
    f.close()
    output_list = []
    for rootname in rootname_list:
        rootname = rootname.strip()
        fullfilepath = ""
        # Find the best available version of the file Assume flc preferred to flt, local preferred to network.
        for fits_ext in ["flc", "flt"]:
            for file_path in [os.getcwd(), os.getenv("DATA_PATH")]:
                if len(fullfilepath) > 0:
                    continue
                if file_path == os.getcwd():
                    fullfilepath_guess = "{}/{}_{}.fits".format(file_path, rootname, fits_ext)
                else:
                    fullfilepath_guess = "{}/{}/{}/{}_{}.fits".format(os.getenv("DATA_PATH"), rootname[:4],
                                                                      rootname, rootname, fits_ext)
                if os.path.exists(fullfilepath_guess):
                    fullfilepath = fullfilepath_guess
        if len(fullfilepath) > 0:
            print("Rootname {}: Found fits file {}".format(rootname, fullfilepath))
            imgname = fullfilepath.split("/")[-1]
        else:
            # Warn user if no fits file can be located for a given rootname, and skip processing of the file.
            print("WARNING: No fits file found for rootname '{}'. This rootname will be omitted from poller file generation.".format(rootname))
            continue

        linelist = []
        linelist.append(imgname)
        imghdu = fits.open(imgname)
        imghdr = imghdu[0].header
        linelist.append("{}".format(imghdr['proposid']))
        linelist.append(imgname[1:4].upper())
        linelist.append(imghdr['linenum'].split(".")[0])
        linelist.append("{}".format(imghdr['exptime']))
        if imghdr['INSTRUME'].lower() == "acs":
            filter = poller_utils.determine_filter_name("{};{}".format(imghdr['FILTER1'], imghdr['FILTER2']))
        elif imghdr['INSTRUME'].lower() == "wfc3":
            filter = poller_utils.determine_filter_name(imghdr['FILTER'])
        linelist.append(filter.upper())
        linelist.append(imghdr['detector'].upper())
        linelist.append(fullfilepath)

        output_list.append(",".join(linelist) + "\n")

        imghdu.close()

    with open(output_poller_filename, 'w') as f:
        f.writelines(output_list)
    print("wrote {}".format(output_poller_filename))

if __name__ == '__main__':
    # Parse input arguments

    parser = argparse.ArgumentParser(description='Create a HAP SVM or MVM poller file')
    parser.add_argument('poller_file_type', choices=['svm', 'mvm'],
                        help='Type of poller file to be created. "smv" to create a poller file for use with '
                             'the single-visit mosaics pipeline and "mvm" to create  a poller file for use '
                             'with the multiple-visit mosaics pipeline. NOTE: if creating a MVM poller file, '
                             'one must specify the skycell name using the "-s" input argument.')
    parser.add_argument('-i', '--input_list', required=True,
                        help='Name of a file containing a list of rootnames (9 characters, usually ending '
                             'with a "q" to process. The corresponding flc.fits or flt.fits files must '
                             'exist either in the current working direcotry on in the online cache') # TODO: VERIFY NAME
    parser.add_argument('-o', '--output_poller_filename', required=False, default="poller_file.out",
                        help='Name of an output poller file that will be created. If not explicitly '
                             'specified, the poller file will be named "poller_file.out".')
    parser.add_argument('-s', '--skycell_name', required=False, default="None",
                        help='Name of the skycell. NOTE: this input argument is *REQUIRED* only for MVM '
                             'poller file creation. ')
    in_args = parser.parse_args()

    # reformat input args
    if in_args.input_list == 'None':
        in_args.input_list = None
    if in_args.skycell_name == 'None':
        in_args.skycell_name = None

    # logic to make sure user has specified the skycell name if a MVM poller file is to be created.
    if in_args.poller_file_type == "mvm" and in_args.skycell_name is None:
        parser.error("ERROR: To create a MVM poller file, a skycell name must be specified with the '-s' argument.")

    generate_poller_file(poller_file_type=in_args.poller_file_type,input_list=in_args.input_list,
                         skycell=in_args.skycell_name)
