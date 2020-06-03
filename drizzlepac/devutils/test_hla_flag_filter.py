#!/usr/bin/env python
"""script to run hla_flag_filter and compare_sourcelists to test hla_flag_filter parameter updates"""
import argparse
from datetime import datetime
import glob
import json
import logging
import os
import pdb
import pickle
import shutil
import sys

from stsci.tools import logutil

from drizzlepac import hapsequencer
from drizzlepac.devutils.comparison_tools import compare_sourcelists
from drizzlepac.haputils import hla_flag_filter


__taskname__ = 'test_hla_flag_filter'
MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)

# ----------------------------------------------------------------------------------------------------------------------

def hff_parameter_manager(hff_inputs,qc_json_filename):
    """displays and updates (as needed) hla_flag_filter.py parameters stored in hff_inputs['param_dict'].
    NOTE: For this to work, the hff_inputs quality control parameters and the parameters from the user-specified QC
    param file MUST have the exact same overall structure.

    Parameters
    ----------
    hff_inputs : dict
        dictionary containing all necessary inputs to run hla_flag_filter.run_source_list_flagging().

    qc_json_filename : str
        Name of the quality_control .json file that will provide the parameters for the hla_flag_filter run.

    Returns
    -------
    hff_inputs : dict
        if qc_json_filename was specified, an updated version of hff_inputs will be returned with the updated quality
        control parameters. Otherwise, hff_inputs will be returned unchanged
    """
    instrument = hff_inputs['drizzled_image'].split("_")[3]
    detector = hff_inputs['drizzled_image'].split("_")[4]
    full_pars_path = "{}pars/hap_pars/default_parameters/{}/{}/{}_{}_quality_control_all.json".format(os.path.realpath(__file__).split("devutils")[0],instrument,detector,instrument,detector)
    local_pars_path = os.path.basename(full_pars_path).replace(".json","_default.json")
    if not os.path.exists(local_pars_path):
        shutil.copy(full_pars_path,local_pars_path)
        log.info("Created local copy of default quality control .json file {}.".format(local_pars_path))
    hff_params = hff_inputs['param_dict']['quality control']
    if qc_json_filename:
        extra_text_string = " and parameters updates from {}".format(qc_json_filename)
        with open(qc_json_filename) as f_cfg:
            new_params = json.load(f_cfg)
    else:
        extra_text_string = ""
        new_params = hff_params
    log.info("Summary of original hla_flag_filter parameters{}".format(extra_text_string))
    if qc_json_filename:
        log.info("NOTE: updated parameters listed with a exclamation point")
    resursive_print_all_nested_dict_values(hff_params,new_params)

    if qc_json_filename:
        hff_inputs['param_dict']['quality control'] = new_params
    return hff_inputs

# ----------------------------------------------------------------------------------------------------------------------

def preserve_orig_files(hff_inputs,source_path,dest_path,verbose):
    """either move original files already produced by the original run of hla_flag_filter.py to a safe place
    so they don't get clobbered by the new run, or move new files to final directory and move orig. files back

    """
    if source_path == os.getcwd():
        files_to_move = [hff_inputs['catalog_name']]

        found_files = glob.glob("{}*.txt".format(hff_inputs['catalog_name'].replace(".ecsv","")))
        if len(found_files) > 0:
            files_to_move = files_to_move + found_files

        if hff_inputs['catalog_name'].endswith("point-cat.ecsv"):
            hla_classic_phot_type = "dao"
            files_to_move.append(hff_inputs['catalog_name'].replace("_point-cat.ecsv","_matched_sources_only_point-cat.ecsv"))
        else:
            hla_classic_phot_type = "sex"
            files_to_move.append(
            hff_inputs['catalog_name'].replace("_segment-cat.ecsv", "_matched_sources_only_point-cat.ecsv"))
        files_to_move.append("{}_{}phot_corrected.txt".format(hff_inputs['drizzled_image'][:-16],hla_classic_phot_type))
        files_to_move.append("{}_matched_sources_only_{}phot_corrected.txt".format(hff_inputs['drizzled_image'][:-16],hla_classic_phot_type))

        for item in ["_INTERMEDIATE.txt","_ALL_FLT_SAT_FLAG_PIX.txt"]:
            files_to_move.append("{}{}".format( hff_inputs['drizzled_image'][:-5],item))

        for fitsname in hff_inputs['flt_list']:
            found_files = glob.glob("{}_sci?.txt".format(fitsname.replace(".fits","")))
            if len(found_files) > 0:
                files_to_move = files_to_move + found_files

        found_files = glob.glob("*_comparision_plots.pdf")
        if len(found_files) > 0:
            files_to_move = files_to_move + found_files
        if not os.path.exists(dest_path):
            os.mkdir(dest_path)
        for item in files_to_move:
            full_source_path = os.path.join(source_path,item)
            full_dest_path = os.path.join(dest_path,item)
            if os.path.exists(item):
                if verbose:
                    log.info("Move {} \u21E8  {}".format(item,full_dest_path))
                shutil.move(item,full_dest_path) # TODO: change from copy to move
    else:
        for file_name in glob.glob(os.path.join(source_path,"*")):
            full_dest_path = os.path.join(dest_path,file_name.split("/")[1])
            if os.path.exists(file_name):
                if verbose:
                    log.info("Move {} \u21E8  {}".format(file_name,full_dest_path))
                shutil.move(file_name,full_dest_path) # TODO: change from copy to move
        log.info("Removing temp dir {}".format(source_path))
        os.rmdir(source_path)



# ----------------------------------------------------------------------------------------------------------------------

def resursive_print_all_nested_dict_values(old_dict,new_dict,recursion_level=0,recursion_limit=20):
    """recursively print all elemnets of dictonary and highlight any changes.
    NOTE: This is a recursive subroutine.

    Parameters
    ----------
    old_dict : dictionary
        origional dictionary whose elements will be printed and compared to corresponding new_dict elements

    new_dict : dictionary
        new dictioanry whose elements will be printed and compared to corresponding old_dict elements

    recursion_level : int, optional
        current recursive depth. if not explicitly specified, the default value is 0.

    recursion_limit : int, optional
        maximum allowed recursive depth. This here to prevent unexpected recursive runaway. if not explicitly
        specified, the default value is 20.

    Returns
    -------
    Nothing.
    """
    if recursion_level == recursion_limit:
        sys.exit("RECURSION LIMIT REACHED!")
    sorted_key_list = list(old_dict.keys())
    for item in sorted(sorted_key_list):
        if isinstance(old_dict[item], dict):
            log.info("  {}{}\u2798".format("     "*recursion_level,item))
            recursion_level+=1
            resursive_print_all_nested_dict_values(old_dict[item],new_dict[item],recursion_level=recursion_level)
            recursion_level-=1
        else:
            if old_dict[item] == new_dict[item]:
                log.info("  {}{}: {}".format("     "*recursion_level,item,new_dict[item]))
            else:
                log.info("! {}{}: {} \u21E8  {}".format("     "*recursion_level,item,old_dict[item],new_dict[item]))

# ----------------------------------------------------------------------------------------------------------------------

def run_compare_sourcelists(hff_inputs, log_level):
    """locate HLA classic image and sourcelists, convert HLA classic sorucelist X, Y, RA and Dec to HAP ref frame for
    apples-to-apples comparision, and run comparision code

    Parameters
    ----------
    hff_inputs : dict
        dictionary containing all necessary inputs to run hla_flag_filter.run_source_list_flagging().

    log_level : int
        desired log level

    Returns
    -------
    Nothing.
    """
    # get HLA classic path details from environment variables
    hla_classic_basepath = os.getenv('HLA_CLASSIC_BASEPATH')
    hla_build_ver = os.getenv("HLA_BUILD_VER")
    parse_imgname = hff_inputs['drizzled_image'].split("_")
    instrument = parse_imgname[3]
    prop_id = parse_imgname[1]
    obset_id = parse_imgname[2]
    filter = parse_imgname[5]
    hla_product_basename = "_".join(parse_imgname[0:6])
    hap_product_basename = "_".join(parse_imgname[0:7])
    if hla_classic_basepath and hla_build_ver and os.path.exists(hla_classic_basepath):
        hla_cassic_basepath = os.path.join(hla_classic_basepath, instrument, hla_build_ver)
        hla_classic_path = os.path.join(hla_cassic_basepath, prop_id, prop_id + "_" + obset_id)  # Generate path to HLA classic products
    elif os.path.exists(os.path.join(os.getcwd(), "hla_classic")):  # For local testing
        hla_classic_basepath = os.path.join(os.getcwd(), "hla_classic")
        hla_classic_path = hla_classic_basepath
    else:
        return  # bail out if HLA classic path can't be found.


    hap_imgname = hff_inputs['drizzled_image']
    hla_imgname = glob.glob("{}/{}_dr*.fits".format(hla_classic_path,hla_product_basename))[0]

    if not os.path.exists(hap_imgname) or not os.path.exists(hla_imgname):  # bail if one or both of the images can't be found
        sys.exit("Error: one or both of the images can't be found")
    hap_sourcelist_name = hff_inputs['catalog_name']
    if hap_sourcelist_name.endswith("point-cat.ecsv"):
        hla_classic_cat_type = "dao"
        plotfile_prefix = hap_product_basename + "_point"
    else:
        hla_classic_cat_type = "sex"
        plotfile_prefix = hap_product_basename + "_segment"

    if hla_classic_basepath and hla_build_ver and os.path.exists(hla_classic_basepath):
        hla_sourcelist_name = "{}/logs/{}_{}phot.txt".format(hla_classic_path,
                                                             hla_product_basename,
                                                             hla_classic_cat_type)
    else:
        hla_sourcelist_name = "{}/{}_{}phot.txt".format(hla_classic_path,
                                                        hla_product_basename,
                                                        hla_classic_cat_type)
    plotfile_prefix += "_testing"
    if not os.path.exists(hap_sourcelist_name) or not os.path.exists(hla_sourcelist_name):  # Skip catalog type if one or both of the catalogs can't be found
        sys.exit() # TODO: Circle back to this

    # convert HLA Classic RA and Dec values to HAP reference frame so the RA and Dec comparisons are correct
    updated_hla_sourcelist_name = hapsequencer.correct_hla_classic_ra_dec(hla_sourcelist_name, hap_imgname,
                                                                          hla_classic_cat_type, log_level)
    log.info("HAP image:                   {}".format(os.path.basename(hap_imgname)))
    log.info("HLA Classic image:           {}".format(os.path.basename(hla_imgname)))
    log.info("HAP catalog:                 {}".format(os.path.basename(hap_sourcelist_name)))
    log.info("HLA Classic catalog:         {}".format(os.path.basename(updated_hla_sourcelist_name)))

    # once all file exist checks are passed, execute sourcelist comparision

    return_status = compare_sourcelists.comparesourcelists([updated_hla_sourcelist_name, hap_sourcelist_name],
                                                           [hla_imgname, hap_imgname], good_flag_sum=255,
                                                           plotGen="file", plotfile_prefix=plotfile_prefix,
                                                           verbose=True, log_level=log_level, debugMode=True)

# ----------------------------------------------------------------------------------------------------------------------

def run_hla_flag_filter(hff_inputs):
    """execute hla_flag_filter.run_source_list_flagging() and write out the catalog to file

    Parameters
    ----------
    hff_inputs : dict
        dictionary containing all necessary inputs to run hla_flag_filter.run_source_list_flagging().

    Returns
    -------
    Nothing.
    """
    source_list = hla_flag_filter.run_source_list_flagging(hff_inputs["drizzled_image"],
                                                           hff_inputs["flt_list"],
                                                           hff_inputs["param_dict"],
                                                           hff_inputs["exptime"],
                                                           hff_inputs["plate_scale"],
                                                           hff_inputs["median_sky"],
                                                           hff_inputs["catalog_name"],
                                                           hff_inputs["catalog_data"],
                                                           hff_inputs["cat_type"],
                                                           hff_inputs["drz_root_dir"],
                                                           hff_inputs["hla_flag_msk"],
                                                           hff_inputs["ci_lookup_file_path"],
                                                           hff_inputs["output_custom_pars_file"],
                                                           hff_inputs["log_level"],
                                                           True)
    # log.info("renamed existing catalog {} -> {}".format(hff_inputs["catalog_name"],
    #                                                     hff_inputs["catalog_name"].replace(".ecsv","_old.ecsv")))
    # os.rename(hff_inputs["catalog_name"],hff_inputs["catalog_name"].replace(".ecsv","_old.ecsv"))
    source_list.write(hff_inputs['catalog_name'], format='ascii.ecsv')
    log.info('Wrote new catalog {}.'.format(hff_inputs["catalog_name"]))

# ----------------------------------------------------------------------------------------------------------------------


def run_stuff(input_pickle_filename, qc_json_filename, verbose):
    timestamp = "{}".format(datetime.now().strftime("D%m_%d_%YT%H_%M_%S"))
    log_level = 10
    log.setLevel(log_level)
    logname = "test_hla_flag_filter_{}.log".format(timestamp)
    logging.basicConfig(filename=logname, format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)
    timestamp = "{}".format(datetime.now().strftime("D%m_%d_%YT%H_%M_%S"))
    hff_inputs = pickle.load(open(input_pickle_filename, "rb"))

    hff_inputs = hff_parameter_manager(hff_inputs,qc_json_filename)

    temp_dirname = "tempdir_{}_orig_hla_flag_files".format(timestamp)
    final_dir = "{}_hla_flag_files".format(timestamp)
    log.info("Moving original files generated by hla_flag_filter.py to temp dir "+temp_dirname)
    preserve_orig_files(hff_inputs,os.getcwd(), temp_dirname, verbose)


    run_hla_flag_filter(hff_inputs)

    run_compare_sourcelists(hff_inputs, log_level)
    log.info("Moving new files generated by hla_flag_filter.py to final path " + final_dir)
    preserve_orig_files(hff_inputs,os.getcwd(),final_dir,verbose)
    if qc_json_filename:
        shutil.copy(qc_json_filename,os.path.join(final_dir,qc_json_filename))

    log.info("Restoring original files from temp dir {} to CWD.".format(temp_dirname))
    preserve_orig_files(hff_inputs,temp_dirname,os.getcwd(),verbose)
    logging.shutdown()
    shutil.move(logname,os.path.join(final_dir,logname))


# =======================================================================================================================
if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description='Test hla_flag_filter.py')
    PARSER.add_argument('input_pickle_filename',nargs=1,help="Name of pickle file containing hla_flag_filter.py input parameters")
    PARSER.add_argument('-p', '--qc_json_filename',required=False,default=None,help="Name of quality control .json file to use instead for hla_flag_filter params")
    PARSER.add_argument('-v', '--verbose',required=False, action='store_true',help="Print detailed info about file moves? Default value is 'False'.")

    ARGS = PARSER.parse_args()
    run_stuff(ARGS.input_pickle_filename[0],ARGS.qc_json_filename,ARGS.verbose)

    # TODO: get logging working so this script produces a log file!