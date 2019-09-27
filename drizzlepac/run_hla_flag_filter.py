#!/usr/bin/env python
"""This script simply calls drizzlepac/hlautils/hla_flag_filter.py for test purposes"""
import json
import os
import pdb

from astropy.table import Table
import drizzlepac
from drizzlepac.hlautils import hla_flag_filter
from drizzlepac.hlautils import config_utils
from drizzlepac.hlautils import poller_utils

def run_hla_flag_filter():






    #   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +
    # All below lines are to get it working, not actual final code.
    # Get parameter values
    for cmd in ['rm -f *.*', 'cp orig/* .']:
        print(cmd)
        os.system(cmd)


    obs_info_dict, total_list = poller_utils.interpret_obset_input("j92c01.out")

    out_pars_file = None
    for total_item in total_list:

        total_item.configobj_pars = config_utils.HapConfig(total_item, output_custom_pars_file=out_pars_file,
                                                           use_defaults=True)
        for filter_item in total_item.fdp_list:
            filter_item.configobj_pars = config_utils.HapConfig(filter_item, output_custom_pars_file=out_pars_file,
                                                                use_defaults=True)
        for expo_item in total_item.edp_list:
            expo_item.configobj_pars = config_utils.HapConfig(expo_item, output_custom_pars_file=out_pars_file,
                                                              use_defaults=True)


    # * * * * hla_flag_filter.run_source_list_flagging inputs for HLA Classic test run* * * *
    # mode = "dao"
    # all_drizzled_filelist = ["hst_10265_01_acs_wfc_f606w_drz.fits"]
    # working_hla_red = os.getcwd()
    # filter_sorted_flt_dict = {"f606w": ["j92c01b4q_flc.fits", "j92c01b5q_flc.fits", "j92c01b7q_flc.fits", "j92c01b9q_flc.fits"]}
    # param_dict = total_list[0].fdp_list[0].configobj_pars.as_single_giant_dict()
    # param_dict['catalog generation']['dao']['bthresh'] = 5.0 #force it to use the value from HLA classic
    # readnoise_dictionary_drzs = {"hst_10265_01_acs_wfc_f606w_drz.fits": 4.97749985}
    # scale_dict_drzs = {"hst_10265_01_acs_wfc_f606w_drz.fits": 0.05}
    # zero_point_AB_dict = {"hst_10265_01_acs_wfc_f606w_drz.fits": 26.5136022236}
    # exp_dictionary_scis = {"hst_10265_01_acs_wfc_f606w_drz.fits": 5060.0}
    # detection_image = "hst_10265_01_acs_wfc_total_drz.fits"
    # dict_newTAB_matched2drz = {"hst_10265_01_acs_wfc_f606w_drz.fits": "hst_10265_01_acs_wfc_f606w_{}phot.txt".format(mode)}
    # phot_table_matched2cat = {all_drizzled_filelist[0]: Table.read(dict_newTAB_matched2drz[all_drizzled_filelist[0]], format='ascii')}
    # proc_type = "{}phot".format(mode)
    # drz_root_dir = os.getcwd()
    # rms_dict = {"hst_10265_01_acs_wfc_f606w_drz.fits": "hst_10265_01_acs_wfc_f606w_rms.fits"}


    mode = "dao"
    all_drizzled_filelist = ["hst_10265_01_acs_wfc_f606w_drz.fits"]
    working_hla_red = os.getcwd()
    filter_sorted_flt_dict = {"f606w": ["j92c01b4q_flc.fits", "j92c01b5q_flc.fits", "j92c01b7q_flc.fits", "j92c01b9q_flc.fits"]}
    param_dict = total_list[0].fdp_list[0].configobj_pars.as_single_giant_dict()
    param_dict['catalog generation']['dao']['bthresh'] = 5.0 #force it to use the value from HLA classic
    readnoise_dictionary_drzs = {"hst_10265_01_acs_wfc_f606w_drz.fits": 4.97749985}
    scale_dict_drzs = {"hst_10265_01_acs_wfc_f606w_drz.fits": 0.05}
    zero_point_AB_dict = {"hst_10265_01_acs_wfc_f606w_drz.fits": 26.5136022236}
    exp_dictionary_scis = {"hst_10265_01_acs_wfc_f606w_drz.fits": 5060.0}
    detection_image = "hst_10265_01_acs_wfc_total_drz.fits"
    dict_newTAB_matched2drz = {"hst_10265_01_acs_wfc_f606w_drz.fits": "hst_10265_01_acs_wfc_f606w_{}phot.txt".format(mode)}
    phot_table_matched2cat = {all_drizzled_filelist[0]: Table.read(dict_newTAB_matched2drz[all_drizzled_filelist[0]], format='ascii')}
    proc_type = "{}phot".format(mode)
    drz_root_dir = os.getcwd()
    rms_dict = {"hst_10265_01_acs_wfc_f606w_drz.fits": "hst_10265_01_acs_wfc_f606w_rms.fits"}
    #   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +
    # Execute hla_flag_filter.run_source_list_flaging
    hla_flag_filter.run_source_list_flaging(all_drizzled_filelist, working_hla_red, filter_sorted_flt_dict,
                                            param_dict, readnoise_dictionary_drzs,
                                            scale_dict_drzs, zero_point_AB_dict, exp_dictionary_scis,
                                            detection_image, dict_newTAB_matched2drz, phot_table_matched2cat,
                                            proc_type, drz_root_dir, rms_dict)

if __name__ == "__main__":
    run_hla_flag_filter()
