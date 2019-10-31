#!/usr/bin/env python
"""This script simply calls drizzlepac/hlautils/hla_flag_filter.py for test purposes"""
import json
import glob
import os
import pdb
import sys

from astropy.table import Table
import drizzlepac

from drizzlepac.hlautils import config_utils
from drizzlepac.hlautils import poller_utils

def run_hla_flag_filter():
    from drizzlepac.hlautils import hla_flag_filter

    #   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +
    # All below lines are to get it working, not actual final code.
    out_file = glob.glob("??????.out")[0]
    # out_file = "j92c01.out" # acs_10265_01
    # #out_file = "j9es06.out" # acs_10595_06
    # Get parameter values
    if os.getcwd().endswith("orig"): sys.exit("Don't run in the orig dir! YOU'LL RUIN EVERYTHING!")
    for cmd in ['rm -f *.*', 'cp orig/* .']:
        print(cmd)
        os.system(cmd)



    obs_info_dict, total_list = poller_utils.interpret_obset_input(out_file)
    out_pars_file = "pars.json"
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
    if out_file == "j92c01.out": # acs_10265_01
        # settings for testing ~/Documents/HLAtransition/runhlaprocessing_testing/acs_10265_01/flag_testing/hla
        mode = "dao"
        drizzled_image = "hst_10265_01_acs_wfc_f606w_drz.fits"
        flt_list = ["j92c01b4q_flc.fits", "j92c01b5q_flc.fits", "j92c01b7q_flc.fits", "j92c01b9q_flc.fits"]
        param_dict = total_list[0].fdp_list[0].configobj_pars.as_single_giant_dict()
        param_dict['quality control']['ci filter']['sourcex_bthresh'] = 5.0  # force it to use the value from HLA classic
        param_dict['quality control']['ci filter']['dao_bthresh'] = 5.0  # force it to use the value from HLA classic
        exptime = 5060.0
        catalog_name = "hst_10265_01_acs_wfc_f606w_{}phot.txt".format(mode)
        catalog_data = Table.read(catalog_name, format='ascii')
        proc_type = "{}phot".format(mode)
        drz_root_dir = os.getcwd()


        # for filt_key in filter_sorted_flt_dict.keys(): flt_list = filter_sorted_flt_dict[filt_key]
        # os.remove("hst_10265_01_acs_wfc_f606w_msk.fits")
        # from devutils import make_mask_file
        # make_mask_file.make_mask_file_old(all_drizzled_filelist[0].replace("drz.fits","wht.fits"))

        comp_cmd = "python /Users/dulude/Documents/Code/HLATransition/drizzlepac/drizzlepac/devutils/comparison_tools/compare_sourcelists.py orig/hst_10265_01_acs_wfc_f606w_{}phot_orig.txt hst_10265_01_acs_wfc_f606w_{}phot.txt -i hst_10265_01_acs_wfc_f606w_drz.fits hst_10265_01_acs_wfc_f606w_drz.fits -m absolute -p none".format(mode,mode)


    if out_file == "j9es06.out": # acs_10595_06
        # settings for testing ~/Documents/HLAtransition/runhlaprocessing_testing/acs_10595_06_flag_testing/
        mode = "sex"
        drizzled_image = "hst_10595_06_acs_wfc_f435w_drz.fits"
        flt_list = ["j9es06rbq_flc.fits", "j9es06rcq_flc.fits", "j9es06req_flc.fits", "j9es06rgq_flc.fits"]
        param_dict = total_list[0].fdp_list[0].configobj_pars.as_single_giant_dict()
        param_dict['quality control']['ci filter']['sourcex_bthresh'] = 5.0 #force it to use the value from HLA classic
        param_dict['quality control']['ci filter']['dao_bthresh'] = 5.0  # force it to use the value from HLA classic
        exptime = 710.0
        catalog_data = Table.read(catalog_name, format='ascii')
        catalog_data = Table.read(dict_newTAB_matched2drz[all_drizzled_filelist[0]], format='ascii')
        proc_type = "{}phot".format(mode)
        drz_root_dir = os.getcwd()

        # os.remove("hst_10595_06_acs_wfc_f435w_msk.fits")
        # from devutils import make_mask_file
        # make_mask_file.make_mask_file("hst_10595_06_acs_wfc_f435w_wht.fits")
        comp_cmd = "python /Users/dulude/Documents/Code/HLATransition/drizzlepac/drizzlepac/devutils/comparison_tools/compare_sourcelists.py orig_cats/hst_10595_06_acs_wfc_f435w_{}phot.txt hst_10595_06_acs_wfc_f435w_{}phot.txt -i hst_10595_06_acs_wfc_f435w_drz.fits hst_10595_06_acs_wfc_f435w_drz.fits -m absolute -p none".format(mode,mode)
    #   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +
    # Execute hla_flag_filter.run_source_list_flaging
    catalog_data = hla_flag_filter.run_source_list_flaging(drizzled_image, flt_list,
                                            param_dict, exptime,
                                            catalog_name, catalog_data,
                                            proc_type, drz_root_dir, debug = True)


    catalog_data.write(catalog_name, delimiter=",",format='ascii',overwrite=True)
    print("Wrote {}".format(catalog_name))
    try:
        os.system(comp_cmd)
    except:
        print("skipping automatic comparision run")

#=======================================================================================================================
def run_hla_flag_filter_HLAClassic():
    from drizzlepac.hlautils import hla_flag_filter_HLAClassic

    #   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +
    # All below lines are to get it working, not actual final code.
    out_file = glob.glob("??????.out")[0]
    # out_file = "j92c01.out" # acs_10265_01
    # #out_file = "j9es06.out" # acs_10595_06
    # Get parameter values
    if os.getcwd().endswith("orig"): sys.exit("Don't run in the orig dir! YOU'LL RUIN EVERYTHING!")
    for cmd in ['rm -f *.*', 'cp orig/* .']:
        print(cmd)
        os.system(cmd)



    obs_info_dict, total_list = poller_utils.interpret_obset_input(out_file)
    out_pars_file = "pars.json"
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
    if out_file == "j92c01.out": # acs_10265_01
        # settings for testing ~/Documents/HLAtransition/runhlaprocessing_testing/acs_10265_01/flag_testing/hla
        mode = "dao"
        all_drizzled_filelist = ["hst_10265_01_acs_wfc_f606w_drz.fits"]
        working_hla_red = os.getcwd()
        filter_sorted_flt_dict = {"f606w": ["j92c01b4q_flc.fits", "j92c01b5q_flc.fits", "j92c01b7q_flc.fits", "j92c01b9q_flc.fits"]}
        param_dict = total_list[0].fdp_list[0].configobj_pars.as_single_giant_dict()
        param_dict['quality control']['ci filter']['sourcex_bthresh'] = 5.0  # force it to use the value from HLA classic
        param_dict['quality control']['ci filter']['dao_bthresh'] = 5.0  # force it to use the value from HLA classic
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

        # for filt_key in filter_sorted_flt_dict.keys(): flt_list = filter_sorted_flt_dict[filt_key]
        # os.remove("hst_10265_01_acs_wfc_f606w_msk.fits")
        # from devutils import make_mask_file
        # make_mask_file.make_mask_file_old(all_drizzled_filelist[0].replace("drz.fits","wht.fits"))

        comp_cmd = "python /Users/dulude/Documents/Code/HLATransition/drizzlepac/drizzlepac/devutils/comparison_tools/compare_sourcelists.py orig/hst_10265_01_acs_wfc_f606w_{}phot_orig.txt hst_10265_01_acs_wfc_f606w_{}phot.txt -i hst_10265_01_acs_wfc_f606w_drz.fits hst_10265_01_acs_wfc_f606w_drz.fits -m absolute -p none".format(mode,mode)


    if out_file == "j9es06.out": # acs_10595_06
        # settings for testing ~/Documents/HLAtransition/runhlaprocessing_testing/acs_10595_06_flag_testing/
        mode = "sex"
        all_drizzled_filelist = ["hst_10595_06_acs_wfc_f435w_drz.fits"]
        working_hla_red = os.getcwd()
        filter_sorted_flt_dict = {"f435w": ["j9es06rbq_flc.fits", "j9es06rcq_flc.fits", "j9es06req_flc.fits", "j9es06rgq_flc.fits"]}
        param_dict = total_list[0].fdp_list[0].configobj_pars.as_single_giant_dict()
        param_dict['quality control']['ci filter']['sourcex_bthresh'] = 5.0 #force it to use the value from HLA classic
        param_dict['quality control']['ci filter']['dao_bthresh'] = 5.0  # force it to use the value from HLA classic
        readnoise_dictionary_drzs = {"hst_10595_06_acs_wfc_f435w_drz.fits": 5.247499925}
        scale_dict_drzs = {"hst_10595_06_acs_wfc_f435w_drz.fits": 0.05}
        zero_point_AB_dict = {"hst_10595_06_acs_wfc_f435w_drz.fits": 25.6888167958}
        exp_dictionary_scis = {"hst_10595_06_acs_wfc_f435w_drz.fits": 710.0}
        detection_image = "hst_10595_06_acs_wfc_total_drz.fits"
        dict_newTAB_matched2drz = {"hst_10595_06_acs_wfc_f435w_drz.fits": "hst_10595_06_acs_wfc_f435w_{}phot.txt".format(mode)}
        phot_table_matched2cat = {all_drizzled_filelist[0]: Table.read(dict_newTAB_matched2drz[all_drizzled_filelist[0]], format='ascii')}
        proc_type = "{}phot".format(mode)
        drz_root_dir = os.getcwd()
        rms_dict = {"hst_10595_06_acs_wfc_f435w_drz.fits": "hst_10595_06_acs_wfc_f435w_rms.fits"}
        # os.remove("hst_10595_06_acs_wfc_f435w_msk.fits")
        # from devutils import make_mask_file
        # make_mask_file.make_mask_file("hst_10595_06_acs_wfc_f435w_wht.fits")
        comp_cmd = "python /Users/dulude/Documents/Code/HLATransition/drizzlepac/drizzlepac/devutils/comparison_tools/compare_sourcelists.py orig_cats/hst_10595_06_acs_wfc_f435w_{}phot.txt hst_10595_06_acs_wfc_f435w_{}phot.txt -i hst_10595_06_acs_wfc_f435w_drz.fits hst_10595_06_acs_wfc_f435w_drz.fits -m absolute -p none".format(mode,mode)
    #   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +   +
    # Execute hla_flag_filter.run_source_list_flaging
    catalog_data = hla_flag_filter_HLAClassic.run_source_list_flaging(all_drizzled_filelist, filter_sorted_flt_dict,
                                            param_dict, exp_dictionary_scis,
                                            dict_newTAB_matched2drz, phot_table_matched2cat,
                                            proc_type, drz_root_dir, debug = True)

    catalog_name = dict_newTAB_matched2drz[all_drizzled_filelist[0]]
    catalog_data.write(catalog_name, delimiter=",",format='ascii',overwrite=True)
    print("Wrote {}".format(catalog_name))
    try:
        os.system(comp_cmd)
    except:
        print("skipping automatic comparision run")

if __name__ == "__main__":
    run_hla_flag_filter_HLAClassic()
