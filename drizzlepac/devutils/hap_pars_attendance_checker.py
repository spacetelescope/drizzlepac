#!/usr/bin/env python
"""This script checks to make sure all necessary HAP json parameter files are present"""
import glob
import os
from pathlib import Path
import pdb
import json
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def json2list(json_filename):
    with open(json_filename, 'r') as jfin:
        json_data = json.load(jfin)

    foo = flatten_dict(json_data, '', {})

    pdb.set_trace()
    out_list=[]
    return out_list

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def flatten_dict(current, key, result):
    """Flatten nested dictionaries into a non-nested dictionary. Assumes that there are no non-unique keys.
    Code credit: https://stackoverflow.com/questions/24448543/how-would-i-flatten-a-nested-dictionary-in-python-3
    Solution submitted by user 'Matthew Franglen'.
    """
    if isinstance(current, dict):
        for k in current:
            new_key = "{1}".format(key, k) if len(key) > 0 else k
            flatten_dict(current[k], new_key, result)
    else:
        result[key] = current
    return result
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def json_content_checker(json_file_dict, unique_list,base_pars_path):

    json_filename = "/Users/dulude/Documents/Code/HLAtransition/drizzlepac/drizzlepac/pars/hap_pars/default_parameters/any/any_astrodrizzle_filter_hap_basic.json"
    json2list(json_filename)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def attendance_checker(base_pars_path):
    full_json_list = Path(base_pars_path).glob("**/*.json")

    param_dict = {'default_parameters': [], 'user_parameters':[]}
    inst_det_list = ['acs_hrc', 'acs_sbc', 'acs_wfc', 'wfc3_ir', 'wfc3_uvis']
    step_list = ['alignment', 'astrodrizzle', 'catalog_generation', 'quality_control']

    print("DEFAULT USER     FILE")
    # check for the case where all files are missing for given instrument/detector/step
    for inst_det in inst_det_list:
        inst = inst_det.split("_")[0]
        det = inst_det.split("_")[1]
        for step in step_list:
            status = ""
            for param_set in param_dict.keys():
                search_string = "{}/{}/{}/{}/{}_{}*.json".format(base_pars_path, param_set, inst,det, inst_det, step)
                files_found = glob.glob(search_string)

                if param_set == "default_parameters":
                    d_status = "OK     "
                    if not files_found:
                        d_status = "MISSING"
                    status="{}{}".format(status,d_status)
                if param_set == "user_parameters":
                    u_status = "OK     "
                    if not files_found:
                        u_status = "MISSING"
                    short_search_string = "{}/{}/{}_{}*.json".format(inst, det, inst_det, step)
                    status = "{} {}  {}".format(status, d_status, short_search_string)
            if d_status == "MISSING" or u_status =="MISSING":
                print(status)

    # Check that all json files in default_parameters/ are also in user_parameters
    for json_file in full_json_list:
        json_file = str(json_file)
        for param_file_type in param_dict.keys():
            json_path = "{}/{}/".format(base_pars_path, param_file_type)
            if json_file.startswith(json_path):
                param_dict[param_file_type].append(json_file.split(json_path)[1])

    unique_combo_list = list(set(param_dict["default_parameters"] + param_dict["user_parameters"]))

    for json_file in unique_combo_list:
        d_status = "OK     "
        u_status = "OK     "
        if not json_file in param_dict["default_parameters"]:
            d_status = "MISSING"
        if not json_file in param_dict["user_parameters"]:
            u_status = "MISSING"

        print("{} {}  {}".format(d_status, u_status, json_file))
    return param_dict, unique_combo_list

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if __name__ == '__main__':
    code_dir = os.path.abspath(__file__)
    base_dir = os.path.dirname(os.path.dirname(code_dir))
    base_pars_path = os.path.join(base_dir, "pars/hap_pars")
    json_file_dict, unique_list = attendance_checker(base_pars_path)
    json_content_checker(json_file_dict, unique_list,base_pars_path)
