#!/usr/bin/env python
"""This script checks to make sure all necessary HAP json parameter files are present"""
import glob
import os
from pathlib import Path


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def attendance_checker(base_pars_path):
    """check JSON parameter file attendance

    parameters
    ----------
    base_pars_path : string
        path that contains 'user_parameters' and 'default_parameters' subdirectories

    returns
    -------
    nothing.
    """
    full_json_list = Path(base_pars_path).glob("**/*.json")

    param_dict = {'default_parameters': [], 'user_parameters': []}
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
                search_string = "{}/{}/{}/{}/{}_{}*.json".format(base_pars_path, param_set, inst, det, inst_det, step)
                files_found = glob.glob(search_string)

                if param_set == "default_parameters":
                    d_status = "OK     "
                    if not files_found:
                        d_status = "MISSING"
                    status = "{}{}".format(status, d_status)
                if param_set == "user_parameters":
                    u_status = "OK     "
                    if not files_found:
                        u_status = "MISSING"
                    short_search_string = "{}/{}/{}_{}*.json".format(inst, det, inst_det, step)
                    status = "{} {}  {}".format(status, d_status, short_search_string)
            if d_status == "MISSING" or u_status == "MISSING":
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
        if json_file not in param_dict["default_parameters"]:
            d_status = "MISSING"
        if json_file not in param_dict["user_parameters"]:
            u_status = "MISSING"

        print("{} {}  {}".format(d_status, u_status, json_file))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if __name__ == '__main__':
    code_dir = os.path.abspath(__file__)
    base_dir = os.path.dirname(os.path.dirname(code_dir))
    base_pars_path = os.path.join(base_dir, "pars/hap_pars")
    attendance_checker(base_pars_path)
