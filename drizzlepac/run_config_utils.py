#!/usr/bin/env python
import pdb
import sys


import drizzlepac
from drizzlepac import hapsequencer
from drizzlepac.hlautils import config_utils
from drizzlepac.hlautils import poller_utils



# input_filename = sys.argv[1]

obs_info_dict, total_list = poller_utils.interpret_obset_input('j92c01.out')
product_list = total_list.copy()

for total_item in total_list:
    product_list += total_item.fdp_list
    product_list += total_item.edp_list

param_filename = "superparamfile.json"

for item in product_list:
    item.pars = config_utils.HapConfig(item,output_custom_pars_file=param_filename,use_defaults=False)
#
print(" ")
print(product_list[1].pars.get_pars("alignment"))


# rv = hapsequencer.run_hla_processing(input_filename)
# pdb.set_trace()

# config_utils.reformat_json_file("superparamfile.json","new_superparamfile.json")

# config_utils.batch_run_cfg2json()











