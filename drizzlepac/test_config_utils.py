#!/usr/bin/env python
import pdb
import sys


import drizzlepac
from drizzlepac import hapsequencer
from drizzlepac.hlautils import config_utils
from drizzlepac.hlautils import poller_utils



input_filename = sys.argv[1]

# obs_info_dict, expo_list, filt_list, total_list = poller_utils.interpret_obset_input('j92c01.out')
#
# param_filename = "superparamfile.json"
#
# for item in expo_list+filt_list+total_list:
#     item.pars = config_utils.HapConfig(item,input_custom_pars_file=param_filename,use_defaults=False)
# #
# print(" ")
# print(expo_list[0].pars.get_pars("alignment"))
#

rv = hapsequencer.run_hla_processing(input_filename)
# pdb.set_trace()

# config_utils.reformat_json_file("superparamfile.json","new_superparamfile.json")

# config_utils.batch_run_cfg2json()











