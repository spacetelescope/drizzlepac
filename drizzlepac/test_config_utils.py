#!/usr/bin/env python
import pdb
import sys


import drizzlepac
from drizzlepac.hlautils import config_utils
from drizzlepac.hlautils import poller_utils


#input_filename = sys.argv[1]

obs_info_dict, expo_list, filt_list, total_list = poller_utils.interpret_obset_input('j92c01.out')

param_filename = "superparamfile.json"
# for item in expo_list+filt_list+total_list:
#     print(hasattr(item,"edp_list"),hasattr(item,"fdp_list"))
#     item.pars =config_utils.hap_config(item)


total_list[0].pars = config_utils.hap_config(total_list[0],input_custom_pars_file = param_filename,output_custom_pars_file="afdadf")
# print(total_list[0].pars.get_pars("astrodrizzle"))

filt_list[0].pars = config_utils.hap_config(filt_list[0],input_custom_pars_file = param_filename)
print(filt_list[0].pars.get_pars("astrodrizzle"))





# pdb.set_trace()
