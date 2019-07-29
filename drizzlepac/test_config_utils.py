#!/usr/bin/env python
import pdb
import sys


import drizzlepac
from drizzlepac.hlautils import config_utils
from drizzlepac.hlautils import poller_utils


input_filename = sys.argv[1]
obs_info_dict, expo_list, filt_list, total_list = poller_utils.interpret_obset_input(input_filename)


for item in expo_list+filt_list+total_list:
    print(hasattr(item,"edp_list"),hasattr(item,"fdp_list"))
    item.pars =config_utils.hap_config(item)
# total_list[0].pars = config_utils.hap_config(total_list[0])
# print(total_list[0].pars.get_pars("catalog generation"))
# pdb.set_trace()


# pdb.set_trace()
