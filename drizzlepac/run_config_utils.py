#!/usr/bin/env python
import pdb
import sys


import drizzlepac
from drizzlepac import hapsequencer
from drizzlepac.hlautils import config_utils
from drizzlepac.hlautils import poller_utils



input_filename = sys.argv[1]

obs_info_dict, total_list = poller_utils.interpret_obset_input(input_filename)
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

# #
# print(" ")
# print(expo_list[0].configobj_pars.get_pars("alignment"))
#

# rv = hapsequencer.run_hap_processing(input_filename)
pdb.set_trace()

# config_utils.reformat_json_file("superparamfile.json","new_superparamfile.json")

# config_utils.batch_run_cfg2json()











