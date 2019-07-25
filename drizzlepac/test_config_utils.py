#!/usr/bin/env python
import pdb
import sys


import drizzlepac
from drizzlepac.hlautils import config_utils
from drizzlepac.hlautils import poller_utils


obs_info_dict, expo_list, filt_list, total_list = poller_utils.interpret_obset_input(sys.argv[1])

foo = config_utils.hap_config("acs", "wfc", use_defaults=False)
blarg = foo.get_pars("catalog generation")

print(blarg)

pdb.set_trace()