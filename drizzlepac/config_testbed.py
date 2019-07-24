import pdb
import os
import json
from stsci.tools import teal
# __taskname__ = 'config_testbed'
#
# configobj = teal.load(__taskname__)

inst_det = "acs_wfc"
code_dir = os.path.abspath(__file__)
base_dir = os.path.dirname(os.path.dirname(code_dir))
pars_dir = os.path.join(base_dir,"drizzlepac","pars")
cfg_index_fileanme = inst_det+"_index.cfg"
cfg_index_filename = os.path.join(pars_dir,cfg_index_fileanme)


with open(cfg_index_filename) as jsonFile:
    jsonString = jsonFile.read()
    full_cfg_index = json.loads(jsonString)

cfgfile = os.path.join(pars_dir,full_cfg_index['catalog generation']['default'])



configobj = teal.load(cfgfile)
print(configobj)
pdb.set_trace()