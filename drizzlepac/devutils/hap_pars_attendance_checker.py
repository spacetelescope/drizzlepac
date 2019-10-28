#!/usr/bin/env python
"""This script checks to make sure all nessessary HAP json parameter files are present"""
import glob
import os
from pathlib import Path
import pdb

code_dir = os.path.abspath(__file__)
base_dir = os.path.dirname(os.path.dirname(code_dir))
full_pars_path = os.path.join(base_dir, "pars/hap_pars")
json_list = Path(full_pars_path).glob("**/*.json")
pdb.set_trace()

param_list = ['default_parameters', 'user_parameters']
inst_det_list = ['acs_hrc', 'acs_sbc', 'acs_wfc', 'wfc3_ir', 'wfc3_uvis']