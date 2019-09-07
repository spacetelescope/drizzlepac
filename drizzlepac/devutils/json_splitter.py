#!/usr/bin/env python
"""searches through hap_pars subdirectory tree and splits the json file 'default_values' dictinary into one
file, and the 'parameter_values' into another."""

import json
import glob
import os
from pathlib import Path
import pdb

# build hap_path full path and recursivly find all .json pars files
code_dir = os.path.abspath(__file__)
base_dir = os.path.dirname(os.path.dirname(code_dir))
full_pars_path = os.path.join(base_dir, "pars/hap_pars")
json_list = Path(full_pars_path).glob("*/*/*_test.json")


for json_file in json_list:
    json_file = str(json_file)
    print(json_file)
    with open(json_file, 'r') as jfin:
        json_data = json.load(jfin)
    out_file = json_file.replace(".json", "_user.json")
    with open(out_file, 'w') as f:
        json.dump(json_data['parameters'], f, indent=4)
    print('Wrote {}'.format(out_file))

    out_file = json_file.replace('.json', '_default.json')
    with open(out_file, 'w') as f:
        json.dump(json_data['default_values'], f, indent=4)
    print('Wrote {}'.format(out_file))

    os.remove(json_file)
    print('removed {}'.format(json_file))
    print("\n")
