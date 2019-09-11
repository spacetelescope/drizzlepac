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
json_list = Path(full_pars_path).glob("**/*.json")


for json_file in json_list:
    json_file = str(json_file)
    if not json_file.endswith("_index.json"):
        print(json_file)
        # open and read in current unsplit json file
        with open(json_file, 'r') as jfin:
            json_data = json.load(jfin)

        # write out user-editable parameter set to file
        # generate output filename
        out_file = json_file.replace('hap_pars/', 'hap_pars/user_parameters/')
        # Generate output path if it doesn't already exist
        out_path = os.path.dirname(os.path.abspath(out_file))
        if not os.path.exists(out_path):
            cmd = "mkdir -p {}".format(out_path)
            print(cmd)
            os.system(cmd)
        # write output file
        with open(out_file, 'w') as f:
            json.dump(json_data['parameters'], f, indent=4)
        print('Wrote {}'.format(out_file))

        # write out user-editable parameter set to file
        # generate output fiename
        out_file = json_file.replace('hap_pars/','hap_pars/default_parameters/')
        # Generate output path if it doesn't already exist
        out_path = os.path.dirname(os.path.abspath(out_file))
        if not os.path.exists(out_path):
            cmd = "mkdir -p {}".format(out_path)
            print(cmd)
            os.system(cmd)
        # write output file
        with open(out_file, 'w') as f:
            json.dump(json_data['default_values'], f, indent=4)
        print('Wrote {}'.format(out_file))

        # remove existing unsplit file
        # print("removing existing unsplit file {}".format("json_file"))
        # os.remove(json_file)

        print("\n")
