#!/usr/bin/env python

"""This script removes the 'input' and 'output' fields and sets the 'build' and
['STATE OF INPUT FILES']['clean'] fields to 'True' in all astrodrizzle parameter .json files."""

import json
import glob
from pathlib import Path

code_dir = os.path.abspath(__file__)
base_dir = os.path.dirname(os.path.dirname(code_dir))
full_pars_path = os.path.join(base_dir, "pars/hap_pars")
json_list = Path('full_pars_path').glob('**/*astrodrizzle*.json')
for json_file in json_list:
    json_file = str(json_file)
    with open(json_file) as jfin:
        json_data = json.load(jfin)
    update_flag=False
    for dkey in json_data.keys():
        if "input" in json_data[dkey]:
            print("{}: [{}]['input'] removed.".format(json_file, dkey))
            del json_data[dkey]['input']
            update_flag = True

        if "output" in json_data[dkey]:
            print("{}: [{}]['output'] removed.".format(json_file, dkey))
            del json_data[dkey]['output']
            update_flag = True

        if json_data[dkey]['build'] == False:
            print("{}: [{}]['build'] reset from False to True.".format(json_file,dkey))
            json_data[dkey]['build'] = True
            update_flag = True

        if json_data[dkey]['STATE OF INPUT FILES']['clean'] == False:
            print("{}: [{}]['STATE OF INPUT FILES']['clean'] reset from False to True.".format(json_file,dkey))
            json_data[dkey]['clean'] = True
            update_flag = True
    if update_flag:
        with open(json_file, 'w') as f:
            json.dump(json_data, f, indent=4)
        print("Updated custom pars file {}\n".format(json_file))
    else:
        print("{}: NO CHANGES\n".format(json_file))