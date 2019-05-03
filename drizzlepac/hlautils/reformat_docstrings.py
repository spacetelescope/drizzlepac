#!/usr/bin/env python
import pdb
import collections
import sys
# python_filename = sys.argv[1]
# subroutine_name = sys.argv[2]
python_filename = '/Users/dulude/Documents/Code/HLApipeline/HLApipe/scripts/hla_reduction.py'
subroutine_name = 'pixarea'
with open(python_filename, 'r') as infile:
    file_lines = infile.readlines()
found_the_subroutine = False
ds_begin = None
param_dict = collections.OrderedDict()
for file_line in file_lines:
    file_line = file_line.strip()
    if file_line.startswith("def {}(".format(subroutine_name)):
        found_the_subroutine = True
    if ds_begin == True and file_line.startswith('"""'):
        ds_begin = False
        break
    if ds_begin:
        print(ds_begin,file_line)
        last_line = ""
        if file_line.startswith(":param "):
            param_name = file_line.replace(":param ","").split(":")[0]
            param_descrip = file_line.replace(":param {}:".format(param_name),"")
            param_dict[param_name] = {"descrip": param_descrip, "type": 'n/a'}
            print(">>>>>",file_line,">> ",param_name, param_descrip)
            last_line = file_line
        if file_line.startswith(":type "):
            param_name = file_line.replace(":type ", "").split(":")[0]
            param_dict[param_name]["type"] = file_line.replace(":type {}:".format(param_name),"")

        if last_line != "":
            if last_line.startswith(":param ") and file_line.startswith(":param ") == False and file_line.startswith(":type ")==False:
                print("*********",file_line)

        # if line

    if found_the_subroutine and file_line.startswith('"""') and ds_begin == None:
        ds_begin = True

if ds_begin == False:
    for item in param_dict.keys():
        print(">>>>>-----",item,param_dict[item])
    print("\tParameters")
    print("\t----------")
    for key in param_dict.keys():
        print("\t{} : {}".format(key,param_dict[key]['type']))
        print("\t\t{}\n".format(param_dict[key]['descrip']))


