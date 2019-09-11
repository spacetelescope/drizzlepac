import os
import drizzlepac
from stsci.tools import teal
import json

def cfg2json(cfgfilename, outpath=None):
    """Convert config files to json format

    Parameters
    ----------
    cfgfilename : str
        Input .cfg file to be converted to json format.

    outpath : str, optional
        Destination path of the json file.

    Returns
    -------
    Nothing
    """
    # open cfg file and load up the output dictionary
    cfg_data = teal.load(cfgfilename, strict=False)
    del cfg_data['_task_name_']
    del cfg_data['_RULES_']

    out_dict = {"parameters": cfg_data, "default_values": cfg_data}

    # build output json filename
    json_filename = cfgfilename.split("/")[-1].replace(".cfg", ".json")

    if not outpath:
        code_dir = os.path.abspath(__file__)
        base_dir = os.path.dirname(os.path.dirname(code_dir))
        out_dir = os.path.join(base_dir, "pars/hap_pars")
        det = json_filename.split("_")[0]
        json_filename = json_filename.replace(det, det+"_astrodrizzle")
        if det == "any":
            json_filename = os.path.join(out_dir, det, json_filename)
        else:
            if det in ["hrc", "sbc", "wfc"]:
                inst = "acs"
            if det in ["ir", "uvis"]:
                inst = "wfc3"
            json_filename = "{}_{}".format(inst, json_filename)
            json_filename = os.path.join(out_dir, inst, det, json_filename)
    else:
        json_filename = os.path.join(outpath, "any_"+json_filename)
        json_filename = json_filename.replace("hap.json", "hap_basic.json")

    # write out data.
    if os.path.exists(json_filename):
        os.remove(json_filename)
    with open(json_filename, 'w') as fp:
        json.dump(out_dict, fp, indent=4)
    print("Wrote {}".format(json_filename))


# ----------------------------------------------------------------------------------------------------------------------


def batch_run_cfg2json():
    """run cfg2json() on a predefined list of .cfg files"""
    cfg_path = "/user/mack/hla_cfgs/"
    cfg_list = ['any_n1.cfg',
                'ir_grism_n2.cfg',
                'ir_grism_n4.cfg',
                'ir_any_n2.cfg',
                'ir_any_n4.cfg',
                'uvis_any_n2.cfg',
                'uvis_any_n4.cfg',
                'uvis_any_n6.cfg',
                'uvis_any_pre2012_n2.cfg',
                'uvis_any_pre2012_n4.cfg',
                'uvis_any_pre2012_n6.cfg',
                'wfc_any_n2.cfg',
                'wfc_any_n4.cfg',
                'wfc_any_n6.cfg',
                'sbc_blue_n2.cfg',
                'sbc_blue_n6.cfg',
                'sbc_any_n2.cfg',
                'sbc_any_n6.cfg',
                'hrc_any_n2.cfg',
                'hrc_any_n4.cfg',
                'hrc_any_n6.cfg']
    for cfgfile in cfg_list:
        cfgfile = os.path.join(cfg_path, cfgfile)
        cfg2json(cfgfile)

    cfg_path = "/Users/dulude/Documents/Code/HLAtransition/drizzlepac/drizzlepac/pars/"
    out_path = "/Users/dulude/Documents/Code/HLAtransition/drizzlepac/drizzlepac/pars/hap_pars/any/"
    cfg_list = ["astrodrizzle_filter_hap.cfg", "astrodrizzle_single_hap.cfg", "astrodrizzle_total_hap.cfg"]
    for cfgfile in cfg_list:
        cfgfile = os.path.join(cfg_path, cfgfile)
        cfg2json(cfgfile, outpath=out_path)


# ----------------------------------------------------------------------------------------------------------------------


def reformat_json_file(infilename, outfilename, clobber=False):
    """Reformat user-specifed input json file to use standard (indent = 4) format

    Parameters
    ----------
    infilename: str
        name of json file to reformat

    outfilename : str
        name of output reformatted json file

    clobber : bool, optional
        if file named in outfilename already exists, should it be overwritten? Default value is False.

    Returns
    -------
    Nothing.
    """
    # open input json file
    with open(infilename) as json_file:
        json_string = json_file.read()
        json_data = json.loads(json_string)

    # see if output file already exists and determine course of action
    if os.path.exists(outfilename):
        if clobber:
            os.remove(outfilename)
        else:
            sys.exit("Error: output file {} already exists and would be overwritten!".format(outfilename))

    # write out reformatted json file
    with open(outfilename, 'w') as f:
        json.dump(json_data, f, indent=4)


# ----------------------------------------------------------------------------------------------------------------------