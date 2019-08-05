#!/usr/bin/env python

"""This script contains code to create the complete set of configuration parameters required to run XXXXX given the
specified observation conditions and instrument/detector used in the observations"""

import collections
import json
import os
import pdb
import sys

from astropy.time import Time

# TODO: Does this module require logging?
# ======================================================================================================================


class HapConfig(object):
    def __init__(self, prod_obj, use_defaults=True, input_custom_pars_file=None, output_custom_pars_file=None):
        """
        A set of routines to generate appropriate set of configuration parameters.

        Parameters
        ----------
        prod_obj : drizzlepac.hlautils.Product.TotalProduct, drizzlepac.hlautils.Product.FilterProduct, or
        drizzlepac.hlautils.Product.ExposureProduct, depending on input
            Product to get configuration values for.

        use_defaults : bool, optional
            Use default configuration parameters? Default value is True.

        input_custom_pars_file: str, optional
            Name of the full configuration file (with full path) to use for ALL input params. WARNING: Specifying a
            file will turn off automatic parameter determination.

        output_custom_pars_file: str, optional
            Name of the full configuration file (with full path) that all parameters will be written to.

        Returns
        -------
        Nothing.
        """
        if input_custom_pars_file and input_custom_pars_file and input_custom_pars_file == output_custom_pars_file:
            sys.exit("ERROR: Input and output parameter files must have unique names!")
        self.label = "hap_config"
        self.description = "A set of routines to generate appropriate set of configuration parameters"
        self.instrument = prod_obj.instrument
        self.detector = prod_obj.detector
        self.inst_det = "{}_{}".format(prod_obj.instrument, prod_obj.detector).lower()
        self.use_defaults = use_defaults
        self.input_custom_pars_file = input_custom_pars_file
        self.output_custom_pars_file = output_custom_pars_file
        self._determine_conditions(prod_obj)
        self._get_cfg_index()

        # Instantiate the parameter set
        self.pars = {}

        # open input parameter file if specified by user
        if self.input_custom_pars_file:
            with open(self.input_custom_pars_file) as f_cfg:
                self.input_cfg_json_data = json.load(f_cfg)[prod_obj.product_basename]
        else:
            self.input_cfg_json_data = None

        # generate parameter sets for each pipeline step
        step_name_list = [AlignmentPars, AstrodrizzlePars, CatalogGenerationPars, QualityControlPars]
        step_title_list = ['alignment', 'astrodrizzle', 'catalog generation', 'quality control']
        # step_name_list = [AstrodrizzlePars, CatalogGenerationPars, QualityControlPars]
        # step_title_list = ['astrodrizzle', 'catalog generation', 'quality control']
        for step_title, step_name in zip(step_title_list, step_name_list):
            cfg_index = self.full_cfg_index[step_title]
            self.pars[step_title] = step_name(cfg_index,
                                              self.conditions,
                                              self.pars_dir,
                                              step_title,
                                              self.use_defaults,
                                              self.input_cfg_json_data)

        # write out all parameters to file if specified by user
        if output_custom_pars_file:
            self.write_pars(prod_obj)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def _determine_conditions(self, prod_obj):
        """Determine observing condition or conditions present for a given step

        Parameters
        ----------
        prod_obj : drizzlepac.hlautils.Product.TotalProduct, drizzlepac.hlautils.Product.FilterProduct, or
        drizzlepac.hlautils.Product.ExposureProduct, depending on input
            Product to get configuration values for.

        Returns
        -------
        Nothing.

        """

        # determine product type, initialize and build conditions list
        if hasattr(prod_obj, "edp_list") and hasattr(prod_obj, "fdp_list"):  # For total products
            self.conditions = ["total_basic"]
        elif hasattr(prod_obj, "edp_list") and not hasattr(prod_obj, "fdp_list"):  # For filter products
            self.conditions = ["filter_basic"]
            n_exp = len(prod_obj.edp_list)
            if n_exp == 1:
                self.conditions.append("any_n1")
            else:
                if self.instrument == "acs":
                    if self.detector == "hrc":
                        if n_exp in [2, 3]:
                            self.conditions.append("acs_hrc_any_n2")
                        if n_exp in [4, 5]:
                            self.conditions.append("acs_hrc_any_n4")
                        if n_exp >= 6:
                            self.conditions.append("acs_hrc_any_n6")
                    elif self.detector == "sbc":
                        if self.filters.lower() in ["f115lp", "f122m"]:
                            if n_exp in [2, 3, 4, 5]:
                                self.conditions.append("acs_sbc_blue_n2")
                            if n_exp >= 6:
                                self.conditions.append("acs_sbc_blue_n6")
                        else:
                            if n_exp in [2, 3, 4, 5]:
                                self.conditions.append("acs_sbc_any_n2")
                            if n_exp >= 6:
                                self.conditions.append("acs_sbc_any_n6")
                    elif self.detector == "wfc":
                        if n_exp in [2, 3]:
                            self.conditions.append("acs_wfc_any_n2")
                        if n_exp in [4, 5]:
                            self.conditions.append("acs_wfc_any_n4")
                        if n_exp >= 6:
                            self.conditions.append("acs_wfc_any_n6")
                    else:
                        sys.exit("INVALID ACS DETECTOR!")
                elif self.instrument == "wfc3":
                    if self.detector == "ir":
                        if self.filters.lower() in ["g102", "g141"]:
                            if n_exp in [2, 3]:
                                self.conditions.append("wfc3_ir_grism_n2")
                            if n_exp >= 4:
                                self.conditions.append("wfc3_ir_grism_n4")
                        else:
                            if n_exp in [2, 3]:
                                self.conditions.append("wfc3_ir_any_n2")
                            if n_exp >= 4:
                                self.conditions.append("wfc3_ir_any_n4")
                    elif self.detector == "uvis":
                        thresh_time = Time("2012-11-08T02:59:15", format='isot', scale='utc').mjd
                        if self.mjd >= thresh_time:
                            if n_exp in [2, 3]:
                                self.conditions.append("wfc3_uvis_any_post_n2")
                            if n_exp in [4, 5]:
                                self.conditions.append("wfc3_uvis_any_post_n4")
                            if n_exp >= 6:
                                self.conditions.append("wfc3_uvis_any_post_n6")
                        else:
                            if n_exp in [2, 3]:
                                self.conditions.append("wfc3_uvis_any_pre_n2")
                            if n_exp in [4, 5]:
                                self.conditions.append("wfc3_uvis_any_pre_n4")
                            if n_exp >= 6:
                                self.conditions.append("wfc3_uvis_any_pre_n6")
                    else:
                        sys.exit("INVALID WFC3 DETECTOR!")
                else:
                    sys.exit("INVALID HST INSTRUMENT!")
        else:  # For single-exposure products
            self.conditions = ["single_basic"]
            self.conditions.append("any_n1")  # TODO: verify that single-exposure products should use nexp=1 cfg file.


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def _get_cfg_index(self):
        """return the contents of the appropriate index cfg file."""
        code_dir = os.path.abspath(__file__)
        base_dir = os.path.dirname(os.path.dirname(code_dir))
        self.pars_dir = os.path.join(base_dir, "pars/hap_pars")
        cfg_index_fileanme = self.inst_det + "_index.json"
        cfg_index_filename = os.path.join(self.pars_dir, cfg_index_fileanme)

        with open(cfg_index_filename) as jsonFile:
            json_string = jsonFile.read()
            self.full_cfg_index = json.loads(json_string)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def get_pars(self, step_name):
        """This method returns the parameter set for a specified step (alignment, astrodrizzle, etc.)

        Parameters
        ----------
        step_name : str
            Name of the step for which parameters will be returned.

        Returns
        -------
        Dictionary of parameter values keyed by parameter name.
        """
        step_list = ['alignment', 'astrodrizzle', 'catalog generation', 'quality control']
        if step_name in step_list:
            return self.pars[step_name].outpars
        else:
            print("ERROR! '{}' is not a recognized step name.".format(step_name))
            print("Recognized step names: \n{}".format(str(step_list)[2:-2].replace("', '", "\n")))
            sys.exit(1)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def write_pars(self, prod_obj):
        """This method writes the current parameter set to the specified file.

        Parameters
        ----------
        prod_obj : drizzlepac.hlautils.Product.TotalProduct, drizzlepac.hlautils.Product.FilterProduct, or
        drizzlepac.hlautils.Product.ExposureProduct, depending on input
            Product to write configuration values for.

        Returns
        -------
        Nothing.
        """
        new_json_data = {}
        for stepname in self.pars.keys():
            new_json_data[stepname] = self.pars[stepname].outpars
        new_json_data = {prod_obj.product_basename: {"parameters": new_json_data, "default_values": new_json_data}}

        if os.path.exists(self.output_custom_pars_file):
            with open(self.output_custom_pars_file) as f:
                json_data = json.load(f)

            json_data.update(new_json_data)

            with open(self.output_custom_pars_file, 'w') as f:
                json.dump(json_data, f, indent=4)
            print("Updated custom pars file {}".format(self.output_custom_pars_file))
        else:
            with open(self.output_custom_pars_file, 'w') as f:
                json.dump(new_json_data, f, indent=4)
            print("Wrote custom pars file {}".format(self.output_custom_pars_file))


# ----------------------------------------------------------------------------------------------------------------------


class Par():
    def __init__(self, cfg_index, conditions, pars_dir, step_title, use_defaults, input_cfg_json_data):
        """Parent class for alignment_pars, astrodrizzle_pars, catalog_generation_pars, and quality_control_pars

        Parameters
        ----------
        cfg_index : dictionary
            portion of the index config file returned for a specific step

        conditions : list
            list of observing conditions that will be used to build the final composite parameter set.

        pars_dir : str
            full path of the directory that contains the config files

        step_title : str
            name of the specified step (alignment, astrodrizzle, catalog generation, or quality control)

        use_defaults : bool
            Use default values for all configuration parameters?

        input_cfg_json_data : dict
            dictionary containing custom parameter settings from input custom parameter file. NOTE: if no custom input
            parameter file was specified, *input_cfg_json_data* is set to None.

        Returns
        -------
        Nothing.

        """
        self.cfg_index = cfg_index
        self.conditions = conditions
        self.pars_dir = pars_dir
        self.step_title = step_title
        self.use_defaults = use_defaults
        self.input_cfg_json_data = input_cfg_json_data
        if not self.input_cfg_json_data:
            self._get_params()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def _combine_conditions(self):
        """Combine parameters from multiple conditions into a single parameter set.
        """
        self.outpars = {}
        for cfg_key in self.pars_multidict.keys():
            self.outpars = self._dict_merge(self.outpars, self.pars_multidict[cfg_key])


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def _dict_merge(self, dct, merge_dct, add_keys=True):
        """ Recursive dict merge. Inspired by :meth:``dict.update()``, instead of
        updating only top-level keys, dict_merge recurses down into dicts nested
        to an arbitrary depth, updating keys. The ``merge_dct`` is merged into
        ``dct``.

        This version will return a copy of the dictionary and leave the original
        arguments untouched.

        The optional argument ``add_keys``, determines whether keys which are
        present in ``merge_dict`` but not ``dct`` should be included in the
        new dict.

        NOTE: this method was found on stack overflow
        (URL: https://gist.github.com/angstwad/bf22d1822c38a92ec0a9)
        It was written by user "DomWeldon".

        Parameters
        ----------
        dct : dictionary
            dictionary onto which the merge is executed

        merge_dct : dictionary
            dictionary merged into dct

        add_keys : Boolean, optional
            whether to add new keys if they don't exist in dct. Default value is True

        Returns
        -------
        dict: dictionary
            updated dictionary
        """
        dct = dct.copy()
        if not add_keys:
            merge_dct = {
                k: merge_dct[k]
                for k in set(dct).intersection(set(merge_dct))
            }

        for k, v in merge_dct.items():
            if (k in dct and isinstance(dct[k], dict)
                    and isinstance(merge_dct[k], collections.Mapping)):
                dct[k] = self._dict_merge(dct[k], merge_dct[k], add_keys=add_keys)
            else:
                dct[k] = merge_dct[k]

        return dct

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def _get_params(self):
        """read in params from config files based on instrument, detector, and condition(s), and return a ordered
        dictionary of these values."""
        self.pars_multidict = collections.OrderedDict()
        found_cfg = False
        for condition in self.conditions:
            if condition in self.cfg_index.keys():
                found_cfg = True
                self.subcfgfilename = os.path.join(self.pars_dir, self.cfg_index[condition])
                self.pars_multidict[condition] = self._read_json_file()
        # if no specific cfg files can be found for the specified conditions, use the generic cfg file.
        if not found_cfg:
            self.subcfgfilename = os.path.join(self.pars_dir, self.cfg_index["all"])
            self.pars_multidict["all"] = self._read_json_file()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def _read_custom_pars(self):
        """process parameters from user-specified input pars file"""
        if self.use_defaults:
            param_set = "default_values"
        else:
            param_set = "parameters"
        self.outpars = self.input_cfg_json_data[param_set][self.step_title]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def _read_json_file(self):
        """read in json file, return read-in information as dictionary of values

        Parameters
        ----------
        None.

        Returns
        -------
        json_data : dictionary
            information from json file
        """
        with open(self.subcfgfilename) as json_file:
            json_string = json_file.read()
            json_data = json.loads(json_string)
            if self.use_defaults:
                return json_data['default_values']
            else:
                return json_data['parameters']


# ----------------------------------------------------------------------------------------------------------------------


class AlignmentPars(Par):
    def __init__(self, cfg_index, conditions, pars_dir, step_title, use_defaults, input_cfg_json_data):
        """Configuration parameters for the image alignment step. See Par.__init__() for input argument definitions."""
        super().__init__(cfg_index, conditions, pars_dir, step_title, use_defaults, input_cfg_json_data)
        self.set_name = "alignment"
        if input_cfg_json_data:
            self._read_custom_pars()
        else:
            self._combine_conditions()


# ----------------------------------------------------------------------------------------------------------------------


class AstrodrizzlePars(Par):
    def __init__(self, cfg_index, conditions, pars_dir, step_title, use_defaults, input_cfg_json_data):
        """Configuration parameters for the AstroDrizzle step. See Par.__init__() for input argument definitions."""
        super().__init__(cfg_index, conditions, pars_dir, step_title, use_defaults, input_cfg_json_data)
        if input_cfg_json_data:
            self._read_custom_pars()
        else:
            self._combine_conditions()


# ----------------------------------------------------------------------------------------------------------------------


class CatalogGenerationPars(Par):
    def __init__(self, cfg_index, conditions, pars_dir, step_title, use_defaults, input_cfg_json_data):
        """Configuration parameters for the photometric catalog generation step. See Par.__init__() for input argument
        definitions."""
        super().__init__(cfg_index, conditions, pars_dir, step_title, use_defaults, input_cfg_json_data)
        if input_cfg_json_data:
            self._read_custom_pars()
        else:
            self._combine_conditions()


# ----------------------------------------------------------------------------------------------------------------------


class QualityControlPars(Par):
    def __init__(self, cfg_index, conditions, pars_dir, step_title, use_defaults, input_cfg_json_data):
        """Configuration parameters for the quality control step. See Par.__init__() for input argument definitions."""
        super().__init__(cfg_index, conditions, pars_dir, step_title, use_defaults, input_cfg_json_data)
        if input_cfg_json_data:
            self._read_custom_pars()
        else:
            self._combine_conditions()


# ======================================================================================================================


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
    import drizzlepac
    from stsci.tools import teal

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
