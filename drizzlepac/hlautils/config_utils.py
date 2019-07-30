#!/usr/bin/env python

"""This script contains code to create the complete set of configuration parameters required to run XXXXX given the
specified observation conditions and instrument/detector used in the observations"""

import collections
import json
import os
import pdb
import sys

from astropy.time import Time


# ======================================================================================================================


class hap_config(object):
    def __init__(self,prod_obj,use_defaults=False,cfg_index_file=None):
        """
        A set of routines to generate appropriate set of configuration parameters

        Parameters
        ----------
        instrument : str
            instrument name

        detector : str
            detector name

        use_defaults : bool, optional
            Use default values for all configuration parameters? Default value is False.

        cfg_index_file : str, optional
            Name of the full configuration file (with full path) to use for ALL input params. WARNING: Specifying a
            file will turn off automatic parameter determination.

        Returns
        -------
        Nothing.
        """
        self.label = "hap_config"
        self.description = "A set of routines to generate appropriate set of configuration parameters"
        self.instrument = prod_obj.instrument
        self.detector = prod_obj.detector
        self.inst_det = "{}_{}".format(prod_obj.instrument,prod_obj.detector).lower()
        self.cfg_index_file = cfg_index_file
        self.use_defaults = use_defaults
        self._determine_conditions(prod_obj)
        self._get_cfg_index()

        # Instantiate the parameter set
        self.pars = {}
        #step_list = [alignment_pars,astrodrizzle_pars,catalog_generation_pars,quality_control_pars] # TODO: uncomment when everything is working
        step_list = [astrodrizzle_pars,catalog_generation_pars] # TODO: Just a placeholder until we add complexity!

        for step_name in step_list:
            step_title = step_name.__name__.replace("_pars","").replace("_"," ")
            cfg_index = self.full_cfg_index[step_title]
            self.pars[step_title] = step_name(cfg_index,self.conditions,self.pars_dir,step_title,self.use_defaults)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def _determine_conditions(self,prod_obj):
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

        #determine product type, initialize and build conditions list
        if (hasattr(prod_obj,"edp_list") and hasattr(prod_obj,"fdp_list")): # For total products
            self.conditions = ["total_basic"]
        elif (hasattr(prod_obj,"edp_list") and not hasattr(prod_obj,"fdp_list")): # For filter products
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
                            if n_exp in [2,3,4,5]:
                                self.conditions.append("acs_sbc_blue_n2")
                            if n_exp >= 6:
                                self.conditions.append("acs_sbc_blue_n6")
                        else:
                            if n_exp in [2, 3, 4, 5]:
                                self.conditions.append("acs_sbc_any_n2")
                            if n_exp >= 6:
                                self.conditions.append("acs_sbc_any_n6")
                    elif self.detector == "wfc":
                        if n_exp in [2,3]:
                            self.conditions.append("acs_wfc_any_n2")
                        if n_exp in [4,5]:
                            self.conditions.append("acs_wfc_any_n4")
                        if n_exp >=6:
                            self.conditions.append("acs_wfc_any_n6")
                    else:
                        sys.exit("INVALID ACS DETECTOR!")
                elif self.instrument == "wfc3":
                    if self.detector == "ir":
                        if self.filters.lower() in ["g102","g141"]: # grisms
                            if n_exp in [2,3]:
                                self.conditions.append("wfc3_ir_grism_n2")
                            if n_exp >= 4:
                                self.conditions.append("wfc3_ir_grism_n4")
                        else: # everything else that's not a grism
                            if n_exp in [2,3]:
                                self.conditions.append("wfc3_ir_any_n2")
                            if n_exp >= 4:
                                self.conditions.append("wfc3_ir_any_n4")
                    elif self.detector == "uvis":
                        thresh_time = Time("2012-11-08T02:59:15", format='isot', scale='utc').mjd # TODO: Verify that this time is UTC!
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
        else: # For single-exposure products
            self.conditions = ["single_basic"]
            self.conditions.append("any_n1") # TODO: Double check that single-exposure products should use Jen's nexp=1 cfg file.


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def _get_cfg_index(self):
        """return the contents of the appropriate index cfg file."""
        code_dir = os.path.abspath(__file__)
        base_dir = os.path.dirname(os.path.dirname(code_dir))
        self.pars_dir = os.path.join(base_dir, "pars/hap_pars")
        cfg_index_fileanme = self.inst_det + "_index.json"
        cfg_index_filename = os.path.join(self.pars_dir, cfg_index_fileanme)

        with open(cfg_index_filename) as jsonFile:
            jsonString = jsonFile.read()
            self.full_cfg_index = json.loads(jsonString)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def get_pars(self,step_name):
        """This method returns the parameter set for a specified step (alignment, astrodrizzle, etc.)

        Parameters
        ----------
        step_name : str
            Name of the step for which parameters will be returned.

        Returns
        -------
        Dictionary of parameter values keyed by parameter name.
        """
        step_list = ['alignment', 'astrodrizzle','catalog generation','quality control']
        if step_name in step_list:
            return(self.pars[step_name].outpars)
        else:
            print("ERROR! '{}' is not a recognized step name.".format(step_name))
            print("Recognized step names: \n{}".format(str(step_list)[2:-2].replace("', '","\n")))
            sys.exit(1)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def write_pars(self,out_filename):
        """This method writes the current parameter set to the specified file."""

        pass


#-----------------------------------------------------------------------------------------------------------------------


class par():
    def __init__(self,cfg_index,conditions,pars_dir,step_title,use_defaults):
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

        Returns
        -------
        Nothing.

        """
        self.cfg_index = cfg_index
        self.conditions = conditions
        self.pars_dir = pars_dir
        self.step_title = step_title
        self.use_defaults = use_defaults
        self._get_params()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def _combine_conditions(self):
        """Combine parameters from multiple conditions into a single parameter set.
        """
        self.outpars = {}
        for cfg_key in self.pars_multidict.keys():
            self.outpars = self._dict_merge(self.outpars, self.pars_multidict[cfg_key])


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def _dict_merge(self,dct, merge_dct, add_keys=True):
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
        import collections
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
        if not found_cfg: # if no specific cfg files can be found for the specified conditions, use the generic cfg file.
            self.subcfgfilename = os.path.join(self.pars_dir, self.cfg_index["all"])
            self.pars_multidict["all"] = self._read_json_file()


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
        with open(self.subcfgfilename) as jsonFile:
            jsonString = jsonFile.read()
            json_data = json.loads(jsonString)
            if self.use_defaults:
                return(json_data['default_values'])
            else:
                return(json_data['parameters'])


#-----------------------------------------------------------------------------------------------------------------------


class alignment_pars(par):
    def __init__(self,cfg_index,conditions,pars_dir,step_title,use_defaults):
        """Configuration parameters for the image alignment step"""
        super().__init__(cfg_index,conditions,pars_dir,step_title,use_defaults)
        self._combine_conditions()


#-----------------------------------------------------------------------------------------------------------------------


class astrodrizzle_pars(par):
    def __init__(self,cfg_index,conditions,pars_dir,step_title,use_defaults):
        """Configuration parameters for the AstroDrizzle step"""
        super().__init__(cfg_index,conditions,pars_dir,step_title,use_defaults)
        self._combine_conditions()


#-----------------------------------------------------------------------------------------------------------------------


class catalog_generation_pars(par):
    def __init__(self,cfg_index,conditions,pars_dir,step_title,use_defaults):
        """Configuration parameters for the photometric catalog generation step"""
        super().__init__(cfg_index,conditions,pars_dir,step_title,use_defaults)
        self._combine_conditions()




#-----------------------------------------------------------------------------------------------------------------------


class quality_control_pars(par):
    def __init__(self,cfg_index,conditions,pars_dir,step_title,use_defaults):
        """Configuration parameters for the quality control step"""
        super().__init__(cfg_index,conditions,pars_dir,step_title,use_defaults)
        self._combine_conditions()



# ======================================================================================================================


def cfg2json(cfgfilename,outpath=None):
    """Convert config files to json format

    Parameters
    ----------
    cfgfilename : str
        Input .cfg file to be converted to json format.

    outpath : str
        Destination path of the json file.

    Returns
    -------
    Nothing
    """
    import drizzlepac
    from stsci.tools import teal

    #open cfg file and load up the output dictionary
    cfg_data = teal.load(cfgfilename,strict=False)
    del cfg_data['_task_name_']
    del cfg_data['_RULES_']
    out_dict = {"parameters": cfg_data, "default_values": cfg_data}

    #build output json filename
    json_filename = cfgfilename.split("/")[-1].replace(".cfg", ".json")

    if not outpath:
        code_dir = os.path.abspath(__file__)
        base_dir = os.path.dirname(os.path.dirname(code_dir))
        out_dir = os.path.join(base_dir, "pars/hap_pars")
        det=json_filename.split("_")[0]
        json_filename=json_filename.replace(det,det+"_astrodrizzle")
        if det == "any":
            json_filename = os.path.join(out_dir,det,json_filename)
        else:
            if det in ["hrc","sbc","wfc"]:
                inst = "acs"
            if det in ["ir","uvis"]:
                inst = "wfc3"
            json_filename = "{}_{}".format(inst,json_filename)
            json_filename = os.path.join(out_dir, inst, det, json_filename)
    else:
        json_filename = os.path.join(outpath,json_filename)

    # write out data.
    if not os.path.exists(json_filename):
        with open(json_filename, 'w') as fp:
            json.dump(out_dict, fp)
        print("Wrote {}".format(json_filename))
    else:
        print("Skipped writing {}. File already exists.".format(json_filename))


#-----------------------------------------------------------------------------------------------------------------------


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
    # cfg_path = "/Users/dulude/Documents/Code/HLAtransition/drizzlepac/drizzlepac/pars/"
    # out_path = "/Users/dulude/Documents/Code/HLAtransition/drizzlepac/drizzlepac/pars/hap_pars/any/"
    # cfg_list = ["astrodrizzle_filter_hap.cfg", "astrodrizzle_single_hap.cfg","astrodrizzle_total_hap.cfg"]
    for cfgfile in cfg_list:
        cfgfile = os.path.join(cfg_path,cfgfile)
        cfg2json(cfgfile)


#\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\


if __name__ == '__main__':
    """Super simple testing interface for the above code."""

    import pdb

    # foo = hap_config("acs", "wfc",use_defaults=False)
    # blarg = foo.get_pars("catalog generation")
    #
    # print(blarg)

    # batch_run_cfg2json()
