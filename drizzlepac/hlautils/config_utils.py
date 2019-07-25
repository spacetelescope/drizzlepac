#!/usr/bin/env python

"""This script contains code to create the complete set of configuration parameters required to run XXXXX given the
specified observation conditions and instrument/detector used in the observations"""

import collections
import json
import os
import sys


from stsci.tools import teal
# ======================================================================================================================


class hap_config(object):
    def __init__(self, instrument,detector,use_defaults=False,cfg_index_file=None):
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
        """
        self.label = "hap_config"
        self.description = "A set of routines to generate appropriate set of configuration parameters"
        self.instrument = instrument
        self.detector = detector
        self.inst_det = "{}_{}".format(instrument,detector).lower()
        self.cfg_index_file = cfg_index_file
        self.use_defaults = use_defaults
        self._determine_conditions()
        self._get_cfg_index()

        # Instantiate the parameter set
        self.pars = {}
        #step_list = [alignment_pars,astrodrizzle_pars,catalog_generation_pars,quality_control_pars] # TODO: uncomment when everything is working
        step_list = [catalog_generation_pars] # TODO: Just a placeholder until we add complexity!
        for step_name in step_list:
            step_title = step_name.__name__.replace("_pars","").replace("_"," ")
            cfg_index = self.full_cfg_index[step_title]
            self.pars[step_title] = step_name(cfg_index,self.conditions,self.pars_dir,step_title,self.use_defaults)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def _determine_conditions(self):
        """Determine observing condition or conditions present for a given step
        List of possible conditions:
        * n_exp in mosaic >= 4
        * n_exp in mosaic < 4
        * total drizzle product
        * filter drizzle product
        * single drizzle product
        * Long/short

        """
        self.conditions = ["nexpGTE4"]  # TODO: Just a placeholder until we add complexity!


        pass

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def _get_cfg_index(self):
        """return the contents of the appropriate index cfg file."""
        code_dir = os.path.abspath(__file__)
        base_dir = os.path.dirname(os.path.dirname(code_dir))
        self.pars_dir = os.path.join(base_dir, "pars")
        cfg_index_fileanme = self.inst_det + "_index.config"
        cfg_index_filename = os.path.join(self.pars_dir, cfg_index_fileanme)

        with open(cfg_index_filename) as jsonFile:
            jsonString = jsonFile.read()
            self.full_cfg_index = json.loads(jsonString)
        pass


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

        """
        self.cfg_index = cfg_index
        self.conditions = conditions
        self.pars_dir = pars_dir
        self.step_title = step_title
        self.use_defaults = use_defaults
        self._get_params()



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def _get_params(self):
        """read in params from config files based on instrument, detector, and condition(s), and return a ordered
        dictionary of these values."""
        self.pars_multidict = collections.OrderedDict()
        found_cfg = False
        for condition in self.conditions:
            if condition in self.cfg_index.keys():
                found_cfg = True
                subcfgfilename = os.path.join(self.pars_dir, self.cfg_index[condition])
                configobj = teal.load(subcfgfilename)
                del configobj['_task_name_']
                self.pars_multidict[condition] = configobj
        if not found_cfg: # if no specific cfg files can be found for the specified conditions, use the generic cfg file.
            subcfgfilename = os.path.join(self.pars_dir, self.cfg_index["all"])
            configobj = teal.load(subcfgfilename,defaults=self.use_defaults,strict=False) # TODO: set strict=True for deployment
            del configobj['_task_name_']
            self.pars_multidict["all"] = configobj



#-----------------------------------------------------------------------------------------------------------------------


class alignment_pars(par):
    def __init__(self,cfg_index,conditions,pars_dir,step_title,use_defaults):
        """Configuration parameters for the image alignment step"""
        super().__init__(cfg_index,conditions,pars_dir,step_title,use_defaults)
        pass


#-----------------------------------------------------------------------------------------------------------------------


class astrodrizzle_pars(par):
    def __init__(self,cfg_index,conditions,pars_dir,step_title,use_defaults):
        """Configuration parameters for the AstroDrizzle step"""
        super().__init__(cfg_index,conditions,pars_dir,step_title,use_defaults)
        pass


#-----------------------------------------------------------------------------------------------------------------------


class catalog_generation_pars(par):
    def __init__(self,cfg_index,conditions,pars_dir,step_title,use_defaults):
        """Configuration parameters for the photometric catalog generation step"""
        super().__init__(cfg_index,conditions,pars_dir,step_title,use_defaults)
        self.outpars = {}
        if len(self.pars_multidict.keys()) == 1:
            for mdkey in self.pars_multidict.keys():
                for key in self.pars_multidict[mdkey].keys():
                    self.outpars[key] = self.pars_multidict[mdkey][key]




#-----------------------------------------------------------------------------------------------------------------------


class quality_control_pars(par):
    def __init__(self,cfg_index,conditions,pars_dir,step_title,use_defaults):
        """Configuration parameters for the quality control step"""
        super().__init__(cfg_index,conditions,pars_dir,step_title,use_defaults)
        pass



# ======================================================================================================================



if __name__ == '__main__':
    """Super simple testing interface for the above code."""

    import pdb

    foo = hap_config("acs", "wfc",use_defaults=False)
    blarg = foo.get_pars("catalog generation")

    print(blarg)


