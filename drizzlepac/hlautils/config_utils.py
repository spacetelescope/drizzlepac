#!/usr/bin/env python

"""This script contains code to create the complete set of configuration parameters required to run XXXXX given the
specified observation conditions and instrument/detector used in the observations"""

import collections
import json
import os


from stsci.tools import teal
# ======================================================================================================================


class hap_config(object):
    def __init__(self, instrument,detector,cfg_index_file=None):
        """
        A set of routines to generate appropriate set of configuration parameters

        Parameters
        ----------
        instrument : str
            instrument name

        detector : str
            detector name

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
        self._determine_conditions()
        self._get_cfg_index()

        # Instantiate the parameter set
        self.pars = {}
        #step_list = [alignment_pars,astrodrizzle_pars,catalog_generation_pars,quality_control_pars] # TODO: uncomment when everything is working
        step_list = [catalog_generation_pars] # TODO: Just a placeholder until we add complexity!
        for step_name in step_list:
            step_title = step_name.__name__.replace("_pars","").replace("_"," ")
            self.pars[step_title] = step_name(self.conditions,self.inst_det,self.full_cfg_index[step_title],self.pars_dir,step_title)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def _determine_conditions(self):
        """Determine observing condition or conditions present for a given step"""
        self.conditions = ["default"]  # TODO: Just a placeholder until we add complexity!
        pass

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def _get_cfg_index(self):
        """return the contents of the appropriate index cfg file."""
        code_dir = os.path.abspath(__file__)
        base_dir = os.path.dirname(os.path.dirname(code_dir))
        self.pars_dir = os.path.join(base_dir, "pars")
        cfg_index_fileanme = self.inst_det + "_index.cfg"
        cfg_index_filename = os.path.join(self.pars_dir, cfg_index_fileanme)

        with open(cfg_index_filename) as jsonFile:
            jsonString = jsonFile.read()
            self.full_cfg_index = json.loads(jsonString)
        pass


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def get_pars(self,step_name):
        """This method returns the parameter set for a specified step (alignment, drizzle, etc.)"""

        return(self.pars[step_name].outpars)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def write_pars(self,out_filename):
        """This method writes the current parameter set to the specified file."""

        pass


#-----------------------------------------------------------------------------------------------------------------------


class par():
    def __init__(self,conditions,inst_det,cfg_index,pars_dir,step_title):
        """INSIGHTFUL SUMMARY HERE"""
        print("---->>>>>",conditions,inst_det,cfg_index,pars_dir,step_title)
        self.conditions = conditions
        self.inst_det = inst_det
        self.cfg_index = cfg_index
        self.pars_dir = pars_dir
        self.step_title = step_title
        self._get_params()



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def _get_params(self):
        """read in params from config files based on instrument, detector, and condition(s), and return a ordered
        dictionary of these values."""
        self.pars_multidict = collections.OrderedDict()
        for condition in self.conditions:
            subcfgfilename = os.path.join(self.pars_dir, self.cfg_index[condition])
            print(condition, subcfgfilename)
            configobj = teal.load(subcfgfilename)
            del configobj['_task_name_']
            self.pars_multidict[condition] = configobj


#-----------------------------------------------------------------------------------------------------------------------


class alignment_pars(par):
    def __init__(self,conditions,inst_det,cfg_index,pars_dir,step_title):
        """INSIGHTFUL SUMMARY HERE"""
        super().__init__(conditions,inst_det,cfg_index,pars_dir,step_title)
        pass


#-----------------------------------------------------------------------------------------------------------------------


class astrodrizzle_pars(par):
    def __init__(self,conditions,inst_det,cfg_index,pars_dir,step_title):
        """INSIGHTFUL SUMMARY HERE"""
        super().__init__(conditions,inst_det,cfg_index,pars_dir,step_title)
        pass


#-----------------------------------------------------------------------------------------------------------------------


class catalog_generation_pars(par):
    def __init__(self,conditions,inst_det,cfg_index,pars_dir,step_title):
        """INSIGHTFUL SUMMARY HERE"""
        super().__init__(conditions,inst_det,cfg_index,pars_dir,step_title)
        if len(self.pars_multidict.keys()) == 1:
            self.outpars = {}
            for mdkey in self.pars_multidict.keys():
                for key in self.pars_multidict[mdkey].keys():
                    self.outpars[key] = self.pars_multidict[mdkey][key]



#-----------------------------------------------------------------------------------------------------------------------


class quality_control_pars(par):
    def __init__(self,conditions,inst_det,cfg_index,pars_dir,step_title):
        """INSIGHTFUL SUMMARY HERE"""
        super().__init__(conditions,inst_det,cfg_index,pars_dir,step_title)
        pass



# ======================================================================================================================



if __name__ == '__main__':
    """Super simple testing interface for the above code."""

    import pdb

    foo = hap_config("acs", "wfc")
    blarg = foo.get_pars("catalog generation")

    pdb.set_trace()
