#!/usr/bin/env python

"""This script contains code to create the complete set of configuration parameters required to run XXXXX given the
specified observation conditions and instrument/detector used in the observations"""

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
        pars_dir = os.path.join(base_dir, "pars")
        cfg_index_fileanme = self.inst_det + "_index.cfg"
        cfg_index_filename = os.path.join(pars_dir, cfg_index_fileanme)

        with open(cfg_index_filename) as jsonFile:
            jsonString = jsonFile.read()
            self.full_cfg_index = json.loads(jsonString)
        pass

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def _determine_conditions(self):
        """Determine observing condition or conditions present for a given step"""
        self.conditions = ["default"] #TODO: Just a placeholder until we add complexity!
        pass



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def get_pars(self):
        """This method returns the parameter set for a specified step (alignment, drizzle, etc.)"""

        pass


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def write_pars(self,out_filename):
        """This method writes the current parameter set to the specified file."""

        pass


#-----------------------------------------------------------------------------------------------------------------------


class par(object):
    def __init__(self):
        """INSIGHTFUL SUMMARY HERE"""
        pass


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def get_params(self):
        pass



#-----------------------------------------------------------------------------------------------------------------------


class alignment_pars(par):
    def __init__(self):
        """INSIGHTFUL SUMMARY HERE"""
        super().__init__()
        pass


#-----------------------------------------------------------------------------------------------------------------------


class astrodrizzle_pars(par):
    def __init__(self):
        """INSIGHTFUL SUMMARY HERE"""
        super().__init__()
        pass


#-----------------------------------------------------------------------------------------------------------------------


class catalog_generation_pars(par):
    def __init__(self):
        """INSIGHTFUL SUMMARY HERE"""
        super().__init__()
        pass


#-----------------------------------------------------------------------------------------------------------------------


class quality_control_pars(par):
    def __init__(self):
        """INSIGHTFUL SUMMARY HERE"""
        super().__init__()
        pass


# ======================================================================================================================



if __name__ == '__main__':
    """Super simple testing interface for the above code."""

    import pdb

    foo = hap_config("acs", "wfc")
    print(foo.full_cfg_index)

    pdb.set_trace()
