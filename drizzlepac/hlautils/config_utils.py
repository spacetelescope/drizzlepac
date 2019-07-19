#!/usr/bin/env python

"""This script contains code to create the complete set of configuration parameters required to run XXXXX given the
specified observation conditions and instrument/detector used in the observations"""


# ======================================================================================================================


class hap_config(object):
    def __init__(self, condition, inst_det):
        """
        A set of routines to generate appropriate set of configuration parameters

        Parameters
        ----------
        condition : string
            observation condition

        inst_det : string
            instrument/detector combined together as a single space-separated lower-case string (e.g. wfc3 uvis)
        """
        self.label = "hap_config"
        self.description = "A set of routines to generate appropriate set of configuration parameters"
        self.condition = condition
        self.inst_det = inst_det


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def _astrodrizzle_param_setup(self):
        """Merge settings from the mdriztab, condition-specific and instrument/detector-specific parameter sets, and
        any additional user inputs into a single combined set of configuration parameters for Astrodrizzle"""
        pass


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def _get_condition_specific_params(self):
        """Get observation condition-specific HAP parameters"""
        pass


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def _get_inst_det_specific_params(self):
        """Get instrument/detector-specific HAP parameters"""
        pass


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



    def _get_user_params(self):
        """Add user-specified HAP parameters to the co"""
        pass


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def get_hap_parameters(self):
        """returns parameters for a user-specified section or section/subsection of the HAP config parameter set"""
        pass

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    def hap_config(self):
        """Merge all instrument/detector-specific, condition-specific and user-specified parameters into a the final
        set of input parameters used by XXXX."""

        pass




