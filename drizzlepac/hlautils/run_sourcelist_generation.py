""" script that directly calls sourcelist_generation.create_sourcelists() using an obs_info_dict hardcoded to
acs_10265_01 values. This will speed development of sourcelist_generation.py becasue it allows the user to just run it
without having to run runhlaprocessing.py first.

"""
import pickle
import sys

import sourcelist_generation

pickle_in = open(sys.argv[1], "rb")
obs_info_dict = pickle.load(pickle_in)
pickle_in.close()

sourcelist_generation.create_sourcelists(obs_info_dict)
