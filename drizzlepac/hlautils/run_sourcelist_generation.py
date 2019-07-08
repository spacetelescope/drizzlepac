""" script that directly calls sourcelist_generation.create_sourcelists() using an obs_info_dict hardcoded to
acs_10265_01 values. This will speed development of sourcelist_generation.py becasue it allows the user to just run it
without having to run runhlaprocessing.py first.

Test dir: /Users/dulude/Documents/HLAtransition/runhlaprocessing_testing/acs_10265_01

"""
import os
import pickle
import sys
import traceback
import sourcelist_generation


# ------------------------------------------------------------------------------
def confirm_execution():
    """
    This subroutine prevents accidental execution by requiring the user to type in a randomly generated 6-character
    confirmation string. If the string is typed in incorrectly, the script will simply exit to the command line.

    :return: Nothing
    """
    import string
    import random
    confirm_string = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(4))
    foo = input("Confirm execution by entering the following randomized text string: %s \n" % (confirm_string))
    if foo != confirm_string: sys.exit("Execution aborted.")
    if foo == confirm_string: print("Execution confirmed.")


# ------------------------------------------------------------------------------
print("\n"*25)
print("Current path: ",os.getcwd())
confirm_execution()
try:
    cmd = "rm -f *.*"
    print(cmd)
    os.system(cmd)

    cmd = "cp sl_gen_orig/* ."
    print(cmd)
    os.system(cmd)

    pickle_in = open(sys.argv[1], "rb")
    [obs_info_dict,param_dict] = pickle.load(pickle_in)
    pickle_in.close()

    #copy Warren's custom-created f606W and total drizzle images in for testing!
    # cmdlist = ["cp -f artif_orig/hst_10265_01_acs_wfc_f606w_j92c01_drc.fits ."]
    # cmdlist.append("cp -f artif_orig/hst_10265_01_acs_wfc_total_j92c01_drc.fits .")
    # for cmd in cmdlist:
    #     print(cmd)
    #     os.system(cmd)

    param_dict["ACS WFC"]["sourcex"]["fwhm"] = 0.13
    param_dict["ACS WFC"]["sourcex"]["source_box"] = 5
    param_dict["ACS WFC"]["sourcex"]["thresh"] = None
    sourcelist_generation.run_create_sourcelists(obs_info_dict,param_dict)
except Exception:
    print("\a\a\a")
    exc_type, exc_value, exc_tb = sys.exc_info()
    traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)
    print("Error! exiting...")
print("\a")