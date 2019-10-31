import os
import sys
from drizzlepac import hapsequencer
from drizzlepac.devutils.confirm_execution import confirm_execution
#----------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    print("Current working directory: "+os.getcwd())
    confirm_execution()

    cmd_list = ['rm -f *.*','cp orig/* .']
    for cmd in cmd_list:
        print(cmd)
        os.system(cmd)
    out_pars_file = sys.argv[1].replace(".out", "_config.json")
    x = hapsequencer.run_hap_processing(sys.argv[1], debug=True, output_custom_pars_file=out_pars_file, phot_mode='aperture')
    print("\a\a\a")