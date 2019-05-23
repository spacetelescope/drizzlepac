"""
return CI limits as a function of detector and filter using table



Dependencies
------------
ci_ap_cor_table_ap_20_2016.txt


"""
import os, sys
import pdb

from stsci.tools import logutil

__taskname__ = 'ci_table'

log = logutil.create_logger(__name__, level=logutil.logging.INFO, stream=sys.stdout)

ci_table = None

def read_ci_apcorr_file(debug=False, infile_name='ci_ap_cor_table_ap_20_2016.txt'):

    """
    Read the table of CI values
    
    :param debug: print CI info to screen (True/False)? Default value = False
    :param infile_name: Name of the CI file to use. If not specified, default value = ci_ap_cor_table_ap_20_2016.txt
    :type debug: Boolean
    :type infile_name: string
    :return: CI information
    """

    # assume data file is in same directory with this module
    if not os.path.exists(infile_name):
        localdir = os.path.split(__file__)[0]
        tfile = os.path.join(localdir, infile_name)
        if not os.path.exists(tfile):
            raise ValueError("Cannot find CI table %s" % infile_name)
        infile_name = tfile
    ci_lines = open(infile_name).readlines()

    ci_table={}
    for ci_line in ci_lines:
        if ci_line.startswith("#") == False:
            ci_line=ci_line.strip()
            parse_cil=ci_line.split()
            if len(parse_cil) < 6:
                raise ValueError("Illegal line in %s (too few fields):\n%s" % (infile_name, ci_line))
            obs_config = parse_cil[0].upper()
            try:
                eff_wave = float(parse_cil[1])
                ci_lower = float(parse_cil[2])
                ci_peak = float(parse_cil[3])
                ci_upper = float(parse_cil[4])
                ap_corr = float(parse_cil[5])
            except ValueError as e:
                raise ValueError("Illegal line in %s (bad value):\n%s" % (infile_name, ci_line))
            ci_table[obs_config] = (eff_wave, ci_lower, ci_peak, ci_upper, ap_corr)
    if debug:
        for key in sorted(ci_table.keys()):
            print(key,ci_table[key])
    return ci_table


def get_ci_info(inst, detect, filt, debug=False,
        eff_wave=-999, ci_lower=-999.0, ci_peak=-999.0, ci_upper=-999.0, ap_corr=-999.0):

    """
    Return dictionary with CI info for given detector/filter
    
    :param inst: Instrument name
    :param detect: Detector name
    :param filt: Filter name
    :param debug: Print information to screen (True/False)? If not specified, default value = False
    :param eff_wave: Effective wavelength. If not specified, default value = -999.0
    :param ci_lower: Lower CI limit. If not specified, default value = -999.0
    :param ci_peak: CI peak. If not specified, default value = -999.0
    :param ci_upper: Upper CI limit. If not specified, default value = -999.0
    :param ap_corr: Aperture correction. If not specified, default value = -999.0
    :type inst: string
    :type detect: string
    :type filt: string
    :type debug: Boolean 
    :type eff_wave: float
    :type ci_lower: float
    :type ci_peak: float
    :type ci_upper: float
    :type ap_corr: float
    :return: A dictionary with CI info for given detector/filter
    """

    global ci_table
    if ci_table is None:
        ci_table = read_ci_apcorr_file(debug=debug)

    obs_config=("%s_%s_%s"%(inst,detect,filt)).upper()
    if debug: print("obs_config: ",obs_config)
    if obs_config in ci_table:
        if debug: print("CI values found.")
        eff_wave, ci_lower, ci_peak, ci_upper, ap_corr = ci_table[obs_config]
    else:
        print("get_ci_info: CI values not found for %s, using default values" % obs_config)
    return_dict={'eff_wave':eff_wave,'ci_lower_limit':ci_lower,'ci_peak':ci_peak,'ci_upper_limit':ci_upper,'ap_corr':ap_corr}
    return return_dict


def parse_file(drzfile):
    """
    Parse drizzled file name and return inst, detect, filt
    
    :param drzfile: name of drizzled fits file
    :type drzfile: string
    :return: instrument, detector, filter
    """

    dir, fname = os.path.split(drzfile)
    f = fname.split('_')
    if len(f) < 6:
        raise ValueError("Cannot parse inst/detect/filt from %s" % drzfile)
    inst = f[3].upper()
    if inst == 'WFPC2':
        detect = f[4].upper()[0:2]
    if inst != 'WFPC2':
        detect = f[4].upper()
    filt = f[5].upper()
    # if inst == 'WFPC2':
    #     # wfpc2 detector may be at end instead of after inst
    #     # this should work for old names or new names
    #     if filt in ("WF","PC"):
    #         tmp = filt
    #         filt = detect
    #         detect = tmp
    #     if detect not in ("WF","PC"):
    #         raise ValueError("Cannot parse WFPC2 info from %s" % drzfile)
    return (inst, detect, filt)


def get_ci_from_file(drzfile, **kw):

    """
    Return dictionary with CI info for given filename
    
    :param drzfile: name of drizzled fits file
    :param kw: input parameters for **get_ci_info**
    :type drzfile: string
    :type kw: various
    :return: a dictionary with CI info for given filename
    """

    inst, detect, filt = parse_file(drzfile)
    return get_ci_info(inst, detect, filt, **kw)

if __name__=='__main__':
    inst, detect, filt = ('wfc3','uvis','f606w')
    ciap_dict = get_ci_info(inst,detect,filt)
    print("configuration:  ",inst,detect,filt)
    print("eff_wave:       ",ciap_dict['eff_wave'])
    print("ci_lower:       ",ciap_dict['ci_lower_limit'])
    print("ci_peak:        ",ciap_dict['ci_peak'])
    print("ci_upper:       ",ciap_dict['ci_upper_limit'])
    print("ap_corr:        ",ciap_dict['ap_corr'])

    drizzled_image='hst_12311_03_wfc3_uvis_f275w_drz.fits'
    inst, detect, filt = parse_file(drizzled_image)
    ciap_dict = get_ci_from_file(drizzled_image)
    print("drizzled image: ",drizzled_image)
    print("configuration:  ",inst,detect,filt)
    print("eff_wave:       ",ciap_dict['eff_wave'])
    print("ci_lower:       ",ciap_dict['ci_lower_limit'])
    print("ci_peak:        ",ciap_dict['ci_peak'])
    print("ci_upper:       ",ciap_dict['ci_upper_limit'])
    print("ap_corr:        ",ciap_dict['ap_corr'])

    drizzled_image='hst_11969_04_wfpc2_f606w_pc_drz.fits'
    inst, detect, filt = parse_file(drizzled_image)
    ciap_dict = get_ci_from_file(drizzled_image)
    print("drizzled image: ",drizzled_image)
    print("configuration:  ",inst,detect,filt)
    print("eff_wave:       ",ciap_dict['eff_wave'])
    print("ci_lower:       ",ciap_dict['ci_lower_limit'])
    print("ci_peak:        ",ciap_dict['ci_peak'])
    print("ci_upper:       ",ciap_dict['ci_upper_limit'])
    print("ap_corr:        ",ciap_dict['ap_corr'])
