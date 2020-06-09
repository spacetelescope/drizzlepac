"""
return CI limits as a function of detector and filter using table



Dependencies
------------
ci_ap_cor_table_ap_20_2016.txt


"""
import os
import sys

from stsci.tools import logutil

__taskname__ = 'ci_table'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)

ci_table = None


def read_ci_apcorr_file(ci_lookup_file_path, diagnostic_mode=False, infile_name='ci_ap_cor_table_ap_20_2016.txt'):

    """Read the table of CI values

    Parameters
    ----------
    ci_lookup_file_path : string
        final path elements of the concentration index lookup file

    diagnostic_mode : Boolean, optional
        print CI info to screen (True/False)? Default value = False

    infile_name : string, optional
        Name of the CI file to use. If not specified, default value = ci_ap_cor_table_ap_20_2016.txt

    Returns
    -------
    ci_table : dicitonary
        CI information
    """
    code_dir = os.path.abspath(__file__)
    base_dir = os.path.dirname(os.path.dirname(code_dir))
    pars_dir = os.path.join(base_dir, "pars/hap_pars", ci_lookup_file_path)
    infile_name = os.path.join(pars_dir, infile_name)
    log.info("CI lookup table: {}".format(infile_name))
    ci_lines = open(infile_name).readlines()
    ci_table = {}
    for ci_line in ci_lines:
        if not ci_line.startswith("#"):
            ci_line = ci_line.strip()
            parse_cil = ci_line.split()
            if len(parse_cil) < 6:
                log.warning("Illegal line in {} (too few fields):\n{}".format(infile_name, ci_line))
                raise ValueError("Illegal line in %s (too few fields):\n%s" % (infile_name, ci_line))
            obs_config = parse_cil[0].upper()
            try:
                eff_wave = float(parse_cil[1])
                ci_lower = float(parse_cil[2])
                ci_peak = float(parse_cil[3])
                ci_upper = float(parse_cil[4])
                ap_corr = float(parse_cil[5])
            except ValueError as e:
                log.warning("Illegal line in {} (bad value):\n{}".format(infile_name, ci_line))
                raise ValueError("Illegal line in %s (bad value):\n%s" % (infile_name, ci_line))
            ci_table[obs_config] = (eff_wave, ci_lower, ci_peak, ci_upper, ap_corr)
    if diagnostic_mode:
        for key in sorted(ci_table.keys()):
            log.debug("{} {}".format(key, ci_table[key]))
    return ci_table


def get_ci_info(inst, detect, filt, ci_lookup_file_path, diagnostic_mode=False, eff_wave=-999, ci_lower=-999.0,
                ci_peak=-999.0, ci_upper=-999.0, ap_corr=-999.0):

    """Return dictionary with CI info for given detector/filter

    Parameters
    ----------
    inst : string
        Instrument name

    detect : string
        Detector name

    filt : string
        Filter name

    ci_lookup_file_path : string
        final path elements of the concentration index lookup file

    diagnostic_mode : Boolean
        Print information to screen (True/False)? If not specified, default value = False

    eff_wave : float
        Effective wavelength. If not specified, default value = -999.0

    ci_lower : float
        Lower CI limit. If not specified, default value = -999.0

    ci_peak : float
        CI peak. If not specified, default value = -999.0

    ci_upper : float
        Upper CI limit. If not specified, default value = -999.0

    ap_corr : float
        Aperture correction. If not specified, default value = -999.0

    Returns
    -------
    return_dict : dictionary
        A dictionary with CI info for given detector/filter
    """

    global ci_table
    if ci_table is None:
        ci_table = read_ci_apcorr_file(ci_lookup_file_path, diagnostic_mode=diagnostic_mode)

    obs_config = (("{}_{}_{}").format(inst, detect, filt)).upper()
    if diagnostic_mode:
        log.info("obs_config: {}".format(obs_config))
    if obs_config in ci_table:
        if diagnostic_mode:
            log.info("CI values found.")
        eff_wave, ci_lower, ci_peak, ci_upper, ap_corr = ci_table[obs_config]
    else:
        log.info("get_ci_info: CI values not found for %s, using default values" % obs_config)
    return_dict = {'eff_wave': eff_wave,
                   'ci_lower_limit': ci_lower,
                   'ci_peak': ci_peak,
                   'ci_upper_limit': ci_upper,
                   'ap_corr': ap_corr}
    return return_dict


def parse_file(drzfile):
    """Parse drizzled file name and return inst, detect, filt

    Parameters
    ----------
    drzfile : string
        name of drizzled fits file

    Returns
    -------
    instrument : string
        instrument

    detector : string
        detector

    filter : string
        filter
    """

    dir, fname = os.path.split(drzfile)
    f = fname.split('_')
    if len(f) < 6:
        log.error("Cannot parse inst/detect/filt from {}".format(drzfile))
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


def get_ci_from_file(drzfile, ci_lookup_file_path, log_level, **kw):

    """Return dictionary with CI info for given filename

    Parameters
    ----------
    drzfile : string
        name of drizzled fits file

    ci_lookup_file_path : string
        final path elements of the concentration index lookup file

    log_level : int
        The desired level of verboseness in the log statements displayed on the screen and written to the .log file.

    Returns
    -------
    a dictionary with CI info for given filename
    """
    # set logging level to user-specified level
    log.setLevel(log_level)

    inst, detect, filt = parse_file(drzfile)
    return get_ci_info(inst, detect, filt, ci_lookup_file_path, **kw)


if __name__ == '__main__':
    inst, detect, filt = ('wfc3', 'uvis', 'f606w')
    ciap_dict = get_ci_info(inst, detect, filt)
    log.info("configuration:  {}".format(inst, detect, filt))
    log.info("eff_wave:       {}".format(ciap_dict['eff_wave']))
    log.info("ci_lower:       {}".format(ciap_dict['ci_lower_limit']))
    log.info("ci_peak:        {}".format(ciap_dict['ci_peak']))
    log.info("ci_upper:       {}".format(ciap_dict['ci_upper_limit']))
    log.info("ap_corr:        {}".format(ciap_dict['ap_corr']))

    drizzled_image = 'hst_12311_03_wfc3_uvis_f275w_drz.fits'
    inst, detect, filt = parse_file(drizzled_image)
    ciap_dict = get_ci_from_file(drizzled_image)
    log.info("drizzled image: {}".format(drizzled_image))
    log.info("configuration:  {}".format(inst, detect, filt))
    log.info("eff_wave:       {}".format(ciap_dict['eff_wave']))
    log.info("ci_lower:       {}".format(ciap_dict['ci_lower_limit']))
    log.info("ci_peak:        {}".format(ciap_dict['ci_peak']))
    log.info("ci_upper:       {}".format(ciap_dict['ci_upper_limit']))
    log.info("ap_corr:        {}".format(ciap_dict['ap_corr']))

    drizzled_image = 'hst_11969_04_wfpc2_f606w_pc_drz.fits'
    inst, detect, filt = parse_file(drizzled_image)
    ciap_dict = get_ci_from_file(drizzled_image)
    log.info("drizzled image: {}".format(drizzled_image))
    log.info("configuration:  {}".format(inst, detect, filt))
    log.info("eff_wave:       {}".format(ciap_dict['eff_wave']))
    log.info("ci_lower:       {}".format(ciap_dict['ci_lower_limit']))
    log.info("ci_peak:        {}".format(ciap_dict['ci_peak']))
    log.info("ci_upper:       {}".format(ciap_dict['ci_upper_limit']))
    log.info("ap_corr:        {}".format(ciap_dict['ap_corr']))
