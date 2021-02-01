import pexpect
import os
import sys
import re
import time

acs_filter_dict = {'WFC': ['f435w', 'f475w', 'f555w', 'f606w', 'f625w',
                          'f775w', 'f814w', 'f850lp', 'f550m',
                          'f502n', 'f658n', 'f660n', 'f892n'],
                  'HRC': ['f220w', 'f250w', 'f330w', 'f435w', 'f550m',
                          'f555w', 'f606w', 'f625w', 'f775w', 'f814w',
                          'f344n', 'f658n', 'f660n', 'f892n', 'f502n'],
                  'SBC': ['f115lp', 'f122m', 'f125lp', 'f140lp', 'f150lp',
                                 'f165lp']}

wfc3_filter_dict = {'UVIS': ['f200lp', 'f218w', 'f225w', 'f275w', 'f336w',
                            'f350lp', 'f390m', 'f410m', 'f438w', 'f467m',
                            'f475w', 'f547m', 'f555w', 'f600lp', 'f606w',
                            'f621m', 'f625w', 'f689m', 'f763m', 'f775w',
                            'f814w', 'f845m', 'f850lp', 'fq422m',
                            'f280n', 'f343n', 'f373n', 'f395n',
                            'f469n', 'f487n', 'f502n', 'f631n',
                            'f645n', 'f656n', 'f657n', 'f658n',
                            'f665n', 'f673n', 'f680n', 'f953n',
                            'fq232n', 'fq243n', 'fq378n', 'fq387n',
                            'fq436n', 'fq437n', 'fq492n', 'fq508n',
                            'fq575n', 'fq619n', 'fq634n', 'fq672n',
                            'fq674n', 'fq727n', 'fq750n', 'fq889n',
                            'fq906n', 'fq924n', 'fq937n'],
                    'IR': ['f098m', 'f105w', 'f110w', 'f125w', 'f127m',
                          'f139m', 'f140w', 'f153m', 'f160w',
                          'f126n', 'f128n', 'f130n', 'f164n', 'f167n']}

FILTER_DICT = {'ACS': acs_filter_dict, 'WFC3': wfc3_filter_dict}

SPECTRA = {'short': '10000', 'long': '3000'}
DET_CHIP = 2
SPECTRA_TYPE = 2
FILTER_BREAK = 5


class BasePSF(object):
    pars = {}
    psf_size = 6.0
    prompts = ['Choice',
               'Position',
               'Filter',
               'Choose',
               'Temperature',
               'diameter',
               'Focus',
               'Rootname']
    instrument = None
    detector = None

    def __init__(self, filter_name, logfile=None, verbose=False):
        self.filter_name = filter_name.lower()
        self.parname = '_'.join([self.instrument, self.detector, filter_name]).lower()
        self.parfile = '{}.par'.format(self.parname)

        self.verbose = verbose
        self.logfile = logfile

        # interpret filter name
        self.filter_type = 'narrow' if filter_name.endswith('n') else 'wide'

        filter_indx = 2 if self.filter_name[1] == 'q' else 1
        if self.detector != 'IR':
            self.filter_range = 'short' if int(self.filter_name[filter_indx]) < FILTER_BREAK else 'long'
        else:
            self.filter_range = 'long'

        tiny_exe = pexpect.which('tiny1')
        if tiny_exe is None:
            raise FileNotFoundError("TinyTim executables NOT found! Please install first..")

        # Initialize output PSF names
        self.model_name = None
        self.psf_name = None


    def build_prompts(self):
        self.tiny_prompts = self.prompts.copy()
        if self.filter_type == 'narrow':
            del self.tiny_prompts[3]
            del self.tiny_prompts[3]


    def build_responses(self):
        """Generate answers to provide to tiny1 based on inputs"""
        self.responses = []
        self.custom_vals = {'Filter': self.filter_name,
                            'Choose': SPECTRA_TYPE,
                            'Temperature': SPECTRA[self.filter_range],
                            'diameter': self.psf_size,
                            'Rootname': self.parname}

        for p in self.tiny_prompts:
            val = None
            # Look for common responses first
            for parstr, par in self.pars.items():
                if parstr in p:
                    val = str(par)
                    break
            if val is None:
                for parstr, par in self.custom_vals.items():
                    if parstr in p:
                        val = str(par)
                        break
            self.responses.append(val)


    def tiny1(self, filename=None):
        """Call 'tiny1' based on defined inputs"""

        # Call 'tiny1' with provided inputs
        self.parname = filename if filename else self.parname

        # Start executable
        tiny_cmd = 'tiny1 {}'.format(self.parfile)
        print(tiny_cmd)
        print(self.responses)
        c = pexpect.spawn(tiny_cmd, maxread=4096)
        if self.verbose:
            c.logfile = sys.stdout.buffer

        # Now interact with the 'tiny1'
        for r in self.responses:
            c.expect(self.tiny_prompts)
            c.sendline(r)
            time.sleep(0.1)

    def tiny2(self):
        """Calls 'tiny2' and defines the name of the model PSF created"""
        tiny2_cmd = 'tiny2 {}'.format(self.parfile)
        if self.verbose:
            print(tiny2_cmd)
        pexpect.run(tiny2_cmd, encoding='utf-8')
        self.model_name = '{}00_psf.fits'.format(self.parfile)


    def tiny3(self):
        """ Calls 'tiny3' and returns the name of the PSF that was written out"""
        tiny3_cmd = 'tiny3 {}'.format(self.parfile)
        if self.verbose:
            print(tiny3_cmd)
        pexpect.run(tiny3_cmd, encoding='utf-8')
        psf_name = '{}00.fits'.format(self.parname)
        self.psf_name = psf_name.replace('00.fits', '.fits')
        os.rename(psf_name, self.psf_name)


class WFCPSF(BasePSF):
    pars = {'detector': 2, 'Choice': 15, 'Focus': -2.0, 'Position': "2048 1024"}
    psf_size = 25.0
    prompts = ['Choice',
               'detector',
               'Position',
               'Filter',
               'Choose',
               'Temperature',
               'diameter',
               'Focus',
               'Rootname']
    instrument = 'ACS'
    detector = 'WFC'

    def __init__(self, filter_name, logfile=None):
        super().__init__(filter_name, logfile=logfile)
        self.build_prompts()
        self.build_responses()


    def build_prompts(self):
        self.tiny_prompts = self.prompts.copy()
        if self.filter_type == 'narrow':
            del self.tiny_prompts[4]
            del self.tiny_prompts[4]

    def build_responses(self):
        """Generate answers to provide to tiny1 based on inputs"""
        self.responses = []
        self.custom_vals = {'Filter': self.filter_name,
                            'Enter detector': DET_CHIP,
                            'Choose': SPECTRA_TYPE,
                            'Temperature': SPECTRA[self.filter_range],
                            'diameter': self.psf_size,
                            'Rootname': self.parname}

        for p in self.tiny_prompts:
            val = None
            # Look for common responses first
            for parstr, par in self.pars.items():
                if parstr in p:
                    val = str(par)
                    break
            if val is None:
                for parstr, par in self.custom_vals.items():
                    if parstr in p:
                        val = str(par)
                        break
            self.responses.append(val)


class HRCPSF(BasePSF):
    pars = {'Choice': 16, 'Focus': -2.0, 'Position': "512 512"}
    psf_size = 6.0
    instrument = 'ACS'
    detector = 'HRC'

    def __init__(self, filter_name, logfile=None):
        super().__init__(filter_name, logfile=logfile)
        self.build_prompts()
        self.build_responses()


class SBCPSF(BasePSF):
    pars = {'Choice': 18, 'Focus': -2.0, 'Position': "512 512"}
    psf_size = 6.0
    instrument = 'ACS'
    detector = 'SBC'

    def __init__(self, filter_name, logfile=None):
        super().__init__(filter_name, logfile=logfile)
        self.build_prompts()
        self.build_responses()


class UVISPSF(BasePSF):
    pars = {'detector': 2, 'Choice': 22, 'Focus': -2.0, 'Position': "2048 1024"}
    psf_size = 25.0

    prompts = ['Choice',
               'detector',
               'Position',
               'Filter',
               'Choose',
               'Temperature',
               'diameter',
               'Focus',
               'Rootname']
    instrument = 'WFC3'
    detector = 'UVIS'

    def __init__(self, filter_name, logfile=None):
        super().__init__(filter_name, logfile=logfile)
        self.build_prompts()
        self.build_responses()


    def build_prompts(self):
        self.tiny_prompts = self.prompts.copy()
        if self.filter_type == 'narrow':
            del self.tiny_prompts[4]
            del self.tiny_prompts[4]

    def build_responses(self):
        """Generate answers to provide to tiny1 based on inputs"""
        self.responses = []
        self.custom_vals = {'Filter': self.filter_name,
                            'detector': DET_CHIP,
                            'Choose': SPECTRA_TYPE,
                            'Temperature': SPECTRA[self.filter_range],
                            'diameter': self.psf_size,
                            'Rootname': self.parname}

        for p in self.tiny_prompts:
            val = None
            # Look for common responses first
            for parstr, par in self.pars.items():
                if parstr in p:
                    val = str(par)
                    break
            if val is None:
                for parstr, par in self.custom_vals.items():
                    if parstr in p:
                        val = str(par)
                        break
            self.responses.append(val)


class IRPSF(BasePSF):
    pars = {'Choice': 23, 'Focus': -2.0, 'Position': "512 512"}
    psf_size = 25.0
    instrument = 'WFC3'
    detector = 'IR'

    def __init__(self, filter_name, logfile=None):
        super().__init__(filter_name, logfile=logfile)
        self.build_prompts()
        self.build_responses()


def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)


def run(parent_dir=None, logfile=None, verbose=False):
    """Create TinyTim PSFs using standardized parameters

    This function will automatically generate ALL PSFs listed in the
    global filter dictionaries in the module in a tree with a directory for each
    detector.

    """
    psf_classes = {('ACS', 'WFC'): WFCPSF,
                   ('ACS', 'HRC'): HRCPSF,
                   ('ACS', 'SBC'): SBCPSF,
                   ('WFC3', 'UVIS'): UVISPSF,
                   ('WFC3', 'IR'): IRPSF}

    if pexpect.which('tiny1') is None:
        raise FileNotFoundError("TinyTim executables NOT found! Please install first..")

    if parent_dir is None:
        parent_dir = os.getcwd()

    try:
        for instr, detectors in FILTER_DICT.items():
            instr_path = os.path.join(parent_dir, instr.lower())
            if not os.path.exists(instr_path):
                os.makedirs(instr_path)
            # move into that directory
            os.chdir(instr_path)
            # Now, let's start looping over the detectors
            for det, filter_names in detectors.items():
                det_path = os.path.join(instr_path, det.lower())
                if not os.path.exists(det_path):
                    os.makedirs(det_path)
                # move into that directory
                os.chdir(det_path)
                # Now, start creating the PSFs listed in the filter_dict for this detector
                for filter_name in filter_names:
                    # Create class for this PSF
                    psfobj = psf_classes[(instr, det)](filter_name,
                                                       logfile=logfile,
                                                       verbose=verbose)
                    psfobj.tiny1()
                    psfobj.tiny2()
                    psfobj.tiny3()
                # We are finished with creating all the PSFs for this detector...
                # Return to parent instrument directory
                os.chdir(instr_path)

            # We are finished with all the PSFs for this instrument...
            # return to the parent directory
            os.chdir(parent_dir)
    finally:
        os.chdir(parent_dir)
