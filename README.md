# Drizzlepac

[![Build Status](https://dev.azure.com/spacetelescope/drizzlepac/_apis/build/status/spacetelescope.drizzlepac?branchName=master)](https://dev.azure.com/spacetelescope/drizzlepac/_build/latest?definitionId=2&branchName=master)
[![Build Status](https://ssbjenkins.stsci.edu/job/STScI/job/drizzlepac/job/master/badge/icon)](https://ssbjenkins.stsci.edu/job/STScI/job/drizzlepac/job/master/)
[![Documentation Status](https://readthedocs.org/projects/drizzlepac/badge/?version=latest)](http://drizzlepac.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/spacetelescope/drizzlepac/branch/master/graph/badge.svg)](https://codecov.io/gh/spacetelescope/drizzlepac)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3743274.svg)](https://doi.org/10.5281/zenodo.3743274)

Nightly regression test results are available only from within the STScI network at this time.
https://plwishmaster.stsci.edu:8081/job/RT/job/drizzlepac/

The use of this software on HST data is described at:

    http://drizzlepac.stsci.edu/

A complete description of the documented interfaces in the code itself 
can be found at:

    http://drizzlepac.readthedocs.io


# Installation

## Conda (Recommended)

```bash
$ conda config --add channels http://ssb.stsci.edu/astroconda
$ conda create -n astroconda stsci
```

## From Source

### Clone this repository
```bash
$ git clone https://github.com/spacetelescope/drizzlepac
$ cd drizzlepac
```

### Build the documentation

*Note:* If you intend to use `drizzlepac`'s embedded help feature from within
an interactive `python` or `ipython` session, we recommend you do not skip
this step.

```bash
$ python setup.py build_sphinx
```

### Install drizzlepac

```bash
$ python setup.py install
```
